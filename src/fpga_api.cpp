#include "fpga_api.h"
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <cstring>
#include <cmath>


#define min(x, y) (((x) < (y)) ? (x) : (y))

FPGA::FPGA(off_t data_addr, off_t output_addr, int m_size, int v_size)
{
  m_size_ = m_size;
  v_size_ = v_size;
  data_size_ = (m_size_ + 1) * v_size_; // fpga bram data size

  qvec_ = new char[v_size_];
  qmat_ = new char[m_size_*v_size_];

  m1_size_ = v_size * v_size;
  m2_size_ = v_size * v_size;
  data_size_M = (v_size_+v_size_)*v_size_;
  
  qm1_ = new char[v_size_*v_size_];
  qm2_ = new char[v_size_*v_size_];
  
  qout_ = new int[m_size_];
  qout_M = new int[v_size_*v_size_];

  output_ = new unsigned int[m_size_]; // use output_ as tempolar output
  output_M = new unsigned int[v_size_*v_size_]; // use output_M as tempolar output

  data_ = new float[data_size_];
  data_M = new float[data_size_M];

  fd_ = open("/dev/mem", O_RDWR);

  qdata_ = new char[data_size_];
  qdata_M = static_cast<char *>(mmap(NULL, data_size_M, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, data_addr));
  
  output_ = static_cast<unsigned int *>(mmap(NULL, sizeof(unsigned int), PROT_READ | PROT_WRITE, MAP_SHARED, fd_, output_addr));
  num_block_call_ = 0;
}

FPGA::~FPGA()
{
  munmap(qdata_M, data_size_);
  munmap(output_, sizeof(unsigned int));
  close(fd_);

  delete[] output_;
  delete[] data_;
  delete[] qvec_;
  delete[] qmat_;
  delete[] qout_;
}

char *FPGA::qmatrix(void)
{
  return qdata_ + v_size_;
}

char *FPGA::qvector(void)
{
  return qdata_;
}

char *FPGA::qmatrix_M1(void)
{
  return qdata_M;
}

char *FPGA::qmatrix_M2(void)
{
  return qdata_M + m1_size_;
}

float *FPGA::matrix(void)
{
  return data_ + v_size_;
}

float *FPGA::vector(void)
{
  return data_;
}

float *FPGA::matrix_M1(void)
{
  return data_M;
}

float *FPGA::matrix_M2(void)
{
  return data_M + m1_size_;
}

void FPGA::reset(void)
{
  num_block_call_ = 0;
}

int FPGA::num_block_call(void)
{
  return num_block_call_;
}

void quantize(const float* input, char* quantized, int num_input, int bits_min, int bits_max, char offset, float scale)
{
  for(int i = 0; i < num_input; i++)
  {
    quantized[i] = (char)ceil(input[i]/scale); // TODO: convert floating point to quantized value
  }
}

void dequantize(const int* quantized, float* output, int num_output, int offset, float scale)
{
  for(int i = 0; i < num_output; i++)
  {
    output[i] = scale * (quantized[i]); // TODO: convert quantized value to floating point
  }
}

const int *__attribute__((optimize("O0"))) FPGA::qblockMM(Compute* comp)
{
  num_block_call_ += 1;

  // fpga version
  *output_ = 0x5555;
  while (*output_ == 0x5555)
    ;

  return reinterpret_cast<int *>(qdata_);
}

// const float* FPGA::blockMM(Compute* comp)
// {
  // num_block_call_ += 1;

  // // cpu version
  // int* m1 = this->qmatrix_M1();
  // int* m2 = this->qmatrix_M2();
  // float* out  = reinterpret_cast<float*>(output_M);  

  // if(comp->quantized)
  // {
  //   char act_bits_min = 0;
  //   char act_bits_max = (1<<(comp->act_bits-1))-1;

  //   float act_scale = 0; // TODO calculate the scale factor
  //   char act_offset = 0; // TODO calculate the zero-offset
  //   quantize(); // TODO complete quantize function

  //   char weight_bits_min = 0;
  //   char weight_bits_max = (1<<(comp->weight_bits-1))-1;

  //   float weight_scale = 0; // TODO calculate the scale factor
  //   char weight_offset = 0; // TODO calculate the zero-offset
  //   quantize(); // TODO complete quantize function

  //   for(int i = 0; i < v_size_; ++i)
  //   {
  //     for(int j = 0; j < v_size_; ++j){    
  //       qout_M[v_size_*i+j] = 0;
  //       for(int k = 0; k < v_size_; ++k){
  //         qout_M[v_size_*i+j] += m1[v_size_*i+k] * m2[v_size_*k + j];
  //       }
  //     }
  //   }
  //   dequantize(); // TODO complete dequantize function

  // }
  // else{
  //   for(int i = 0; i < v_size_; ++i)
  //   {
  //     for(int j = 0; j < v_size_; ++j){    
  //       out[v_size_*i+j] = 0;
  //       for(int k = 0; k < v_size_; ++k){
  //         out[v_size_*i+j] += m1[v_size_*i+k] * m2[v_size_*k + j];
  //       }
  //     }
  //   }
  // }

  // for(int i = 0; i < m1_size_; ++i)
  //   data_M[i] = out[i];

  // return data_M;    
// }

const float *FPGA::blockMV(Compute* comp)
{
  num_block_call_ += 1;

  // cpu version
  float *vec = this->vector();
  float *mat = this->matrix();
  float *out = reinterpret_cast<float *>(output_M);

  if(comp->quantized)
  {
    char act_bits_min = 0;
    char act_bits_max = (1<<(comp->act_bits-1))-1;

    float act_scale = (comp->act_max - comp->act_min) / act_bits_max; // TODO calculate the scale factor
    char act_offset = (char)ceil(-comp->act_min/act_scale); // TODO calculate the zero-offset
    quantize(vec, qvec_, v_size_, act_bits_min, act_bits_max, act_offset, act_scale); // TODO complete quantize function

    char weight_bits_min = 0;
    char weight_bits_max = (1<<(comp->weight_bits-1))-1;


    float weight_scale = (comp->weight_max - comp->weight_min) / weight_bits_max; // TODO calculate the scale factor
    char weight_offset = (char)ceil(-comp->weight_min/weight_scale); // TODO calculate the zero-offset
    quantize(mat, qmat_, m_size_*v_size_, weight_bits_min, weight_bits_max, weight_offset, weight_scale); // TODO complete quantize function

    for (int i = 0; i < m_size_; ++i)
    {
      qout_[i] = 0;
      for (int j = 0; j < v_size_; ++j)
        qout_[i] += (qvec_[j]) * (qmat_[v_size_ * i + j]);
        // qout_[i] += (qvec_[j]-act_offset) * (qmat_[v_size_ * i + j]-weight_offset);
    }
    int out_offset = 0;
    float out_scale = act_scale * weight_scale;
    dequantize(qout_, out, m_size_, out_offset, out_scale); // TODO complete dequantize function
  }
  else
  {
    for (int i = 0; i < m_size_; ++i)
    {
      out[i] = 0;
      for (int j = 0; j < v_size_; ++j)
        out[i] += vec[j] * mat[v_size_ * i + j];
    }
  }

  for (int i = 0; i < m_size_; ++i)
    data_[i] = out[i];

  return data_;
}

void FPGA::largeMM(const float* weight_mat, const float* input_mat, float* output, int num_input, int num_output, int num_matrix2, Compute* comp)
{
  float* m1 = this->matrix_M1();
  float* m2 = this->matrix_M2();

  char* qm1 = this->qmatrix_M1();
  char* qm2 = this->qmatrix_M2();

  char act_bits_min = 0;
  char act_bits_max = (1<<(comp->act_bits-1))-1;

  float act_scale = (comp->act_max - comp->act_min) / act_bits_max; // TODO calculate the scale factor
  char act_offset = (char)ceil(-comp->act_min/act_scale); // TODO calculate the zero-offset

  char weight_bits_min = 0;
  char weight_bits_max = (1<<(comp->weight_bits-1))-1;

  float weight_scale = (comp->weight_max - comp->weight_min) / weight_bits_max; // TODO calculate the scale factor
  char weight_offset = (char)ceil(-comp->weight_min/weight_scale); // TODO calculate the zero-offset

  // 0) Initialize output vector		
  for(int i = 0; i < num_output*num_matrix2; ++i)
    output[i] = 0;

  for(int i = 0; i < num_output; i += v_size_)
  {
    for(int j = 0; j < num_input; j += v_size_)
    {			
      for(int k = 0; k < num_matrix2; k += v_size_)
      {
        // 0) Initialize input vector
        int block_row = min(v_size_, num_output-i);
        int block_col_1 = min(v_size_, num_input-j);
        int block_col_2 = min(v_size_, num_matrix2-k);

        // 1) Assign a m1
        // IMPLEMENT THIS
        printf("Asdfasdfadsfa\n");
        memset(m1, 0, sizeof(float) * v_size_ * v_size_);
        printf("Asdfasdfadsfa\n");
        for (int row = 0; row < block_row; row++)
        {
          memcpy(m1 + row * v_size_, weight_mat + (i + row) * num_input + j, sizeof(float)*block_col_1);
        }
        printf("Asdfasdfadsfa\n");

        // 2) Assign a m2
        // IMPLEMENT THIS
        printf("Asdfasdfadsfa\n");
        memset(m2, 0, sizeof(float) * v_size_ * v_size_);
        printf("Asdfasdfadsfa\n");
        for (int row = 0; row < block_col_1; row++)
        {
          memcpy(m2 + row * v_size_, input_mat + (j + row) * num_matrix2 + k, sizeof(float)*block_col_2);
        }
        printf("Asdfasdfadsfa\n");

        quantize(m2, qm2, m2_size_, act_bits_min, act_bits_max, act_offset, act_scale);
        printf("Asdfasdfadsfa\n");
        quantize(m1, qm1, m1_size_, weight_bits_min, weight_bits_max, weight_offset, weight_scale);

        printf("Asdfasdfadsfa\n");
        for (int x=0; x<8; x++) {
          for(int y=0; y<8; y++) {
            printf("%-5d", qm1[x*8+y]);
          }
          printf("\n");
        }
        printf("\n");
        for (int x=0; x<8; x++) {
          for(int y=0; y<8; y++) {
            printf("%-5d", qm2[x*8+y]);
          }
          printf("\n");
        }
        printf("\n");

        // 3) Call a function `blockMM() to execute Matrix matrix multiplication
        const int* ret = this->qblockMM(comp);

        for (int x=0; x<8; x++) {
          for(int y=0; y<8; y++) {
            printf("%-5d", ret[x*8+y]);
          }
          printf("\n");
        }
        printf("\n");

        dequantize(ret, data_M, v_size_*v_size_, 0, weight_scale*act_scale);

        // 4) Accumulate intermediate results
        for(int n = 0; n<block_row; ++n)
        {
          for(int m = 0; m<block_col_2; ++m)
          {
            output[(i + n) + (k + m)*num_output] += data_M[n*v_size_ + m];
          }
        }
      }
    } 
  }
}

void FPGA::largeMV(const float *large_mat, const float *input, float *output, int num_input, int num_output, Compute* comp)
{
  float *vec = this->vector();
  float *mat = this->matrix();

  // 0) Initialize output vector
  for (int i = 0; i < num_output; ++i)
    output[i] = 0;

  for (int i = 0; i < num_output; i += m_size_)
  {
    for (int j = 0; j < num_input; j += v_size_)
    {
      // 0) Initialize input vector
      int block_row = min(m_size_, num_output - i);
      int block_col = min(v_size_, num_input - j);

      // 1) Assign a vector
      // IMPLEMENT THIS
      memcpy(vec, &input[j], sizeof(float) * block_col);
      if (v_size_ - block_col > 0)
        memset(vec + block_col, 0, sizeof(float) * (v_size_ - block_col));

      // 2) Assign a matrix
      // IMPLEMENT THIS
      memset(mat, 0, sizeof(float) * v_size_ * m_size_);
      for (int row = 0; row < block_row; row++)
      {
        memcpy(&mat[row * v_size_], &large_mat[(i + row) * num_input + j], sizeof(float) * block_col);
      }

      // 3) Call a function `blockMV() to execute MV multiplication
      const float* ret = this->blockMV(comp);

      // 4) Accumulate intermediate results
      for (int row = 0; row < block_row; ++row)
        output[i + row] += ret[row];
    }
  }
}

