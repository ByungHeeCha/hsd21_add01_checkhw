#include "fpga_api.h"
#include <cstdio>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#define DATA_SIZE 2*SIZE*(SIZE)*sizeof(float) // fpga bram data size

#define min(x,y) (((x)<(y))?(x):(y))

FPGA::FPGA(off_t data_addr, off_t api_addr)
{
    fd_ = open("/dev/mem", O_RDWR);
    data_ = static_cast<float*>(mmap(NULL, DATA_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, fd_, data_addr));
    api_ = static_cast<unsigned int*>(mmap(NULL, sizeof(unsigned int), PROT_READ|PROT_WRITE, MAP_SHARED,fd_, api_addr));
}

FPGA::~FPGA()
{
    munmap(data_, DATA_SIZE );
    munmap(api_, sizeof(unsigned int));
    close(fd_);
}

float* FPGA::matrix_M1(void)
{
	return data_ ;
}

float* FPGA::matrix_M2(void)
{
	return data_ + SIZE * SIZE;
}

const float* __attribute__((optimize("O0"))) FPGA::run()
{
    *api_ = 0x5555;
    while(*api_ == 0x5555);

    return data_;    
}

// const float* FPGA::real()
// {
//   float* real_ = new float[8*8];
//   float* m1 = this->matrix_M1();
// 	float* m2 = this->matrix_M2();
//   for(int aaa=0; aaa<8; aaa++) {
//     for(int bbb=0; bbb<8; bbb++) {
//       real_[aaa*8+bbb] = 0;
//     }
//   }

//   for(int aaa=0; aaa<8;aaa++) {
//     for(int bbb=0; bbb<8; bbb++) {
//       for(int kkk=0; kkk<8; kkk++) {
//         real_[aaa*8+bbb] += m1[aaa*8+kkk] * m2[kkk*8+bbb];
//       }
//     }
//   }
//   return real_;
// }

// Test code for bitstream
void FPGA::largeMM(const float* weight_mat, const float* input_mat, float* output, 
							int num_input, int num_output, int num_matrix2)
{
	float* m1 = this->matrix_M1();
	float* m2 = this->matrix_M2();
	for(int i = 0; i < num_output*num_matrix2; ++i)
    output[i] = 0;
  
  for(int i=0; i < num_output; i++) {
    for(int j=0; j<num_input; j++) {
      printf("%f ", weight_mat[i*num_input+j]);
    }
    printf("\n");
  }
  printf("\n");

  for(int i=0; i<num_input; i++) {
    for(int j=0; j<num_matrix2; j++) {
      printf("%f ", input_mat[i*num_matrix2+j]);
    }
    printf("\n");
  }
  printf("\n");

  for(int i = 0; i < num_output; i += SIZE)
  {
    for(int j = 0; j < num_input; j += SIZE)
    {			
      for(int k = 0; k < num_matrix2; k += SIZE)
      {
        // 0) Initialize input vector
        int block_row = min(SIZE, num_output-i);
        int block_col_1 = min(SIZE, num_input-j);
        int block_col_2 = min(SIZE, num_matrix2-k);

        // 1) Assign a m1
        // Implement This
        memset(m1, 0, sizeof(float) * SIZE * SIZE);
        for (int row = 0; row < block_row; row++)
        {
          memcpy(m1 + row * SIZE, weight_mat + (i + row) * num_input + j, block_col_1*sizeof(float));
        }

        for(int aaa=0; aaa<8; aaa++) {
          for(int bbb=0; bbb<8; bbb++) {
            printf("%f ", m1[aaa*8+bbb]);
          }
          printf("\n");
        }
        printf("\n");

        // 2) Assign a m2
        // IMPLEMENT THIS
        memset(m2, 0, sizeof(float) * SIZE * SIZE);
        for (int row = 0; row < block_col_1; row++)
        {
          memcpy(m2 + row * SIZE, input_mat + (j + row) * num_matrix2 + k, block_col_2*sizeof(float));
        }
        for(int aaa=0; aaa<8; aaa++) {
          for(int bbb=0; bbb<8; bbb++) {
            printf("%f ", m2[aaa*8+bbb]);
          }
          printf("\n");
        }
        printf("\n");


        // 2) Assign a m2
        // Implement This


		// 3) Call a function `blockMM() to execute Matrix matrix multiplication
		const float* rst = this->run();

    // const float* real_ = this->real();

    // 4) Accumulate intermediate results
    // It is slightly different from the code for the project.
		for(int n = 0; n<block_row; ++n)
        {
          for(int m = 0; m<block_col_2; ++m)
          {
            output[n*SIZE + m] += rst[n*SIZE + m];
            printf("%f ", output[n*SIZE + m]);
          }
          printf("\n");
        }
        printf("\n");
    // for(int n = 0; n<block_row; n++) {
    //   for(int m=0; m<block_col_2; m++) {
    //     printf("%f ", real_[n*8+m]);
    //   }
    //   printf("\n");
    // }
		// 4) Accumulate intermediate results
 	  } 
    
	}
  }
}
