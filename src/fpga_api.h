#ifndef _FPGA_API_H_
#define _FPGA_API_H_
#include <sys/types.h>
#include <vector>
#include "compute.h"

// matrix vector multiplicator
// matrix M: M_SIZE by V_SIZE
// vector V: V_SIZE
// output = M * V

class FPGA
{

private:
  int fd_;
  float *data_;
  float *data_M;
  unsigned int *output_;
  unsigned int *output_M;
  unsigned int *output_MV;

  char *qvec_;
  char *qmat_;
  char *qm1_;
  char *qm2_;
  int *qout_;
  int *qout_M;
  int *qout_MV;
  
  char *qdata_;
  char *qdata_M;

  int m_size_;
  int v_size_;
  int m1_size_;
  int m2_size_;

  int data_size_;
  int data_size_M;
  int num_block_call_;

public:
  FPGA(off_t data_addr, off_t output_addr, int m_size, int v_size);
  ~FPGA();

  // return internal pointer for the data
  float *matrix(void);
  float *vector(void);
  float *matrix_M1(void);
  float *matrix_M2(void);

  char *qmatrix(void);
  char *qvector(void);
  char *qmatrix_M1(void);
  char *qmatrix_M2(void);

  void reset(void);
  int num_block_call(void);

  // perform matrix multiplication and return output array pointer
  const float *blockMV(Compute* comp);
  const int *qblockMV(Compute* comp);
  const int *blockMM(Compute* comp);
  const int *qblockMM(Compute* comp);

  // Input vector size: num_input
  // Matrix size: num_output * num_input
  // Output vector size: num_output
  // O = M * I
  void largeMV(const float *mat, const float *input, float *output, int num_input, int num_output, Compute* comp);
  void largeMM(const float *mat, const float *input, float *output, int num_input, int num_output, int num_matrix, Compute* comp);
};
#endif
