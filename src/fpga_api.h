#ifndef _FPGA_API_H_
#define _FPGA_API_H_

#include <sys/types.h>

// matrix vector multiplicator
// matrix M: SIZE by SIZE
// vector V: SIZE
// output = M * V

#define SIZE 8 // input, weight, output Matrix = 8 x 8

class FPGA
{
private:
    int fd_;	
    int* data_;
    unsigned int* api_;
	// float* real_;
	
public:
    FPGA(off_t data_addr, off_t api_addr);
    ~FPGA();
	
	// return internal pointer for the data
	char* matrix_M1(void);
	char* matrix_M2(void);
	
	// perform matrix multiplication and return output array pointer
    const int* run();	
	// const float* real();
	
	// input vector size: M
	// matrix size: N by M
	// output vector size: N
	// O = M * I d
	void largeMM(const char* weight_mat, const char* input_mat, int* output, int num_input, int num_output, int num_matrix2);
};

#endif

