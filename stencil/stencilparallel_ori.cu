#include<cuda_runtime.h>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
using namespace std;

#define REAL double

#define BX 128
#define BY 2
#define BZ 1
#define GZ 1

const REAL cc = 0.4;
const REAL ce = 0.1;
const REAL cw = 0.1;
const REAL cs = 0.1;
const REAL cn = 0.1;
const REAL ct = 0.1;
const REAL cb = 0.1;

//Must be re-written, including all the parameters
int stencil(REAL *A, REAL *B, int nx, int ny, int nz, int steps)
{
  int i, j, k, s;
#define IDX(i,j,k) ((i)*ny*nz+(j)*nz+(k))
  for(s = 0; s < steps; s ++) {
    for(i = 0; i < nx; i ++) {
      for(j = 0; j < ny; j ++) {
        for(k = 0; k < nz; k ++) {
          REAL r = 0.4*A[IDX(i,j,k)];
          if(k !=  0)   r += 0.1*A[IDX(i,j,k-1)];
          else          r += 0.1*A[IDX(i,j,k)];
          if(k != nz-1) r += 0.1*A[IDX(i,j,k+1)];
          else          r += 0.1*A[IDX(i,j,k)];
          if(j !=  0)   r += 0.1*A[IDX(i,j-1,k)];
          else          r += 0.1*A[IDX(i,j,k)];
          if(j != ny-1) r += 0.1*A[IDX(i,j+1,k)];
          else          r += 0.1*A[IDX(i,j,k)];
          if(i !=  0)   r += 0.1*A[IDX(i-1,j,k)];
          else          r += 0.1*A[IDX(i,j,k)];
          if(i != nx-1) r += 0.1*A[IDX(i+1,j,k)];
          else          r += 0.1*A[IDX(i,j,k)];
          B[IDX(i,j,k)] = r;
        }
      }
    }
     REAL *tmp = NULL;
    tmp = A, A = B, B = tmp;
  }
  return 0;
}
void check(REAL *a, REAL *b, int nx, int ny, int nz) {
    int slice = nx * ny;
    for (int z = 1; z < nz-1; ++z) {
        for (int y = 1; y < ny-1; ++y) {
            for (int x = 1; x < nx-1; ++x) {
                int idx = z * slice + y * nx + x;
                //cout << abs(a[idx]-b[idx]) << endl;
                if (abs(a[idx]-b[idx]) > 1e-5) {
                    printf("%d\n", idx);
                    printf("Wrong!!!!!!!!\n");
                    return;
                }
            }
        }
    }
    printf("Right!!!!!!\n");
    return;
}
__global__ void baseline(REAL* A, REAL* B, int nx, int ny, int nz)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int kb = nz/gridDim.z*blockIdx.z;
	int slice = nx*ny;

    int k = kb;
	//int k = kb > 0? kb: 1;
	int ke = (kb+nz/gridDim.z<nz)? kb+nz/gridDim.z : nz;
	int c = i + j*nx + k*slice;
//#pragma unroll
	for (; k < ke; k++){
        int w = (i==0)?c:c-1;
        int e = (i==nx-1)?c:c+1;
        int n = (j==0)?c:c-nx;
        int s = (j==ny-1)?c:c+nx;
        int b = (k==0)?c:c-slice;
        int t = (k==nz-1)?c:c+slice;
        B[c] = ce*A[e] + cw*A[w] + cs*A[s] + cn*A[n]
            +ct*A[t] + cb*A[b] + cc*A[c];
        c += slice;
		//if (k > 0 && k < nz-1 && i > 0 && i < nx-1 && j > 0 && j < ny-1){
		//	B[idx] = ce*A[idx+1] + cw*A[idx-1] + cs*A[idx+nx] + cn*A[idx-nx]
		//			+ct*A[idx+slice] + cb*A[idx-slice] + cc*A[idx];
		//	idx += slice;
	}
}

__global__ void baseopt(REAL* A, REAL* B, int nx, int ny, int nz)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int kb = nz/gridDim.z*blockIdx.z;
	int slice = nx*ny;

	//int k = kb > 0? kb: 1;
    int k = kb;
	int ke = (kb+nz/gridDim.z<nz)? kb+nz/gridDim.z : nz;
	int c = i + j*nx + k*slice;
    int b = (k==0)?c:c-slice;
    int w = (i==0)?c:c-1;
    int e = (i==nx-1)?c:c+1;
    int n = (j==0)?c:c-nx;
    int s = (j==ny-1)?c:c+nx;
    int t;
    double b_b = A[b];
    double b_c = A[c];
    double b_t;
#pragma unroll
    for (; k < ke; k++){
        t = (k==nz-1)?c:c+slice;
        b_t = A[t];
        B[c] = ce*A[e] + cw*A[w] + cs*A[s] + cn*A[n]
            +ct*b_t + cb*b_b + cc*b_c;
        b_b = b_c;
        b_c = b_t;
        c += slice;
        //b_t = B[idx+slice];
        ////A[idx] = ce*B[idx+1] + cw*B[idx-1] + cs*B[idx+nx] + cn*B[idx-nx]
        ////		+ct*B[idx+slice] + cb*B[idx-slice] + cc*B[idx];
        //A[idx] = ce*B[idx+1] + cw*B[idx-1] + cs*B[idx+nx] + cn*B[idx-nx]
        //		+ct*b_t + cb*b_b + cc*b_c;
        //b_b = b_c;
        //b_c = b_t;
        //idx += slice;
    }
	return;
}

__global__ void roc(const REAL* __restrict__ A, REAL* B, int nx, int ny, int nz)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	int j = threadIdx.y + blockDim.y*blockIdx.y;
	int kb = nz/gridDim.z*blockIdx.z;
	int slice = nx*ny;

	//int k = kb > 0? kb: 1;
    int k = kb;
	int ke = (kb+nz/gridDim.z<nz)? kb+nz/gridDim.z : nz;
	int c = i + j*nx + k*slice;
    int b = (k==0)?c:c-slice;
    int w = (i==0)?c:c-1;
    int e = (i==nx-1)?c:c+1;
    int n = (j==0)?c:c-nx;
    int s = (j==ny-1)?c:c+nx;
    int t;
    double b_b = A[b];
    double b_c = A[c];
    double b_t;
#pragma unroll
    for (; k < ke; k++){
        t = (k==nz-1)?c:c+slice;
        b_t = A[t];
        B[c] = ce*A[e] + cw*A[w] + cs*A[s] + cn*A[n]
            +ct*b_t + cb*b_b + cc*b_c;
        b_b = b_c;
        b_c = b_t;
        c += slice;
    }
	return;
}

int main(int argc, char** argv){

    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    int NZ = atoi(argv[3]);
    int T = atoi(argv[4]);
	int size = sizeof(REAL)*NX*NY*NZ;
	REAL* host_A = (REAL*)malloc(size);
	REAL* host_B = (REAL*)malloc(size);
	REAL* cpu_A = (REAL*)malloc(size);
	REAL* cpu_B = (REAL*)malloc(size);

    for (int k = 0; k < NZ; k++)
		for (int j = 0; j < NY; j++)
			for (int i = 0; i < NX; i++) {
				host_A[k*NY*NX+j*NX+i] = 1.0;	
				host_B[k*NY*NX+j*NX+i] = 1.0;	
				cpu_A[k*NY*NX+j*NX+i] = 1.0;	
				cpu_B[k*NY*NX+j*NX+i] = 1.0;	
            }

	REAL *dev_A, *dev_B;
	cudaMalloc(&dev_A, size);
	cudaMalloc(&dev_B, size);
	cudaMemcpy(dev_A, host_A, size, cudaMemcpyHostToDevice);
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	float elapsed_time;
	double flops;

	dim3 threadPerBlock(BX, BY, BZ); //128,1,1
	dim3 blockPerGrid((NX+BX-1)/BX, (NY+BY-1)/BY, GZ); //512/128,512/1,1 = 4,512,1

	///////////////////////////////////////////////////////////////
	//baseline
	cudaEventRecord(start, 0);
	for (int t = 0; t < T; t++){
		baseline<<<blockPerGrid, threadPerBlock>>>(dev_A, dev_B, NX, NY, NZ);		
		REAL* tmp = dev_A;
		dev_A = dev_B;
		dev_B = tmp;
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	if (cudaGetLastError() != cudaSuccess)
		printf("baseline: wrong!!!\n");
    cudaMemcpy(host_A, dev_A, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_B, dev_B, size, cudaMemcpyDeviceToHost);
	cudaEventElapsedTime(&elapsed_time, start, stop);

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    //stencil(cpu_A, cpu_B, NX, NY, NZ, T);
    gettimeofday(&t2, NULL);
    float cpubaseline_time = (t2.tv_sec-t1.tv_sec)*1e3 + (t2.tv_usec-t1.tv_usec)*1e-3;
    cout << "CPU time:" << cpubaseline_time/T << " ms" << endl;
    //check(cpu_A, host_A, NX, NY, NZ);


	printf("baseline: elapsed time = %f ms\n", elapsed_time/T);
	flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
	flops /= elapsed_time;
	//printf("baseline: Gflops = %lf\n", flops);
	///////////////////////////////////////////////////////////////

/*
	///////////////////////////////////////////////////////////////
	//baseopt
	cudaEventRecord(start, 0);
	for (int t = 0; t < T; t++){
		baseopt<<<blockPerGrid, threadPerBlock>>>(dev_A, dev_B, NX, NY, NZ);		
		REAL* tmp = dev_A;
		dev_A = dev_B;
		dev_B = tmp;
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	if (cudaGetLastError() != cudaSuccess)
		printf("baseopt: wrong!!!\n");
	cudaEventElapsedTime(&elapsed_time, start, stop);

	printf("baseopt: elapsed time = %f ms\n", elapsed_time/T);
	flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
	flops /= elapsed_time;
	//printf("baseopt: Gflops = %lf\n", flops);
	///////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////
	//read-only data cache
	cudaEventRecord(start, 0);
	for (int t = 0; t < T; t++){
		roc<<<blockPerGrid, threadPerBlock>>>(dev_A, dev_B, NX, NY, NZ);		
		REAL* tmp = dev_A;
		dev_A = dev_B;
		dev_B = tmp;
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	if (cudaGetLastError() != cudaSuccess)
		printf("read-only data cache: wrong!!!\n");
	cudaEventElapsedTime(&elapsed_time, start, stop);

	printf("read-only data cache: elapsed time = %f ms\n", elapsed_time/T);
	flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
	flops /= elapsed_time;
    //printf("read-only data cache: Gflops = %lf\n", flops);
	///////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////
	//share memory raw
	cudaEventRecord(start, 0);
	for (int t = 0; t < T; t++){
		shm_raw<<<blockPerGrid, threadPerBlock>>>(dev_A, dev_B, NX, NY, NZ);		
		REAL* tmp = dev_A;
		dev_A = dev_B;
		dev_B = tmp;
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	if (cudaGetLastError() != cudaSuccess)
		printf("share memory raw: wrong!!!\n");
	cudaEventElapsedTime(&elapsed_time, start, stop);

	printf("share memory raw: elapsed time = %f ms\n", elapsed_time/T);
	flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
	flops /= elapsed_time;
	//printf("share memory raw: Gflops = %lf\n", flops);
	///////////////////////////////////////////////////////////////

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
*/
	return 0;
}
