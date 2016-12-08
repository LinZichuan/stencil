#include<cuda_runtime.h>
#include<iostream>
#include<stdio.h>
#include<sys/time.h>
#include<assert.h>
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
            for (int x = 1; x < nz-1; ++x) {
                int idx = z * slice + y * nx + x;
                if (abs(a[idx]-b[idx]) > 1e-5) {
                    cout << a[idx] << " " << b[idx] << endl;
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

int main(int argc, char **argv){

    int NX = atoi(argv[2]);
    int NY = atoi(argv[3]);
    int NZ = atoi(argv[4]);
    int T = atoi(argv[5]);
    cout << NX << " " << NY << " " << NZ << " " << T << endl;
    int num = 3;
    int NZ_ = NZ/num+2;
	int size = sizeof(REAL)*NX*NY*NZ;
	int partsize = sizeof(REAL)*NX*NY*NZ_;
    REAL **host_A = new REAL*[num];
    REAL **host_B = new REAL*[num];
    int size_ = NZ_*NY*NX;
    for (int i = 0; i < num; ++i) {
        host_A[i] = new REAL[size_];
        host_B[i] = new REAL[size_];
        //host_A[i] = (REAL*)malloc(partsize);
        //host_B[i] = (REAL*)malloc(partsize);
    }
	REAL* cpu_A = new REAL[NX*NY*NZ];
    REAL* result_A = new REAL[NX*NY*NZ];
	REAL* cpu_B = new REAL[NX*NY*NZ];

    for (int part = 0; part < num; part++)
        for (int k = 0; k < NZ_; k++)
            for (int j = 0; j < NY; j++)
                for (int i = 0; i < NX; i++) {
                    host_A[part][k*NY*NX+j*NX+i] = 1.0;	
                    host_B[part][k*NY*NX+j*NX+i] = 1.0;
                }

    for (int k = 0; k < NZ; k++)
		for (int j = 0; j < NY; j++)
			for (int i = 0; i < NX; i++) {
                //cout << k*NY*NX + j*NX + i << endl;
				cpu_A[k*NY*NX+j*NX+i] = 1.0;	
				cpu_B[k*NY*NX+j*NX+i] = 1.0;
				result_A[k*NY*NX+j*NX+i] = 1.0;
            }

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    float elapsed_time;
    double flops;
    for (int i = 0; i < num; ++i) {
        REAL *dev_A, *dev_B;
        cudaMalloc(&dev_A, partsize);
        cudaMalloc(&dev_B, partsize);
        cudaMemcpy(dev_A, host_A[i], partsize, cudaMemcpyHostToDevice);

        dim3 threadPerBlock(BX, BY, BZ); //128,1,1
        dim3 blockPerGrid((NX+BX-1)/BX, (NY+BY-1)/BY, GZ); //512/128,512/1,1 = 4,512,1

        ///////////////////////////////////////////////////////////////
        //baseline
        for (int t = 0; t < T; t++){
            baseline<<<blockPerGrid, threadPerBlock>>>(dev_A, dev_B, NX, NY, NZ_);		
            REAL* tmp = dev_A;
            dev_A = dev_B;
            dev_B = tmp;
        }
        ///////////////////////////////////////////////////////////////
        if (cudaGetLastError() != cudaSuccess)
            printf("baseline: wrong!!!\n");

        cudaMemcpy(host_A[i], dev_A, partsize, cudaMemcpyDeviceToHost);
        cudaFree(dev_A);
        cudaFree(dev_B);
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    stencil(cpu_A, cpu_B, NX, NY, NZ, T);
    gettimeofday(&t2, NULL);
    float cpubaseline_time = (t2.tv_sec-t1.tv_sec)*1e3 + (t2.tv_usec-t1.tv_usec)*1e-3;
    cout << "CPU time:" << cpubaseline_time/T << " ms" << endl;

    int begin, end;
    int smallsize = NZ/num * NY * NX;
    int i=0, z=0, y=0, x=0;
    for (i = 0; i < num; ++i) {
        begin = 1;
        end = NZ_-1;
        if (i == 0) { begin=0; end=NZ_-2; } 
        else if (i == num-1) { begin=2; end=NZ_; }
        int index = i*smallsize;
        for (z = begin; z < end; ++z)
            for (y = 0; y < NY; ++y)
                for (x = 0; x < NX; ++x) {
                    result_A[index] = host_A[i][NY*NX*z + y*NX + x];
                    //assert(abs(host_A[i][NY*NX*z + y*NX + x] - 1.0) < 1e-5);
                    //if (i == 2)
                    //    cout <<  host_A[i][NY*NX_*z + y*NX_ + x] << endl;
                    index++;
                }
    }

    check(cpu_A, result_A, NX, NY, NZ);
    //printf("baseline: Gflops = %lf\n", flops);

    printf("baseline: elapsed time = %f ms\n", elapsed_time/T);
    flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
    flops /= elapsed_time;

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
