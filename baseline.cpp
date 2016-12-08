#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sys/time.h>

#define REAL float
#define NX 512
#define NY 512
#define NZ 512
#define T 10

const REAL cc = 0.1;
const REAL ce = 0.2;
const REAL cw = 0.3;
const REAL cs = 0.4;
const REAL cn = 0.5;
const REAL ct = 0.6;
const REAL cb = 0.7;

void cpubaseline(REAL* A, REAL* B)
{
	int slice = NX*NY;
    for (int i = 1; i < NX-1; ++i){
        for (int j = 1; j < NY-1; ++j){
            for (int k = 1; k < NZ-1; ++k){
                int idx = i + j*NX + k*slice;
                A[idx] = ce*B[idx+1] + cw*B[idx-1] + cs*B[idx+NX] + cn*B[idx-NX]
                        +ct*B[idx+slice] + cb*B[idx-slice] + cc*B[idx];
            }
        }
    }
    REAL* tmp;
    tmp = A; A = B; B = tmp;
    return;
}

int main() {
	int size = sizeof(REAL)*NX*NY*NZ;
	REAL* host_A = (REAL*)malloc(size);
	REAL* host_B = (REAL*)malloc(size);

	for (int k = 0; k < NX; k++)
		for (int j = 0; j < NY; j++)
			for (int i = 0; i < NZ; i++)
				host_B[k*NX*NX+j*NX+i] = i - j + 1.0/k;	

    //cpu baseline
    struct timeval t1,t2;
    gettimeofday(&t1, NULL);
    for (int t = 0; t < T; t++){
        cpubaseline(host_A, host_B);
    }
    gettimeofday(&t2, NULL);
    printf("%d\n", t2.tv_sec);
    printf("%d\n", t1.tv_sec);
    printf("%d\n", t2.tv_usec);
    printf("%d\n", t1.tv_usec);
    float elapsed_time = (t2.tv_sec-t1.tv_sec)*1.e+3f + (t2.tv_usec-t1.tv_usec)*1.e-3f;
	printf("cpubaseline: elapsed time = %f ms\n", elapsed_time);
	float flops = 1.0*13*(NX-2)*(NY-2)*(NZ-2)*T/1.e+6;
	flops /= elapsed_time;
	printf("cpubaseline: Gflops = %lf\n", flops);
    return 0;
}

