#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include <sstream>
using namespace std;

void Aproduct(double *A, double *x,double *Ax, int nx, int ny, int nz) {
    int size = nx * ny * nz;
    int slice = ny * nz;
   // double *Ax = new double[size];
#define IDXA(i,j,k,n) (i)*ny*nz*19+(j)*nz*19+(k)*19+(n)
#define IDXx(i,j,k  ) (i+1)*(ny+2)*(nz+2)+(j+1)*(nz+2)+(k+1)
#define IDXAx(i,j,k ) (i)*ny*nz + (j)*nz + (k)
//#pragma omp parallel for
    for (int i = 0; i < nx; i ++ ) {
        for (int j = 0; j < ny; j ++ ) {
            for (int k = 0; k < nz; k ++ ) {
                Ax[IDXAx(i,j,k)] =   A[IDXA(i,j,k,0 )]  * x[IDXx(i  ,j  ,k  )]
                                   + A[IDXA(i,j,k,1 )]  * x[IDXx(i-1,j  ,k  )]
                                   + A[IDXA(i,j,k,2 )]  * x[IDXx(i+1,j  ,k  )]
                                   + A[IDXA(i,j,k,3 )]  * x[IDXx(i  ,j-1,k  )]
                                   + A[IDXA(i,j,k,4 )]  * x[IDXx(i  ,j+1,k  )]
                                   + A[IDXA(i,j,k,5 )]  * x[IDXx(i+1,j+1,k  )]
                                   + A[IDXA(i,j,k,6 )]  * x[IDXx(i+1,j-1,k  )]
                                   + A[IDXA(i,j,k,7 )]  * x[IDXx(i-1,j-1,k  )]
                                   + A[IDXA(i,j,k,8 )]  * x[IDXx(i-1,j+1,k  )]
                                   + A[IDXA(i,j,k,9 )]  * x[IDXx(i  ,j  ,k-1)]
                                   + A[IDXA(i,j,k,10)]  * x[IDXx(i-1,j  ,k-1)]
                                   + A[IDXA(i,j,k,11)]  * x[IDXx(i+1,j  ,k-1)]
                                   + A[IDXA(i,j,k,12)]  * x[IDXx(i  ,j-1,k-1)]
                                   + A[IDXA(i,j,k,13)]  * x[IDXx(i  ,j+1,k-1)]
                                   + A[IDXA(i,j,k,14)]  * x[IDXx(i  ,j  ,k+1)]
                                   + A[IDXA(i,j,k,15)]  * x[IDXx(i-1,j  ,k+1)]
                                   + A[IDXA(i,j,k,16)]  * x[IDXx(i+1,j  ,k+1)]
                                   + A[IDXA(i,j,k,17)]  * x[IDXx(i  ,j-1,k+1)]
                                   + A[IDXA(i,j,k,18)]  * x[IDXx(i  ,j+1,k+1)] ;
            }
        }
    }
   // return Ax;
}
bool judge(const pair<int,int> a, const pair<int,int> b) {
    return a.first < b.first;
}
int getcsrA(double *A, int *Ai, int *Aj, double *value, int* diagpos, int nx, int ny, int nz, int xs, int xe, int ys, int ye, int NX, int NY, int NZ) {
    //TODO: nan, mpi wrong
    int slice = ny*nz;
    int nonzero = 0;
    vector<int> a(19); //only a is global
    int nxslice = (nx-1)*slice;
#define IN(i,j,k)  (((i+NX)%NX>=xs && (i+NX)%NX<xe && (j)>=ys && (j)<ye) ? 1:0)
#define xIDX(i,j,k) (IN(i,j,k)?(((i+NX)%NX-xs)*(ny)*(nz) + (j-ys)*(nz) + (k)):-1) //local
    for (int ii = 0; ii < nx; ++ii) {
        for (int jj = 0; jj < ny; ++jj) {
            for (int kk = 0; kk < nz; ++kk) {
                int i = xs + ii;
                int j = ys + jj;
                int k = kk;
                int idx = (ii*slice + jj*nz + kk)*19; //local idx
                int idxAi = idx/19;
                int a0 = xIDX(i,j,k);  //i*slice + j*nz + k;  //idx/19;
                a[0]  = xIDX(i,j,k);  //i*slice + j*nz + k;
                //i-1,j,k
                a[1] = xIDX(i-1,j,k);
                //i+1,j,k
                a[2] = xIDX(i+1,j,k);
                //i,j-1,k
                if(j == 0){
                    a[3] = xIDX((i+NX/2)%NX, j, k);
                }else{
                    a[3] = xIDX(i, j-1, k);
                }
                //i,j+1,k
                if(j == NY-1){
                    a[4] = xIDX((i+NX/2)%NX, j, k);
                }else{
                    a[4] = xIDX(i, j+1, k);
                }
                //i+1,j+1,k
                if(j == NY-1){
                    a[5] = xIDX((i+NX/2+1)%NX, j, k);
                }else{
                    a[5] = xIDX(i+1, j+1, k);
                }
                //i+1,j-1,k
                if (j == 0) {
                    a[6] = xIDX((i+NX/2+1)%NX, j, k);
                } else {
                    a[6] = xIDX(i+1, j-1, k);
                }
                //i-1,j-1,k
                if (j == 0) {
                    a[7] = xIDX((i+NX/2-1)%NX, j, k);
                } else {
                    a[7] = xIDX(i-1, j-1, k);
                }
                //i-1,j+1,k
                if (j == NY-1) {
                    a[8] = xIDX((i+NX/2-1)%NX, j, k);
                } else {
                    a[8] = xIDX(i-1, j+1, k);
                }
                //i,j,k-1
                a[9]  = xIDX(i, j, k-1);
                //i-1,j,k-1
                a[10] = xIDX(i-1, j, k-1);
                //i+1,j,k-1
                a[11] = xIDX(i+1, j, k-1);
                //i,j-1,k-1
                if (j == 0) {
                    a[12] = xIDX((i+NX/2)%NX, j, k-1);
                } else {
                    a[12] = xIDX(i, j-1, k-1);
                }
                //i,j+1,k-1
                if (j == NY-1) { //TODO
                    a[13] = xIDX((i+NX/2)%NX, j, k-1);
                } else {
                    a[13] = xIDX(i, j+1, k-1);
                }
                //i,j,k+1
                a[14] = xIDX(i,j,k+1);
                //i-1,j,k+1
                a[15] = xIDX(i-1, j, k+1);
                //i+1,j,k+1
                a[16] = xIDX(i+1, j, k+1);
                //i,j-1,k+1
                if (j == 0) {
                    a[17] = xIDX((i+NX/2)%NX, j, k+1);
                } else {
                    a[17] = xIDX(i, j-1, k+1);
                }
                //i,j+1,k+1
                if (j == NY-1) {
                    a[18] = xIDX((i+NX/2)%NX, j, k+1);
                } else {
                    a[18] = xIDX(i, j+1, k+1);
                }
                if (k == 0) a[9] = a[10] = a[11] = a[12] = a[13] = -1;
                else if (k == NZ-1) a[14] = a[15] = a[16] = a[17] = a[18] = -1;
              
                //sort b
                vector< pair<int,int> > b(19);
                for (int t = 0; t < 19; ++t) {
                    b[t] = make_pair<int,int>(a[t], t);
                    //cout << A[idx+t] << endl;
                }
                sort(b.begin(), b.end(), judge);
                //sort(a.begin(), a.end());
                int nonzero_num = 0;
                for (int t = 0; t < 19; ++t) {
                    if (b[t].first == -1) continue;
                    //int xx = b[t].first / slice;
                    //int yy = (b[t].first % slice) / nz;
                 //if ((i+NX)%NX < xs || (i+NX)%NX >= xe || j < ys || j >= ye) continue;
                    //if (xx < xs || xx >= xe || yy < ys || yy >= ye) continue; //TODO
                    value[nonzero] = A[idx+b[t].second]; ///BUG!!
                    Aj[nonzero] = b[t].first;
                    if (b[t].first == a0) {
                        diagpos[idxAi] = nonzero;  //local
                    }
                    //cout << "a[t]:" << a[t] << endl;
                    nonzero++;
                    nonzero_num++;
                }
                Ai[idxAi+1] = Ai[idxAi] + nonzero_num ;  //local
                //cout << Ai[a0+1] << endl;
            }
        }
    }
    return nonzero;
}
double* Matrixminus(double *b, double *Ax, int nx, int ny, int nz) {
    int size = nx * ny * nz;
    double *R = new double[size];
    for (int idx = 0; idx < size; ++idx) {
        R[idx] = b[idx] - Ax[idx];
        //cout << "Ax:" << Ax[idx] << endl;
        //cout << "R:" << R[idx] << endl;
    }
    return R;
}
double dotproduct(double *a, double *b, int n) {
    double res = 0.0;
    double res_sum = 0.0; 
    for (int i = 0; i < n; ++i) {
        res += a[i] * b[i];
        //cout << "a:" << a[i] << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&res,&res_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    return res_sum;
}
void saxpy(double *x, double *p, double alpha, int n) {
    for (int i = 0; i < n; ++i) {
        x[i] = x[i] + alpha*p[i];
    }
}
double* add(double *a, double *b, int n) {
    double *res = new double[n];
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}
double* sigma(double *beta, vector<double*> &vec, int jbegin, int i, int size) {
    double *res = new double[size];
    for (int t = 0; t < size; ++t) {
        res[t] = 0.0;
    }
    for (int j = jbegin; j < i; ++j) {
        for (int t = 0; t < size; ++t) {
            res[t] += beta[j-jbegin] * vec[j][t];
        }
    }
    return res;
}
void preprocess(int *Ai, int *Aj, double *value, int *diagpos, double *R, double *R_, int size, int nx, int ny, int nz, int xs, int xe, int ys, int ye) {
    //R=M*R_ , solve R_
    //diag of L = 1
    //Ly = R
    double *y = new double[size];
    for (int i = 0; i < size; ++i) {
        int begin = Ai[i];
        int end = diagpos[i]; //Ai[i+1];
        double tmp = 0.0;
        for (int j = begin; j < end; ++j) {
            //int xidx = Aj[j] / (ny*nz);
            //int yidx = (Aj[j] % (ny*nz)) / nz;
            //if (xidx<xs || xidx>=xe || yidx<ys || yidx>=ye) continue;
            tmp += value[j] * y[Aj[j]];
        }
        y[i] = (R[i]-tmp) / 1.0;
    }
    //UR_ = y
    for (int i = size-1; i >= 0; --i) {
        int begin = diagpos[i]+1; //Ai[i];
        int end = Ai[i+1];
        double tmp = 0.0;
        for (int j = begin; j < end; ++j) {
            //int xidx = Aj[j] / (ny*nz);
            //int yidx = (Aj[j] % (ny*nz)) / nz;
            //if (xidx<xs || xidx>=xe || yidx<ys || yidx>=ye) continue;
            tmp += value[j] * R_[Aj[j]];
        }
        R_[i] = (y[i]-tmp) / value[diagpos[i]];
    }
}
void communicate(double *A, int nx, int ny, int nz, int px, int py, int PX, int PY) {
    //cout << nx << " " << ny << " " << nz << " " << size << endl;
    int ysize = (ny) * (nz);
    int xsize = (nx) * (nz);
    int csize = nz+2;
    double yhalo[2][ysize];
    double xhalo[2][xsize];
    double chalo[4][csize];
    double s_yhalo[2][ysize];
    double s_xhalo[2][xsize];
    double s_chalo[4][csize];
    int source[9];
    int dest[9];
    int tag[9];
#define IDX(x,y) ((x)*PY+(y))
    //west
    dest[1] = IDX((px-1+PX)%PX, py);
    tag[1] = 13;
    //east
    dest[3] = IDX((px+1+PX)%PX, py);
    tag[3] = 31;
    //north
    if (py == PY-1) {
        dest[2] = IDX((px+PX/2)%PX, py);
        tag[2] = 22;
    } else {
        dest[2] = IDX(px, py+1);
        tag[2] = 24;
    }
    //south
    if (py == 0) {
        dest[4] = IDX((px+PX/2)%PX, py);
        tag[4] = 44;
    } else {
        dest[4] = IDX(px, py-1);
        tag[4] = 42;
    }
    //west-south
    if (py == 0) {
        dest[5] = IDX((px-1+PX/2)%PX, py);
        tag[5] = 58;
    } else {
        dest[5] = IDX((px-1+PX)%PX, py-1);
        tag[5] = 57;
    }
    //west-north
    if (py == PY-1) {
        dest[6] = IDX((px-1+PX/2)%PX, py);
        tag[6] = 67;
    } else {
        dest[6] = IDX((px-1+PX)%PX, py+1);
        tag[6] = 68;
    }
    //east-north
    if (py == PY-1) {
        dest[7] = IDX((px+1+PX/2)%PX, py);
        tag[7] = 76;
    } else {
        dest[7] = IDX((px+1)%PX, py+1);
        tag[7] = 75;
    }
    //east-south
    if (py == 0) {
        dest[8] = IDX((px+1+PX/2)%PX, py);
        tag[8] = 85;
    } else {
        dest[8] = IDX((px+1)%PX, py-1);
        tag[8] = 86;
    }
    for (int i = 1; i < 9; ++i) source[i] = dest[i];
    //assign A,b
#define AIDX(x,y,z) ((x)*(ny+2)*(nz+2)+(y)*(nz+2)+(z))
    int tmp = 0;
    for (int y = 1; y < ny+1; ++y) {
        for (int z = 1; z < nz+1; ++z) {
            //A[AIDX(0, y, z)] = yhalo[0][tmp];
            s_yhalo[0][tmp] = A[AIDX(1, y, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int y = 1; y < ny+1; ++y) {
        for (int z = 1; z < nz+1; ++z) {
            //A[AIDX(nx+1, y, z)] = yhalo[1][tmp];
            s_yhalo[1][tmp] = A[AIDX(nx, y, z)];
            tmp++;
        }
    }
    //return;
    tmp = 0;
    for (int x = 1; x < nx+1; ++x) {
        for (int z = 1; z < nz+1; ++z) {
            //A[AIDX(x, ny+1, z)] = xhalo[0][tmp];
            s_xhalo[0][tmp] = A[AIDX(x, ny, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int x = 1; x < nx+1; ++x) {
        for (int z = 1; z < nz+1; ++z) {
            //A[AIDX(x, 0, z)] = xhalo[1][tmp];
            s_xhalo[1][tmp] = A[AIDX(x, 1, z)]; 
            tmp++;
        }
    }
    tmp = 0;
    for (int z = 1; z < nz+1; ++z) {
       // A[AIDX(0,0,z)] = chalo[0][tmp];
        //A[AIDX(0,ny+1,z)] = chalo[1][tmp];
        //A[AIDX(nx+1,ny+1,z)] = chalo[2][tmp];
        //A[AIDX(nx+1,0,z)] = chalo[3][tmp];
        s_chalo[0][tmp] = A[AIDX(1,1,z)];
        s_chalo[1][tmp] = A[AIDX(1,ny,z)];
        s_chalo[2][tmp] = A[AIDX(nx,ny,z)];
        s_chalo[3][tmp] = A[AIDX(nx,1,z)];
        tmp++;
    }

    MPI_Send(s_yhalo[0], ysize, MPI_DOUBLE, dest[1], tag[1], MPI_COMM_WORLD);
    MPI_Send(s_yhalo[1], ysize, MPI_DOUBLE, dest[3], tag[3], MPI_COMM_WORLD);
    MPI_Send(s_xhalo[0], xsize, MPI_DOUBLE, dest[2], tag[2], MPI_COMM_WORLD);
    MPI_Send(s_xhalo[1], xsize, MPI_DOUBLE, dest[4], tag[4], MPI_COMM_WORLD);
    MPI_Send(s_chalo[0], csize, MPI_DOUBLE, dest[5], tag[5], MPI_COMM_WORLD);
    MPI_Send(s_chalo[1], csize, MPI_DOUBLE, dest[6], tag[6], MPI_COMM_WORLD);
    MPI_Send(s_chalo[2], csize, MPI_DOUBLE, dest[7], tag[7], MPI_COMM_WORLD);
    MPI_Send(s_chalo[3], csize, MPI_DOUBLE, dest[8], tag[8], MPI_COMM_WORLD);

    MPI_Status status;
#define INV(num) num%10*10 + num/10
    MPI_Recv(yhalo[0], ysize, MPI_DOUBLE, source[1], INV(tag[1]), MPI_COMM_WORLD, &status);
    MPI_Recv(yhalo[1], ysize, MPI_DOUBLE, source[3], INV(tag[3]), MPI_COMM_WORLD, &status);
    MPI_Recv(xhalo[0], xsize, MPI_DOUBLE, source[2], INV(tag[2]), MPI_COMM_WORLD, &status);
    MPI_Recv(xhalo[1], xsize, MPI_DOUBLE, source[4], INV(tag[4]), MPI_COMM_WORLD, &status);
    MPI_Recv(chalo[0], csize, MPI_DOUBLE, source[5], INV(tag[5]), MPI_COMM_WORLD, &status);
    MPI_Recv(chalo[1], csize, MPI_DOUBLE, source[6], INV(tag[6]), MPI_COMM_WORLD, &status);
    MPI_Recv(chalo[2], csize, MPI_DOUBLE, source[7], INV(tag[7]), MPI_COMM_WORLD, &status);
    MPI_Recv(chalo[3], csize, MPI_DOUBLE, source[8], INV(tag[8]), MPI_COMM_WORLD, &status);

    tmp = 0;
    for (int y = 1; y < ny+1; ++y) {
        for (int z = 1; z < nz+1; ++z) {
            A[AIDX(0, y, z)] = yhalo[0][tmp];
            //s_yhalo[0][tmp] = A[AIDX(1, y, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int y = 1; y < ny+1; ++y) {
        for (int z = 1; z < nz+1; ++z) {
            A[AIDX(nx+1, y, z)] = yhalo[1][tmp];
            //s_yhalo[1][tmp] = A[AIDX(nx, y, z)];
            tmp++;
        }
    }
    //return;
    tmp = 0;
    for (int x = 1; x < nx+1; ++x) {
        for (int z = 1; z < nz+1; ++z) {
            A[AIDX(x, ny+1, z)] = xhalo[0][tmp];
            //s_xhalo[0][tmp] = A[AIDX(x, ny, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int x = 1; x < nx+1; ++x) {
        for (int z = 1; z < nz+1; ++z) {
            A[AIDX(x, 0, z)] = xhalo[1][tmp];
            //s_xhalo[1][tmp] = A[AIDX(x, 1, z)]; 
            tmp++;
        }
    }
    tmp = 0;
    for (int z = 1; z < nz+1; ++z) {
        A[AIDX(0,0,z)] = chalo[0][tmp];
        A[AIDX(0,ny+1,z)] = chalo[1][tmp];
        A[AIDX(nx+1,ny+1,z)] = chalo[2][tmp];
        A[AIDX(nx+1,0,z)] = chalo[3][tmp];
        //s_chalo[0][tmp] = A[AIDX(1,1,z)];
        //s_chalo[1][tmp] = A[AIDX(1,ny,z)];
        //s_chalo[2][tmp] = A[AIDX(nx,ny,z)];
        //s_chalo[3][tmp] = A[AIDX(nx,1,z)];
        tmp++;
    }
}
void spmv(double *A, double *x, double *Ax, int nx, int ny, int nz, int px, int py, int PX, int PY) {
    int size = nx * ny * nz;
    int sizex = (nx+2)*(ny+2)*(nz+2);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get Process ID
    double *xhalo = new double[sizex];
    for (int i = 0; i < sizex; ++i)
        xhalo[i] = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int xidx = (i+1)*(ny+2)*(nz+2) + (j+1)*(nz+2) + k+1;
                int idx = (i)*(ny)*(nz) + (j)*(nz) + k;
                xhalo[xidx] = x[idx];
            }
        }
    }
    //for (int i = 0; i < nx; ++i) {
    //    int xidx = (i)*(ny+2)*(nz+2) ;//+ (j+1)*(nz+2) + k+1;
    //    cout << xhalo[xidx] << endl;
    //}
    communicate(xhalo, nx, ny, nz, px, py, PX, PY);
    //return;
    Aproduct(A, xhalo,Ax, nx, ny, nz);
}
double check_sum(double *value, int nonzero) {
    double sum = 0.0;
    double res_sum = 0.0;
    for (int i = 0; i < nonzero; ++i) {
        //cout << value[i] << endl;
        sum += value[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&sum,&res_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return res_sum;
}
void gcr(double *x, double *A, double *b, int nx, int ny, int nz, int px, int py, int NX, int NY, int NZ, int PX, int PY) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get Process ID
    //sync x
    int size = nx*ny*nz;
    //int sizex = (nx+2)*(ny+2)*(nz+2);
    //if (rank == 0) {
    //    for (int i = 0; i < nx+2; ++i) {
    //        //cout << x[i*ny*nz] << endl;
    //    }
    //}
    //communicate(x, nx, ny, nz, px, py, PX, PY);
    //if (rank == 0) {
    //    for (int i = 0; i < nx+2; ++i) {
    //        cout << x[i*(ny+2)*(nz+2)] << endl;
    //    }
    //}

    double *Ax = new double[size];
    spmv(A, x, Ax, nx, ny, nz, px, py, PX, PY); //Aproduct(A, x, nx, ny, nz);
    //return;
    //return;
    double *R = Matrixminus(b, Ax, nx, ny, nz);
    double Rnorm2 = dotproduct(R, R, size);
    if (rank == 0) 
        cout << "init rnorm2=" << Rnorm2 << endl; //1e-5
    double *R_ = new double[size];
    double *p0 = new double[size];
    
    int *Ai = new int[size+1];
    Ai[0] = 0;
    int *Aj = new int[size*19+1];
    double *value = new double[size*19+1];
    //cout << 1 << endl;
    int *diagpos = new int[size];
    //TODO
    int xs = px*NX/PX;
    for (int i = 0; i < px; ++i) {
        xs += ((px*2)%PX<(NX%PX)?1:0);
    }
    int xe = xs+nx;
    int ys = py*NY/PY + (py<(NY%PY)?py:(NY%PY));
    int ye = ys+ny;
    //for (int i = 0; i < size*19; ++i) {
    //    cout << A[i] << endl;
    //}
    int nonzero = getcsrA(A, Ai, Aj, value, diagpos, nx, ny, nz, xs, xe, ys, ye, NX, NY, NZ);
    //cout << 2 << endl;
    //ILU
    for (int i = 0; i < size; ++i) {
        int begin = Ai[i];
        int end = Ai[i+1];
        for (int j = begin; j < diagpos[i]; ++j) {
            //TODO: judge, in ns~ne
            //int xidx = Aj[j] / (ny*nz);
            //int yidx = Aj[j] % (ny*nz) / nz;
            //if (xidx<xs || xidx>=xe || yidx<ys || yidx>=ye) continue;
            int jj = diagpos[Aj[j]];
           // if (jj >= size || j >= size) continue;
            value[j] /= value[jj];
            int jbegin = j+1;
            int jend = end;
            int kbegin = jj+1;
            //xidx = (Aj[j]+1) / (ny*nz);
            //yidx = (Aj[j]+1) % (ny*nz) / nz;
            //if (xidx<xs || xidx>=xe || yidx<ys || yidx>=ye) continue;
            int kend = Ai[Aj[j]+1];
            int idxj=jbegin, idxk=kbegin;
            while (idxj < jend && idxk < kend) {
                if (Aj[idxj] == Aj[idxk]) {
                    value[idxj] -= (value[j] * value[idxk]); //Aik = Aik - Aij * Ajk;
                    idxj++;
                    idxk++;
                }
                else if (Aj[idxj] < Aj[idxk]) {
                    idxj += 1;
                } else {
                    idxk += 1;
                }
            }
        }
    }
    //checksum
    int checksum = check_sum(value,nonzero);
    if(rank == 0 ){
        cout<<"check_sum:"<<checksum<<endl;
    }
    //preprocess
    preprocess(Ai, Aj, value, diagpos, R, R_, size, nx, ny, nz, xs, xe, ys, ye);
    //return;
    //cout << 4 << endl;
   // for (int i = 0; i < size; ++i) {
      //  R_[i] = R[i];
     //   p0[i] = R_[i]; //p = R_
   // }
    vector<double*> vecAp;
    vector<double*> vecp;
    //TODO:少一次通信 
    //communicate(p0, nx, ny, nz, px, py, PX, PY);
    //double *Ap0 = Aproduct(A, p0, nx, ny, nz);
    double *Ap0 = new double[size];
    spmv(A, R_, Ap0, nx, ny, nz, px, py, PX, PY); //Aproduct(A, x, nx, ny, nz);
    vecp.push_back(R_);
    vecAp.push_back(Ap0);
    int k = 5;
    double vecAp_2[150];
    for (int i = 1; i < 29; ++i) {
        vecAp_2[i-1] = dotproduct(vecAp[i-1], vecAp[i-1], size);
        double alpha = dotproduct(R, vecAp[i-1], size) / vecAp_2[i-1];
        //cout << i << endl;
        saxpy(x, vecp[i-1], alpha, size);  //x = x + alpha * p;
        saxpy(R, vecAp[i-1], -alpha, size);  //R = R - alpha * Ap;
        double Rnorm2 = dotproduct(R, R, size);
        //TODO:算一下中间值 
		//preprocess
        preprocess(Ai, Aj, value, diagpos, R, R_, size, nx, ny, nz, xs, xe, ys, ye);
       /* for (int j = 0; j < size; ++j) {
            R_[j] = R[j];
        }*/
        //sync R_
        int jbegin = int((i-1)/k)*k;
        double *beta = new double[i-jbegin];
        //communicate(R_, nx, ny, nz, px, py, PX, PY);
        //double *AR_ = Aproduct(A, R_, nx, ny, nz);
        double *AR_ = new double[size];
        spmv(A, R_, AR_, nx, ny, nz, px, py, PX, PY); //Aproduct(A, x, nx, ny, nz);
        for (int j = jbegin; j < i; ++j) {
            beta[j-jbegin] = - dotproduct(AR_, vecAp[j], size) / vecAp_2[j];
        }
        //TODO
        double *p = add(R_, sigma(beta, vecp, jbegin, i, size), size);
        double *Ap = add(AR_, sigma(beta, vecAp, jbegin, i, size), size);
        vecp.push_back(p);
        vecAp.push_back(Ap);
        if (rank == 0)
            cout << "Rnorm2=" << Rnorm2 << endl; //1e-5
    }
    //communicate(x, nx, ny, nz, px, py, PX, PY);
    //double *pp = Aproduct(A, x, nx, ny, nz);
    //double *diff = Matrixminus(pp, b, nx, ny, nz);
    //double norm2_ = dotproduct(diff, diff, size);
    //if (rank == 0)
    //    cout << norm2_ << endl;

    //if(rank == 0){
    //    double *Ax_last = Aproduct(A, x, nx, ny, nz);
    //    double *R_last = Matrixminus(b, Ax_last, nx, ny, nz);
    //    double Rnorm2_last = dotproduct(R_last, R_last, size);
    //    cout << "Rnorm2_last=" << Rnorm2_last << endl; //1e-5
    //}
}
bool load_data(double *A, double *b, double *x, int NX, int NY, int NZ, int size) {
    FILE *fp, *fpb, *fpx;;
    string nxstr, nystr, nzstr;
    stringstream ss1, ss2, ss3;
    ss1 << NX; ss1 >> nxstr;
    ss2 << NY; ss2 >> nystr;
    ss3 << NZ; ss3 >> nzstr;
    string binfile = "case_"+nxstr+"x"+nystr+"x"+nzstr+"/A.bin";
    string binfileb = "case_"+nxstr+"x"+nystr+"x"+nzstr+"/b.bin";
    string binfilex = "case_"+nxstr+"x"+nystr+"x"+nzstr+"/x0.bin";
    if ((fp = fopen(binfile.c_str(), "rb")) == NULL) {
        cout << "cannot open A.bin" << endl;
        return 0;
    }
    if ((fpb = fopen(binfileb.c_str(), "rb")) == NULL) {
        cout << "cannot open b.bin" << endl;
        return 0;
    }
    if ((fpx = fopen(binfilex.c_str(), "rb")) == NULL) {
        cout << "cannot open x0.bin" << endl;
        return 0;
    }
    if (fread(A, sizeof(double), size*19, fp) != size*19) {
        cout << "file A.bin read error" << endl;
        return 0;
    }
    if (fread(b, sizeof(double), size, fpb) != size) {
        cout << "file b.bin read error" << endl;
        return 0;
    }
    if (fread(x, sizeof(double), size, fpx) != size) {
        cout << "file x0.bin read error" << endl;
        return 0;
    }
    fclose(fp);
    fclose(fpb);
    fclose(fpx);
    return 1;
}
int main(int argc,char**argv) {

    int NX=720/2, NY=360/2, NZ=38;
    int slice = NY*NZ;
    int size = NX*NY*NZ;
    double *A = new double[size*19];
    double *b = new double[size];
    double *x = new double[size];
    if (!load_data(A, b, x, NX, NY, NZ, size)) {
        cout << "Load Data Failed!" << endl;
        return 0;
    }

    timeval t1, t2;
    gettimeofday(&t1, NULL);

    int PX = 4, PY = 6;

    MPI_Init(&argc, &argv);
    int rank, nprocs;
    char hostname[100];
    gethostname(hostname, 100);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get Process ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Get # of Processes
    //cout<<"rank number:"<<rank<<endl;
    //parallel
    int px = rank/PY;
    int py = rank%PY;
    int nx_ = NX/PX + (((px*2)%PX<(NX%PX))?1:0);
    int ny_ = NY/PY + (py<(NY%PY)?1:0);
    int nz_ = NZ;

    int size_ = nx_*ny_*nz_;
    //int sizex = (nx_+2) * (ny_+2) * (nz_+2);
    double *A_ = new double[size_*19];
    double *b_ = new double[size_];
    double *x_ = new double[size_];
    for (int i = 0; i < size_*19; ++i)
        A_[i] = 0;
    for (int i = 0; i < size_; ++i)
        b_[i] = 0;
    for (int i = 0; i < size_; ++i)
        x_[i] = 0;
    //assign A,b,x
    int xs = px*NX/PX;
    for (int i = 0; i < px; ++i) {
        xs += ((px*2)%PX<(NX%PX)?1:0);
    }
    int xe = xs+nx_;
    int ys = py*NY/PY + (py<(NY%PY)?py:(NY%PY));
    int ye = ys+ny_;
    //    cout << "rank:" << rank << endl;
    //    cout << xs << " " << xe << " " << ys <<  " " << ye << endl;
#define local_idx(i1,j1,k1) (((i1)*ny_*nz_)+((j1)*nz_)+(k1))
#define local_idxx(i1,j1,k1) (((i1)*(ny_)*(nz_))+((j1)*(nz_))+(k1))
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            for (int k = 0; k < nz_; ++k) {
                int global_idx = (i+xs)*NY*NZ + (j+ys)*NZ + k;
                for (int t = 0; t < 19; ++t) {
                    A_[local_idx(i,j,k)*19 + t] = A[global_idx*19 + t];
                }
                b_[local_idx(i,j,k)] = b[global_idx];
                //cout << local_idx(i+1,j+1,k+1) << endl;
                x_[local_idx(i,j,k)] = x[global_idx];
            }
        }
    }
   /* for (int i = 1; i < nx_+1; ++i) {
        for (int j = 1; j < ny_+1; ++j) {
            for (int k = 1; k < nz_+1; ++k) {
                int global_idx = (i+xs)*NY*NZ + (j+ys)*NZ + k;
                x_[(((i)*ny_*nz_)+((j)*nz_)+(k))] = x[global_idx];
            }
        }
    }*/
   /* if (rank == 0) {
        for (int i = 0; i < nx_+2; ++i) {
            cout << x_[(i)*(ny_+2)*(nz_+2)] <<endl;
        
        }
    }
   */
    //return 0;
    MPI_Barrier(MPI_COMM_WORLD);
    gcr(x_, A_, b_, nx_, ny_, nz_, px, py, NX, NY, NZ, PX, PY);
    MPI_Barrier(MPI_COMM_WORLD);

    //double *Ax = Aproduct(A, x, NX, NY, NZ);
    //double *diff = Matrixminus(Ax, b, NX, NY, NZ);
    //double norm2 = dotproduct(diff, diff, size);
    //cout << norm2 << endl;
    gettimeofday(&t2, NULL);
    if (rank == 0) {
        float time = (t2.tv_sec - t1.tv_sec)*1e3 + (t2.tv_usec-t1.tv_usec)*1e-3;
        cout << time << "ms" << endl;
    }

    MPI_Finalize();

    return 0;

}
