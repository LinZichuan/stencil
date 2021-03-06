#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
using namespace std;

#define mp(a, b) make_pair<int,int>(a,b)
#define BACK(idx, N) (idx-1+N)%N
#define FORW(idx, N) (idx+1)%N

double* Aproduct(double *A, double *x, int nx, int ny, int nz) {
    int size = nx * ny * nz;
    int slice = ny * nz;
    int a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18;
    double *Ax = new double[size];
    int nxslice = (nx-1)*slice;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int idx = (i*slice + j*nz + k)*19;
                int a0 = idx/19;
                //i-1,j,k
                 if (i == 0) {
                    a1 = a0 + nxslice; 
                } else {
                    a1 = a0 - slice;
                }
                //i+1,j,k
                if(i == nx-1){
                    a2 = a0 - nxslice;
                }
                else{
                    a2 = a0 + slice;
                }
                //i,j-1,k
                if(j == 0){
                    a3 = a0 + ((i+nx/2)%nx-i)*slice;
                }
                else{
                    a3 = a0 - nz;
                }
                //i,j+1,k
                if(j == ny-1){
                    a4 = a0 + ((i+nx/2)%nx-i)*slice;
                }
                else{
                    a4 = a0 + nz;
                }
                //i+1,j+1,k
               if( i == nx-1){
                    if(j == ny-1){
                         a5 = a0  + ((i+nx/2+1)%nx-i)*slice;
                    }else{
                        a5 = a0 - nxslice + nz;
                    }
                }else{
                    if(j == ny-1){
                        a5 = a0 + ((i+nx/2+1)%nx-i)*slice;
                    }else{
                        a5 = a0 + slice + nz;
                    }
                }
                //i+1,j-1,k
                if (i == nx-1) {
                    if (j == 0) {
                        a6 = a0 + ((i+nx/2+1)%nx-i)*slice;
                    } else {
                        a6 = a0 - nxslice - nz;
                    }
                } else {
                    if (j == 0) {
                        a6 = a0 + ((i+nx/2+1)%nx-i)*slice;
                    } else {
                        a6 = a0 + slice - nz;
                    }
                }
                //i-1,j-1,k
                if (i == 0) {
                    if (j == 0) {
                        a7 = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a7 = a0 + nxslice - nz;
                    }
                } else {
                    if (j == 0) {
                        a7 = a0  + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a7 = a0 - slice - nz;
                    }
                }
                //i-1,j+1,k
                if (i == 0) {
                    if (j == ny-1) {
                        a8 = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a8 = a0 + nxslice + nz;
                    }
                } else {
                    if (j == ny-1) {
                        a8 = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a8 = a0 - slice + nz;
                    }
                }
                //i,j,k-1
                a9  = a0 - 1;
                //i-1,j,k-1
                if (i == 0) {
                    a10 = a0 + nxslice - 1;
                } else {
                    a10 = a0 - slice - 1;
                }
                //i+1,j,k-1
                if (i == nx-1) {
                    a11 = a0 - nxslice - 1;
                } else {
                    a11 = a0 + slice - 1;
                }
                //i,j-1,k-1
                if (j == 0) {
                    a12 = a0 + ((i+nx/2)%nx-i)*slice - 1;
                } else {
                    a12 = a0 - nz - 1;
                }
                //i,j+1,k-1
                if (j == ny-1) {
                    a13 = a0 + ((i+nx/2)%nx-i)*slice - 1;
                } else {
                    a13 = a0 + nz - 1;
                }
                //i,j,k+1
                a14 = a0 + 1;
                //i-1,j,k+1
                if (i == 0) {
                    a15 = a0 + nxslice + 1;
                } else {
                    a15 = a0 - slice + 1;
                }
                //i+1,j,k+1
                if (i == nx-1) {
                    a16 = a0 - nxslice + 1;
                } else {
                    a16 = a0 + slice + 1;
                }
                //i,j-1,k+1
                if (j == 0) {
                    a17 = a0 + ((i+nx/2)%nx-i)*slice + 1;
                } else {
                    a17 = a0 - nz + 1;
                }
                //i,j+1,k+1
                if (j == ny-1) {
                    a18 = a0 + ((i+nx/2)%nx-i)*slice + 1;
                } else {
                    a18 = a0 + nz + 1;
                }
                if (k == 0) a9 = a10 = a11 = a12 = a13 = 0;
                else if (k == nz-1) a14 = a15 = a16 = a17 = a18 = 0;
                Ax[a0] = A[idx]*x[a0] 
                        + A[idx+1]*x[a1] 
                        + A[idx+2]*x[a2] 
                        + A[idx+3]*x[a3]
                        + A[idx+4]*x[a4]
                        + A[idx+5]*x[a5]
                        + A[idx+6]*x[a6]
                        + A[idx+7]*x[a7]
                        + A[idx+8]*x[a8]
                        + A[idx+9]*x[a9]
                        + A[idx+10]*x[a10]
                        + A[idx+11]*x[a11]
                        + A[idx+12]*x[a12]
                        + A[idx+13]*x[a13]
                        + A[idx+14]*x[a14]
                        + A[idx+15]*x[a15]
                        + A[idx+16]*x[a16]
                        + A[idx+17]*x[a17]
                        + A[idx+18]*x[a18];
            }
        }
    }
    return Ax;
}
bool judge(const pair<int,int> a, const pair<int,int> b) {
    return a.first < b.first;
}
int getcsrA(double *A, int *Ai, int *Aj, double *value, int* diagpos, int nx, int ny, int nz, int xs, int xe, int ys, int ye) {
    int slice = ny*nz;
    int nonzero = 0;
    vector<int> a(19);
    //double a[19];
    int nxslice = (nx-1)*slice;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int idx = (i*slice + j*nz + k)*19;
#define BACK(idx, N) (idx-1+N)%N
#define FORW(idx, N) (idx+1)%N
                int a0 = idx/19;
                a[0]  = i*slice + j*nz + k;
                //i-1,j,k
                 if (i == 0) {
                    a[1] = a0 + nxslice; 
                } else {
                    a[1] = a0 - slice;
                }
                //a[1]  = BACK(i,nx)*slice + j*nz + k;
                //i+1,j,k
                if(i == nx-1){
                    a[2] = a0 - nxslice;
                }
                else{
                    a[2] = a0 + slice;
                }
                // a[2]  = FORW(i,nx)*slice + j*nz + k;
                //i,j-1,k
                if(j == 0){
                    a[3] = a0 + ((i+nx/2)%nx-i)*slice;
                }
                else{
                    a[3] = a0 - nz;
                }
               // a[3]  = i*slice + (j-1)*nz + k;
               // if (j == 0) a[3] = (i+nx/2)%nx*slice + j*nz + k;
                //i,j+1,k
                if(j == ny-1){
                    a[4] = a0 + ((i+nx/2)%nx-i)*slice;
                }
                else{
                    a[4] = a0 + nz;
                }
               // a[4]  = i*slice + (j+1)*nz + k;
               // if (j == ny-1) a[4] = (i+nx/2)%nx*slice + j*nz + k;
                //i+1,j+1,k
               if( i == nx-1){
                    if(j == ny-1){
                         a[5] = a0  + ((i+nx/2+1)%nx-i)*slice;
                    }else{
                        a[5] = a0 - nxslice + nz;
                    }
                }else{
                    if(j == ny-1){
                        a[5] = a0 + ((i+nx/2+1)%nx-i)*slice;
                    }else{
                        a[5] = a0 + slice + nz;
                    }
                }
                //i+1,j-1,k
                if (i == nx-1) {
                    if (j == 0) {
                        a[6] = a0 + ((i+nx/2+1)%nx-i)*slice;
                    } else {
                        a[6] = a0 - nxslice - nz;
                    }
                } else {
                    if (j == 0) {
                        a[6] = a0 + ((i+nx/2+1)%nx-i)*slice;
                    } else {
                        a[6] = a0 + slice - nz;
                    }
                }
                //i-1,j-1,k
                if (i == 0) {
                    if (j == 0) {
                        a[7] = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a[7] = a0 + nxslice - nz;
                    }
                } else {
                    if (j == 0) {
                        a[7] = a0  + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a[7] = a0 - slice - nz;
                    }
                }
                //i-1,j+1,k
                if (i == 0) {
                    if (j == ny-1) {
                        a[8] = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a[8] = a0 + nxslice + nz;
                    }
                } else {
                    if (j == ny-1) {
                        a[8] = a0 + ((i+nx/2-1)%nx-i)*slice;
                    } else {
                        a[8] = a0 - slice + nz;
                    }
                }
                //i,j,k-1
                a[9]  = a0 - 1;
                //i-1,j,k-1
                if (i == 0) {
                    a[10] = a0 + nxslice - 1;
                } else {
                    a[10] = a0 - slice - 1;
                }
                //i+1,j,k-1
                if (i == nx-1) {
                    a[11] = a0 - nxslice - 1;
                } else {
                    a[11] = a0 + slice - 1;
                }
                //i,j-1,k-1
                if (j == 0) {
                    a[12] = a0 + ((i+nx/2)%nx-i)*slice - 1;
                } else {
                    a[12] = a0 - nz - 1;
                }
                //i,j+1,k-1
                if (j == ny-1) {
                    a[13] = a0 + ((i+nx/2)%nx-i)*slice - 1;
                } else {
                    a[13] = a0 + nz - 1;
                }
                //i,j,k+1
                a[14] = a0 + 1;
                //i-1,j,k+1
                if (i == 0) {
                    a[15] = a0 + nxslice + 1;
                } else {
                    a[15] = a0 - slice + 1;
                }
                //i+1,j,k+1
                if (i == nx-1) {
                    a[16] = a0 - nxslice + 1;
                } else {
                    a[16] = a0 + slice + 1;
                }
                //i,j-1,k+1
                if (j == 0) {
                    a[17] = a0 + ((i+nx/2)%nx-i)*slice + 1;
                } else {
                    a[17] = a0 - nz + 1;
                }
                //i,j+1,k+1
                if (j == ny-1) {
                    a[18] = a0 + ((i+nx/2)%nx-i)*slice + 1;
                } else {
                    a[18] = a0 + nz + 1;
                }
                if (k == 0) a[9] = a[10] = a[11] = a[12] = a[13] = -1;
                else if (k == nz-1) a[14] = a[15] = a[16] = a[17] = a[18] = -1;
              
                //sort b
                vector< pair<int,int> > b(19);
                for (int t = 0; t < 19; ++t) {
                    b[t] = make_pair<int,int>(a[t], t);
                }
                sort(b.begin(), b.end(), judge);
                //sort(a.begin(), a.end());
                int nonzero_num = 0;
                //diagpos[a0] = nonzero;
                for (int t = 0; t < 19; ++t) {
                    if (b[t].first == -1) continue;
                    int xx = b[t].first / slice;
                    int yy = b[t].first % slice / nz;
                    if (xx < xs || xx >= xe || yy < ys || yy >= ye) continue;
                    value[nonzero] = A[idx+b[t].second]; ///BUG!!
                    Aj[nonzero] = b[t].first;
                    if (b[t].first == a0) {
                        diagpos[a0] = nonzero;
                    }
                    //cout << "a[t]:" << a[t] << endl;
                    nonzero++;
                    nonzero_num++;
                }
                Ai[a0+1] = Ai[a0] + nonzero_num ;
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
    }
    return R;
}
double dotproduct(double *a, double *b, int n) {
    double res = 0.0;
    double res_sum = 0.0; 
    for (int i = 0; i < n; ++i) {
        res += a[i] * b[i];
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
void preprocess(int *Ai, int *Aj, double *value, int *diagpos, double *R, double *R_, int size) {
    //R=M*R_ , solve R_
    //diag of L = 1
    //Ly = R
    double *y = new double[size];
    for (int i = 0; i < size; ++i) {
        int begin = Ai[i];
        int end = diagpos[i]; //Ai[i+1];
        double tmp = 0.0;
        for (int j = begin; j < end; ++j) {
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
            tmp += value[j] * R_[Aj[j]];
        }
        //R_[i] = (y[i]-tmp) / iluA[mp(i,i)];
        R_[i] = (y[i]-tmp) / value[diagpos[i]];
    }
}
void communicate(double *A, int nx, int ny, int nz, int size, int px, int py, int PX, int PY) {
    //cout << nx << " " << ny << " " << nz << " " << size << endl;
    int ysize = (ny-2) * (nz-2);
    int xsize = (nx-2) * (nz-2);
    int csize = nz-2;
    double yhalo[2][ysize];
    double xhalo[2][xsize];
    double chalo[4][csize];
    double s_yhalo[2][ysize];
    double s_xhalo[2][xsize];
    double s_chalo[4][csize];
    int source[9];
    int dest[9];
    int tag[9];
#define IDX(x,y) x*PY+y
    //west
    dest[1] = IDX((px-1+PX)%PX, py);
    tag[1] = 13;
    //east
    dest[3] = IDX((px+1+PX)%PX, py);
    tag[3] = 31;
    //north
    if (py == 0) {
        dest[2] = IDX((px+PX/2)%PX, py);
        tag[2] = 22;
    } else {
        dest[2] = IDX(px, py-1);
        tag[2] = 24;
    }
    //south
    if (py == PY-1) {
        dest[4] = IDX((px+PX/2)%PX, py);
        tag[4] = 44;
    } else {
        dest[4] = IDX(px, py+1);
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
#define AIDX(x,y,z) (x)*ny*nz+(y)*nz+(z)
    int tmp = 0;
    for (int y = 1; y < ny-1; ++y) {
        for (int z = 1; z < nz-1; ++z) {
            A[AIDX(0, y, z)] = yhalo[0][tmp];
            s_yhalo[0][tmp] = A[AIDX(1, y, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int y = 1; y < ny-1; ++y) {
        for (int z = 1; z < nz-1; ++z) {
            A[AIDX(nx-1, y, z)] = yhalo[1][tmp];
            s_yhalo[1][tmp] = A[AIDX(nx-2, y, z)];
            tmp++;
        }
    }
    //return;
    tmp = 0;
    for (int x = 1; x < nx-1; ++x) {
        for (int z = 1; z < nz-1; ++z) {
            A[AIDX(x, ny-1, z)] = xhalo[0][tmp];
            s_xhalo[0][tmp] = A[AIDX(x, ny-2, z)];
            tmp++;
        }
    }
    tmp = 0;
    for (int x = 1; x < nx-1; ++x) {
        for (int z = 1; z < nz-1; ++z) {
            A[AIDX(x, 0, z)] = xhalo[1][tmp];
            s_xhalo[1][tmp] = A[AIDX(x, 1, z)]; 
            tmp++;
        }
    }
    tmp = 0;
    for (int z = 1; z < nz-1; ++z) {
        A[AIDX(0,0,z)] = chalo[0][tmp];
        A[AIDX(0,ny-1,z)] = chalo[1][tmp];
        A[AIDX(nx-1,ny-1,z)] = chalo[2][tmp];
        A[AIDX(nx-1,0,z)] = chalo[3][tmp];
        s_chalo[0][tmp] = A[AIDX(1,1,z)];
        s_chalo[1][tmp] = A[AIDX(1,ny-2,z)];
        s_chalo[2][tmp] = A[AIDX(nx-2,ny-2,z)];
        s_chalo[3][tmp] = A[AIDX(nx-2,1,z)];
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

}
double check_sum(double *value, int nonzero) {
    double sum = 0.0;
    double res_sum = 0.0;
    for (int i = 0; i < nonzero; ++i) {
        sum += value[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&sum,&res_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return res_sum;
}
void gcr(double *x, double *A, double *b, int nx, int ny, int nz, int size, int px, int py, int NX, int NY, int PX, int PY) {
    //sync x
    communicate(x, nx, ny, nz, size, px, py, PX, PY);
    double *Ax = Aproduct(A, x, nx, ny, nz);
    double *R = Matrixminus(b, Ax, nx, ny, nz);
    double Rnorm2 = dotproduct(R, R, size);
    cout << Rnorm2 << endl; //1e-5
    double *R_ = new double[size];
    double *p0 = new double[size];
    
    int *Ai = new int[size+1];
    Ai[0] = 0;
    int *Aj = new int[size*19+1];
    double *value = new double[size*19+1];
    cout << 1 << endl;
    int *diagpos = new int[size];
    //TODO
    int xs = px*NX/PX;
    for (int i = 0; i < px; ++i) {
        xs += ((px*2)%PX<(NX%PX)?1:0);
    }
    int xe = xs+nx;
    int ys = py*NY/PY + (py<(NY%PY)?py:(NY%PY));
    int ye = ys+ny;
    int nonzero = getcsrA(A, Ai, Aj, value, diagpos, nx, ny, nz, xs, xe, ys, ye);
    cout << 2 << endl;
    //ILU
    for (int i = 0; i < size; ++i) {
        int begin = Ai[i];
        int end = Ai[i+1];
        for (int j = begin; j < diagpos[i]; ++j) {
            //TODO: judge, in ns~ne
            int jj = diagpos[Aj[j]];
            value[j] /= value[jj];
            int jbegin = j+1;
            int jend = end;
            int kbegin = jj+1;
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
    double res_sum = check_sum(value, nonzero);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get Process ID
    if (rank == 0) {
        cout << "CHECKSUM:" << res_sum << endl;
    }
    cout << 3 << endl;
    //preprocess
    //TODO???
    preprocess(Ai, Aj, value, diagpos, R, R_, size);
    cout << 4 << endl;
    for (int i = 0; i < size; ++i) {
        p0[i] = R_[i]; //p = R_
    }
    vector<double*> vecAp;
    vector<double*> vecp;
    double *Ap0 = Aproduct(A, p0, nx, ny, nz);
    vecp.push_back(p0);
    vecAp.push_back(Ap0);
    int k = 5;
    double vecAp_2[150];
    for (int i = 1; i < 21; ++i) {
        //sync R_
        communicate(R_, nx, ny, nz, size, px, py, PX, PY);
        vecAp_2[i-1] = dotproduct(vecAp[i-1], vecAp[i-1], size);
        double alpha = dotproduct(R, vecAp[i-1], size) / vecAp_2[i-1];
        cout << i << endl;
        saxpy(x, vecp[i-1], alpha, size);  //x = x + alpha * p;
        saxpy(R, vecAp[i-1], -alpha, size);  //R = R - alpha * Ap;
        //preprocess
        preprocess(Ai, Aj, value, diagpos, R, R_, size);
        int jbegin = int((i-1)/k)*k;
        double *beta = new double[i-jbegin];
        double *AR_ = Aproduct(A, R_, nx, ny, nz);
        for (int j = jbegin; j < i; ++j) {
            beta[j-jbegin] = - dotproduct(AR_, vecAp[j], size) / vecAp_2[j];
        }
        double *p = add(R_, sigma(beta, vecp, jbegin, i, size), size);
        double *Ap = add(AR_, sigma(beta, vecAp, jbegin, i, size), size);
        vecp.push_back(p);
        vecAp.push_back(Ap);
    }
}
bool load_data(double *A, double *b, double *x, int size) {
    FILE *fp, *fpb, *fpx;;
    string binfile = "case_360x180x38/A.bin";
    string binfileb = "case_360x180x38/b.bin";
    string binfilex = "case_360x180x38/x0.bin";
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
    timeval t1, t2;
    gettimeofday(&t1, NULL);

    int NX=360, NY=180, NZ=38;
    int slice = NY*NZ;
    int size = 360*180*38;
    double *A = new double[size*19];
    double *b = new double[size];
    double *x = new double[size];
    if (!load_data(A, b, x, size)) {
        cout << "Load Data Failed!" << endl;
        return 0;
    }

    int PX = 4, PY = 3;

    MPI_Init(&argc, &argv);
    int rank, nprocs;
    char hostname[100];
    gethostname(hostname, 100);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get Process ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Get # of Processes
    cout<<"rank number:"<<rank<<endl;
    //parallel
    int px = rank/PY;
    int py = rank%PY;
    int nx_ = NX/PX + (((px*2)%PX<(NX%PX))?1:0);
    int ny_ = NY/PY + (py<(NY%PY)?1:0);
    int nz_ = NZ;

    int size_ = (nx_+2) * (ny_+2) * (nz_+2);
    double *A_ = new double[size_*19];
    double *b_ = new double[size_];
    double *x_ = new double[size_];

    MPI_Barrier(MPI_COMM_WORLD);
    gcr(x_, A_, b_, nx_, ny_, nz_, size_, px, py, NX, NY, PX, PY);
    MPI_Barrier(MPI_COMM_WORLD);

    double *Ax = Aproduct(A, x, NX, NY, NZ);
    double *diff = Matrixminus(Ax, b, NX, NY, NZ);
    double norm2 = dotproduct(diff, diff, size);
    cout << norm2 << endl;
    gettimeofday(&t2, NULL);
    int time = t2.tv_sec - t1.tv_sec;
    cout << time << "s" << endl;

    MPI_Finalize();

    return 0;

}
