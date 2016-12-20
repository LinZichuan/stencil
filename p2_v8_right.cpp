#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
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
int getiluA(double *A, int *Ai, int *Aj, double *value, int* diagpos, int nx, int ny, int nz) {
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
                    //if ((t>=9 && t<=13 && k==0) || (t>=14 && k==nz-1)) {
                    //    continue;
                    //}
                    //if (a[t] == -1) continue;
                    //value[nonzero] = A[idx+t]; ///BUG!!
                    //Aj[nonzero] = a[t];
                    if (b[t].first == -1) continue;
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
    for (int i = 0; i < n; ++i) {
        res += a[i] * b[i];
    }
    return res;
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
            //tmp += iluA[mp(i,j)] * y[j];
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
            //tmp += iluA[mp(i,j)] * R_[j];
        }
        //R_[i] = (y[i]-tmp) / iluA[mp(i,i)];
        R_[i] = (y[i]-tmp) / value[diagpos[i]];
    }
}
//Ap = add(AR_, sigma(beta, vecAp, jbegin, i));
int main() {
    timeval t1, t2;
    gettimeofday(&t1, NULL);
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
    int nx=360, ny=180, nz=38;
    int slice = ny*nz;
    int size = 360*180*38;
    double *A = new double[size*19];
    double *b = new double[size];
    double *x = new double[size];
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
    //cout.precision(27);
    //for (int i = 0; i < size; ++i) {
    //    cout << A[i] << endl;
    //    cout << b[i] << endl;
    //    cout << x[i] << endl;
    //    break;
    //}
    fclose(fp);
    fclose(fpb);
    fclose(fpx);

    //double *Ax = new double[size];
    //double *R = new double[size];
    double *Ax = Aproduct(A, x, nx, ny, nz);
    double *R = Matrixminus(b, Ax, nx, ny, nz);
    double Rnorm2 = dotproduct(R, R, size);
    cout << Rnorm2 << endl; //1e-5
    double *R_ = new double[size];
    //double *Ap = new double[size];
    double *p0 = new double[size];
    
    //ILU
    int *Ai = new int[size+1];
    Ai[0] = 0;
    int *Aj = new int[size*19+1];
    double *value = new double[size*19+1];
    cout << 1 << endl;
    int *diagpos = new int[size];
    int nonzero = getiluA(A, Ai, Aj, value, diagpos, nx, ny, nz);
    //for (int i = 0; i < size; ++i) {
    //    cout << value[diagpos[i]] << endl;
    //}
    //return 0;
    cout << 2 << endl;
    for (int i = 0; i < size; ++i) {
        int begin = Ai[i];
        int end = Ai[i+1];
        //int jbegin = Ai[ Aj[begin+1] ];
        //int jend = Ai[ Aj[begin+1]+1 ];
        //cout << i << endl;
        //cout << begin << "," << end << " " << jbegin << "," << jend << endl;
        //cout << begin << " " << diagpos[i] << endl;
        for (int j = begin; j < diagpos[i]; ++j) {
            int jj = diagpos[Aj[j]];
            value[j] /= value[jj];
            //cout << value[j] << endl;
            int jbegin = j+1;
            int jend = end;
            int kbegin = jj+1;
            int kend = Ai[Aj[j]+1];
            //int kbegin = Ai[Aj[j+1]];
            //int kend = Ai[Aj[j]+1];
            //int idxj=begin, idxk=jbegin;
            int idxj=jbegin, idxk=kbegin;
            //cout << begin << "," << end << " " << jbegin << "," << jend << endl;
            while (idxj < jend && idxk < kend) {
                //if (Aj[idxj] < i) {
                //    value[idxj] = value[idxj] / value[diagpos[Aj[idxj]]]; //TODO: sort, save diagpos for diag element in getiluA
                //}
                if (Aj[idxj] == Aj[idxk]) {
                    //cout << idxk << " " << idxj << endl;
                    //if (i == Aj[idxj]) {
                    //    cout << "aaaaaaaaaaa" << endl;
                    //}
                    //cout << value[j] << endl;
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
    double sum = 0.0;
    for (int i = 0; i < nonzero; ++i) {
        sum += value[i];
    }
    cout << sum << endl;
    cout << 3 << endl;
    //preprocess
    //preprocess(iluA, R, R_, Apos, size);
    preprocess(Ai, Aj, value, diagpos, R, R_, size);
    cout << 4 << endl;

    for (int i = 0; i < size; ++i) {
        //R_[i] = R[i]; //R_ = M^-1 * R //TODO
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
        vecAp_2[i-1] = dotproduct(vecAp[i-1], vecAp[i-1], size);
        double alpha = dotproduct(R, vecAp[i-1], size) / vecAp_2[i-1];
        cout << i << endl;
        saxpy(x, vecp[i-1], alpha, size);  //x = x + alpha * p;
        saxpy(R, vecAp[i-1], -alpha, size);  //R = R - alpha * Ap;
        //TODO R_ = Aproduct(inv(M), R, nx, ny, nz);  //R_ = M^-1 * R;
        //preprocess(iluA, R, R_, Apos, size);
        preprocess(Ai, Aj, value, diagpos, R, R_, size);
        //for (int t = 0; t < size; ++t) {
        //    R_[t] = R[t]; //R_ = M^-1 * R //TODO
        //}
        int jbegin = int((i-1)/k)*k;
        double *beta = new double[i-jbegin];
        double *AR_ = Aproduct(A, R_, nx, ny, nz);
        for (int j = jbegin; j < i; ++j) {
            //beta[j-jbegin] = - dotproduct(AR_, vecAp[j], size) / dotproduct(vecAp[j], vecAp[j], size);
            beta[j-jbegin] = - dotproduct(AR_, vecAp[j], size) / vecAp_2[j];
        }
        double *p = add(R_, sigma(beta, vecp, jbegin, i, size), size);
        double *Ap = add(AR_, sigma(beta, vecAp, jbegin, i, size), size);
        vecp.push_back(p);
        vecAp.push_back(Ap);
    }
    Ax = Aproduct(A, x, nx, ny, nz);
    double *diff = Matrixminus(Ax, b, nx, ny, nz);
    //for (int i = 0; i < size; ++i)
    //    cout << diff[i] << " ";
    //cout << endl;
    double norm2 = dotproduct(diff, diff, size);
    cout << norm2 << endl;
    gettimeofday(&t2, NULL);
    int time = t2.tv_sec - t1.tv_sec;
    cout << time << "s" << endl;

    return 0;

}
