#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include<omp.h>
#include<sys/time.h>
#define nx 360
#define ny 180
#define nz 38
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

#define OMP
#define num_thread 24


void csr_init(double A[nx][ny][nz][19]);
void ilufact();
void solve(double xx[nx][ny][nz], double bb[nx][ny][nz] );
//void LUfact(double A[nx][ny][nz][19], double L[nx][ny][nz][19], double U[nx][ny][nz][19]);
void printM(double x[nx][ny][nz][19])
{
  int i,j,k,d;
  for(i=0;i<1;i++){
    for(j=0;j<1;j++){
      for(k=0;k<2;k++){
        for(d=0;d<19;d++){
          printf("%.30lf\n",x[i][j][k][d]);
        }
      }
    }
  }
}

double DotProduct(double x1[nx][ny][nz], double x2[nx][ny][nz])
{
  int i, j, k;
  double sigma = 0;
#ifdef OMP
#pragma omp parallel for reduction(+:sigma) private(j,k) num_threads(num_thread)
#endif
  for ( i = 0; i < nx; i ++ ) { 
    for ( j = 0; j < ny; j ++ ) {
      for ( k = 0; k < nz; k ++ ) {
        sigma += x1[i][j][k] * x2[i][j][k];
      }
    }
  }

  return sigma;
}

double norm2(double Ax[nx][ny][nz], double b[nx][ny][nz])
{
  int i, j, k;
  double residual = 0;
  struct timeval t[2];

  gettimeofday(&t[0], NULL);
#ifdef OMP
#pragma omp parallel for reduction(+:residual) private(j,k) num_threads(num_thread) 
#endif
  for ( i = 0; i < nx; i ++ ) {
    for ( j = 0; j < ny; j ++ ) {
      for ( k = 0; k < nz; k ++ ) {
        residual = residual + (Ax[i][j][k] - b[i][j][k]) * (Ax[i][j][k] - b[i][j][k]);
      }
    }
  }
  gettimeofday(&t[1], NULL);
  printf("norm-Time:%.6lfs  norm-2:%.30lf\n",TIME(t[0],t[1]),residual);
  return residual;;
}
void Vec_Copy(double x1[nx][ny][nz], double x2[nx][ny][nz])
{
  int i, j, k;

#ifdef OMP
#pragma omp parallel for private(j,k) num_threads(num_thread)
#endif
  for ( i = 0; i < nx; i ++ ) {
    for ( j = 0; j < ny; j ++ ) {
      for ( k = 0; k < nz; k ++ ) {
        x1[i][j][k] = x2[i][j][k];
      }
    }
  }
}
void Vec_Add(double x[nx][ny][nz], double x1[nx][ny][nz], double x2[nx][ny][nz], double a, double b)
{
  int i, j, k;
#ifdef OMP
#pragma omp parallel for private(j,k) num_threads(num_thread)
#endif
  for ( i = 0; i < nx; i ++ ) {
    for ( j = 0; j < ny; j ++ ) {
      for ( k = 0; k < nz; k ++ ) {
        x[i][j][k] = a * x1[i][j][k] + b * x2[i][j][k];
      }
    }
  }
}

void form_Ax(double A[nx][ny][nz][19], double x[nx][ny][nz], double Ax[nx][ny][nz])
{
  int i, j, k;
  double x0[19];
#ifdef OMP
#pragma omp parallel for private(j,k,x0) num_threads(num_thread)
#endif
  for ( i = 0; i < nx; i ++ ) {
    for ( j = 0; j < ny; j ++ ) {
      for ( k = 0; k < nz ; k ++ ) {

        x0[0]  =  x[i][j][k]    ;
        x0[1]  =  x[i-1][j][k]  ; 
        x0[2]  =  x[i+1][j][k]  ; 
        x0[3]  =  x[i][j-1][k]  ; 
        x0[4]  =  x[i][j+1][k]  ; 
        x0[5]  =  x[i+1][j+1][k]; 
        x0[6]  =  x[i+1][j-1][k]; 
        x0[7]  =  x[i-1][j-1][k]; 
        x0[8]  =  x[i-1][j+1][k]; 
        x0[9]  =  x[i][j][k-1]  ; 
        x0[10] =  x[i-1][j][k-1]; 
        x0[11] =  x[i+1][j][k-1]; 
        x0[12] =  x[i][j-1][k-1]; 
        x0[13] =  x[i][j+1][k-1]; 
        x0[14] =  x[i][j][k+1]  ; 
        x0[15] =  x[i-1][j][k+1]; 
        x0[16] =  x[i+1][j][k+1]; 
        x0[17] =  x[i][j-1][k+1]; 
        x0[18] =  x[i][j+1][k+1]; 


        if ( j == 0 )
        {
          x0[3] = x[(i+nx/2)%nx][j][k];
          x0[12] = x[(i+nx/2)%nx][j][k-1];
          x0[17] = x[(i+nx/2)%nx][j][k+1];

          x0[6] = x[(i+1+nx/2)%nx][j][k];
          x0[7] = x[(i-1+nx/2)%nx][j][k];
          if ( i == 0 )
          {
            x0[1] = x[nx-1][j][k];
            x0[10] = x[nx-1][j][k-1];
            x0[15] = x[nx-1][j][k+1];
            x0[8] = x[nx-1][j+1][k];
          }
          else if ( i == nx - 1)
          {
            x0[2] = x[0][j][k];
            x0[11] = x[0][j][k-1];
            x0[16] = x[0][j][k+1];
            x0[5] = x[0][j+1][k];
          }
        }
        else if (j == ny - 1)
        {
          x0[4] = x[(i+nx/2)%nx][j][k];
          x0[13] = x[(i+nx/2)%nx][j][k-1];
          x0[18] = x[(i+nx/2)%nx][j][k+1];

          x0[5] = x[(i+1+nx/2)%nx][j][k];
          x0[8] = x[(i-1+nx/2)%nx][j][k];
          if ( i == 0 )
          {
            x0[1] = x[nx-1][j][k];
            x0[10] = x[nx-1][j][k-1];
            x0[15] = x[nx-1][j][k+1];
            x0[7] = x[nx-1][j-1][k];
          }
          else if ( i == nx - 1)
          {
            x0[2] = x[0][j][k];
            x0[11] = x[0][j][k-1];
            x0[16] = x[0][j][k+1];
            x0[6] = x[0][j-1][k];
          }
        }
        else if ( i == 0 )
        {
          x0[1] = x[nx-1][j][k];
          x0[7] = x[nx-1][j-1][k];
          x0[8] = x[nx-1][j+1][k];
          x0[10] = x[nx-1][j][k-1];
          x0[15] = x[nx-1][j][k+1];
        }
        else if ( i == nx - 1 )
        {
          x0[2] = x[0][j][k];
          x0[5] = x[0][j+1][k];
          x0[6] = x[0][j-1][k];
          x0[11] = x[0][j][k-1];
          x0[16] = x[0][j][k+1];
        }
        if ( k == 0 ) x0[9] = x0[10] = x0[11] = x0[12] = x0[13] = 0;
        else if ( k == nz - 1 ) x0[14] = x0[15] = x0[16] = x0[17] = x0[18] = 0;

        Ax[i][j][k] =   A[i][j][k][0]  * x0[0] \
                        + A[i][j][k][1]  * x0[1] \
                        + A[i][j][k][2]  * x0[2] \
                        + A[i][j][k][3]  * x0[3] \
                        + A[i][j][k][4]  * x0[4] \
                        + A[i][j][k][5]  * x0[5] \
                        + A[i][j][k][6]  * x0[6] \
                        + A[i][j][k][7]  * x0[7] \
                        + A[i][j][k][8]  * x0[8] \
                        + A[i][j][k][9]  * x0[9] \
                        + A[i][j][k][10] * x0[10]\
                        + A[i][j][k][11] * x0[11]\
                        + A[i][j][k][12] * x0[12]\
                        + A[i][j][k][13] * x0[13]\
                        + A[i][j][k][14] * x0[14]\
                        + A[i][j][k][15] * x0[15]\
                        + A[i][j][k][16] * x0[16]\
                        + A[i][j][k][17] * x0[17]\
                        + A[i][j][k][18] * x0[18];

        //        Ax[i][j][k] = A[i][j][k][0]  * x[i][j][k]     \
        +A[i][j][k][1]  * x[i-1][j][k]   \
          +A[i][j][k][2]  * x[i+1][j][k]   \
          +A[i][j][k][3]  * x[i][j-1][k]   \
          +A[i][j][k][4]  * x[i][j+1][k]   \
          +A[i][j][k][5]  * x[i+1][j+1][k] \
          +A[i][j][k][6]  * x[i+1][j-1][k] \
          +A[i][j][k][7]  * x[i-1][j-1][k] \
          +A[i][j][k][8]  * x[i-1][j+1][k] \
          +A[i][j][k][9]  * x[i][j][k-1]   \
          +A[i][j][k][10] * x[i-1][j][k-1] \
          +A[i][j][k][11] * x[i+1][j][k-1] \
          +A[i][j][k][12] * x[i][j-1][k-1] \
          +A[i][j][k][13] * x[i][j+1][k-1] \
          +A[i][j][k][14] * x[i][j][k+1]   \
          +A[i][j][k][15] * x[i-1][j][k+1] \
          +A[i][j][k][16] * x[i+1][j][k+1] \
          +A[i][j][k][17] * x[i][j-1][k+1] \
          +A[i][j][k][18] * x[i][j+1][k+1] ;
      }
    }
  }

}
int main()
{
  int i, j, k, d;
  int r;
  FILE *fp;
  double (*A)[ny][nz][19], (*x0)[ny][nz], (*b)[ny][nz];
  double (*Ax)[ny][nz], (*AR_)[ny][nz];
  double (*x)[nx][ny][nz];
  double (*Ap)[nx][ny][nz];
  double (*R)[nx][ny][nz];
  double (*R_)[nx][ny][nz];

  double (*p)[nx][ny][nz];

  double beta[150];

  struct timeval t[10];


  gettimeofday(&t[0], NULL);

  A  = malloc(nx*ny*nz*19*sizeof(double));
  x0 = malloc(nx*ny*nz*sizeof(double));
  b  = malloc(nx*ny*nz*sizeof(double));

  Ax  = malloc(nx*ny*nz*sizeof(double));
  AR_ = malloc(nx*ny*nz*sizeof(double));

  Ap = malloc(150*nx*ny*nz*sizeof(double));
  R = malloc(150*nx*ny*nz*sizeof(double));
  R_ = malloc(150*nx*ny*nz*sizeof(double));

  x = malloc(151*nx*ny*nz*sizeof(double));
  p = malloc(151*nx*ny*nz*sizeof(double));

  // Initialization all malloc-array for 0 
  //  memset(A, 0, nx*ny*nz*19*sizeof(double));
  //  memset(x0, 0, nx*ny*nz*sizeof(double));
  //  memset(b, 0, nx*ny*nz*sizeof(double));
  //  memset(Ax, 0, nx*ny*nz*sizeof(double));
  //  memset(AR_, 0, nx*ny*nz*sizeof(double));
  //  memset(x, 0, nx*ny*nz*150*sizeof(double));
  //  memset(Ap, 0, nx*ny*nz*150*sizeof(double));
  //  memset(R, 0, nx*ny*nz*150*sizeof(double));
  //  memset(R_, 0, nx*ny*nz*150*sizeof(double));
  //  memset(p, 0, nx*ny*nz*150*sizeof(double));

  // Read ASCII file and store data to Binary file
#if 0
  // ASCII To Binary For A.txt
  if((fp = fopen("A.txt","r")) == NULL)
    printf("fail to open A.txt\n");
  else 
  {
    for ( i = 0; i < nx; i ++ ) {
      for ( j = 0; j < ny; j ++ ) {
        for ( k = 0; k < nz; k ++ ) {
          for ( d = 0; d < 19; d ++ ) {
            fscanf(fp,"%d %d %d %d %lf",&r,&r,&r,&r,&A[i][j][k][d]);
          }
        }
      }
    }
  }
  if((fp = fopen("Ab.txt","wb")) == NULL)
    printf("fail to open Ab.txt\n");
  else
  {
    fwrite(A,8,nx*ny*nz*19,fp);
  }

  // ASCII To Binary For b.txt
  if((fp = fopen("b.txt","r")) == NULL)
    printf("fail to open b.txt\n");
  else 
  {
    for ( i = 0; i < nx; i ++ ) {
      for ( j = 0; j < ny; j ++ ) {
        for ( k = 0; k < nz; k ++ ) {
          fscanf(fp,"%d %d %d %lf",&r,&r,&r,&b[i][j][k]);
        }
      }
    }
  }
  if((fp = fopen("bb.txt","wb")) == NULL)
    printf("fail to open bb.txt\n");
  else
  {
    fwrite(b,8,nx*ny*nz,fp);
  }

  // ASCII To Binary For x0.txt
  if((fp = fopen("x0.txt","r")) == NULL)
    printf("fail to open x0.txt\n");
  else 
  {
    for ( i = 0; i < nx; i ++ ) {
      for ( j = 0; j < ny; j ++ ) {
        for ( k = 0; k < nz; k ++ ) {
          fscanf(fp,"%d %d %d %lf",&r,&r,&r,&x0[i][j][k]);
        }
      }
    }
  }
  if((fp = fopen("x0b.txt","wb")) == NULL)
    printf("fail to open x0b.txt\n");
  else
  {
    fwrite(x0,8,nx*ny*nz,fp);
  }
#endif

  gettimeofday(&t[1], NULL);
  // Read Binary file
  if((fp = fopen("Ab.txt","rb")) == NULL)
    printf("fail to open Ab.txt\n");
  else 
    fread(A,8,nx*ny*nz*19,fp);

  if((fp = fopen("bb.txt","rb")) == NULL)
    printf("fail to open bb.txt\n");
  else 
    fread(b,8,nx*ny*nz,fp);

  if((fp = fopen("x0b.txt","rb")) == NULL)
    printf("fail to open x0b.txt\n");
  else 
    fread(x0,8,nx*ny*nz,fp);

  gettimeofday(&t[2], NULL);
//  {
//    // check sum of A
//    double sum = 0;
//    for ( i = 0; i < nx; i ++ ) {
//      for ( j = 0; j < ny; j ++ ) {
//        for ( k = 0; k < nz; k ++ ) {
//          for ( d = 0; d < 19; d ++ ) {
//            sum += A[i][j][k][d];
//          }
//        }
//      }
//    }
//    printf("sum of A is %.30lf\n", sum);
//  }
//  {
//    // check sum of b
//    double sum = 0;
//    for ( i = 0; i < nx; i ++ ) {
//      for ( j = 0; j < ny; j ++ ) {
//        for ( k = 0; k < nz; k ++ ) {
//          sum += b[i][j][k];
//        }
//      }
//    }
//    printf("sum of b is %.30lf\n", sum);
//  }
//  {
//    // check sum of x0
//    double sum = 0;
//    for ( i = 0; i < nx; i ++ ) {
//      for ( j = 0; j < ny; j ++ ) {
//        for ( k = 0; k < nz; k ++ ) {
//          sum += x0[i][j][k];
//        }
//      }
//    }
//    printf("sum of x0 is %.30lf\n", sum);
//  }


  // 防止越界
  x = &x[1];
  p = &p[1];

  Vec_Copy(x[0],x0); // 赋值:x[0] = x0
  form_Ax(A,x[0],Ax); // Ax = A*x[0] 
  Vec_Add(R[0],b,Ax,1,-1); // R[0] = b - Ax

  // Precondition
    csr_init(A);
    ilufact2();
    solve(R_[0],R[0]);

//  Vec_Copy(R_[0],R[0]); // 赋值:x[0] = x0

  Vec_Copy(p[0],R_[0]); // p[0] = R^[0]
  form_Ax(A,p[0],Ap[0]); // Ap[0] = A*p[0]

  //  printM(A);

  i = 1;
  k = 5;
  printf("Residual 0__");

  gettimeofday(&t[3], NULL);
  while ( norm2(Ax,b) >= 1.0E-10 && i <= 150 ) {
    double alpha;

    //    gettimeofday(&t[0], NULL);
    alpha = DotProduct(R[i-1],Ap[i-1]) / DotProduct(Ap[i-1],Ap[i-1]);
    //    gettimeofday(&t[1], NULL);
    //    printf("DotProduct-Time:%.6lfs\n",TIME(t[0],t[1]));

    Vec_Add(x[i],x[i-1],p[i-1],1,alpha); 
    Vec_Add(R[i],R[i-1],Ap[i-1],1,-alpha);
    // Precondition
    {
            solve(R_[i],R[i]);
     // Vec_Copy(R_[i],R[i]); // 赋值:x[0] = x0
    }
    ////////////////
    j = i - 1;
    form_Ax(A,R_[i],AR_);
    while ( j >= 0 && j >= i - k ) {
      beta[j] = - DotProduct(AR_,Ap[j]) / DotProduct(Ap[j], Ap[j]); 
      j --;
    }  
    j = i - 1;
    Vec_Copy(p[i],R_[i]);
    while ( j >= 0 && j >= i - k ) {
      Vec_Add(p[i],p[i],p[j],1,beta[j]);
      j --;
    }
    j = i - 1;
    Vec_Copy(Ap[i],AR_);
    while ( j >= 0 && j >= i - k ) {
      Vec_Add(Ap[i],Ap[i],Ap[j],1,beta[j]);
      j --;
    }

    form_Ax(A,x[i],Ax);
    printf("Residual %d__",i);
    i ++;
  }
  gettimeofday(&t[4], NULL);
  printf("Time for Malloc:%.6lfs\n", TIME(t[0],t[1]));
  printf("Time for IO:%.6lfs\n", TIME(t[1],t[2]));
  printf("Time for Initialization:%.6lfs\n", TIME(t[2],t[3]));
  printf("Time for Iteration:%.6lfs\n", TIME(t[3],t[4]));
  printf("Time for GCR:%.6lfs\n", TIME(t[0],t[4]));

  return 0;
}

