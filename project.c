#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {float r; float i;} complex;
static complex ctmp;

#define C_SWAP(a,b) {ctmp=(a);(a)=(b);(b)=ctmp;}

#define N 512

void c_fft1d(complex *r, int n, int isign)
{
   int     m,i,i1,j,k,i2,l,l1,l2;
   float   c1,c2,z;
   complex t, u;

   if (isign == 0) return;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j)
         C_SWAP(r[i], r[j]);
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* m = (int) log2((double)n); */
   for (i=n,m=0; i>1; m++,i/=2);

   /* Compute the FFT */
   c1 = -1.0;
   c2 =  0.0;
   l2 =  1;
   for (l=0;l<m;l++) {
      l1   = l2;
      l2 <<= 1;
      u.r = 1.0;
      u.i = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;

            /* t = u * r[i1] */
            t.r = u.r * r[i1].r - u.i * r[i1].i;
            t.i = u.r * r[i1].i + u.i * r[i1].r;

            /* r[i1] = r[i] - t */
            r[i1].r = r[i].r - t.r;
            r[i1].i = r[i].i - t.i;

            /* r[i] = r[i] + t */
            r[i].r += t.r;
            r[i].i += t.i;
         }
         z =  u.r * c1 - u.i * c2;

         u.i = u.r * c2 + u.i * c1;
         u.r = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (isign == -1) /* FWD FFT */
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for inverse transform */
   if (isign == 1) {       /* IFFT*/
      for (i=0;i<n;i++) {
         r[i].r /= n;
         r[i].i /= n;
      }
   }
}

void getData(char fileName[15], complex data[N][N]){
    FILE *fp = fopen(fileName, "r");

    int i, j;

    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            fscanf(fp,"%f",&data[i][j].r);
            data[i][j].i = 0.00;
        }
    }

    fclose(fp);
}

void fft2d(complex data[N][N], complex transpose[N][N], int isign){

    int i, j;

    for(i=0;i<N;i++){
        c_fft1d(data[i], N, isign);
    }

    //transpose
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            transpose[j][i] = data[i][j];
        }
    }

    for(i=0;i<N;i++){
        c_fft1d(data[i], N, isign);
    }

    //transpose back
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            data[j][i] = transpose[i][j];
        }
    }

}

void fft2d_new(complex data[N][N], int isign){

    int i, j;

    complex* vec;

    vec = (complex *)malloc(N * sizeof(complex));

    for (j=0;j<N;j++) {
      for (i=0;i<N;i++) {
         vec[i] = data[i][j];
      }
      c_fft1d(vec, N, isign);
      for (i=0;i<N;i++) {
         data[i][j] = vec[i];
      }
   }

   free(vec);

   vec = (complex *)malloc(N * sizeof(complex));

   for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
         vec[j] = data[i][j];
      }
      c_fft1d(vec, N, isign);
      for (j=0;j<N;j++) {
         data[i][j] = vec[j];
      }
   }
}

void mmpoint(complex data1[N][N], complex data2[N][N], complex data3[N][N]){

    int i, j;

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            data3[i][j].r = data1[i][j].r * data2[i][j].r;
        }
    }
}

void printfile(char fileName[15], complex data[N][N]){

    FILE *fp = fopen(fileName, "w");

    int i, j;

    for (i=0;i<N;i++) {
        for (j=0;j<N;j++){
            fprintf(fp,"%6.2f",data[i][j].r);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

int main(){
    complex data1[N][N], data2[N][N], data3[N][N];

    char fileName1[15] = "sample/1_im1";
    char fileName2[15] = "sample/1_im2";
    char fileName3[15] = "out_test";

    getData(fileName1, data1);
    getData(fileName2, data2);

    fft2d(data1, data3, 1);
    fft2d(data2, data3, 1);

    mmpoint(data1, data2, data3);

    fft2d(data3, data1, -1);

    printfile(fileName3, data3);

    return 0;
}
