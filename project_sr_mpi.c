#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

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
            fscanf(fp,"%g",&data[i][j].r);
            data[i][j].i = 0.00;
        }
    }

    fclose(fp);
}

void transpose(complex data[N][N]){
    int i, j;
    for (i = 0; i < N; i++){
      for(j = 0 ; j < N ; j++){
         C_SWAP(data[j][i], data[i][j]);
      }
    }
}

void mmpoint(complex data1[N][N], complex data2[N][N], complex data3[N][N]){

    int i, j;

    float real, imag;

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            data3[i][j].r = (data1[i][j].r * data2[i][j].r) - (data1[i][j].i * data2[i][j].i);
            data3[i][j].i = (data1[i][j].r * data2[i][j].i) + (data1[i][j].i * data2[i][j].r);
        }
    }
}

void printfile(char fileName[15], complex data[N][N]){

    FILE *fp = fopen(fileName, "w");

    int i, j;

    for (i=0;i<N;i++) {
        for (j=0;j<N;j++){
            fprintf(fp,"   %.7e",data[i][j].r);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

int main(int argc, char **argv){
    int my_rank, p, source = 0, dest;

    complex data1[N][N], data2[N][N], data3[N][N];
    complex *vec;

    char fileName1[15] = "sample/1_im1";
    char fileName2[15] = "sample/1_im2";
    char fileName3[15] = "mpi_out_test";

    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &p);


    int i,j;

    double startTime, stopTime;

    //Starting and send rows of data1, data2

    int offset;

    int tag = 345;

    int rows = N/p;

    int lb = my_rank*rows;
    int hb = lb + rows;

    printf("%d have lb = %d and hb = %d\n", my_rank, lb, hb);

    //Starting and send rows of data1, data2

    if(my_rank == 0){
        getData(fileName1, data1);
        getData(fileName2, data2);

        /* Start Clock */
        printf("\nStarting clock.\n");
        startTime = MPI_Wtime();

        for(i=1;i<p;i++){
            offset=i*rows;
                //printf("source sending data1 to %d\n", dest);
                MPI_Send(&data1[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
                //printf("source sending data2 to %d\n", dest);
                MPI_Send(&data2[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
        }
    }else{
        //printf("%d reading data1\n", my_rank);
        MPI_Recv(data1[lb], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
        //printf("%d reading data2\n", my_rank);
        MPI_Recv(data2[lb], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);

    }

    //Doing fft1d forward for data1 and data2 rows

    vec = (complex *)malloc(N * sizeof(complex));

    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            vec[j] = data1[i][j];
        }
        printf("vec1[%d] = %f\n",i, vec[0].r);
        c_fft1d(vec, N, -1);
        for (j=0;j<N;j++) {
            data1[i][j] = vec[j];
        }
    }

    free(vec);

    vec = (complex *)malloc(N * sizeof(complex));

    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            vec[j] = data2[i][j];
        }
        c_fft1d(vec, N, -1);
        for (j=0;j<N;j++) {
            data2[i][j] = vec[j];
        }
    }

    free(vec);


    //Receving rows of data1, data2

    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            //printf("%d reading data1\n", my_rank);
            MPI_Recv(&data1[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
            //printf("%d reading data2\n", my_rank);
            MPI_Recv(&data2[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
        }
    }else{
        //printf("source sending data1 to %d\n", dest);
        MPI_Send(&data1[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
        //printf("source sending data2 to %d\n", dest);
        MPI_Send(&data2[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);

    }

    //Starting and send columns of data1, data2

    if(my_rank == 0){
        //transpose(data1);
        //transpose(data2);

        for(i=1;i<p;i++){
            offset=i*rows;
                //printf("source sending data1 to %d\n", dest);
                MPI_Send(&data1[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
                //printf("source sending data2 to %d\n", dest);
                MPI_Send(&data2[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
        }
    }else{
        //printf("%d reading data1\n", my_rank);
        MPI_Recv(&data1[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
        //printf("%d reading data2\n", my_rank);
        MPI_Recv(&data2[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);

    }

    //Doing fft1d forward for data1 and data2 columns

    vec = (complex *)malloc(N * sizeof(complex));

    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            vec[j] = data1[i][j];
        }
        c_fft1d(vec, N, -1);
        for (j=0;j<N;j++) {
            data1[i][j] = vec[j];
        }
    }

    free(vec);

    vec = (complex *)malloc(N * sizeof(complex));

    for (j=0;j<N;j++) {
        for (i=lb;i<hb;i++) {
            vec[j] = data2[i][j];
        }
        c_fft1d(vec, N, -1);
        for (i=lb;i<hb;i++) {
            data2[i][j] = vec[j];
        }
    }

    free(vec);

    //Receving columns of data1, data2

    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            //printf("%d reading data1\n", my_rank);
            MPI_Recv(&data1[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
            //printf("%d reading data2\n", my_rank);
            MPI_Recv(&data2[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
        }
    }else{
        //printf("source sending data1 to %d\n", dest);
        MPI_Send(&data1[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
        //printf("source sending data2 to %d\n", dest);
        MPI_Send(&data2[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);

    }


    if(my_rank == 0){
        //transpose(data1);
        //transpose(data2);
        mmpoint(data1, data2, data3);
    }

    //Starting and send rows of data1, data2

    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
                //printf("source sending data1 to %d\n", dest);
                MPI_Send(&data3[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
        }
    }else{
        //printf("%d reading data1\n", my_rank);
        MPI_Recv(&data3[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);

    }

    //Doing fft1d forward for data1 and data2 rows

    vec = (complex *)malloc(N * sizeof(complex));

    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            vec[j] = data3[i][j];
        }
        c_fft1d(vec, N, 1);

        for (j=0;j<N;j++) {
            data3[i][j] = vec[j];
        }
    }

    free(vec);


    //Receving rows of data1, data2

    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            //printf("%d reading data1\n", my_rank);
            MPI_Recv(&data3[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
        }
    }else{
        //printf("source sending data1 to %d\n", dest);
        MPI_Send(&data3[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);

    }

    //Starting and send columns of data1, data2

    if(my_rank == 0){
        //transpose(data3);

        for(i=1;i<p;i++){
            offset=i*rows;
                //printf("source sending data1 to %d\n", dest);
                MPI_Send(&data3[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
        }
    }else{
        //printf("%d reading data1\n", my_rank);
        MPI_Recv(&data3[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);

    }

    //Doing fft1d forward for data1 and data2 columns

    vec = (complex *)malloc(N * sizeof(complex));

    for (i=lb;i<hb;i++) {
        for (j=0;j<N;j++) {
            vec[j] = data3[i][j];
        }
        c_fft1d(vec, N, 1);
        for (j=0;j<N;j++) {
            data3[i][j] = vec[j];
        }
    }

    free(vec);

    //Receving columns of data1, data2

    if(my_rank == 0){
        for(i=1;i<p;i++){
            offset=i*rows;
            MPI_Recv(&data3[offset][0], rows*N, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
        }
    }else{
        MPI_Send(&data3[lb][0], rows*N, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);

    }

    if(my_rank == 0){
        //transpose(data3);
        /* Stop Clock */
        stopTime = MPI_Wtime();

        printf("\nElapsed time = %lf s.\n",(stopTime - startTime));
        printf("--------------------------------------------\n");
    }

    MPI_Finalize();

    printfile(fileName3, data3);

    return 0;
}
