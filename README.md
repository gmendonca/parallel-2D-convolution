# parallel-2D-convolution

This project entails implementing two-dimensional convolution on parallel computers using
different parallelization techniques and models.

This project uses [MPICH2](https://www.mpich.org/) and [OpenMP](http://openmp.org/wp/)

[Reference](http://www.clarku.edu/~djoyce/complex/mult.html) for pointwise multiplication of a matrix with a struct of complex numbers.

[Reference1](http://stackoverflow.com/questions/4377127/sending-blocks-of-2d-array-rows-using-mpi-in-c) and [Reference2](http://stackoverflow.com/questions/5104847/mpi-bcast-a-dynamic-2d-array) for sending 2d array over using MPI.

[Reference](http://stackoverflow.com/questions/20228772/mpi-c-send-2d-array-of-structures) for creating a MPI data type


## Description of approaches

project.c - Sequential version ([Reference](http://paulbourke.net/miscellaneous/dft/))
project_transpose.c - Sequential version using matrix transpose
project_sr_mpi.c - Parallel version using MPI_Send and MPI_Recv.
project_bcast_mpi.c - Parallel version using MPI_Bcast.
project_openmp_mpi.c - Parallel version using MPI and OpenMP.
project_tasks_mpi.c - Parallel version using 2 processors for each task, has to be run with 8 processors.

## Compile

To compile just run the Makefile:

```bash
$ make
```

To clean:
```bash
make clean
```
