MPICC:=mpicc
CC:=gcc

program: test.c

	$(CC) project.c

	$(CC) project_transpose.c

	$(MPICC) -c project_sr_mpi.c

	$(MPICC) -o project_sr_mpi project_sr_mpi.o

	$(MPICC) -c project_openmp_mpi.c -fopenmp

	$(MPICC) -o project_openmp_mpi project_openmp_mpi.o

	$(MPICC) -c project_bcast_mpi.c

	$(MPICC) -o project_bcast_mpi project_bcast_mpi.o

	$(MPICC) -c project_tasks_mpi.c

	$(MPICC) -o project_tasks_mpi project_tasks_mpi.o

clean:
	rm -rf *.o
