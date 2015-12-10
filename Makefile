MPICC:=mpicc
CC:=gcc

program: project.c project_transpose.c project_sr_mpi.c project_openmp_mpi.c project_bcast_mpi.c project_tasks_mpi.c

	$(CC) project.c -o project.o -lm

	$(CC) project_transpose.c -o project_transpose.o -lm

	$(MPICC) -c project_sr_mpi.c

	$(MPICC) -o project_sr_mpi project_sr_mpi.o -lm

	$(MPICC) -c project_openmp_mpi.c -fopenmp

	$(MPICC) -o project_openmp_mpi project_openmp_mpi.o -lm

	$(MPICC) -c project_bcast_mpi.c

	$(MPICC) -o project_bcast_mpi project_bcast_mpi.o -lm

	$(MPICC) -c project_tasks_mpi.c

	$(MPICC) -o project_tasks_mpi project_tasks_mpi.o -lm

clean:
	rm -rf *.o
	rm -rf *_mpi
