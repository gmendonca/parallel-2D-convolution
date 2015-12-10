MPICC:=mpicc

program: test.c

	$(MPICC) -c project_sr_mpi.c

	$(MPICC) -o project_sr_mpi project_sr_mpi.o

#	 $(MPICC) -c project_openmp_mpi.c -fopenmp

#	 $(MPICC) -o project_openmp_mpi project_openmp_mpi.o

	$(MPICC) -c project_bcast_mpi.c

	$(MPICC) -o project_bcast_mpi project_bcast_mpi.o

clean:
	rm -rf *.o
