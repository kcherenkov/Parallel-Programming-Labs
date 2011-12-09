#include <mpi.h>
#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

#pragma omp parallel
	{
		int iam, np;
		np = omp_get_num_threads();
		iam = omp_get_thread_num();

#pragma omp critical
		cout << "Hello from thread " << iam << " out of " << np << " from process " << rank << " out of " << numprocs << " on " << processor_name << endl;
	}

	MPI_Finalize();

#if _DEBUG
	system("pause");
#endif
	return EXIT_SUCCESS;
}
