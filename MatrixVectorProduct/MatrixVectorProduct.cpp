/* Product of a vector by a matrix */

#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	int rank, numprocs, M, N, count, remainder, myRowsSize;
	int* matrix = NULL;
	int* vector;
	int* result = NULL;
	int* sendcounts = NULL;
	int* senddispls = NULL;
	int* recvcounts = NULL;
	int* recvdispls = NULL;

	MPI_Init (&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (0 == rank)
	{
		cout << "Product of a vector by a matrix" << endl << endl;
		cout << "Enter the number of matrix rows:" << endl;
		cin >> M;
		if (M < 1)	return EXIT_FAILURE;
		cout << "Enter the number of matrix columns:" << endl;
		cin >> N;
		if (N < 1)	return EXIT_FAILURE;

		cout << endl;
		matrix = new int[M * N];
		// generate matrix
		cout << "matrix:" << endl;
		for (int i = 0; i < M; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				matrix[N * i + j] = rand() % 100;
				cout << matrix[N * i + j] << '\t';
			}
			cout << endl;
		}
		cout << endl;

		vector = new int[N];
		// generate vector
		cout << "vector:" << endl;
		for (int i = 0; i < N; ++i)
		{
			vector[i] = rand() & 100;
			cout << vector[i] << ' ';
		}
		cout << endl << endl;

		sendcounts = new int[numprocs];
		senddispls = new int[numprocs];
		recvcounts = new int[numprocs];
		recvdispls = new int[numprocs];
		
		count = M / numprocs;
		remainder = M - count * numprocs;
		int prefixSum = 0;
		for (int i = 0; i < numprocs; ++i)
		{
			recvcounts[i] = (i < remainder) ? count + 1 : count;
			sendcounts[i] = recvcounts[i] * N;
			recvdispls[i] = prefixSum;
			senddispls[i] = prefixSum * N;
			prefixSum += recvcounts[i];
		}
	}

	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 != rank)
		vector = new int[N];
	MPI_Bcast(vector, N, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 != rank)
	{
		count = M / numprocs;
		remainder = M - count * numprocs;
	}
	myRowsSize = rank < remainder ? count + 1 : count;
	int* matrixPart = new int[myRowsSize * N];
	MPI_Scatterv(matrix, sendcounts, senddispls, MPI_INT, matrixPart, myRowsSize * N, MPI_INT, 0, MPI_COMM_WORLD);
	if (0 == rank)
	{
		delete[] sendcounts;
		delete[] senddispls;
		delete[] matrix;
	}

	int* resultPart = new int[myRowsSize];
#pragma omp parallel for
	for (int i = 0; i < myRowsSize; ++i)
	{
		resultPart[i] = 0;
		for (int j = 0; j < N; ++j)
		{
			resultPart[i] += matrixPart[i * N + j] * vector[j];
		}
	}
	delete[] matrixPart;
	delete[] vector;

	if (0 == rank)
		result = new int[M];
	MPI_Gatherv(resultPart, myRowsSize, MPI_INT, result, recvcounts, recvdispls, MPI_INT, 0, MPI_COMM_WORLD); 
	delete[] resultPart;
	if (0 == rank)
	{
		delete[] recvcounts;
		delete[] recvdispls;
	}

	MPI_Finalize();

	if (0 == rank)
	{
		cout << "result:" << endl;
		for (int i = 0; i < M; ++i)
			cout << result[i] << ' ';
		cout << endl;
		delete[] result;
#if _DEBUG
		system("pause");
#endif
	}
	return EXIT_SUCCESS;
}