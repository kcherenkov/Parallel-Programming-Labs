/* Matrix multiplication */

#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	int rank, numprocs, M1, N1, M2, N2, count, remainder, myRowsSize;
	int* matrix1 = NULL;
	int* matrix2;
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
		cout << "Matrix multiplication" << endl << endl;
		cout << "Enter the number of first matrix rows:" << endl;
		cin >> M1;
		if (M1 < 1)	return EXIT_FAILURE;
		cout << "Enter the number of first matrix columns:" << endl;
		cin >> N1;
		if (N1 < 1)	return EXIT_FAILURE;
		M2 = N1;
		cout << "Enter the number of second matrix columns:" << endl;
		cin >> N2;
		if (N2 < 1)	return EXIT_FAILURE;

		cout << endl;
		matrix1 = new int[M1 * N1];
		// generate matrix
		cout << "matrix 1:" << endl;
		for (int i = 0; i < M1; ++i)
		{
			for (int j = 0; j < N1; ++j)
			{
				matrix1[N1 * i + j] = rand() % 100;
				cout << matrix1[N1 * i + j] << '\t';
			}
			cout << endl;
		}
		cout << endl;

		matrix2 = new int[M2 * N2];
		cout << "matrix 2:" << endl;
		for (int i = 0; i < M2; ++i)
		{
			for (int j = 0; j < N2; ++j)
			{
				matrix2[N2 * i + j] = rand() % 100;
				cout << matrix2[N2 * i + j] << '\t';
			}
			cout << endl;
		}
		cout << endl;


		sendcounts = new int[numprocs];
		senddispls = new int[numprocs];
		recvcounts = new int[numprocs];
		recvdispls = new int[numprocs];

		count = M1 / numprocs;
		remainder = M1 - count * numprocs;
		int prefixSum = 0;
		for (int i = 0; i < numprocs; ++i)
		{
			int t1 = (i < remainder) ? count + 1 : count;
			recvcounts[i] = t1 * N2;
			sendcounts[i] = t1 * N1;
			int t2 = prefixSum;
			recvdispls[i] = t2 * N2;
			senddispls[i] = t2 * N1;
			prefixSum += t1;
		}
	}

	MPI_Bcast(&M1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N2, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 != rank)
		matrix2 = new int[M2 * N2];
	MPI_Bcast(matrix2, M2 * N2, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 != rank)
	{
		count = M1 / numprocs;
		remainder = M1 - count * numprocs;
	}
	myRowsSize = rank < remainder ? count + 1 : count;
	int* matrixPart = new int[myRowsSize * N1];
	MPI_Scatterv(matrix1, sendcounts, senddispls, MPI_INT, matrixPart, myRowsSize * N1, MPI_INT, 0, MPI_COMM_WORLD);
	if (0 == rank)
	{
		delete[] sendcounts;
		delete[] senddispls;
		delete[] matrix1;
	}

	int* resultPart = new int[myRowsSize * N2];
#pragma omp parallel for
	for (int i = 0; i < myRowsSize; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			resultPart[i * N2 + j] = 0;
			for (int k = 0; k < N1; ++k)
			{
				resultPart[i * N2 + j] += matrixPart[i * N1 + k] * matrix2[k * N2 + j];
			}
		}
	}
	delete[] matrixPart;
	delete[] matrix2;

	if (0 == rank)
		result = new int[M1 * N2];
	MPI_Gatherv(resultPart, myRowsSize * N2, MPI_INT, result, recvcounts, recvdispls, MPI_INT, 0, MPI_COMM_WORLD); 
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
		for (int i = 0; i < M1; ++i)
		{
			for (int j = 0; j < N2; ++j)
				cout << result[i * N2 + j] << '\t';
			cout << endl;
		}
		delete[] result;
#if _DEBUG
		system("pause");
#endif
	}
	return EXIT_SUCCESS;
}