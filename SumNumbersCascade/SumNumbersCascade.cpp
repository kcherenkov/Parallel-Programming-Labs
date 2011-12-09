#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
	int myid, numprocs, n, count, remainder, myBlockSize;
	int* data = NULL;
	int* sendcounts = NULL;
	int* displs = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (0 == myid)
	{
		cout << "Summation of numbers with advanced cascade algorithm" << endl << endl;
		cout << "Enter the quantity of elements for summation:" << endl;
		cin >> n;
		data = new int[n];

		int sum = 0;
		for (int i = 0; i < n; ++i)
		{
			data[i] = rand() % 100;
			cout << data[i] << ' ';
			sum += data[i];
		}
		cout << endl;
		cout << "Exact sum: " << sum << endl;

		sendcounts = new int[numprocs];
		displs = new int[numprocs];
		
		count = n / numprocs;
		remainder = n - count * numprocs;
		int prefixSum = 0;
		for (int i = 0; i < numprocs; ++i)
		{
			sendcounts[i] = (i < remainder) ? count + 1 : count;
			displs[i] = prefixSum;
			prefixSum += sendcounts[i];
		}
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 != myid)
	{
		count = n / numprocs;
		remainder = n - count * numprocs;
	}
	myBlockSize = myid < remainder ? count + 1 : count;

	int* res = new int[myBlockSize];

	MPI_Scatterv(data, sendcounts, displs, MPI_INT, res, myBlockSize, MPI_INT, 0, MPI_COMM_WORLD);

	if (0 == myid)
	{
		delete[] sendcounts;
		delete[] displs;
		delete[] data;
	}

	int total = 0;
#pragma omp parallel for reduction(+:total)
	for (int i = 0; i < myBlockSize; ++i)
		total += res[i];
	delete[] res;

	int temp;
	for (int i = 1; i < numprocs; i *= 2)
	{
		if (myid % (i*2) != 0)
		{
			MPI_Send(&total, 1, MPI_INT, myid - i, 0, MPI_COMM_WORLD);
			break;
		}
		else if (myid + i < numprocs)
		{
			MPI_Recv(&temp, 1, MPI_INT, myid + i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total += temp;
		}
	}

	MPI_Finalize();

	if (myid == 0)
	{
		cout << "Computed sum: " << total << endl;
#if _DEBUG
		system("pause");
#endif
	}
	return EXIT_SUCCESS;
}