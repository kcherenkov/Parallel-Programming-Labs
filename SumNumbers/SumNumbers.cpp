/* This program sums all rows in an array using MPI parallelism.
    * The root process acts as a master and sends a portion of the
    * array to each child process.  Master and child processes then
    * all calculate a partial sum of the portion of the array assigned
    * to them, and the child processes send their partial sums to 
    * the master, who calculates a grand total.
    **/

#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
	int myid, numprocs, n, count, remainder, myBlockSize, result;
	int* data = NULL;
	int* sendcounts = NULL;
	int* displs = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (0 == myid)
	{
		cout << "Summation of numbers with simple 2-step algorithm" << endl << endl;
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

	MPI_Reduce(&total, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Finalize();

	if (myid == 0)
	{
		cout << "Computed sum: " << result << endl;
#if _DEBUG
		system("pause");
#endif
	}
	return EXIT_SUCCESS;
}