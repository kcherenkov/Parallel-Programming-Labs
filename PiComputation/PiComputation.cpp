/*
This exercise presents a simple program to determine the value of pi. The algorithm suggested here is chosen for its simplicity. The method evaluates the integral of 4/(1+x*x) between -1/2 and 1/2. The method is simple: the integral is approximated by a sum of n intervals; the approximation to the integral in each interval is (1/n)*4/(1+x*x). The master process (rank 0) asks the user for the number of intervals; the master should then broadcast this number to all of the other processes. Each process then adds up every n'th interval (x = -1/2+rank/n, -1/2+rank/n+size/n,...). Finally, the sums computed by each process are added together using a reduction.
*/

#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) 
{ 
	int n, rank, size, i; 
	double PI25DT = 3.141592653589793238462643; 
	double mypi, pi, h, sum;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		cout << "Determining the value of pi with integral method" << endl << endl;
		cout << "Enter the number of intervals:" << endl; 
		cin >> n; 
	} 

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (n < 1)
		return EXIT_FAILURE;

	h = 1.0 / (double) n;
	sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (i = rank + 1; i <= n; i += size)
	{
		double x = h * ((double)i - 0.5);
		sum += (4.0 / (1.0 + x*x));
	} 
	mypi = h * sum;

	MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Finalize();

	if (rank == 0)
	{
		cout << "pi is approximately " << pi << ", Error is " << fabs(pi - PI25DT) << endl; 
#if _DEBUG
		system("pause");
#endif
	}
	return EXIT_SUCCESS; 
} 