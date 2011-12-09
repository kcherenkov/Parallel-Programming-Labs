#Parallel Programming Labs: MPI and OpenMP Examples#
This labs will help you to understand C++ parallel programming with [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) and [OpenMP](https://en.wikipedia.org/wiki/OpenMP).
Visual Studio 2010 solution, [Microsoft MPI](https://en.wikipedia.org/wiki/Microsoft_Messaging_Passing_Interface), Intel Compiler with /Qopenmp.
[Compiled binaries and libraries](https://github.com/downloads/kcherenkov/Parallel-Programming-Labs/ParallelProgrammingLabs.zip)

##HelloWorld##
See Hello message from every thread and every process.

##SumNumbers##
This program sums all rows in an array using parallelism.
The root process acts as a master and sends a portion of the
array to each child process.  Master and child processes then
all calculate a partial sum of the portion of the array assigned
to them, and the child processes send their partial sums to
the master, who calculates a grand total.

##SumNumbersCascade##
Sums all rows in an array. Cascade algorithm.

##PiComputation##
This exercise presents a simple program to determine the value of pi.
The algorithm suggested here is chosen for its simplicity.
The method evaluates the integral of 4/(1+x*x) between 0 and 1.
The method is simple: the integral is approximated by a sum of n intervals;
the approximation to the integral in each interval is (1/n)*4/(1+x*x).
The master process (rank 0) asks the user for the number of intervals;
the master should then broadcast this number to all of the other processes.
Each process then adds up every n'th interval (x = rank/n, rank/n+size/n,...).
Finally, the sums computed by each process are added together using a reduction.

##MatrixVectorProduct##
MPI code for Matrix Vector Multiplication with Row wise block striped decomposition.
The main idea behind implementation of said phenomena can be understood by the following steps:

1. Each process contains a copy of the complete vector and one or more rows of the matrix.
2. Each process multiplies it’s row elements with the vector and saves it’s portion of solution vector.
3. Root process gathers and combines the portion of solution vector from every process and presents it as the output.

##MatrixMatrixProduct##
This example is a simple matrix multiplication program. AxB=C.

1. Matrix A is divided into blocks and distributed among processors.
2. Matrix B is copied to every processor.
3. The data is distributed among the workers who perform the actual multiplication in smaller blocks and send back their results to the master.