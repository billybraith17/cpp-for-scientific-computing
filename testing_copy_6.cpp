
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

// Prototype of function f for  Δu = f
double f(double x, double y);
double u_exact(double x, double y);
double boundary(int i, int N_x, int N_y, double h_x, double h_y, double x_0, double y_0);


// Solves Δu(x,y) = -sin(x) - cos(y)
// with u(0,y) = cos(y), u(1,y) = sin(1) + cos(y), u(x,0) = sin(x) + 1, u(x,1) = sin(x) + cos(1)
//
// exact sol: u(x) = sin(x) + cos(y);



int main(int argc, char* argv[])
{

//  Boundaries
//	We use x_N and y_N to refer to the points on the far boundary but  really these are x_{N-1} and y_{N-1}
	double x_0=2*M_PI, x_N=4*M_PI;
	double y_0=2*M_PI, y_N=4*M_PI;



//	no. of iterations with different no. of gridpoints
	int n_its = 4;  // use 5 when getting plots for special topic

//	Lists to store things
	int* N_list;
	double* h_list;
	double* L2_list;

	N_list = new int [n_its];
	h_list = new double [n_its];
	L2_list = new double [n_its];



//	Looping over different numbers of grid points
	for(int j=1; j<n_its+1; j++)
	{
		// no. of grid points
		N_list[j-1]  = 1 + pow(2,1+j);
		int N_x = N_list[j-1];
		int N_y = N_list[j-1];



	//	step sizes
		double h_x = (x_N-x_0) / ((double) N_x-1 );
		double h_y = (y_N-y_0) / ((double) N_y-1 );



	//	Defining the matrix and vector st Ax = y
		Matrix A(N_x*N_y, N_x*N_y);
		Vector y(N_x*N_y);

	//	remember to ensure data at the boundaries is smooth
		for(int i=0; i<N_x*N_y; i++)
		{

	//		when x = x_0 (ignoring last entry as included later)
			if (i<N_y-1)
			{
				A[i][i] = 1;
				y[i] = boundary(i, N_x, N_y, h_x, h_y, x_0, y_0);
			}

	//		when x = x_N (ignoring first entry as included later)
			else if (i> N_x*N_y - N_y)
			{
				A[i][i] = 1;
				y[i] = boundary(i, N_x, N_y, h_x, h_y, x_0, y_0);
			}

	//		when y = y_0
			else if ((i+1)%N_y == 0)
			{
				A[i][i] = 1;
				y[i] = boundary(i, N_x, N_y, h_x, h_y, x_0, y_0);
			}

	//		when y = y_N
			else if (i%N_y == 0)
			{
				A[i][i] = 1;
				y[i] = boundary(i, N_x, N_y, h_x, h_y, x_0, y_0);
			}

	//		all other entries
			else
			{
				A[i][i] = -2*(1/pow(h_x,2) + 1/pow(h_y,2));
				A[i][i-1] = 1/ pow(h_y,2);
				A[i][i+1] = 1/ pow(h_y,2);
				A[i][i-N_y] = 1/ pow(h_x,2);
				A[i][i+N_y] = 1/ pow(h_x,2);
				y[i] = f(x_0 + floor(i/N_y) * h_x, y_0 + (i%N_y)*h_y);
			}
		}


	//	for testing matrices have been defined correctly (change boundary values)
//		std::cout << A << "\n";
//		std::cout << y << "\n \n";
//		std::cout.flush();


	//	Finding our solution
		Vector u_pred(N_x*N_y);
		try
		{
			u_pred = A/y;
		}
		catch (Exception& err)
		{
			err.DebugPrint();
		}
	//	std::cout << u_pred << "\n";


		// Checking A*u_pred = y
	//	It does therefore the solver is doing its job correctly
	//	Matrix u_pred_mat(N_x*N_y,1);
	//	for(int i=0; i<N_x*N_y; i++)
	//	{
	//		u_pred_mat[i][0] = u_pred[i];
	//	}
	//	std::cout << y << " " << A*u_pred_mat << "\n \n";


	//	Defining the exact solution and error at the gridpoints
		Vector exact_vec(N_x*N_y);
		Vector errors(N_x*N_y);
		for(int i=0; i<N_x*N_y; i++)
		{
			exact_vec[i] = u_exact(x_0 + floor(i/N_y) * h_x, y_0 + (i%N_y)*h_y);

	//		Checking the exact solution being calculated correctly
	//		std::cout << x_0 + floor(i/N) * h_x << " " << y_0 + (i%N)*h_y << "\n";
		}

	//	Calculating errors
		errors = exact_vec - u_pred;
		L2_list[j-1] = norm(errors) / ((double) N_x);
		std::cout << exact_vec << "\n";
		std::cout << u_pred << "\n";
//		std::cout << errors << "\n";
		std::cout << L2_list[j-1] << "\n";
	}

//	Write the list of N and L2 to a file so we can plot using matlab
	std::ofstream out("output.dat");
	assert(out.is_open());
	out.setf(std::ios::scientific);
	for(int k=0; k<n_its; k++)
	{
		out<< N_list[k] << " " << L2_list[k] <<"\n";
	}
	out.close();
}









// Function for f in  Δu = f
double f(double x, double y)
{
	return -sin(x) - cos(y);
}


// Exact solution of the problem
double u_exact(double x, double y)
{
	return sin(x) + cos(y);
}


// Gives boundary conditions given an index i
double boundary(int i, int N_x, int N_y, double h_x, double h_y, double x_0, double y_0)
{
	return u_exact(x_0 + floor(i/N_y) * h_x, y_0 + (i%N_y)*h_y);
}









