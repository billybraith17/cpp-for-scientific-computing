
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

// Prototype of function f for  Δu = f
double f(double x, double y, double z);
double u_exact(double x, double y, double z);
double boundary(int i, int N_x, int N_y, int N_z, double h_x, double h_y, double h_z, double x_0, double y_0, double z_0);
double xfromindex(int i, double x_0, int N_y, double N_z, double h_x);
double yfromindex(int i, double y_0, int N_y, int N_z, double h_y);
double zfromindex(int i, double z_0, int N_z, double h_z);


// Solves Δu(x,y) = -sin(x) - cos(y) - sin(z),     on [0,1]^3
//
// exact sol: u(x,y,z) = sin(x) + cos(y) + sin(z);


int main(int argc, char* argv[])
{

//  Boundaries
//	We use x_N and y_N to refer to the points on the far boundary but  really these are x_{N-1} and y_{N-1}
	double x_0=0, x_N=2*M_PI;
	double y_0=0, y_N=2*M_PI;
	double z_0=0, z_N=2*M_PI;

///	no. of iterations with different no. of gridpoints
	int n_its = 34;

//	Lists to store things
	int* k_list;
	double* L2_list;
	double* inf_list;

	k_list = new int [n_its];
	L2_list = new double [n_its];
	inf_list = new double [n_its];

	// no. of grid points
	int N_x = 17;
	int N_y = 17;
	int N_z = 17;

//	step sizes
	double h_x = (x_N-x_0) / ((double) N_x-1 );
	double h_y = (y_N-y_0) / ((double) N_y-1 );
	double h_z = (z_N-z_0) / ((double) N_z-1 );



//	Looping over different k
	for(int j=1; j<n_its+1; j++)
	{
		int k = 5*j;
		k_list[j-1] = k;

	//	Defining the matrix and vector st Ax = y
		Matrix A(N_x*N_y*N_z,N_x*N_y*N_z);
		Vector y_vec(N_x*N_y*N_z);

	//	remember to ensure data at the boundaries is smooth
		for(int i=0; i<N_x*N_y*N_z; i++)
		{
			double x = xfromindex(i, x_0, N_y, N_z, h_x);
			double y = yfromindex(i, y_0, N_y, N_z, h_y);
			double z = zfromindex(i, z_0, N_z, h_z);

			if ((x!=x_0 && x!=x_N) && (y!=y_0 && y!=y_N) && (z!=z_0 && z!=z_N))
			{
				A[i][i] = -2*(1/pow(h_x,2) + 1/pow(h_y,2) + 1/pow(h_z,2));
				A[i][i-1] = 1/ pow(h_z,2);
				A[i][i+1] = 1/ pow(h_z,2);
				A[i][i-N_y] = 1/ pow(h_y,2);
				A[i][i+N_y] = 1/ pow(h_y,2);
				A[i][i-N_y*N_z] = 1/ pow(h_x,2);
				A[i][i+N_y*N_z] = 1/ pow(h_x,2);
				y_vec[i] = f(x, y, z);
			}

			else
			{
				A[i][i] = 1;
				y_vec[i] = boundary(i, N_x, N_y, N_z, h_x, h_y, h_z, x_0, y_0, z_0);
			}
		}


	//	for testing matrices have been defined correctly (change boundary values)
	//	std::cout << A << "\n";
	//	std::cout << y_vec;


	//	Finding our solution
//		Vector u_pred(N_x*N_y*N_z);
		Vector u_GMRES(N_x*N_y*N_z);
		try
		{
//			u_pred = A/y_vec;
			u_GMRES = GMRES(A, y_vec, k);
		}
		catch (Exception& err)
		{
			err.DebugPrint();
		}
	//	std::cout << u_pred << "\n";


	//	Defining the exact solution and error at the gridpoints
		Vector exact_vec(N_x*N_y*N_z);
		Vector errors(N_x*N_y*N_z);
//		Vector errors2(N_x*N_y*N_z);
		for(int i=0; i<N_x*N_y*N_z; i++)
		{
			exact_vec[i] = u_exact(xfromindex(i, x_0, N_y, N_z, h_x), yfromindex(i, y_0, N_y, N_z, h_y), zfromindex(i, z_0, N_z, h_z));

	//		Checking the exact solution being calculated correctly
	//		std::cout << x_0 + floor(i/N) * h_x << " " << y_0 + (i%N)*h_y << "\n";
		}

	//	Calculating errors
//		errors = exact_vec - u_pred;
		errors = exact_vec - u_GMRES;
		L2_list[j-1] = norm(errors) / pow((double) N_x, 1.5);
		inf_list[j-1] = max(errors);
		std::cout << exact_vec << "\n";
//		std::cout << u_pred << "\n";
//		std::cout << u_GMRES << "\n";
		std::cout << L2_list[j-1] << "\n";
		std::cout << inf_list[j-1] << "\n";
	}
	

//	Write the list of N and L2 to a file so we can plot using matlab
	std::ofstream out("output2.dat");
	assert(out.is_open());
	out.setf(std::ios::scientific);
	for(int k=0; k<n_its; k++)
	{
		out<< k_list[k] << " " << L2_list[k] << " " << inf_list[k] <<"\n";
	}
	out.close();

}




// Function for f in  Δu = f
double f(double x, double y, double z)
{
	return -sin(x) - cos(y) - sin(z);
}

double u_exact(double x, double y, double z)
{
	return sin(x) + cos(y) + sin(z);
}


// Gives boundary conditions given an index i
double boundary(int i, int N_x, int N_y, int N_z, double h_x, double h_y, double h_z, double x_0, double y_0, double z_0)
{
	return u_exact(xfromindex(i, x_0, N_y, N_z, h_x), yfromindex(i, y_0, N_y, N_z, h_y), zfromindex(i, z_0, N_z, h_z));
}


//Gives the x position given an index i
double xfromindex(int i, double x_0, int N_y, double N_z, double h_x)
{
	return x_0 + floor(i/(N_y*N_z)) * h_x;
}

//Gives the y position given an index i
double yfromindex(int i, double y_0, int N_y, int N_z, double h_y)
{
	return y_0 + ((i/N_z)%( N_y))*h_y;
}

//Gives the z position given an index i
double zfromindex(int i, double z_0, int N_z, double h_z)
{
	return z_0 + (i%N_z)*h_z;
}
