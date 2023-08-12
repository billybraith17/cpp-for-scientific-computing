
#include <stdlib.h>
#include <iostream>
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


// Solves Δu(x,y) = -sin(x) - cos(y) - sin(z),     on [2pi,4pi]^3
//
// exact sol: u(x,y,z) = sin(x) + cos(y) + sin(z);


int main(int argc, char* argv[])
{

//  Boundaries
//	We use x_N and y_N to refer to the points on the far boundary but  really these are x_{N-1} and y_{N-1}
	double x_0=2*M_PI, x_N=4*M_PI;
	double y_0=2*M_PI, y_N=4*M_PI;
	double z_0=2*M_PI, z_N=4*M_PI;

//	no. of grid points per row
	int N_x = 17;
	int N_y = 17;
	int N_z = 17;

//	step sizes
	double h_x = (x_N-x_0) / ((double) N_x-1 );
	double h_y = (y_N-y_0) / ((double) N_y-1 );
	double h_z = (z_N-z_0) / ((double) N_z-1 );


//	checking the index to real number function works correctly
//	Tested with  N_x = 2, N_y = 4, N_z = 5; so we have different no. of points for each
//	for(int i=0; i<N_x*N_y*N_z; i++)
//	{
//		std::cout << xfromindex(i, x_0, N_y, N_z, h_x) << " " << yfromindex(i, y_0, N_y, N_z, h_y) << " " << zfromindex(i, z_0, N_z, h_z) << " " << "\n";
//	}



//	Defining the matrix and vector st A*p = y_vec
	Matrix A(N_x*N_y*N_z,N_x*N_y*N_z);
	Vector y_vec(N_x*N_y*N_z);

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


//	for testing matrices have been defined correctly
//	std::cout << A << "\n";
//	std::cout << y_vec;


//	Finding our solution
//	Vector u_pred(N_x*N_y*N_z);
	Vector u_GMRES(N_x*N_y*N_z);
	try
	{
//		u_pred = A/y_vec;
		u_GMRES = GMRES(A, y_vec, 120);
	}
	catch (Exception& err)
	{
		err.DebugPrint();
	}
//	std::cout << u_pred << "\n";


//	Defining the exact solution and error at the gridpoints
	Vector exact_vec(N_x*N_y*N_z);
	Vector errors(N_x*N_y*N_z);
	for(int i=0; i<N_x*N_y*N_z; i++)
	{
		exact_vec[i] = u_exact(xfromindex(i, x_0, N_y, N_z, h_x), yfromindex(i, y_0, N_y, N_z, h_y), zfromindex(i, z_0, N_z, h_z));

//		Checking the exact solution being calculated correctly
//		std::cout << x_0 + floor(i/N) * h_x << " " << y_0 + (i%N)*h_y << "\n";
	}

//	Calculating errors
	errors = exact_vec - u_GMRES;
	double inf = max(errors);
	double L2 = norm(errors) / pow((double) N_x, 1.5);
	std::cout << "exact_vec = " << exact_vec << "\n";
//	std::cout << u_pred << "\n";
	std::cout << "GMRES = " << u_GMRES << "\n";
	std::cout << "errors = " << errors << "\n";
	std::cout << "L2 = " << L2 << "\n";
	std::cout << "infinity = " << inf << "\n";
	


//	std::cout << GMRES(A, y_vec, 2) << "\n";



//	Testing the norm function for matrices
	int m=20;
	Matrix B(m,m);
//	norm(B);
//	std::cout << norm(B) << "\n";
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<m; j++)
		{
			B[i][i] = i+j;
		}
	}
	B[1][1] = 4;
	B[2][2] = 8;
//	B = eye(m);
	B[1][0] = 4;
	B[2][1] = 7;
//	B[3][2] = 6;

	Vector a_vec(m);
	for(int j=0; j<m; j++)
	{
		a_vec[j] = j;
	}
	
	Vector output_exact(m);
	Vector output_GMRES(m);
	Vector errors2(m);
	output_exact = B / a_vec;
	output_GMRES = GMRES(B, a_vec, m-10);
	errors2 = output_exact - output_GMRES;

	std::cout << "GMRES = \n" << output_GMRES << "\n \n";
	std::cout << "exact sol = \n" << output_exact << "\n \n";
	std::cout << "errors = \n" << output_exact - output_GMRES << "\n \n";
	std::cout << "MSE = \n" << norm(errors2)/m << "\n \n";
//	std::cout << norm(B) << "\n";


////	B = B/3;
//	B = eye(4);
//	std::cout << B << "\n";

//	Matrix Q_T(4,4);
//	Q_T = eye(4);
//	std::cout << Q_T << "\n";
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
