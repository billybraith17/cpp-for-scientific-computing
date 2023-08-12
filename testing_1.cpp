
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

//	Testing the norm function for matrices
	int m=20;
	Matrix B(m,m);
//	norm(B);
//	std::cout << norm(B) << "\n";
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<m; j++)
		{
			B[i][i] = i+1;
		}
	}
//	B[1][1] = 4;
//	B[2][2] = 8;
//////	B = eye(m);
//	B[1][0] = 4;
//	B[2][1] = 7;
//	B[3][2] = 6;

	Vector a_vec(m);
//	for(int j=0; j<m; j++)
//	{
//		a_vec[j] = j+1;
//	}
	a_vec[0] = 1;
	
	Vector output_exact(m);
	Vector output_GMRES(m);
	Vector errors2(m);
	try
	{
		output_exact = B / a_vec;
//		output_GMRES = GMRES(B, a_vec, 2);
	}
		catch (Exception& err)
		{
			err.DebugPrint();
		}
	//	std::cout << u_pred << "\n";
//	output_exact = B / a_vec;
	output_GMRES = GMRES(B, a_vec, 19);
	errors2 = output_exact - output_GMRES;

	std::cout << "GMRES = \n" << output_GMRES << "\n \n";
	std::cout << "exact sol = \n" << output_exact << "\n \n";
	std::cout << "errors = \n" << output_exact - output_GMRES << "\n \n";
	std::cout << "MSE = \n" << norm(errors2)/m << "\n \n";
	std::cout << "supremum = \n" << max(errors2) << "\n \n";
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
