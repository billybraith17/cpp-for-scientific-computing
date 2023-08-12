#include <iostream>
#include "Matrix.hpp"




// constructor that creates matrix of given size with
// double precision entries all initially set to zero
Matrix::Matrix(int rows, int cols)
{
	mData = new double* [rows];

	for(int i = 0; i<rows; i++)
	{
		mData[i] = new double [cols];
	}

	// defining matrix
	for(int i = 0; i<rows; i++)
	{
		for(int j = 0; j<cols; j++)
		{
			mData[i][j] = 0;
		}
	}

	mSize = new int[2];
	mSize[0] = rows;
	mSize[1] = cols;
}




// copy constructor - creates matrix with the same entries as m1
Matrix::Matrix(const Matrix& m1)
{
  mSize = m1.mSize;
  mData = new double* [mSize[0]];

  for(int i = 0; i<mSize[0]; i++)
  {
	  mData[i] = new double [mSize[1]];
  }

  for(int i = 0; i<mSize[0]; i++)
  {
	  for(int j = 0; j<mSize[1]; j++)
	  {
	      mData[i][j] = m1.mData[i][j];
	  }
  }


}




// destructor - deletes pointer
Matrix::~Matrix()
{
  for (int i=0; i<mSize[0]; i++)
  {
	  delete[] mData[i];
  }
  delete[] mData;
}







// definition of << for matrices
std::ostream& operator<<(std::ostream& output, const Matrix& m) {
  for (int i=0; i<m.mSize[0]; i++)
  {
	  output << "[";
	  for (int j=0; j<m.mSize[1]; j++)
	  {
	      output <<  m.mData[i][j];
	      if (j != m.mSize[1]-1)
	    	  output  << ", ";
	      else
	    	  output  << "]" << "\n";
	  }
    }
  return output;  // for multiple << operators.
}










// definition of + between two matrices
Matrix operator+(const Matrix& v1,
                        const Matrix& v2)
{
	assert ((v1.mSize[0] == v2.mSize[0]) && (v1.mSize[1] == v2.mSize[1]));

	Matrix w(v1.mSize[0], v1.mSize[1]);

    for (int i=0; i<v1.mSize[0]; i++)
	{
    	for (int j=0; j<v1.mSize[1]; j++)
    	{
    		w.mData[i][j] = v1.mData[i][j] + v2.mData[i][j];
    	}
	}

  return w;
}






// definition of - between two matrices
Matrix operator-(const Matrix& v1,
                        const Matrix& v2)
{
	assert ((v1.mSize[0] == v2.mSize[0]) && (v1.mSize[1] == v2.mSize[1]));

	Matrix w(v1.mSize[0], v1.mSize[1]);

    for (int i=0; i<v1.mSize[0]; i++)
	{
    	for (int j=0; j<v1.mSize[1]; j++)
    	{
    		w.mData[i][j] = v1.mData[i][j] - v2.mData[i][j];
    	}
	}

  return w;
}













// definition of matrix operator =
Matrix& Matrix::operator=(const Matrix& v)
{

//  check both matrices have same size
	assert ((v.mSize[0] == mSize[0]) && (v.mSize[1] == mSize[1]));

    for (int i=0; i<mSize[0]; i++)
	{
    	for (int j=0; j<mSize[1]; j++)
    	{
    		mData[i][j] = v.mData[i][j];
    	}
	}

  return *this;
}











// definition of matrix operator []
// allows v.mData[i][j] to be written as v[i][j]
// v[i] returns a pointer to the first element of row i
double* Matrix::operator[](int i)
{

  if (i < 0)
    {
      throw Exception("out of range",
		  "accessing matrix through () - index too small");
    }
  else if (i > mSize[0]-1)
    {
      throw Exception("length mismatch",
		  "accessing matrix through () - index too high");
    }


  return mData[i];

}

















// definition of the unary operator -
Matrix operator-(const Matrix& v)
{

//  create a matrix w with entries equal to -v

  Matrix w(v.mSize[0], v.mSize[1]);

  for (int i=0; i<v.mSize[0]; i++)
  {
	  for (int j=0; j<v.mSize[1]; j++)
	  {
		  w.mData[i][j] = -v.mData[i][j];
	  }
  }

  return w;
}












// definition of multiplication between a matrix and a scalar
Matrix operator*(const Matrix& v, const double& a)
{

//  create a matrix of the same length as v with entries equal to a*v

  Matrix w(v.mSize[0], v.mSize[1]);

  for (int i=0; i<v.mSize[0]; i++)
  {
	  for (int j=0; j<v.mSize[1]; j++)
	  {
		  w.mData[i][j] = a * v.mData[i][j];
	  }
  }

  return w;
}




// definition of multiplication between a scalar and a matrix
Matrix operator*(const double& a, const Matrix& v)
{

//  create a matrix of the same length as v with entries equal to a*v

  Matrix w(v.mSize[0], v.mSize[1]);

  for (int i=0; i<v.mSize[0]; i++)
  {
	  for (int j=0; j<v.mSize[1]; j++)
	  {
		  w.mData[i][j] = a * v.mData[i][j];
	  }
  }

  return w;
}











// definition of scalar product between two matrices
Matrix operator*(const Matrix& v1, const Matrix& v2)
{

  assert ((v1.mSize[1] == v2.mSize[0]));

  Matrix w(v1.mSize[0], v2.mSize[1]);


  for(int i = 0; i<v1.mSize[0]; i++)
	{
		for(int j = 0; j<v2.mSize[1]; j++)
		{
			w.mData[i][j] = 0;
			for(int k=0; k< v1.mSize[1]; k++)
			{
				w.mData[i][j] += v1.mData[i][k]*v2.mData[k][j];
			}
		}
	}

  return w;
}










// definition of / between a matrix and a scalar
Matrix operator/(const Matrix& v, const double& a)
{
	  if (a == 0.0)
	  {
	     throw Exception("div 0", "Attempt to divide by zero");
	  }
	//  create a matrix of the same length as v with entries equal to v/a

	  Matrix w(v.mSize[0], v.mSize[1]);

	  for (int i=0; i<v.mSize[0]; i++)
	  {
		  for (int j=0; j<v.mSize[1]; j++)
			{
			  w.mData[i][j] = v.mData[i][j] / a;
			}
	  }

	  return w;
}









// Solving a linear system using Gaussian elimination (for matrix and vector)
Vector operator/(const Matrix& A, const Vector& b)
{
	// check matrix and vector of correct dimensions
	assert ((size(A)[0] == length(b)) && (size(A)[1] == length(b)));

	int i=0,j=0;
	int n = size(A)[0];

	//	Create matrix and vector to work with
	Matrix A_copy(A);
	Vector b_copy(b);

	// Use j to decide when we finish. If non-singular should be when j=n
	while (j<n)
	{
		// make A[i][j] nonzero
		int homerow = i;
		while (A_copy[i][j] == 0)
		{
			i += 1;
		}

		// swapping row i with row pivot row
		double* swap;
		swap = new double [n];
		for(int k=j; k<n; k++)
		{
			swap[k] = A_copy[i][k];
			A_copy[i][k] = A_copy[homerow][k];
			A_copy[homerow][k] = swap[k];
		}
		swap[0] = b_copy[i];
		b_copy[i] = b_copy[homerow];
		b_copy[homerow] = swap[0];


		// set i back to the top row
		i=homerow;

		// set first entry of row to 1
		if (A_copy[i][j] != 1)
		{
			// set first1 to divide through the other elements as A[i][j] would change
			double first1 = A_copy[i][j];
			for(int k=0; k<n; k++)
			{
				A_copy[i][k] /= first1;
			}
			b_copy[i] /= first1;
		}


		// subtract top row from lower rows
		for(int k=i+1; k<n; k++)
		{
			double first2 = A_copy[k][j];
			for(int l=0; l<n; l++)
			{
				A_copy[k][l] -= A_copy[i][l]*first2;
			}
			b_copy[k] -= b_copy[i]*first2;
		}

		// subtract top row from higher rows
		for(int k=0; k<i; k++)
		{
			double first = A_copy[k][j];
			for(int l=0; l<n; l++)
			{
				A_copy[k][l] -= A_copy[i][l]*first;
			}
			b_copy[k] -= b_copy[i]*first;
		}

		// move pivot diagonally
		i += 1;
		j += 1;
	}
	return b_copy;
}












// Solving a linear system using Gaussian elimination (for matrix and vector)
Vector operator/(const Matrix& A, const Matrix& b)
{
	// check vectors of correct dimensions
	assert (((size(A)[1] == size(b)[0]) && (size(A)[0] == size(b)[0])) && (size(b)[1] == 1));

	int i=0,j=0;
	int n = size(A)[0];

	//	Create matrix and vector to work with
	Matrix A_copy(A);
	Matrix b_mat(b);
	Vector b_copy(n);

	for(int i=0; i<n; i++)
	{
		b_copy[i] = b_mat[i][0];
	}

	// Use j to decide when we finish. As non-singular should be when j=n
	while (j<n)
	{
		// make A[i][j] nonzero
		int homerow = i;
		while (A_copy[i][j] == 0)
		{
			i += 1;
		}

		// swapping row i with row pivot row
		double* swap;
		swap = new double [n];
		for(int k=j; k<n; k++)
		{
			swap[k] = A_copy[i][k];
			A_copy[i][k] = A_copy[homerow][k];
			A_copy[homerow][k] = swap[k];
		}
		swap[0] = b_copy[i];
		b_copy[i] = b_copy[homerow];
		b_copy[homerow] = swap[0];


		// set i back to the top row
		i=homerow;

		// set first entry of row to 1
		if (A_copy[i][j] != 1)
		{
			// set first1 to divide through the other elements as A[i][j] would change
			double first1 = A_copy[i][j];
			for(int k=0; k<n; k++)
			{
				A_copy[i][k] /= first1;
			}
			b_copy[i] /= first1;
		}


		// subtract top row from lower rows
		for(int k=i+1; k<n; k++)
		{
			double first2 = A_copy[k][j];
			for(int l=0; l<n; l++)
			{
				A_copy[k][l] -= A_copy[i][l]*first2;
			}
			b_copy[k] -= b_copy[i]*first2;
		}

		// subtract top row from higher rows
		for(int k=0; k<i; k++)
		{
			double first = A_copy[k][j];
			for(int l=0; l<n; l++)
			{
				A_copy[k][l] -= A_copy[i][l]*first;
			}
			b_copy[k] -= b_copy[i]*first;
		}

		// move pivot diagonally
		i += 1;
		j += 1;
	}
	return b_copy;
}

















// Solving a linear system using GMRES
Vector GMRES(const Matrix& A, const Vector& b, int k)
{
	// check vectors of correct dimensions
	assert ((size(A)[0] == length(b)) && (size(A)[1] == length(b)));

//	Check k is viable
	assert (k<=size(A)[0]);

//	Make copy of A and b so can keep type const
	Matrix A_copy(A);
	Vector b_copy(b);

//	Setting the size of the vectors to n
	int n = size(A)[1];

////////////////////////////////////////
//////// Arnoldi iteration /////////////
////////////////////////////////////////

//	creating the matrix Q_{k+1} but we just call it Q_k for brevity
	Matrix Q_k(n, k+1);


//	Defining the vector we will use for q_i (as a matrix so can do operations more easily)
	Matrix q(n,1);
	for(int i = 0; i<n; i++)
	{
		q[i][0] = b_copy[i];
	}
	q = q/norm(q);

//	Setting first column of Q_k to q
	for(int l=0; l<n; l++)
	{
		Q_k[l][0] = q[l][0];
	}


//	creating the matrix H which is (k+1) x k
	Matrix H(k+1, k);


//	Calculating the Arnoldi iteration k times
//	Need to check how many times to run this through and
	for(int i=0; i<k; i++)
	{
		Matrix v(n,1);
		v = A*q;


		for(int j=0; j<i+1; j++)		// always will do at least one iteration so the indexing is correct
		{
//			calculating dot product of Q_j and v. Can do like this as all matrix data initialised to 0
			for(int l=0; l<n; l++)
			{
				H[j][i] += Q_k[l][j]*v[l][0];
			}

//			Calculating the new vector v
			for(int l=0; l<n; l++)
			{
				v[l][0] = v[l][0] - H[j][i]*Q_k[l][j];
			}
		}
		H[i+1][i] = norm(v);

//		setting new column of Q_k
//		if  (i<n-1)
//		{
			for(int l=0; l<n; l++)
			{
				Q_k[l][i+1] = v[l][0] / H[i+1][i];
				q[l][0] = Q_k[l][i+1];
			}
//		}
	}
//	checking the output of the arnoldi iteration
	std::cout << "Q_k = \n" << Q_k << "\n";
	std::cout << "H = \n" << H << "\n";

//	checking the columns of Q_k are orthogonal
	std::cout << "Q_k^T * Q_k = \n" << transpose(Q_k)*Q_k << "\n";


/////////////////////////////////////////////
//////// QR by Givens rotations of H ////////
/////////////////////////////////////////////

//	define matrices for R and Q_T for the QR factoriation of H
//	Q_T stands for Q^T with Q from the QR factorisation
	Matrix R(H);
	Matrix Q_T(k+1,k+1);
	Q_T = eye(k+1);

//	We need to apply k Givens rotations
	for(int i=0; i<k; i++)
	{
//		Defining each Givens rotation
		Matrix G(k+1,k+1);
		for(int j=0; j<k+1; j++)
		{
			G[j][j] = 1;
		}
//		Maybe have to change this part if H[i][i] = 0
		double c = R[i][i] / pow( pow(R[i][i],2) + pow(R[i+1][i],2), 0.5);
		double s = R[i+1][i] / pow( pow(R[i][i],2) + pow(R[i+1][i],2), 0.5);
		G[i][i] = c;
		G[i+1][i] = -s;
		G[i][i+1] = s;
		G[i+1][i+1] = c;

		R = G*R;
		Q_T = G*Q_T;

//		testing the Givens rotations work
//		std::cout<< "R = " << "\n" << R << "\n";
//		std::cout.flush();
	}


/////////////////////////////////////////////
/////// Solving Ry = Q^T Q_{k+1}^T b ////////
/////////////////////////////////////////////

//	Write b as a matrix so we can multipy by it
	Matrix b_mat(n,1);
	for(int i = 0; i<n; i++)
	{
		b_mat[i][0] = b_copy[i];
	}

//	To solve the Hessenberg least squares problem we solve Ry = Q_T Q_{k+1}^T b
//	but ignoring the bottom row of the matrices on either side
//	where x = Q_k*y
	Vector y(k);

//	Creating Vector and matrix without the bottom row
	Matrix R2(k,k);
	Matrix RHS(k+1,1);
	Matrix RHS2(k,1);

	RHS = Q_T*(transpose(Q_k)*b_mat);

	for(int i=0; i<k; i++)
	{
		for(int j=0; j<k; j++)
		{
			R2[i][j] = R[i][j];
		}
		RHS2[i][0] = RHS[i][0];
	}

	std::cout<< "RHS = " << "\n" << RHS2 << "\n";
	std::cout.flush();


	y = R2 / RHS2;

//	Writing y as a matrix so we can multipy by it
	Matrix y_mat(k,1);
	for(int i = 0; i<k; i++)
	{
		y_mat[i][0] = y[i];
	}


// 	creating Q_k so we can find x from y
	Matrix Q_k2(n,k);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<k; j++)
		{
			Q_k2[i][j] = Q_k[i][j];
		}
	}

////	checking the gaussian elimination solver is working correctly
//	std::cout<< R2*y_mat << "\n";
//	std::cout<< RHS2 << "\n";
//	std::cout.flush();


	Matrix output_mat(n,1);
	Vector output(n);

//	calculating x = Q_k*y
//	changing from y back to x
	output_mat = Q_k2*y_mat;

//	converting to a Vector
	for(int i=0; i<n; i++)
	{
			output[i] = output_mat[i][0];
	}

	return output;
}




















// definition of transposing a matrix
Matrix transpose(Matrix& v)
{
	  Matrix w(v.mSize[1], v.mSize[0]);

	  for(int i = 0; i<v.mSize[1]; i++)
		{
			for(int j = 0; j<v.mSize[0]; j++)
			{
				w.mData[i][j] = v.mData[j][i];
			}
		}

	return w;
}









// return size of a matrix (as an array)
int* size(const Matrix& v)
{
	int* pointer = new int[2];
	pointer[0] = v.mSize[0];
	pointer[1] = v.mSize[1];

  return pointer;
}






// Returns the k dimensional identity matrix
Matrix eye(int k)
{
	Matrix w(k, k);

	for(int i = 0; i<k; i++)
	{
		w.mData[i][i] = 1;
	}

	return w;
}







// Definition of the diag function from Matlab
Matrix diag(Vector& v, int k)
{
	Matrix w(length(v)+abs(k), length(v)+abs(k));

	if (k>=0)
	{
		for(int i = 0; i<length(v); i++)
		{
			w[i][i+k] = v[i];
		}
	}
	else
	{
		for(int i = 0; i<length(v); i++)
				{
					w[i+abs(k)][i] = v[i];
				}
	}

	return w;
}
















// calculate p-norm of a matrix v
// default value for p is 2
double norm(Matrix& v, int p)
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=0; i<size(v)[0]; i++)
  {
	  for (int j=0; j<size(v)[1]; j++)
		{
		  temp = fabs(v[i][j]);
		  norm_val += pow(temp, p);
		}
  }

  return pow(norm_val, 1.0/((double) (p)));
}








