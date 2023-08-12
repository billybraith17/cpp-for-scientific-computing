#ifndef MATRIX_HPP_
#define MATRIX_HPP_



//  **********************
//  *  Class of matrices  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written


#include <cmath>
#include "Exception.hpp"//  This class throws errors using the class "error"
#include "Vector.hpp"
// Not including thr vector class in this so they can be run independently (and the makefile is less tricky)

class Matrix
{
private:
   // member variables
   double** mData;   // data stored in vector
   int* mSize;      // size of vector

public:
   // constructors
   // No default constructor
   // overridden copy constructor
   Matrix(const Matrix& m1);
   // construct vector of given size
   Matrix(int rows, int cols);

   // destructor
   ~Matrix();


   // All "friend" external operators and functions are declared as friend inside the class (here)
   // but their actual prototype definitions occur outside the class.
   // Binary operators
   friend Matrix operator+(const Matrix& v1, const Matrix& v2);
   friend Matrix operator-(const Matrix& v1, const Matrix& v2);
   friend Matrix operator*(const Matrix& v, const double& a);
   friend Matrix operator*(const double& a, const Matrix& v);
   friend Matrix operator*(const Matrix& v1, const Matrix& v2);
   friend int* size(const Matrix& v);
   friend Matrix operator/(const Matrix& v, const double& a);
   friend Vector operator/(const Matrix& A, const Vector& b);
   friend Vector operator/(const Matrix& A, const Matrix& b);

   // Unary operator
   friend Matrix operator-(const Matrix& v);

   //other operators
   //assignment
   Matrix& operator=(const Matrix& v);

   //indexing
   double* operator[](int i);

   // Output
   friend std::ostream& operator<<(std::ostream& output, const Matrix& m);

   // functions that are friends
   friend Matrix eye(int k);
   friend Matrix diag(Vector& v);
   friend Vector GMRES(const Matrix& A, const Vector& b, int k);
   friend Matrix transpose(Matrix& v);
   friend double norm(Matrix& v, int p);

};


// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Matrix operator+(const Matrix& v1, const Matrix& v2);
Matrix operator-(const Matrix& v1, const Matrix& v2);
Matrix operator*(const Matrix& v, const double& a);
Matrix operator*(const double& a, const Matrix& v);
Matrix operator*(const Matrix& v1, const Matrix& v2);
Matrix operator/(const Matrix& v, const double& a);
Vector operator/(const Matrix& A, const Vector& b);
Vector operator/(const Matrix& A, const Matrix& b);

// Unary operator
Matrix operator-(const Matrix& v);

// Transpose
Matrix transpose(Matrix& v);

// Prototype signature of size() friend function
int* size(const Matrix& v);

//prototype for eye, diag, GMRES and norm
Matrix eye(int k);
Matrix diag(Vector& v, int k=0);
Vector GMRES(const Matrix& A, const Vector& b, int k);
double norm(Matrix& v, int p=2);




#endif /* MATRIX_HPP_ */
