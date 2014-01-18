#include "matrix.h"

#include "singularvaluedecomposition.h"
#include "ludecomposition.h"
#include "qrdecomposition.h"
#include "choleskydecomposition.h"
#include "eigenvaluedecomposition.h"

#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace std;

using namespace CPAMA;

Matrix::Matrix (int m, int n) {
    this->m = m;
    this->n = n;
    A = new double*[m];
    for (int i = 0; i < m; i++) {
        A[i] = new double[n];
    }
}

Matrix::Matrix (int m, int n, double s) {
   this->m = m;
   this->n = n;
    A = new double*[m];
    for (int i = 0; i < m; i++) {
        A[i] = new double[n];
    }
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         A[i][j] = s;
      }
   }
}

Matrix::Matrix (double** A, int m, int n) {
    this->m = m;
    this->n = n;
    this->A = new double*[m];
    for (int i = 0; i < m; i++) {
        this->A[i] = new double[n];
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            this->A[i][j] = A[i][j];
        }
    }
}

Matrix::Matrix (double* vals, int valslen, int m) {
   this->m = m;
   n = (m != 0 ? valslen/m : 0);
   if (m*n != valslen) {
      throw invalid_argument("Array length must be a multiple of m.");
   }

   A = new double*[m];
   for (int i = 0; i < m; i++) {
       A[i] = new double[n];
       for (int j = 0; j < n; j++) {
          A[i][j] = vals[i+j*m];
       }
   }
}

Matrix::Matrix(const Matrix &other) {
    m = other.m;
    n = other.n;

    A = new double*[m];
    for (int i = 0; i < m; i++) {
        A[i] = new double[n];
        for (int j = 0; j < n; j++) {
            A[i][j] = other.A[i][j];
        }
    }


}

Matrix::~Matrix() {
    delete[] A;
}

Matrix Matrix::constructWithCopy(double** A, int m, int n) {
   Matrix X = Matrix(m,n);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[i][j] = A[i][j];
      }
   }
   return X;
}

Matrix Matrix::copy () {
   Matrix X = Matrix(m,n);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[i][j] = A[i][j];
      }
   }
   return X;
}

double** Matrix::getArray () {
   return A;
}

double** Matrix::getArrayCopy () {
   double** C = new double*[m];
   for (int i = 0; i < m; i++) {
       C[i] = new double[n];
       for (int j = 0; j < n; j++) {
           C[i][j] = A[i][j];
       }
   }

   return C;
}

double* Matrix::getColumnPackedCopy () {
   double* vals = new double[m*n];
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         vals[i+j*m] = A[i][j];
      }
   }
   return vals;
}

double* Matrix::getRowPackedCopy () {
   double* vals = new double[m*n];
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         vals[i*n+j] = A[i][j];
      }
   }
   return vals;
}

int Matrix::getRowDimension () const {
   return m;
}

int Matrix::getColumnDimension () const {
   return n;
}

double Matrix::get (int i, int j) {
   return A[i][j];
}

Matrix Matrix::getMatrix (int i0, int i1, int j0, int j1) {
   Matrix X = Matrix(i1-i0+1,j1-j0+1);
   double** B = X.getArray();
   for (int i = i0; i <= i1; i++) {
      for (int j = j0; j <= j1; j++) {
         B[i-i0][j-j0] = A[i][j];
      }
   }

   return X;
}

Matrix Matrix::getMatrix (int* r, int rlen, int* c, int clen) {
   Matrix X = Matrix(rlen,clen);
   double** B = X.getArray();
   for (int i = 0; i < rlen; i++) {
      for (int j = 0; j < clen; j++) {
         B[i][j] = A[r[i]][c[j]];
      }
   }

   return X;
}

Matrix Matrix::getMatrix (int i0, int i1, int* c, int clen) {
   Matrix X = Matrix(i1-i0+1,clen);
   double** B = X.getArray();
   for (int i = i0; i <= i1; i++) {
      for (int j = 0; j < clen; j++) {
         B[i-i0][j] = A[i][c[j]];
      }
   }

   return X;
}

Matrix Matrix::getMatrix (int* r, int rlen, int j0, int j1) {
   Matrix X = Matrix(rlen,j1-j0+1);
   double** B = X.getArray();
   for (int i = 0; i < rlen; i++) {
      for (int j = j0; j <= j1; j++) {
         B[i][j-j0] = A[r[i]][j];
      }
   }

   return X;
}

void Matrix::set (int i, int j, double s) {
   A[i][j] = s;
}

void Matrix::setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
   for (int i = i0; i <= i1; i++) {
      for (int j = j0; j <= j1; j++) {
         A[i][j] = X.get(i-i0,j-j0);
      }
   }
}

void Matrix::setMatrix (int* r, int rlen, int* c, int clen, Matrix X) {
    for (int i = 0; i < rlen; i++) {
        for (int j = 0; j < clen; j++) {
            A[r[i]][c[j]] = X.get(i,j);
        }
    }
}

void Matrix::setMatrix (int* r, int rlen, int j0, int j1, Matrix X) {
   for (int i = 0; i < rlen; i++) {
      for (int j = j0; j <= j1; j++) {
         A[r[i]][j] = X.get(i,j-j0);
      }
   }
}

void Matrix::setMatrix (int i0, int i1, int* c, int clen, Matrix X) {
   for (int i = i0; i <= i1; i++) {
      for (int j = 0; j < clen; j++) {
         A[i][c[j]] = X.get(i-i0,j);
      }
   }
}

Matrix Matrix::transpose () {
   Matrix X = Matrix(n,m);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[j][i] = A[i][j];
      }
   }
   return X;
}

double Matrix::norm1 () {
   double f = 0;
   for (int j = 0; j < n; j++) {
      double s = 0;
      for (int i = 0; i < m; i++) {
         s += abs(A[i][j]);
      }
      f = max(f,s);
   }
   return f;
}

double Matrix::norm2 () {
   return SingularValueDecomposition(*this).norm2();
}

double Matrix::normInf () {
   double f = 0;
   for (int i = 0; i < m; i++) {
      double s = 0;
      for (int j = 0; j < n; j++) {
         s += abs(A[i][j]);
      }
      f = max(f,s);
   }
   return f;
}

double Matrix::normF () {
   double f = 0;
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         f = hypot(f,A[i][j]);
      }
   }
   return f;
}

Matrix Matrix::operator -() {
    Matrix X = Matrix(m,n);
    double** C = X.getArray();
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          C[i][j] = -A[i][j];
       }
    }
    return X;
}

Matrix Matrix::operator +(Matrix B) {
    checkMatrixDimensions(B);
    Matrix X = Matrix(m,n);
    double** C = X.getArray();
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          C[i][j] = A[i][j] + B.A[i][j];
       }
    }
    return X;
}

Matrix Matrix::operator +=(Matrix B) {
    checkMatrixDimensions(B);
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          A[i][j] = A[i][j] + B.A[i][j];
       }
    }
    return *this;
}

Matrix Matrix::operator -(Matrix B) {
    checkMatrixDimensions(B);
    Matrix X = Matrix(m,n);
    double** C = X.getArray();
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          C[i][j] = A[i][j] - B.A[i][j];
       }
    }
    return X;
}

Matrix Matrix::operator -=(Matrix B) {
    checkMatrixDimensions(B);
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          A[i][j] = A[i][j] - B.A[i][j];
       }
    }
    return *this;
}

Matrix Matrix::arrayTimes (Matrix B) {
   checkMatrixDimensions(B);
   Matrix X = Matrix(m,n);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[i][j] = A[i][j] * B.A[i][j];
      }
   }
   return X;
}

Matrix Matrix::arrayTimesEquals (Matrix B) {
   checkMatrixDimensions(B);
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         A[i][j] = A[i][j] * B.A[i][j];
      }
   }
   return *this;
}

Matrix Matrix::arrayRightDivide (Matrix B) {
   checkMatrixDimensions(B);
   Matrix X = Matrix(m,n);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[i][j] = A[i][j] / B.A[i][j];
      }
   }
   return X;
}

Matrix Matrix::arrayRightDivideEquals (Matrix B) {
   checkMatrixDimensions(B);
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         A[i][j] = A[i][j] / B.A[i][j];
      }
   }
   return *this;
}

Matrix Matrix::arrayLeftDivide (Matrix B) {
   checkMatrixDimensions(B);
   Matrix X = Matrix(m,n);
   double** C = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         C[i][j] = B.A[i][j] / A[i][j];
      }
   }
   return X;
}

Matrix Matrix::arrayLeftDivideEquals (Matrix B) {
   checkMatrixDimensions(B);
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         A[i][j] = B.A[i][j] / A[i][j];
      }
   }
   return *this;
}

Matrix Matrix::operator *(double s) {
    Matrix X = Matrix(m,n);
    double** C = X.getArray();
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          C[i][j] = s*A[i][j];
       }
    }
    return X;
}

Matrix Matrix::operator *=(double s) {
    for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++) {
          A[i][j] = s*A[i][j];
       }
    }
    return *this;
}

Matrix Matrix::operator *(Matrix B) {
    if (B.m != n) {
       throw invalid_argument("Matrix inner dimensions must agree.");
    }
    Matrix X = Matrix(m,B.n);
    double** C = X.getArray();
    double* Bcolj = new double[n];
    for (int j = 0; j < B.n; j++) {
       for (int k = 0; k < n; k++) {
          Bcolj[k] = B.A[k][j];
       }
       for (int i = 0; i < m; i++) {
          double* Arowi = A[i];
          double s = 0;
          for (int k = 0; k < n; k++) {
             s += Arowi[k]*Bcolj[k];
          }
          C[i][j] = s;
       }
    }

    delete[] Bcolj;
    return X;
}

LUDecomposition Matrix::lu () {
   return LUDecomposition(*this);
}

QRDecomposition Matrix::qr () {
   return QRDecomposition(*this);
}

CholeskyDecomposition Matrix::chol () {
   return CholeskyDecomposition(*this);
}

SingularValueDecomposition Matrix::svd () {
   return SingularValueDecomposition(*this);
}

EigenvalueDecomposition Matrix::eig () {
   return EigenvalueDecomposition(*this);
}

Matrix Matrix::solve (Matrix B) {
   return (m == n ? (LUDecomposition(*this)).solve(B) :
                    (QRDecomposition(*this)).solve(B));
}

Matrix Matrix::solveTranspose (Matrix B) {
   return transpose().solve(B.transpose());
}

Matrix Matrix::inverse () {
   return solve(identity(m,m));
}

double Matrix::det () {
   return LUDecomposition(*this).det();
}

int Matrix::rank () {
   return SingularValueDecomposition(*this).rank();
}

double Matrix::cond () {
   return SingularValueDecomposition(*this).cond();
}

double Matrix::trace () {
   double t = 0;
   for (int i = 0; i < min(m,n); i++) {
      t += A[i][i];
   }
   return t;
}

Matrix Matrix::random (int m, int n) {
   Matrix A = Matrix(m,n);
   double** X = A.getArray();
   srand(time(NULL));
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
          X[i][j] = (double)rand() / RAND_MAX;
      }
   }
   return A;
}

Matrix Matrix::identity (int m, int n) {
   Matrix A = Matrix(m,n);
   double** X = A.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         X[i][j] = (i == j ? 1.0 : 0.0);
      }
   }
   return A;
}

void Matrix::print (int w, int d) {
   cout << endl;  // start on new line.
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
          cout << ' ' << setw(w) << fixed << setprecision(d) << A[i][j];
      }
      cout << endl;
   }
   cout << endl;   // end with blank line.
}

void Matrix::checkMatrixDimensions (Matrix B) {
   if (B.m != m || B.n != n) {
      throw  invalid_argument("Matrix dimensions must agree.");
   }
}
