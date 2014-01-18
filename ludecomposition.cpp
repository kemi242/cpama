#include "ludecomposition.h"

#include "matrix.h"

#include <cmath>
#include <stdexcept>

using namespace std;

using namespace CPAMA;

LUDecomposition::LUDecomposition(Matrix A) {
    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

       LU = A.getArrayCopy();
       m = A.getRowDimension();
       n = A.getColumnDimension();
       piv = new int[m];
       for (int i = 0; i < m; i++) {
          piv[i] = i;
       }
       pivsign = 1;
       double* LUrowi;
       double* LUcolj = new double[m];

       // Outer loop.

       for (int j = 0; j < n; j++) {

          // Make a copy of the j-th column to localize references.

          for (int i = 0; i < m; i++) {
             LUcolj[i] = LU[i][j];
          }

          // Apply previous transformations.

          for (int i = 0; i < m; i++) {
             LUrowi = LU[i];

             // Most of the time is spent in the following dot product.

             int kmax = fmin(i,j);
             double s = 0.0;
             for (int k = 0; k < kmax; k++) {
                s += LUrowi[k]*LUcolj[k];
             }

             LUrowi[j] = LUcolj[i] -= s;
          }

          // Find pivot and exchange if necessary.

          int p = j;
          for (int i = j+1; i < m; i++) {
             if (abs(LUcolj[i]) > abs(LUcolj[p])) {
                p = i;
             }
          }
          if (p != j) {
             for (int k = 0; k < n; k++) {
                double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
             }
             int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
             pivsign = -pivsign;
          }

          // Compute multipliers.

          if (j < m & LU[j][j] != 0.0) {
             for (int i = j+1; i < m; i++) {
                LU[i][j] /= LU[j][j];
             }
          }
       }
       delete[] LUcolj;
}

LUDecomposition::~LUDecomposition() {
    delete[] LU;
    delete[] piv;
}

bool LUDecomposition::isNonsingular () {
   for (int j = 0; j < n; j++) {
      if (LU[j][j] == 0)
         return false;
   }
   return true;
}

Matrix LUDecomposition::getL () {
   Matrix X = Matrix(m,n);
   double** L = X.getArray();
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         if (i > j) {
            L[i][j] = LU[i][j];
         } else if (i == j) {
            L[i][j] = 1.0;
         } else {
            L[i][j] = 0.0;
         }
      }
   }
   return X;
}

Matrix LUDecomposition::getU () {
   Matrix X = Matrix(n,n);
   double** U = X.getArray();
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
         if (i <= j) {
            U[i][j] = LU[i][j];
         } else {
            U[i][j] = 0.0;
         }
      }
   }
   return X;
}

int* LUDecomposition::getPivot () {
   int* p = new int[m];
   for (int i = 0; i < m; i++) {
      p[i] = piv[i];
   }
   return p;
}

double* LUDecomposition::getDoublePivot () {
   double* vals = new double[m];
   for (int i = 0; i < m; i++) {
      vals[i] = (double) piv[i];
   }
   return vals;
}

double LUDecomposition::det () {
   if (m != n) {
      throw invalid_argument("Matrix must be square.");
   }
   double d = (double) pivsign;
   for (int j = 0; j < n; j++) {
      d *= LU[j][j];
   }
   return d;
}

Matrix LUDecomposition::solve (Matrix B) {
   if (B.getRowDimension() != m) {
      throw invalid_argument("Matrix row dimensions must agree.");
   }
   if (!this->isNonsingular()) {
      throw new runtime_error("Matrix is singular.");
   }

   // Copy right hand side with pivoting
   int nx = B.getColumnDimension();
   Matrix Xmat = B.getMatrix(piv,m,0,nx-1);
   double** X = Xmat.getArray();

   // Solve L*Y = B(piv,:)
   for (int k = 0; k < n; k++) {
      for (int i = k+1; i < n; i++) {
         for (int j = 0; j < nx; j++) {
            X[i][j] -= X[k][j]*LU[i][k];
         }
      }
   }
   // Solve U*X = Y;
   for (int k = n-1; k >= 0; k--) {
      for (int j = 0; j < nx; j++) {
         X[k][j] /= LU[k][k];
      }
      for (int i = 0; i < k; i++) {
         for (int j = 0; j < nx; j++) {
            X[i][j] -= X[k][j]*LU[i][k];
         }
      }
   }
   return Xmat;
}
