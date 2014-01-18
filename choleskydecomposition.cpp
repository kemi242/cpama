#include "choleskydecomposition.h"

#include "matrix.h"

#include <cmath>
#include <stdexcept>

using namespace std;

using namespace CPAMA;

CholeskyDecomposition::CholeskyDecomposition(Matrix Arg) {
    // Initialize.
     double** A = Arg.getArray();
     n = Arg.getRowDimension();
     L = new double*[n];
     for (int i = 0; i < n; i++) {
         L[i] = new double[n];
     }
     isspd = (Arg.getColumnDimension() == n);
     // Main loop.
     for (int j = 0; j < n; j++) {
        double* Lrowj = L[j];
        double d = 0.0;
        for (int k = 0; k < j; k++) {
           double* Lrowk = L[k];
           double s = 0.0;
           for (int i = 0; i < k; i++) {
              s += Lrowk[i]*Lrowj[i];
           }
           Lrowj[k] = s = (A[j][k] - s)/L[k][k];
           d = d + s*s;
           isspd = isspd & (A[k][j] == A[j][k]);
        }
        d = A[j][j] - d;
        isspd = isspd & (d > 0.0);
        L[j][j] = sqrt(fmax(d,0.0));
        for (int k = j+1; k < n; k++) {
           L[j][k] = 0.0;
        }
     }
}

bool CholeskyDecomposition::isSPD () const {
   return isspd;
}

Matrix CholeskyDecomposition::getL () {
   return Matrix(L,n,n);
}

Matrix CholeskyDecomposition::solve (Matrix B) {
   if (B.getRowDimension() != n) {
      throw invalid_argument("Matrix row dimensions must agree.");
   }
   if (!isspd) {
      throw runtime_error("Matrix is not symmetric positive definite.");
   }

   // Copy right hand side.
   double** X = B.getArrayCopy();
   int nx = B.getColumnDimension();

       // Solve L*Y = B;
       for (int k = 0; k < n; k++) {
         for (int j = 0; j < nx; j++) {
            for (int i = 0; i < k ; i++) {
                X[k][j] -= X[i][j]*L[k][i];
            }
            X[k][j] /= L[k][k];
         }
       }

       // Solve L'*X = Y;
       for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            for (int i = k+1; i < n ; i++) {
                X[k][j] -= X[i][j]*L[i][k];
            }
            X[k][j] /= L[k][k];
         }
       }


   return Matrix(X,n,nx);
}
