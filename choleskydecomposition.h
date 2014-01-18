#ifndef CHOLESKYDECOMPOSITION_H
#define CHOLESKYDECOMPOSITION_H

namespace CPAMA {

class Matrix;

/** Cholesky Decomposition.
<P>
For a symmetric, positive definite matrix A, the Cholesky decomposition
is an lower triangular matrix L so that A = L*L'.
<P>
If the matrix is not symmetric or positive definite, the constructor
returns a partial decomposition and sets an internal flag that may
be queried by the isSPD() method.
*/

class CholeskyDecomposition {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @brief internal array storage.
    */
    double** L;

    /** Row and column dimension (square matrix).
    @brief matrix dimension.
    */
    int n;

    /** Symmetric and positive definite flag.
    @brief is symmetric and positive definite flag.
    */
    bool isspd;

public:
/* ------------------------
   Constructor
 * ------------------------ */

    /** Cholesky algorithm for symmetric and positive definite matrix.
    @param  A   Square, symmetric matrix.
    @return     Structure to access L and isspd flag.
    */

    CholeskyDecomposition(Matrix Arg);

/* ------------------------
   Destructor
 * ------------------------ */

    ~CholeskyDecomposition();

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Is the matrix symmetric and positive definite?
    @return     true if A is symmetric and positive definite.
    */
    bool isSPD () const;

    /** Return triangular factor.
    @return     L
    */
    Matrix getL ();

    /** Solve A*X = B
    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     X so that L*L'*X = B
    @exception  invalid_argument  Matrix row dimensions must agree.
    @exception  runtime_error  Matrix is not symmetric positive definite.
    */
    Matrix solve (Matrix B);


};

}

#endif // CHOLESKYDECOMPOSITION_H
