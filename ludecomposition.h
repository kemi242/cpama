#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H

namespace CPAMA {

class Matrix;

/** LU Decomposition.
<P>
For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
unit lower triangular matrix L, an n-by-n upper triangular matrix U,
and a permutation vector piv of length m so that A(piv,:) = L*U.
If m < n, then L is m-by-m and U is m-by-n.
<P>
The LU decompostion with pivoting always exists, even if the matrix is
singular, so the constructor will never fail.  The primary use of the
LU decomposition is in the solution of square systems of simultaneous
linear equations.  This will fail if isNonsingular() returns false.
*/

class LUDecomposition {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @brief internal array storage.
    */
    double** LU;

    /** Row and column dimensions, and pivot sign.
    @brief column dimension.
    @brief row dimension.
    @brief pivot sign.
    */
    int m, n, pivsign;

    /** Internal storage of pivot vector.
    @brief pivot vector.
    */
    int* piv;

public:
/* ------------------------
   Constructor
 * ------------------------ */

    /** LU Decomposition
    @param  A   Rectangular matrix
    @return     Structure to access L, U and piv.
    */
    LUDecomposition(Matrix A);

/* ------------------------
   Destructor
 * ------------------------ */

    ~LUDecomposition();


/* ------------------------
   Public Methods
 * ------------------------ */

    /** Is the matrix nonsingular?
    @return     true if U, and hence A, is nonsingular.
    */
    bool isNonsingular ();

    /** Return lower triangular factor
    @return     L
    */
    Matrix getL ();

    /** Return upper triangular factor
    @return     U
    */
    Matrix getU ();

    /** Return pivot permutation vector
    @return     piv
    */
    int* getPivot ();

    /** Return pivot permutation vector as a one-dimensional double array
    @return     (double) piv
    */
    double* getDoublePivot ();

    /** Determinant
    @return     det(A)
    @exception  invalid_argument  Matrix must be square
    */
    double det ();

    /** Solve A*X = B
    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     X so that L*U*X = B(piv,:)
    @exception  invalid_argument Matrix row dimensions must agree.
    @exception  runtime_error  Matrix is singular.
    */
    Matrix solve (Matrix B);



};

}

#endif // LUDECOMPOSITION_H
