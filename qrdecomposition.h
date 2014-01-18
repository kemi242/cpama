#ifndef QRDECOMPOSITION_H
#define QRDECOMPOSITION_H

namespace  CPAMA {

/** QR Decomposition.
<P>
   For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns false.
*/

class Matrix;

class QRDecomposition {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @brief internal array storage.
    */
    double** QR;

    /** Row and column dimensions.
    @brief column dimension.
    @brief row dimension.
    */
    int m, n;

    /** Array for internal storage of diagonal of R.
    @brief diagonal of R.
    */
    double* Rdiag;

public:
/* ------------------------
   Constructor
 * ------------------------ */

    /** QR Decomposition, computed by Householder reflections.
    @param A    Rectangular matrix
    @return     Structure to access R and the Householder vectors and compute Q.
    */
    QRDecomposition(Matrix A);

/* ------------------------
   Destructor
 * ------------------------ */

    ~QRDecomposition();

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Is the matrix full rank?
    @return     true if R, and hence A, has full rank.
    */
    bool isFullRank ();

    /** Return the Householder vectors
    @return     Lower trapezoidal matrix whose columns define the reflections
    */
    Matrix getH ();

    /** Return the upper triangular factor
    @return     R
    */
    Matrix getR ();

    /** Generate and return the (economy-sized) orthogonal factor
    @return     Q
    */
    Matrix getQ ();

    /** Least squares solution of A*X = B
    @param B    A Matrix with as many rows as A and any number of columns.
    @return     X that minimizes the two norm of Q*R*X-B.
    @exception  invalid_argument  Matrix row dimensions must agree.
    @exception  runtime_error  Matrix is rank deficient.
    */
    Matrix solve (Matrix B);


};

}

#endif // QRDECOMPOSITION_H
