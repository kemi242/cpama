#ifndef SINGULARVALUEDECOMPOSITION_H
#define SINGULARVALUEDECOMPOSITION_H

//#include "matrix.h"


namespace CPAMA {

class Matrix;

/** Singular Value Decomposition.
<P>
For an m-by-n matrix A with m >= n, the singular value decomposition is
an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
an n-by-n orthogonal matrix V so that A = U*S*V'.
<P>
The singular values, sigma[k] = S[k][k], are ordered so that
sigma[0] >= sigma[1] >= ... >= sigma[n-1].
<P>
The singular value decompostion always exists, so the constructor will
never fail.  The matrix condition number and the effective numerical
rank can be computed from this decomposition.
*/

class SingularValueDecomposition {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Arrays for internal storage of U and V.
    @brief internal storage of U.
    @brief internal storage of V.
    */
    double **U, **V;

    /** Array for internal storage of singular values.
    @brief internal storage of singular values.
    */
    double* s;

    /**
     * @brief Length of s
     */
    int slen;

    /** Row and column dimensions.
    @brief row dimension.
    @brief column dimension.
    */
    int m, n;



public:
/* ------------------------
   Constructor
 * ------------------------ */

    /** Construct the singular value decomposition
    @param A    Rectangular matrix
    @return     Structure to access U, S and V.
    */

    SingularValueDecomposition(Matrix Arg);

/* ------------------------
   Destructor
 * ------------------------ */
    ~SingularValueDecomposition();

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Return the left singular vectors
    @return     U
    */
    Matrix getU ();

    /** Return the right singular vectors
    @return     V
    */
    Matrix getV ();

    /** Return the one-dimensional array of singular values
    @return     diagonal of S.
    */
    double* getSingularValues ();

    /** Return the diagonal matrix of singular values
    @return     S
    */
    Matrix getS ();

    /** Two norm
    @return     max(S)
    */
    double norm2 ();

    /** Two norm condition number
    @return     max(S)/min(S)
    */
    double cond ();

    /** Effective numerical matrix rank
    @return     Number of nonnegligible singular values.
    */
    int rank ();

};

}

#endif // SINGULARVALUEDECOMPOSITION_H
