#ifndef EIGENVALUEDECOMPOSITION_H
#define EIGENVALUEDECOMPOSITION_H

namespace CPAMA {

class Matrix;

/** Eigenvalues and eigenvectors of a real matrix.
<P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and
    V.times(V.transpose()) equals the identity matrix.
<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
**/

class EigenvalueDecomposition {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Row and column dimension (square matrix).
    @brief matrix dimension.
    */
    int n;

    /** Symmetry flag.
    @brief internal symmetry flag.
    */
    bool issymmetric;

    /** Arrays for internal storage of eigenvalues.
    @brief internal storage of eigenvalues.
    */
    double *d, *e;

    /** Array for internal storage of eigenvectors.
    @brief internal storage of eigenvectors.
    */
    double** V;

    /** Array for internal storage of nonsymmetric Hessenberg form.
    @brief internal storage of nonsymmetric Hessenberg form.
    */
    double** H;

    /** Working storage for nonsymmetric algorithm.
    @brief working storage for nonsymmetric algorithm.
    */
    double* ort;

/* ------------------------
   Private Methods
 * ------------------------ */

    // Symmetric Householder reduction to tridiagonal form.
    void tred2 ();

    // Symmetric tridiagonal QL algorithm.
    void tql2 ();

    // Nonsymmetric reduction to Hessenberg form.
    void orthes ();

    // Complex scalar division.
    double cdivr, cdivi;
    void cdiv(double xr, double xi, double yr, double yi);

    // Nonsymmetric reduction from Hessenberg to real Schur form.
    void hqr2 ();

public:
/* ------------------------
   Constructor
 * ------------------------ */

    /** Check for symmetry, then construct the eigenvalue decomposition
    @param A    Square matrix
    @return     Structure to access D and V.
    */
    EigenvalueDecomposition(Matrix Arg);

/* ------------------------
   Destructor
 * ------------------------ */

    ~EigenvalueDecomposition();

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Return the eigenvector matrix
    @return     V
    */
    Matrix getV ();

    /** Return the real parts of the eigenvalues
    @return     real(diag(D))
    */
    double* getRealEigenvalues () const;

    /** Return the imaginary parts of the eigenvalues
    @return     imag(diag(D))
    */
    double* getImagEigenvalues () const;

    /** Return the block diagonal eigenvalue matrix
    @return     D
    */
    Matrix getD ();

};

}

#endif // EIGENVALUEDECOMPOSITION_H
