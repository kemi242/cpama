#ifndef MATRIX_H
#define MATRIX_H


namespace CPAMA {

class SingularValueDecomposition;
class LUDecomposition;
class QRDecomposition;
class CholeskyDecomposition;
class EigenvalueDecomposition;

/**
   CPAMA = C++ Matrix class. (C++ port of JAMA - Java Matrix Class)
<P>
   The C++ Matrix Class provides the fundamental operations of numerical
   linear algebra.  Various constructors create Matrices from two dimensional
   arrays of double precision floating point numbers.  Various "gets" and
   "sets" provide access to submatrices and matrix elements.  Several methods
   implement basic matrix arithmetic, including matrix addition and
   multiplication, matrix norms, and element-by-element array operations.
   Methods for reading and printing matrices are also included.  All the
   operations in this version of the Matrix Class involve real matrices.
   Complex matrices may be handled in a future version.
<P>
   Five fundamental matrix decompositions, which consist of pairs or triples
   of matrices, permutation vectors, and the like, produce results in five
   decomposition classes.  These decompositions are accessed by the Matrix
   class to compute solutions of simultaneous linear equations, determinants,
   inverses and other matrix functions.  The five decompositions are:
<P><UL>
   <LI>Cholesky Decomposition of symmetric, positive definite matrices.
   <LI>LU Decomposition of rectangular matrices.
   <LI>QR Decomposition of rectangular matrices.
   <LI>Singular Value Decomposition of rectangular matrices.
   <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
</UL>
<DL>
<DT><B>Example of use:</B></DT>
<P>
<DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
<P><PRE>
    double** vals = new double*[3];
    for (int i = 0; i < 3; i++) {
        vals[i] = new double[3];
    }

    vals[0][0] = 1; vals[0][1] = 2; vals[0][2] = 3;
    vals[1][0] = 4; vals[1][1] = 5; vals[1][2] = 6;
    vals[2][0] = 7; vals[2][1] = 8; vals[2][2] = 10;

    Matrix A = Matrix(vals, 3, 3);
    Matrix b = Matrix::random(3,1);
    Matrix x = A.solve(b);
    Matrix r = (A * x) - b;
    double rnorm = r.normInf();
</PRE></DD>
</DL>

@author The MathWorks, Inc. and the National Institute of Standards and Technology. Ported to C++ by kemi.
@version 5 August 1998
*/

class Matrix {
private:
/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of elements.
    @brief internal array storage.
    */
    double** A;

    /** Row and column dimensions.
    @brief row dimension.
    @brief column dimension.
    */
    int m, n;



public:
/* ------------------------
   Constructors
 * ------------------------ */

    /** Construct an m-by-n matrix of zeros.
    @param m    Number of rows.
    @param n    Number of colums.
    */
    Matrix(int m, int n);

    /** Construct an m-by-n constant matrix.
    @param m    Number of rows.
    @param n    Number of colums.
    @param s    Fill the matrix with this scalar value.
    */
    Matrix(int m, int n, double s);

    /** Construct a matrix from a 2-D array.
    @param A    Two-dimensional array of doubles.
    @param m    Number of rows
    @param n    Number of columns
    @exception  invalid_argument All rows must have the same length
    @see        #constructWithCopy
    */
    Matrix(double** A, int m, int n);

    /** Construct a matrix from a one-dimensional packed array
    @param vals     One-dimensional array of doubles, packed by columns (ala Fortran).
    @param m        Number of rows.
    @param valslen  Length of vals
    @exception      invalid_argument Array length must be a multiple of m.
    */
    Matrix(double* vals, int m, int valslen);

    /**
      Copy constructor
      */
    Matrix(const Matrix& other);

/* ------------------------
   Destructor
 * ------------------------ */

    ~Matrix();

/* ------------------------
   Public Methods
 * ------------------------ */


    /** Construct a matrix from a copy of a 2-D array.
    @param A    Two-dimensional array of doubles.
    @exception  invalid_argument All rows must have the same length
    */
    static Matrix constructWithCopy(double** A, int m, int n);

    /** Make a deep copy of a matrix
    */
    Matrix copy ();


    /** Access the internal two-dimensional array.
    @return     Pointer to the two-dimensional array of matrix elements.
    */
    double** getArray();

    /** Copy the internal two-dimensional array.
    @return     Two-dimensional array copy of matrix elements.
    */
    double** getArrayCopy();

    /** Make a one-dimensional column packed copy of the internal array.
    @return     Matrix elements packed in a one-dimensional array by columns.
    */
    double* getColumnPackedCopy();

    /** Make a one-dimensional row packed copy of the internal array.
    @return     Matrix elements packed in a one-dimensional array by rows.
    */
    double* getRowPackedCopy();

    /** Get row dimension.
    @return     m, the number of rows.
    */
    int getRowDimension() const;

    /** Get column dimension.
    @return     n, the number of columns.
    */
    int getColumnDimension() const;

    /** Get a single element.
    @param i    Row index.
    @param j    Column index.
    @return     A(i,j)
    */
    double get(int i, int j);

    /** Get a submatrix.
    @param i0   Initial row index
    @param i1   Final row index
    @param j0   Initial column index
    @param j1   Final column index
    @return     A(i0:i1,j0:j1)
    */
    Matrix getMatrix(int i0, int i1, int j0, int j1);

    /** Get a submatrix.
    @param r    Array of row indices.
    @param rlen Length of r
    @param c    Array of column indices.
    @param clen Length of c
    @return     A(r(:),c(:))
    */
    Matrix getMatrix(int* r, int rlen, int* c, int clen);

    /** Get a submatrix.
    @param i0   Initial row index
    @param i1   Final row index
    @param c    Array of column indices.
    @param clen Length of c
    @return     A(i0:i1,c(:))
    */
    Matrix getMatrix (int i0, int i1, int* c, int clen);

    /** Get a submatrix.
    @param r    Array of row indices.
    @param rlen Length of r
    @param i0   Initial column index
    @param i1   Final column index
    @return     A(r(:),j0:j1)
    */
    Matrix getMatrix (int* r, int rlen, int j0, int j1);

    /** Set a single element.
    @param i    Row index.
    @param j    Column index.
    @param s    A(i,j).
    */
    void set(int i, int j, double s);

    /** Set a submatrix.
    @param i0   Initial row index
    @param i1   Final row index
    @param j0   Initial column index
    @param j1   Final column index
    @param X    A(i0:i1,j0:j1)
    */
    void setMatrix(int i0, int i1, int j0, int j1, Matrix X);

    /** Set a submatrix.
    @param r    Array of row indices.
    @param rlen Length of r
    @param c    Array of column indices.
    @param clen Length of c
    @param X    A(r(:),c(:))
    */
    void setMatrix (int* r, int rlen, int* c, int clen, Matrix X);

    /** Set a submatrix.
    @param r    Array of row indices.
    @param rlen Length of r
    @param j0   Initial column index
    @param j1   Final column index
    @param X    A(r(:),j0:j1)
    */
    void setMatrix (int* r, int rlen, int j0, int j1, Matrix X);

    /** Set a submatrix.
    @param i0   Initial row index
    @param i1   Final row index
    @param c    Array of column indices.
    @param clen Length of c
    @param X    A(i0:i1,c(:))
    */
    void setMatrix (int i0, int i1, int* c, int clen, Matrix X);

    /** Matrix transpose.
    @return    A'
    */
    Matrix transpose();

    /** One norm
    @return    maximum column sum.
    */
    double norm1 ();

    /** Two norm
    @return    maximum singular value.
    */
    double norm2 ();

    /** Infinity norm
    @return    maximum row sum.
    */
    double normInf ();

    /** Frobenius norm
    @return    sqrt of sum of squares of all elements.
    */
    double normF ();

    /**  Unary minus
    @return    -A
    */
    Matrix operator -();

    /** C = A + B
    @param B    another matrix
    @return     A + B
    */
    Matrix operator +(Matrix B);

    /** A = A + B
    @param B    another matrix
    @return     A + B
    */
    Matrix operator +=(Matrix B);

    /** C = A - B
    @param B    another matrix
    @return     A - B
    */
    Matrix operator -(Matrix B);

    /** A = A - B
    @param B    another matrix
    @return     A - B
    */
    Matrix operator -=(Matrix B);

    /** Element-by-element multiplication, C = A.*B
    @param B    another matrix
    @return     A.*B
    */
    Matrix arrayTimes (Matrix B);

    /** Element-by-element multiplication in place, A = A.*B
    @param B    another matrix
    @return     A.*B
    */
    Matrix arrayTimesEquals (Matrix B);

    /** Element-by-element right division, C = A./B
    @param B    another matrix
    @return     A./B
    */
    Matrix arrayRightDivide (Matrix B);

    /** Element-by-element right division in place, A = A./B
    @param B    another matrix
    @return     A./B
    */
    Matrix arrayRightDivideEquals (Matrix B);

    /** Element-by-element left division, C = A.\B
    @param B    another matrix
    @return     A.\B
    */
    Matrix arrayLeftDivide (Matrix B);

    /** Element-by-element left division in place, A = A.\B
    @param B    another matrix
    @return     A.\B
    */
    Matrix arrayLeftDivideEquals (Matrix B);

    /** Multiply a matrix by a scalar, C = s*A
    @param s    scalar
    @return     s*A
    */
    Matrix operator *(double s);

    /** Multiply a matrix by a scalar in place, A = s*A
    @param s    scalar
    @return     replace A by s*A
    */
    Matrix operator *=(double s);

    /** Linear algebraic matrix multiplication, A * B
    @param B    another matrix
    @return     Matrix product, A * B
    @exception  invalid_argument Matrix inner dimensions must agree.
    */
    Matrix operator *(Matrix B);

    /** LU Decomposition
    @return     LUDecomposition
    @see LUDecomposition
    */
    LUDecomposition lu ();

    /** QR Decomposition
    @return     QRDecomposition
    @see QRDecomposition
    */
    QRDecomposition qr ();

    /** Cholesky Decomposition
    @return     CholeskyDecomposition
    @see CholeskyDecomposition
    */
    CholeskyDecomposition chol ();

    /** Singular Value Decomposition
    @return     SingularValueDecomposition
    @see SingularValueDecomposition
    */
    SingularValueDecomposition svd ();

    /** Eigenvalue Decomposition
    @return     EigenvalueDecomposition
    @see EigenvalueDecomposition
    */
    EigenvalueDecomposition eig ();

    /** Solve A*X = B
    @param B    right hand side
    @return     solution if A is square, least squares solution otherwise
    */
    Matrix solve (Matrix B);

    /** Solve X*A = B, which is also A'*X' = B'
    @param B    right hand side
    @return     solution if A is square, least squares solution otherwise.
    */
    Matrix solveTranspose (Matrix B);

    /** Matrix inverse or pseudoinverse
    @return     inverse(A) if A is square, pseudoinverse otherwise.
    */
    Matrix inverse ();

    /** Matrix determinant
    @return     determinant
    */
    double det ();

    /** Matrix rank
    @return     effective numerical rank, obtained from SVD.
    */
    int rank ();

    /** Matrix condition (2 norm)
    @return     ratio of largest to smallest singular value.
    */
    double cond ();

    /** Matrix trace.
    @return     sum of the diagonal elements.
    */
    double trace ();

    /** Generate matrix with random elements
    @param m    Number of rows.
    @param n    Number of colums.
    @return     An m-by-n matrix with uniformly distributed random elements.
    */
    static Matrix random (int m, int n);

    /** Generate identity matrix
    @param m    Number of rows.
    @param n    Number of colums.
    @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
    */
    static Matrix identity (int m, int n);

    /** Print the matrix to stdout.   Line the elements up in columns
      * with a Fortran-like 'Fw.d' style format.
    @param w    Column width.
    @param d    Number of digits after the decimal.
    */
    void print (int w, int d);

private:
/* ------------------------
   Private Methods
 * ------------------------ */

    /** Check if size(A) == size(B) **/
    void checkMatrixDimensions (Matrix B);

};

}

#endif // MATRIX_H
