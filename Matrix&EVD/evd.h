/*****************************************************************************
 *                                    evd.h
 *
 * Class template of eigenvalues and eigenvectors decomposition.
 *
 * For a real matrix A, we have A*V = V*D, where the eigenvalue matrix D is
 * diagonal and the eigenvector matrix V is linear independence. That is the
 * kth diagonal value of D is the eigenvalue and the kth column of V
 * represents the corresponding eigenvector of D[k][k]. If A is symmetric,
 * then V is a orthogonal matrix, which means A = V*D*V', and eigenvalues
 * are all real numbers.
 *
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal
 * with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 * a+i*b, in 2-by-2 blocks [a,b; -b,a]. That is, if the complex eigenvalues
 * look like
 *
 *    u + iv     .        .          .      .    .
 *      .      u - iv     .          .      .    .
 *      .        .      a + ib       .      .    .
 *      .        .        .        a - ib   .    .
 *      .        .        .          .      x    .
 *      .        .        .          .      .    y
 *
 * then D looks like
 *
 *      u        v        .          .      .    .
 *     -v        u        .          .      .    .
 *      .        .        a          b      .    .
 *      .        .       -b          a      .    .
 *      .        .        .          .      x    .
 *      .        .        .          .      .    y
 *
 * This keeps V a real matrix in both symmetric and non-symmetric cases, and
 * A*V = V*D.
 *
 * The matrix V may be badly conditioned, or even singular, so the validity
 * of the equation A=V*D*inverse(V) depends upon the condition number of V.
 *
 * Adapted form Template Numerical Toolkit.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/

#ifndef EVD_H
#define EVD_H


#include "matrix.h"


namespace splab
{

	template <typename Real>
	class EVD
	{

    public:

        EVD();
		~EVD();

        // decomposition
		void dec( const Matrix<Real> &A );

		// the eigenvalues are real or complex
		bool isSymmetric() const;
		bool isComplex( Real tol=Real(EPS) );

        // get eigenvectors
		Matrix<Real> getV() const;
		Matrix<complex<Real> > getCV();

        // get eigenvalues
        Vector<Real> getD() const;
        Vector<complex<Real> > getCD();
//        Matrix<Real> getDM();
//        Matrix<complex<Real> > getCDM();

    private:

		int     n;
		bool    symmetric;
        Real    cdivr,
                cdivi;

		// eigenvalues and its real and image part
		Matrix<Real> V;
		Vector<Real> d;
		Vector<Real> e;

		// temporary storage for internal variables
		Vector<Real> ort;
		Matrix<Real> H;

		void tred2();
		void tql2();
		void cdiv( Real xr, Real xi, Real yr, Real yi );
		void hqr2();
		void others();
		void normalized();

	};
	// class EVD


    #include "evd_impl.h"

}
// namespace splab


#endif
// EVD_H
