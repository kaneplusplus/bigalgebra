#include <string>
#include "bigmemory/BigMatrix.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef REFBLAS
#include "refblas64longlong.h"
#define INT long long
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#define INT int
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);

  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                      SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB,
                      SEXP BETA, SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM,
                      SEXP C_isBM, SEXP C_offset);
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM);
  SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO,
                      SEXP A_isBM);
  SEXP dadd(SEXP N, SEXP ALPHA, SEXP Y, SEXP Y_isBM, SEXP SIGN, SEXP ALPHA_LHS);

#ifdef __cplusplus
}
#endif


/* Pointer utility, returns a double pointer for either a BigMatrix or a
 * standard R matrix.
 */
double *
make_double_ptr (SEXP matrix, SEXP isBigMatrix)
{
  double *matrix_ptr;

  if (LOGICAL_VALUE (isBigMatrix) == (Rboolean) TRUE)   // Big Matrix
    {
      SEXP address = GET_SLOT (matrix, install ("address"));
      BigMatrix *pbm =
        reinterpret_cast < BigMatrix * >(R_ExternalPtrAddr (address));
      if (!pbm)
        return (NULL);

      // Check that have acceptable big.matrix
      if (pbm->row_offset () > 0 && pbm->ncol () > 1)
        {
          std::string errMsg =
            string ("sub.big.matrix objects cannoth have row ") +
            string
            ("offset greater than zero and number of columns greater than 1");
          Rf_error (errMsg.c_str ());
          return (NULL);
        }

      index_type offset = pbm->nrow () * pbm->col_offset ();
      matrix_ptr = reinterpret_cast < double *>(pbm->matrix ()) + offset;
    }
  else                          // Regular R Matrix
    {
      matrix_ptr = NUMERIC_DATA (matrix);
    }

  return (matrix_ptr);
};


/* Wrappers for miscellaneous BLAS and LAPACK routines. */
SEXP
dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
               SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
               SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM,
               SEXP C_offset)
{
  long j = *(DOUBLE_DATA (C_offset));
  double *pA = make_double_ptr (A, A_isBM);
  double *pB = make_double_ptr (B, B_isBM);
  double *pC;
  SEXP ans;
  INT MM = (INT) * (DOUBLE_DATA (M));
  INT NN = (INT) * (DOUBLE_DATA (N));
  INT KK = (INT) * (DOUBLE_DATA (K));
  INT LDAA = (INT) * (DOUBLE_DATA (LDA));
  INT LDBB = (INT) * (DOUBLE_DATA (LDB));
  INT LDCC = (INT) * (DOUBLE_DATA (LDC));
  if(LOGICAL_VALUE(C_isBM) == (Rboolean) TRUE)
  {
/* Return results in a big matrix */
    pC = make_double_ptr (C, C_isBM) + j;
    PROTECT(ans = C);
  } else {
/* Allocate an output R matrix and return results there
   XXX Add check for size of MM and NN XXX 
 */
    PROTECT(ans = allocMatrix(REALSXP, (int)MM, (int)NN));
    pC = NUMERIC_DATA(ans);
  }
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
  dgemm (*((char *) CHARACTER_VALUE (TRANSA)),
         *((char *) CHARACTER_VALUE (TRANSB)),
         MM, NN, KK, *(NUMERIC_DATA (ALPHA)), pA, LDAA, pB,
         LDBB, *(NUMERIC_DATA (BETA)), pC, LDCC);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
  int8_dgemm ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#else
/* Standard Fortran interface from R's blas */
  dgemm_ ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#endif
  unprotect(1);
  return ans;
}



/* Compute A*X + Y for scalar a, vectors X and Y of length N.
 * Y must be a big.matrix, X can be an R vector or big.matrix.
 * The contents of Y are *replaced* by this routine and a reference
 * to Y is returned.
 */
SEXP
daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM)
{
  SEXP ans, Tr;
  double *pY;
  double *pA = DOUBLE_DATA(A);
  double *pX = make_double_ptr (X, X_isBM);
  INT incx = 1;
  INT incy = 1;
  INT NN = (INT) * (DOUBLE_DATA (N));
  PROTECT(ans = Y);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pY = make_double_ptr (Y, Tr);
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
  daxpy_ (NN, pA, pX, incx, pY, incy);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
  int8_daxpy (&NN, pA, pX, &incx, pY, &incy);
#else
/* Standard Fortran interface from R's blas */
  daxpy_ (&NN, pA, pX, &incx, pY, &incy);
#endif
  unprotect(2);
  return ans;
}
  
SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM)
{
  SEXP ans;
  const char *_UPLO = CHARACTER_VALUE(UPLO);
  INT _N = (INT)* (DOUBLE_DATA(N));
  double *_A = make_double_ptr(A, A_isBM);
  INT _LDA = (INT) *(DOUBLE_DATA(LDA));
  INT _INFO = (INT) *(DOUBLE_DATA(INFO));
/* An example of an alternate C-blas interface (e.g., ACML) */
#ifdef CBLAS
//  dpotrf_ (_UPLO, &_N, _A, &_LDA, &_INFO);
#elif REFBLAS
/* Standard Fortran interface without underscoring */
//  int8_dpotrf (_UPLO, &_N, _A, &_LDA, &_INFO);
#else
/* Standard Fortran interface from R's blas */
//  dpotrf_ (_UPLO, &_N, _A, &_LDA, &_INFO);
#endif
  PROTECT(ans = A);
  unprotect(1);
  return ans;
}

// Note that the following two functions should be updated to use
// the proper BLAS functions.
SEXP
dadd(SEXP N, SEXP ALPHA, SEXP Y, SEXP Y_isBM, SEXP SIGN, SEXP ALPHA_LHS) {
  SEXP ans;
  double *pY;
  INT NN = (INT) * (DOUBLE_DATA(N));
  double alpha = *(DOUBLE_DATA(ALPHA));
  int alpha_lhs = *(INTEGER_DATA(ALPHA_LHS));
  double sign = *(DOUBLE_DATA(SIGN));
  pY = make_double_ptr(Y, Y_isBM);
  PROTECT(ans = Y);
  if (alpha_lhs) 
  {
    for (INT i=0; i < NN; ++i)
    {
      pY[i] = alpha + sign * pY[i];
    }
  }
  else
  {
    for (INT i=0; i < NN; ++i)
    {
      pY[i] += sign*alpha;
    }
  }
  unprotect(1);
  return ans;
}
