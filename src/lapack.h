#ifndef _LAPACK_
#define _LAPACK_

// Diagonalize symmetric matrix
// jobz 'N'= eigenvalues only. 'V' = include eigenvectors
// uplo 'U'= upper triangle storrage, 'L'= lower triangle
// n = matrix order
// a = nxn matrix to diagonalize
// lda = n
// w = the eigenvalues
// work = work array of dimension lwork, work(1) is optimal value of lwork
// lwork dimension of work, if lwork=-1 the optimal value is returned in work(1)
// info 0 = successfull run
extern void ssyev_(
              char *jobz,
              char *uplo,
              int *n,
              float *a,
              int *lda,
              float *w,
              float *work,
              int *lwork,
              int *info
              );

// Diagonalize general matrix
// jobvl 'N' = no left eigenvectors. 'V' = include left eigenvectors
// jobvr 'N' = no right eigenvectors. 'V' = include right eigenvectors
// n = matrix order
// a = nxn matrix to diagonalize
// lda = n
// wr = real part of the computed eigenvalues
// wi = imaginary part of the computed eigenvalues
// vl = 
// ldvl = n
// vr 
// ldvr = n
// work = work array of dimension lwork, work(1) is optimal value of lwork
// lwork dimension of work, if lwork=-1 the optimal value is returned in work(1)
// info 0 = successful run
extern void sgeev_(
              char *jobvl,
              char *jobvr,
              int *n,
              float *a,
              int *lda,
              float *wr,
              float *wi,
              float *vl,
              int *ldvl,
              float *vr,
              int *ldvr,
              float *work,
              int *lwork,
              int *info
              );

#endif // LAPACK
