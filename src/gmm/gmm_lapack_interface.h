/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2003-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file gmm_lapack_interface.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date October 7, 2003.
   @brief gmm interface for LAPACK
*/

#ifndef GMM_LAPACK_INTERFACE_H
#define GMM_LAPACK_INTERFACE_H

#include "gmm_blas_interface.h"
#include "gmm_dense_lu.h"
#include "gmm_dense_qr.h"

#if defined(GMM_USES_LAPACK) && !defined(GMM_MATLAB_INTERFACE)

namespace gmm {

  /* ********************************************************************** */
  /* Operations interfaced for T = float, double, std::complex<float>       */
  /*    or std::complex<double> :                                           */
  /*                                                                        */
  /* lu_factor(dense_matrix<T>, std::vector<long>)                          */
  /* lu_solve(dense_matrix<T>, std::vector<T>, std::vector<T>)              */
  /* lu_solve(dense_matrix<T>, std::vector<long>, std::vector<T>,           */
  /*          std::vector<T>)                                               */
  /* lu_solve_transposed(dense_matrix<T>, std::vector<long>, std::vector<T>,*/
  /*          std::vector<T>)                                               */
  /* lu_inverse(dense_matrix<T>)                                            */
  /* lu_inverse(dense_matrix<T>, std::vector<long>, dense_matrix<T>)        */
  /*                                                                        */
  /* qr_factor(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)           */
  /*                                                                        */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>)                 */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>,                 */
  /*                       dense_matrix<T>)                                 */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >)  */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >,  */
  /*                       dense_matrix<T>)                                 */
  /*                                                                        */
  /* geev_interface_right                                                   */
  /* geev_interface_left                                                    */
  /*                                                                        */
  /* schur(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)               */
  /*                                                                        */
  /* svd(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>, std::vector<T>) */
  /* svd(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>,                 */
  /*     std::vector<std::complex<T> >)                                     */
  /*                                                                        */
  /* ********************************************************************** */

  /* ********************************************************************** */
  /* LAPACK functions used.                                                 */
  /* ********************************************************************** */

  extern "C" {
    /* IMPORTANT: every function below MUST have a real, non-variadic       *
     * prototype. A "(...)" declaration uses the variadic ABI, which on     *
     * Darwin/arm64 differs from the non-variadic ABI (all args on stack,   *
     * not in registers x0-x7). This causes EXC_BAD_ACCESS in libLAPACK on  *
     * Apple Silicon. See gmm_blas_interface.h for the full diagnosis.      *
     *                                                                      *
     * Several parameters below are declared as `void *` rather than their  *
     * proper LAPACK type because gmm's interface macros pass buffers whose *
     * C++ types do not match the LAPACK Fortran spec. These are noted      *
     * inline. The mismatches are pre-existing and harmless (the relevant   *
     * outputs are unused by gmm), but cleaning them up requires editing    *
     * the *_interface macros and is out of scope for this ABI fix.         */

    /* LU factorization */
#define GMM_LAPACK_DECL_GETRF(name, T)                                      \
    void name(const BLAS_INT *m, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda,                                    \
              BLAS_INT *ipiv, BLAS_INT *info)
    GMM_LAPACK_DECL_GETRF(sgetrf_, BLAS_S);
    GMM_LAPACK_DECL_GETRF(dgetrf_, BLAS_D);
    GMM_LAPACK_DECL_GETRF(cgetrf_, BLAS_C);
    GMM_LAPACK_DECL_GETRF(zgetrf_, BLAS_Z);
#undef GMM_LAPACK_DECL_GETRF

    /* LU solve */
#define GMM_LAPACK_DECL_GETRS(name, T)                                      \
    void name(const char *trans, const BLAS_INT *n, const BLAS_INT *nrhs,   \
              const T *A, const BLAS_INT *lda, const BLAS_INT *ipiv,        \
              T *B, const BLAS_INT *ldb, BLAS_INT *info)
    GMM_LAPACK_DECL_GETRS(sgetrs_, BLAS_S);
    GMM_LAPACK_DECL_GETRS(dgetrs_, BLAS_D);
    GMM_LAPACK_DECL_GETRS(cgetrs_, BLAS_C);
    GMM_LAPACK_DECL_GETRS(zgetrs_, BLAS_Z);
#undef GMM_LAPACK_DECL_GETRS

    /* Inverse from LU factorization */
#define GMM_LAPACK_DECL_GETRI(name, T)                                      \
    void name(const BLAS_INT *n, T *A, const BLAS_INT *lda,                 \
              const BLAS_INT *ipiv,                                         \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GETRI(sgetri_, BLAS_S);
    GMM_LAPACK_DECL_GETRI(dgetri_, BLAS_D);
    GMM_LAPACK_DECL_GETRI(cgetri_, BLAS_C);
    GMM_LAPACK_DECL_GETRI(zgetri_, BLAS_Z);
#undef GMM_LAPACK_DECL_GETRI

    /* QR factorization */
#define GMM_LAPACK_DECL_GEQRF(name, T)                                      \
    void name(const BLAS_INT *m, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda, T *tau,                            \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GEQRF(sgeqrf_, BLAS_S);
    GMM_LAPACK_DECL_GEQRF(dgeqrf_, BLAS_D);
    GMM_LAPACK_DECL_GEQRF(cgeqrf_, BLAS_C);
    GMM_LAPACK_DECL_GEQRF(zgeqrf_, BLAS_Z);
#undef GMM_LAPACK_DECL_GEQRF

    /* Generate Q from elementary reflectors (real: orgqr, complex: ungqr) */
#define GMM_LAPACK_DECL_ORGQR(name, T)                                      \
    void name(const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *k,      \
              T *A, const BLAS_INT *lda, const T *tau,                      \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_ORGQR(sorgqr_, BLAS_S);
    GMM_LAPACK_DECL_ORGQR(dorgqr_, BLAS_D);
    GMM_LAPACK_DECL_ORGQR(cungqr_, BLAS_C);
    GMM_LAPACK_DECL_ORGQR(zungqr_, BLAS_Z);
#undef GMM_LAPACK_DECL_ORGQR

    /* Apply Q from QR (real: ormqr, complex: unmqr); declared in case      *
     * of future use, not currently called by gmm's macros.                 */
#define GMM_LAPACK_DECL_ORMQR(name, T)                                      \
    void name(const char *side, const char *trans,                          \
              const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *k,      \
              const T *A, const BLAS_INT *lda, const T *tau,                \
              T *C, const BLAS_INT *ldc,                                    \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_ORMQR(sormqr_, BLAS_S);
    GMM_LAPACK_DECL_ORMQR(dormqr_, BLAS_D);
    GMM_LAPACK_DECL_ORMQR(cunmqr_, BLAS_C);
    GMM_LAPACK_DECL_ORMQR(zunmqr_, BLAS_Z);
#undef GMM_LAPACK_DECL_ORMQR

    /* Schur form (used by implicit_qr_algorithm). gmm passes a            *
     * `bool (*)(...)` function pointer for SELECT, which we match here.   *
     * Real and complex variants have different shapes: real has separate  *
     * WR/WI eigenvalue arrays, complex has a single W array plus an       *
     * RWORK workspace. gmm passes std::vector<double>* for the eigenvalue *
     * outputs and BWORK regardless of scalar type, so those positions     *
     * are `void *` to accept whatever gmm hands in.                       */
#define GMM_LAPACK_DECL_GEES_REAL(name, T)                                  \
    void name(const char *jobvs, const char *sort, bool (*select)(...),     \
              const BLAS_INT *n,                                            \
              T *A, const BLAS_INT *lda, BLAS_INT *sdim,                    \
              void *wr, void *wi,  /* gmm: vector<double>* for any T */     \
              T *vs, const BLAS_INT *ldvs,                                  \
              T *work, const BLAS_INT *lwork,                               \
              void *bwork,         /* gmm: vector<double>*; spec: logical* */\
              BLAS_INT *info)
    GMM_LAPACK_DECL_GEES_REAL(sgees_, BLAS_S);
    GMM_LAPACK_DECL_GEES_REAL(dgees_, BLAS_D);
#undef GMM_LAPACK_DECL_GEES_REAL
#define GMM_LAPACK_DECL_GEES_CPLX(name, T)                                  \
    void name(const char *jobvs, const char *sort, bool (*select)(...),     \
              const BLAS_INT *n,                                            \
              T *A, const BLAS_INT *lda, BLAS_INT *sdim,                    \
              void *w,             /* gmm: vector<double>*; spec: complex* */\
              T *vs, const BLAS_INT *ldvs,                                  \
              T *work, const BLAS_INT *lwork,                               \
              void *rwork, void *bwork,  /* gmm: vector<double>* both */    \
              BLAS_INT *info)
    GMM_LAPACK_DECL_GEES_CPLX(cgees_, BLAS_C);
    GMM_LAPACK_DECL_GEES_CPLX(zgees_, BLAS_Z);
#undef GMM_LAPACK_DECL_GEES_CPLX

    /* General eigenvalue problem (used by geev_interface_left/right).     *
     * gmm's call shape matches LAPACK exactly here.                       */
#define GMM_LAPACK_DECL_GEEV_REAL(name, T)                                  \
    void name(const char *jobvl, const char *jobvr, const BLAS_INT *n,      \
              T *A, const BLAS_INT *lda, T *wr, T *wi,                      \
              T *vl, const BLAS_INT *ldvl, T *vr, const BLAS_INT *ldvr,     \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GEEV_REAL(sgeev_, BLAS_S);
    GMM_LAPACK_DECL_GEEV_REAL(dgeev_, BLAS_D);
#undef GMM_LAPACK_DECL_GEEV_REAL
#define GMM_LAPACK_DECL_GEEV_CPLX(name, T, R)                               \
    void name(const char *jobvl, const char *jobvr, const BLAS_INT *n,      \
              T *A, const BLAS_INT *lda, T *w,                              \
              T *vl, const BLAS_INT *ldvl, T *vr, const BLAS_INT *ldvr,     \
              T *work, const BLAS_INT *lwork, R *rwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GEEV_CPLX(cgeev_, BLAS_C, BLAS_S);
    GMM_LAPACK_DECL_GEEV_CPLX(zgeev_, BLAS_Z, BLAS_D);
#undef GMM_LAPACK_DECL_GEEV_CPLX

    /* Extended Schur (used by schur). Two oddities in gmm's call shape:   *
     *  - SELECT receives `&select` where `select` is a stack `bool`,      *
     *    rather than the function pointer LAPACK expects -> `void *`.     *
     *  - In the complex variants RCONDE/RCONDV are passed as complex*     *
     *    where LAPACK expects real* -> `void *`.                          */
#define GMM_LAPACK_DECL_GEESX_REAL(name, T)                                 \
    void name(const char *jobvs, const char *sort,                          \
              void *select,        /* gmm: bool*; spec: function ptr */     \
              const char *sense, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda, BLAS_INT *sdim,                    \
              T *wr, T *wi, T *vs, const BLAS_INT *ldvs,                    \
              T *rconde, T *rcondv,                                         \
              T *work, const BLAS_INT *lwork,                               \
              BLAS_INT *iwork, const BLAS_INT *liwork,                      \
              BLAS_INT *bwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GEESX_REAL(sgeesx_, BLAS_S);
    GMM_LAPACK_DECL_GEESX_REAL(dgeesx_, BLAS_D);
#undef GMM_LAPACK_DECL_GEESX_REAL
#define GMM_LAPACK_DECL_GEESX_CPLX(name, T, R)                              \
    void name(const char *jobvs, const char *sort,                          \
              void *select,        /* gmm: bool*; spec: function ptr */     \
              const char *sense, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda, BLAS_INT *sdim,                    \
              T *w, T *vs, const BLAS_INT *ldvs,                            \
              void *rconde, void *rcondv,  /* gmm: complex*; spec: real* */ \
              T *work, const BLAS_INT *lwork,                               \
              R *rwork, BLAS_INT *bwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GEESX_CPLX(cgeesx_, BLAS_C, BLAS_S);
    GMM_LAPACK_DECL_GEESX_CPLX(zgeesx_, BLAS_Z, BLAS_D);
#undef GMM_LAPACK_DECL_GEESX_CPLX

    /* SVD (used by svd). gmm's call shape matches LAPACK exactly. */
#define GMM_LAPACK_DECL_GESVD_REAL(name, T)                                 \
    void name(const char *jobu, const char *jobvt,                          \
              const BLAS_INT *m, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda, T *S,                              \
              T *U, const BLAS_INT *ldu, T *VT, const BLAS_INT *ldvt,       \
              T *work, const BLAS_INT *lwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GESVD_REAL(sgesvd_, BLAS_S);
    GMM_LAPACK_DECL_GESVD_REAL(dgesvd_, BLAS_D);
#undef GMM_LAPACK_DECL_GESVD_REAL
#define GMM_LAPACK_DECL_GESVD_CPLX(name, T, R)                              \
    void name(const char *jobu, const char *jobvt,                          \
              const BLAS_INT *m, const BLAS_INT *n,                         \
              T *A, const BLAS_INT *lda, R *S,                              \
              T *U, const BLAS_INT *ldu, T *VT, const BLAS_INT *ldvt,       \
              T *work, const BLAS_INT *lwork, R *rwork, BLAS_INT *info)
    GMM_LAPACK_DECL_GESVD_CPLX(cgesvd_, BLAS_C, BLAS_S);
    GMM_LAPACK_DECL_GESVD_CPLX(zgesvd_, BLAS_Z, BLAS_D);
#undef GMM_LAPACK_DECL_GESVD_CPLX
  }

  /* ********************************************************************** */
  /* LU decomposition.                                                      */
  /* ********************************************************************** */

# define getrf_interface(lapack_name, base_type) inline                       \
  size_type lu_factor(dense_matrix<base_type> &A, lapack_ipvt &ipvt) {        \
    GMMLAPACK_TRACE("getrf_interface");                                       \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A)), lda(m);\
    BLAS_INT info(-1);                                                        \
    if (m && n) lapack_name(&m, &n, &A(0,0), &lda, &ipvt[0], &info);          \
    return size_type(abs(info));                                              \
  }

  getrf_interface(sgetrf_, BLAS_S)
  getrf_interface(dgetrf_, BLAS_D)
  getrf_interface(cgetrf_, BLAS_C)
  getrf_interface(zgetrf_, BLAS_Z)

  /* ********************************************************************* */
  /* LU solve.                                                             */
  /* ********************************************************************* */

# define getrs_interface(f_name, trans, lapack_name, base_type) inline     \
  void f_name(const dense_matrix<base_type> &A,                            \
              const lapack_ipvt &ipvt, std::vector<base_type> &x,          \
              const std::vector<base_type> &b) {                           \
    GMMLAPACK_TRACE("getrs_interface");                                    \
    const BLAS_INT n=BLAS_INT(mat_nrows(A)), nrhs(1);                      \
    BLAS_INT info(0); gmm::copy(b, x); trans;                              \
    if (n)                                                                 \
      lapack_name(&t, &n, &nrhs, &A(0,0), &n, &ipvt[0], &x[0], &n, &info); \
  }

# define getrs_trans_n const char t = 'N'
# define getrs_trans_t const char t = 'T'

  getrs_interface(lu_solve, getrs_trans_n, sgetrs_, BLAS_S)
  getrs_interface(lu_solve, getrs_trans_n, dgetrs_, BLAS_D)
  getrs_interface(lu_solve, getrs_trans_n, cgetrs_, BLAS_C)
  getrs_interface(lu_solve, getrs_trans_n, zgetrs_, BLAS_Z)
  getrs_interface(lu_solve_transposed, getrs_trans_t, sgetrs_, BLAS_S)
  getrs_interface(lu_solve_transposed, getrs_trans_t, dgetrs_, BLAS_D)
  getrs_interface(lu_solve_transposed, getrs_trans_t, cgetrs_, BLAS_C)
  getrs_interface(lu_solve_transposed, getrs_trans_t, zgetrs_, BLAS_Z)

  /* ********************************************************************* */
  /* LU inverse.                                                           */
  /* ********************************************************************* */

# define getri_interface(lapack_name, base_type)                           \
  inline void lu_inverse(const dense_matrix<base_type> &LU,                \
                         const lapack_ipvt &ipvt,                          \
                         dense_matrix<base_type> &A) {                     \
    GMMLAPACK_TRACE("getri_interface");                                    \
    const BLAS_INT n=BLAS_INT(mat_nrows(A));                               \
    BLAS_INT info(0), lwork(-1); base_type work1;                          \
    if (n) {                                                               \
      gmm::copy(LU, A);                                                    \
      lapack_name(&n, &A(0,0), &n, &ipvt[0], &work1, &lwork, &info);       \
      lwork = int(gmm::real(work1));                                       \
      std::vector<base_type> work(lwork);                                  \
      lapack_name(&n, &A(0,0), &n, &ipvt[0], &work[0], &lwork, &info);     \
    }                                                                      \
  }

  getri_interface(sgetri_, BLAS_S)
  getri_interface(dgetri_, BLAS_D)
  getri_interface(cgetri_, BLAS_C)
  getri_interface(zgetri_, BLAS_Z)

  /* ********************************************************************** */
  /* QR factorization.                                                      */
  /* ********************************************************************** */

# define geqrf_interface(lapack_name, base_type)                           \
  inline void qr_factor(dense_matrix<base_type> &A) {                      \
    GMMLAPACK_TRACE("geqrf_interface");                                    \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A));     \
    BLAS_INT info(0), lwork(-1); base_type work1;                          \
    if (m && n) {                                                          \
      std::vector<base_type> tau(n);                                       \
      lapack_name(&m, &n, &A(0,0), &m, &tau[0], &work1, &lwork, &info);    \
      lwork = BLAS_INT(gmm::real(work1));                                  \
      std::vector<base_type> work(lwork);                                  \
      lapack_name(&m, &n, &A(0,0), &m, &tau[0], &work[0], &lwork, &info);  \
      GMM_ASSERT1(!info, "QR factorization failed");                       \
    }                                                                      \
  }

  geqrf_interface(sgeqrf_, BLAS_S)
  geqrf_interface(dgeqrf_, BLAS_D)
    // For complex values, housholder vectors are not the same as in
    // gmm::lu_factor. Impossible to interface for the moment.
    //  geqrf_interface(cgeqrf_, BLAS_C)
    //  geqrf_interface(zgeqrf_, BLAS_Z)

# define geqrf_interface2(lapack_name1, lapack_name2, base_type) inline    \
  void qr_factor(const dense_matrix<base_type> &A,                         \
                 dense_matrix<base_type> &Q, dense_matrix<base_type> &R) { \
    GMMLAPACK_TRACE("geqrf_interface2");                                   \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A));     \
    BLAS_INT info(0), lwork(-1); base_type work1;                          \
    if (m && n) {                                                          \
      std::copy(A.begin(), A.end(), Q.begin());                            \
      std::vector<base_type> tau(n);                                       \
      lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work1  , &lwork, &info); \
      lwork = BLAS_INT(gmm::real(work1));                                  \
      std::vector<base_type> work(lwork);                                  \
      lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work[0], &lwork, &info); \
      GMM_ASSERT1(!info, "QR factorization failed");                       \
      base_type *p = &R(0,0), *q = &Q(0,0);                                \
      for (BLAS_INT j = 0; j < n; ++j, q += m-n)                           \
        for (BLAS_INT i = 0; i < n; ++i, ++p, ++q)                         \
          *p = (j < i) ? base_type(0) : *q;                                \
      lapack_name2(&m, &n, &n, &Q(0,0), &m,&tau[0],&work[0],&lwork,&info); \
    }                                                                      \
    else gmm::clear(Q);                                                    \
  }

  geqrf_interface2(sgeqrf_, sorgqr_, BLAS_S)
  geqrf_interface2(dgeqrf_, dorgqr_, BLAS_D)
  geqrf_interface2(cgeqrf_, cungqr_, BLAS_C)
  geqrf_interface2(zgeqrf_, zungqr_, BLAS_Z)

  /* ********************************************************************** */
  /* QR algorithm for eigenvalues search.                                   */
  /* ********************************************************************** */

# define gees_interface(lapack_name, base_type)                            \
  template <typename VECT> inline void implicit_qr_algorithm(              \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q,                                       \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    GMMLAPACK_TRACE("gees_interface");                                     \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    BLAS_INT n=BLAS_INT(mat_nrows(A)), info(0), lwork(-1), sdim;           \
    base_type work1;                                                       \
    if (!n) return;                                                        \
    dense_matrix<base_type> H(n,n); gmm::copy(A, H);                       \
    char jobvs = (compvect ? 'V' : 'N'), sort = 'N';                       \
    std::vector<double> rwork(n), eigv1(n), eigv2(n);                      \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
                &eigv2[0], &Q(0,0), &n, &work1, &lwork, &rwork[0], &info); \
    lwork = BLAS_INT(gmm::real(work1));                                    \
    std::vector<base_type> work(lwork);                                    \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
                &eigv2[0], &Q(0,0), &n, &work[0], &lwork, &rwork[0],&info);\
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    extract_eig(H, eigval_, tol);                                          \
  }

# define gees_interface2(lapack_name, base_type)                           \
  template <typename VECT> inline void implicit_qr_algorithm(              \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q,                                       \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    GMMLAPACK_TRACE("gees_interface2");                                    \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    BLAS_INT n=BLAS_INT(mat_nrows(A)), info(0), lwork(-1), sdim;           \
    base_type work1;                                                       \
    if (!n) return;                                                        \
    dense_matrix<base_type> H(n,n); gmm::copy(A, H);                       \
    char jobvs = (compvect ? 'V' : 'N'), sort = 'N';                       \
    std::vector<double> rwork(n), eigvv(n*2);                              \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work1, &lwork, &rwork[0], &rwork[0], &info); \
    lwork = BLAS_INT(gmm::real(work1));                                    \
    std::vector<base_type> work(lwork);                                    \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work[0], &lwork, &rwork[0], &rwork[0],&info);\
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    extract_eig(H, eigval_, tol);                                          \
  }

  gees_interface(sgees_, BLAS_S)
  gees_interface(dgees_, BLAS_D)
  gees_interface2(cgees_, BLAS_C)
  gees_interface2(zgees_, BLAS_Z)


# define jobv_right char jobvl = 'N', jobvr = 'V';
# define jobv_left char jobvl = 'V', jobvr = 'N';

# define geev_interface(lapack_name, base_type, side)                      \
  template <typename VECT> inline void geev_interface_ ## side(            \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q) {                                     \
    GMMLAPACK_TRACE("geev_interface");                                     \
    BLAS_INT n = BLAS_INT(mat_nrows(A)), info(0), lwork(-1);               \
    base_type work1;                                                       \
    if (!n) return;                                                        \
    dense_matrix<base_type> H(n,n); gmm::copy(A, H);                       \
    jobv_ ## side                                                          \
    std::vector<base_type> eigvr(n), eigvi(n);                             \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigvr[0], &eigvi[0],     \
                &Q(0,0), &n, &Q(0,0), &n, &work1, &lwork, &info);          \
    lwork = BLAS_INT(gmm::real(work1));                                    \
    std::vector<base_type> work(lwork);                                    \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigvr[0], &eigvi[0],     \
                &Q(0,0), &n, &Q(0,0), &n, &work[0], &lwork, &info);        \
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    gmm::copy(eigvr, gmm::real_part(eigval_));                             \
    gmm::copy(eigvi, gmm::imag_part(eigval_));                             \
  }

# define geev_interface2(lapack_name, base_type, side)                     \
  template <typename VECT> inline void geev_interface_ ## side(            \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q) {                                     \
    GMMLAPACK_TRACE("geev_interface");                                     \
    BLAS_INT n = BLAS_INT(mat_nrows(A)), info(0), lwork(-1);               \
    base_type work1;                                                       \
    if (!n) return;                                                        \
    dense_matrix<base_type> H(n,n); gmm::copy(A, H);                       \
    jobv_ ## side                                                          \
    std::vector<base_type::value_type> rwork(2*n);                         \
    std::vector<base_type> eigv(n);                                        \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigv[0], &Q(0,0), &n,    \
                &Q(0,0), &n, &work1, &lwork, &rwork[0], &info);            \
    lwork = BLAS_INT(gmm::real(work1));                                    \
    std::vector<base_type> work(lwork);                                    \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigv[0], &Q(0,0), &n,    \
                &Q(0,0), &n, &work[0], &lwork,  &rwork[0],  &info);        \
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    gmm::copy(eigv, eigval_);                                              \
  }

  geev_interface(sgeev_, BLAS_S, right)
  geev_interface(dgeev_, BLAS_D, right)
  geev_interface2(cgeev_, BLAS_C, right)
  geev_interface2(zgeev_, BLAS_Z, right)

  geev_interface(sgeev_, BLAS_S, left)
  geev_interface(dgeev_, BLAS_D, left)
  geev_interface2(cgeev_, BLAS_C, left)
  geev_interface2(zgeev_, BLAS_Z, left)


  /* ********************************************************************** */
  /* SCHUR algorithm:                                                       */
  /*  A = Q*S*(Q^T), with Q orthogonal and S upper quasi-triangula          */
  /* ********************************************************************** */

# define geesx_interface(lapack_name, base_type)                        \
  inline void schur(dense_matrix<base_type> &A,                         \
                    dense_matrix<base_type> &S,                         \
                    dense_matrix<base_type> &Q) {                       \
    GMMLAPACK_TRACE("geesx_interface");                                 \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A));  \
    GMM_ASSERT1(m == n, "Schur decomposition requires square matrix");  \
    char jobvs = 'V', sort = 'N', sense = 'N';                          \
    bool select = false;                                                \
    BLAS_INT lwork = 8*n, sdim = 0, liwork = 1;                         \
    std::vector<base_type> work(lwork), wr(n), wi(n);                   \
    std::vector<BLAS_INT> iwork(liwork);                                \
    std::vector<BLAS_INT> bwork(1);                                     \
    resize(S, n, n); copy(A, S);                                        \
    resize(Q, n, n);                                                    \
    base_type rconde(0), rcondv(0);                                     \
    BLAS_INT info(0);                                                   \
    lapack_name(&jobvs, &sort, &select, &sense, &n, &S(0,0), &n,        \
                &sdim, &wr[0], &wi[0], &Q(0,0), &n, &rconde, &rcondv,   \
                &work[0], &lwork, &iwork[0], &liwork, &bwork[0], &info);\
    GMM_ASSERT1(!info, "SCHUR algorithm failed");                       \
  }

# define geesx_interface2(lapack_name, base_type)                       \
  inline void schur(dense_matrix<base_type> &A,                         \
                    dense_matrix<base_type> &S,                         \
                    dense_matrix<base_type> &Q) {                       \
    GMMLAPACK_TRACE("geesx_interface");                                 \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A));  \
    GMM_ASSERT1(m == n, "Schur decomposition requires square matrix");  \
    char jobvs = 'V', sort = 'N', sense = 'N';                          \
    bool select = false;                                                \
    BLAS_INT lwork = 8*n, sdim = 0;                                     \
    std::vector<base_type::value_type> rwork(lwork);                    \
    std::vector<base_type> work(lwork), w(n);                           \
    std::vector<BLAS_INT> bwork(1);                                     \
    resize(S, n, n); copy(A, S);                                        \
    resize(Q, n, n);                                                    \
    base_type rconde(0), rcondv(0);                                     \
    BLAS_INT info(0);                                                   \
    lapack_name(&jobvs, &sort, &select, &sense, &n, &S(0,0), &n,        \
                &sdim, &w[0], &Q(0,0), &n, &rconde, &rcondv,            \
                &work[0], &lwork, &rwork[0], &bwork[0], &info);         \
    GMM_ASSERT1(!info, "SCHUR algorithm failed");                       \
  }

  geesx_interface(sgeesx_, BLAS_S)
  geesx_interface(dgeesx_, BLAS_D)
  geesx_interface2(cgeesx_, BLAS_C)
  geesx_interface2(zgeesx_, BLAS_Z)

  template <typename MAT>
  void schur(const MAT &A_, MAT &S, MAT &Q) {
   MAT A(A_);
   schur(A, S, Q);
  }


  /* ********************************************************************** */
  /* Interface to SVD. Does not correspond to a Gmm++ functionnality.       */
  /* Author : Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>        */
  /* ********************************************************************** */

# define gesvd_interface(lapack_name, base_type)                        \
  inline void svd(dense_matrix<base_type> &X,                           \
                  dense_matrix<base_type> &U,                           \
                  dense_matrix<base_type> &Vtransposed,                 \
                  std::vector<base_type> &sigma) {                      \
    GMMLAPACK_TRACE("gesvd_interface");                                 \
    BLAS_INT m = BLAS_INT(mat_nrows(X)), n = BLAS_INT(mat_ncols(X));    \
    BLAS_INT mn_min = m < n ? m : n;                                    \
    sigma.resize(mn_min);                                               \
    std::vector<base_type> work(15 * mn_min);                           \
    BLAS_INT lwork = BLAS_INT(work.size());                             \
    resize(U, m, m);                                                    \
    resize(Vtransposed, n, n);                                          \
    char job = 'A';                                                     \
    BLAS_INT info(0);                                                   \
    lapack_name(&job, &job, &m, &n, &X(0,0), &m, &sigma[0], &U(0,0),    \
                &m, &Vtransposed(0,0), &n, &work[0], &lwork, &info);    \
  }

# define cgesvd_interface(lapack_name, base_type, base_type2)           \
  inline void svd(dense_matrix<base_type> &X,                           \
                  dense_matrix<base_type> &U,                           \
                  dense_matrix<base_type> &Vtransposed,                 \
                  std::vector<base_type2> &sigma) {                     \
    GMMLAPACK_TRACE("gesvd_interface");                                 \
    BLAS_INT m = BLAS_INT(mat_nrows(X)), n = BLAS_INT(mat_ncols(X));    \
    BLAS_INT mn_min = m < n ? m : n;                                    \
    sigma.resize(mn_min);                                               \
    std::vector<base_type> work(15 * mn_min);                           \
    std::vector<base_type2> rwork(5 * mn_min);                          \
    BLAS_INT lwork = BLAS_INT(work.size());                             \
    resize(U, m, m);                                                    \
    resize(Vtransposed, n, n);                                          \
    char job = 'A';                                                     \
    BLAS_INT info(0);                                                   \
    lapack_name(&job, &job, &m, &n, &X(0,0), &m, &sigma[0], &U(0,0),    \
                &m, &Vtransposed(0,0), &n, &work[0], &lwork,            \
                &rwork[0], &info);                                      \
  }

  gesvd_interface(sgesvd_, BLAS_S)
  gesvd_interface(dgesvd_, BLAS_D)
  cgesvd_interface(cgesvd_, BLAS_C, BLAS_S)
  cgesvd_interface(zgesvd_, BLAS_Z, BLAS_D)

  template <typename MAT, typename VEC>
  void svd(const MAT &X_, MAT &U, MAT &Vtransposed, VEC &sigma) {
   MAT X(X_);
   svd(X, U, Vtransposed, sigma);
  }

}

#else

namespace gmm
{
template <typename MAT>
void schur(const MAT &, MAT &, MAT &)
{
  GMM_ASSERT1(false, "Use of function schur(A,S,Q) requires GetFEM "
                     "to be built with Lapack");
}

template <typename BLAS_TYPE>
inline void svd(dense_matrix<BLAS_TYPE> &, dense_matrix<BLAS_TYPE> &,
         dense_matrix<BLAS_TYPE> &, std::vector<BLAS_TYPE> &)
{
  GMM_ASSERT1(false, "Use of function svd(X,U,Vtransposed,sigma) requires GetFEM "
                     "to be built with Lapack");
}

}// namespace gmm

#endif // GMM_USES_LAPACK

#endif // GMM_LAPACK_INTERFACE_H
