
#::::::::::::
#   dqrslm
#::::::::::::

subroutine  dqrslm (x, ldx, n, k, qraux, a, lda, job, info, work)

#  Acronym:  `dqrsl' Matrix version

#  Purpose:  This routine generates the matrix Q^{T}AQ or QAQ^{T}, where 
#      Q is the products of Householder matrix stored in factored form in 
#      the LOWER triangle of `x' and `qraux', and A is assumed to be 
#      symmetric.  This routine is designed to be compatible with LINPACK's 
#      `dqrdc' subroutine.
 
#  References:  1. Dongarra et al. (1979) LINPACK Users' Guide. (chap. 9)
#               2. Golud and Van Loan (1983) Matrix Computation. (pp.276-7)
 
integer           ldx, n, k, lda, job, info
double precision  x(ldx,*), qraux(*), a(lda,*), work(*)

#  On entry:
#      x          output from `dqrdc', of size (ldx,k).
#      ldx        leading dimension of x.
#      n          size of matrix A and Q.
#      k          number of factors in Q.
#      qraux      output from `dqrdc'.
#      a          matrix A (of size (lda,n)), only LOWER triangle refered.
#      lda        leading dimension of a.
#      job        0:  Q^{T} A Q.
#                 1:  Q A Q^{T}.
 
#  On Exit:
#      a          matrix Q^{T}AQ or QAQ^{T} in LOWER triangle.
#      info        0:  normal termination.
#                  1:  `job' is out of scope.
#                 -1:  dimension error.
#      others     unchanged.

#  Work array:
#      work       of size at least (n).

#  Routines called:
#      Blas    -- ddot, daxpy
#      Blas2   -- dsymv, dsyr2

#  Written:  Chong Gu, Statistics, UW-Madison, latest version 8/29/88.

double precision  tmp, alph, ddot
integer           i, j, step

info = 0

#   check input
if ( lda < n | n < k | k < 1 ) {
    info = -1
    return
}
if ( job != 0 & job != 1 ) {
    info = 1
    return
}

#   set operation sequence
if ( job == 0 ) {
    j = 1
    step = 1
}
else {
    j = k
    step = -1
}

#   main process
while ( j >= 1 & j <= k ) {
    if ( qraux(j) == 0.0d0 ) {
        j = j + step
        next
    }

    tmp = x(j,j)
    x(j,j) = qraux(j)

    #   update the columns 1 thru j-1
    for (i=1;i<j;i=i+1) {
        alph = - ddot (n-j+1, x(j,j), 1, a(j,i), 1) / x(j,j)
        call  daxpy (n-j+1, alph, x(j,j), 1, a(j,i), 1)
    }

    #   update the submatrix at bottom-right corner

    #   compute  p = D v / v_{1} 
    alph = 1.d0 / x(j,j)
    call  dsymv ('l', n-j+1, alph, a(j,j), lda, x(j,j), 1, 0.d0, work(j), 1)

    #   compute  w = p - ( p^{T} v / 2 v_{1} ) v 
    alph = - ddot (n-j+1, work(j), 1, x(j,j), 1) / 2.d0 / x(j,j)
    call  daxpy (n-j+1, alph, x(j,j), 1, work(j), 1)

    #   compute  D = D - v w^{T} - w v^{T} 
    call  dsyr2 ('l', n-j+1, -1.d0, x(j,j), 1, work(j), 1, a(j,j), lda)

    x(j,j) = tmp
    j = j + step
}

return
end

#...............................................................................

