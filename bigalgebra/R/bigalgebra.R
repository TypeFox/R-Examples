
is_transposed = function( tcode )
{
  if ( sum(tcode == c('n', 'N')) > 0 )
    return(FALSE)
  if ( sum(tcode == c('T', 't', 'C', 'c')) > 0 )
    return(TRUE)
  stop("Invalid transpose code given") 
}

check_matrix = function(A, classes=c('big.matrix', 'matrix'), 
  types='double')
{
  if (!any( class(A) == classes))
  {
    stop("A is not the correct class type")
  }
  if (!any(typeof(A) == types))
  {
    stop("The matrix type is not correct")
  }
  return( ifelse( class(A) == 'big.matrix', TRUE, FALSE ) )
}

# Do I need a function to add a scalar to each element of a matrix?

# Copy a matrix
# Y := X
dcopy = function(N=NULL, X, INCX=1, Y, INCY=1)
{
  X.is.bm = check_matrix(X)
  Y.is.bm = check_matrix(Y)
  if (is.null(N))
  {
    N = as.double(nrow(X))*as.double(ncol(X))
  }
  .Call('dcopy_wrapper', N, X, as.double(INCX), Y, as.double(INCY),
    X.is.bm, Y.is.bm)
  return(0)
}

# Multiply by a scalar
# Y := ALPHA * Y
dscal = function(N=NULL, ALPHA, Y, INCY=1)
{
  Y.is.bm = check_matrix(Y)
  if (is.null(N))
  {
    N = as.double(nrow(Y))*as.double(ncol(Y))
  } 
  .Call('dscal_wrapper', N, as.double(ALPHA), Y, as.double(INCY), Y.is.bm)
  return(0)
}

# Add two matrices.
# Y := ALPHA * X + Y
daxpy = function(N=NULL, ALPHA=1, X, INCX=1, Y, INCY=1)
{
  X.is.bm = check_matrix(X)
  Y.is.bm = check_matrix(Y)
  if (is.null(N))
  {
    N = as.double(nrow(X))*as.double(ncol(X))
  }
  .Call('daxpy_wrapper', as.double(N), as.double(ALPHA), X, as.double(INCX),
    Y, as.double(INCY), X.is.bm, Y.is.bm)
  return(0)
}

# Matrix Multiply
# C := ALPHA * op(A) * op(B) + BETA * C
# This is function provides dgemm functionality.  The matrix types
# are handled on the C++ side.
dgemm = function(TRANSA='n', TRANSB='n', M=NULL, N=NULL, K=NULL,
  ALPHA=1, A, LDA=NULL, B, LDB=NULL, BETA=0, C, LDC=NULL, COFF=0) 
{
  A.is.bm = check_matrix(A)
  B.is.bm = check_matrix(B)
  C.is.bm = check_matrix(C)
  # The matrices look OK.  Now, if they haven't been specified, let's
  # specify some reasonable dimension information.
  if ( is.null(M) )
  {
    M = ifelse ( is_transposed(TRANSA), ncol(A), nrow(A) )
  }
  if ( is.null(N) ) 
  {
    N = ifelse ( is_transposed(TRANSB), nrow(B), ncol(B) )
  }
  if ( is.null(K) )
  {
    K = ifelse ( is_transposed(TRANSA), nrow(A), ncol(A) )
  }
  if ( is.null(LDA) ) LDA = nrow (A)
  if ( is.null(LDB) ) LDB = nrow (B)
  if ( is.null(LDC) ) LDC = nrow (C)

  .Call('dgemm_wrapper', as.character(TRANSA), as.character(TRANSB),
    as.double(M), as.double(N), as.double(K), as.double(ALPHA), A, 
    as.double(LDA), B, as.double(LDB),
    as.double(BETA), C, as.double(LDC), as.logical(A.is.bm), 
    as.logical(B.is.bm), as.logical(C.is.bm), COFF)
  return(invisible(C))
}

# QR factorization
# return 0 if successful, -i if ith argument has illegal value
dgeqrf=function(M=NULL, N=NULL, A, LDA=NULL, TAU=NULL, WORK=NULL,
  LWORK=NULL)
{
  A.is.bm = check_matrix(A)
  if (is.null(M))
  {
    M = nrow(A)
  }
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  if (is.null(TAU))
  {
    TAU = as.matrix(rep(0.0, min(M,N)))
  }
  if (is.null(LWORK))
  {
    LWORK = max(1, N)
  }
  if (is.null(WORK))
  { 
    WORK = as.matrix(rep(0.0, max(1, LWORK)))
  }
  TAU.is.bm = check_matrix(TAU)
  WORK.is.bm = check_matrix(WORK)
  INFO = 0
  .Call('dgeqrf_wrapper', as.double(M), as.double(N), A, as.double(LDA), 
    TAU, WORK, as.double(LWORK), as.double(INFO), A.is.bm, TAU.is.bm, 
    WORK.is.bm)
  return(INFO) 
}

# Cholesky factorization
# return 0 if successful, <0 if -i-th argument is invalid, > 0 if leading minor
# is not positive definite
dpotrf=function(UPLO='U', N=NULL, A, LDA=NULL)
{
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  A.is.bm = check_matrix(A)
  INFO = 0
  .Call('dpotrf_wrapper', as.character(UPLO), as.double(N), A, as.double(LDA),
    as.double(INFO), A.is.bm)
  return(INFO)
}

# General eigenvalue
# return 0 if successful, <0 i-th argument has illegal value, >0 QR 
# algorithm failed.
# for now, VL and VR have to be matrices but they could be NULL
dgeev=function(JOBVL='V', JOBVR='V', N=NULL, A, LDA=NULL, WR, WI, VL, 
  LDVL=NULL,  VR=NULL, LDVR=NULL, WORK=NULL, LWORK=NULL)
{
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  if (is.null(LDVL) && (JOBVL=='V'))
  {
    LDVL = N
  }
  if (is.null(LDVR) && (JOBVR=='V'))
  {
    LDVR = N
  }
  if (is.null(LWORK))
  {
    LWORK = ifelse( (JOBVL=='V' || JOBVR == 'V'), 4*N, max(1, 3*N) )
  }
  if (is.null(WORK))
  {
    WORK = as.matrix(rep(0.0, max(1, LWORK)))
  }
  # Take car of the case where someone doesn't want to get the 
  # eigen vectors and passed NULL.
  if (is.null(VL))
  {
    VL = matrix(0.0, nrow=1, ncol=1)
  }
  if (is.null(VR))
  {
    VR = matrix(0.0, nrow=1, ncol=1)
  }
  INFO = 0
  A.is.bm = check_matrix(A)
  WR.is.bm = check_matrix(WR)
  WI.is.bm = check_matrix(WI)
  VL.is.bm = check_matrix(VL)
  VR.is.bm = check_matrix(VR)
  WORK.is.bm = check_matrix(WORK)
  INFO=0
  .Call('dgeev_wrapper', as.character(JOBVL), as.character(JOBVR), as.double(N),    A, as.double(LDA), WR, WI, VL, as.double(LDVL), VR, as.double(LDVR),
    WORK, as.double(LWORK), as.double(INFO), A.is.bm, WR.is.bm, WI.is.bm, 
    VL.is.bm, VR.is.bm, WORK.is.bm)
  return(INFO)
}

# Returns: = 0 if successful
#          < 0 if INFO = -i had an illegal value
#          > 0 if DBDSDC did not converge
dgesdd = function( JOBZ='A', M=NULL, N=NULL, A, LDA=NULL, S, U, LDU=NULL, 
  VT, LDVT=NULL, WORK=NULL, LWORK=NULL)
{
  A.is.bm = check_matrix(A)
  S.is.bm = check_matrix(S)
  U.is.bm = check_matrix(U)
  VT.is.bm = check_matrix(VT)
  if (is.null(M))
  {
    M=nrow(A)
  }
  if (is.null(N))
  {
    N=ncol(A)
  }
  if (is.null(LDA))
  {
    LDA=nrow(A)
  }
  if (is.null(LDU))
  {
    LDU=nrow(U)
  }
  if (is.null(LDVT))
  {
    LDVT=nrow(VT)
  }
  if (is.null(LWORK) && is.null(WORK))
  {
    if (JOBZ=='N')
    {
      LWORK = 3 * min(M, N) + max( max(M, N), 7*min(M, N) )
    }
    else if (JOBZ=='O')
    {
      LWORK = 3 * (min(M, N))^2 + max( max(M, N), 5 * (min(M, N))^2
        + 4 * min(M, N) )
    }
    else if (JOBZ == 'S' || JOBZ == 'A')
    {
      LWORK = 3 * (min(M, N))^2 + max( max(M, N), 4 * (min(M, N))^2
        + 4 * min(M, N) )
    }
    else
    {
      stop("Invalid JOBZ argument specified")
    }
    WORK = as.matrix(rep(0.0, max(1, LWORK) ))
  }
  if (is.null(LWORK))
  {
    LWORK = length(WORK)
  }
  if (is.null(WORK))
  {
    WORK = as.matrix(rep(0.0, max(1, LWORK)))
  }
  WORK.is.bm = check_matrix(WORK)
  INFO = 0
  .Call('dgesdd_wrapper', as.character(JOBZ), as.double(M), as.double(N), A, 
    as.double(LDA), S, U, as.double(LDU), VT, as.double(LDVT), WORK, 
    as.double(LWORK), as.double(INFO), A.is.bm, S.is.bm, U.is.bm, VT.is.bm, 
    WORK.is.bm)
  return(INFO)
}
