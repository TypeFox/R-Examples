# XXX TODO
#
# Each routine below may return either a big matrix or an R matrix
# depending on how they are called (or a set of such matrices for
# the decomposition methods). In each case the returned value is
# allocated by the routine.
#
# When big matrices are returned, we need to register a finalizer
# with the big matrix that removes its allocated files when the
# garbage collector is called on the big matrix, to clean up.
# Add this!
#


# Matrix Multiply
# C := ALPHA * op(A) * op(B) + BETA * C
# This is function provides dgemm functionality.
dgemm = function(TRANSA='N', TRANSB='N', M=NULL, N=NULL, K=NULL,
  ALPHA=1, A, LDA=NULL, B, LDB=NULL, BETA=0, C, LDC=NULL, COFF=0) 
{
  A.is.bm = check_matrix(A)
  B.is.bm = check_matrix(B)
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
  if ( is.null(LDA) ) LDA = ifelse (is_transposed(TRANSA), K, M)
  if ( is.null(LDB) ) LDB = ifelse (is_transposed(TRANSB), N, K) 
  if ( is.null(LDC) ) LDC = M

# Default to big matrix output
  if(missing(C)) C = anon_matrix(M, N)
  C.is.bm = "big.matrix" %in% class(C)

  .Call('dgemm_wrapper', as.character(TRANSA), as.character(TRANSB),
    as.double(M), as.double(N), as.double(K), as.double(ALPHA), A, 
    as.double(LDA), B, as.double(LDB),
    as.double(BETA), C, as.double(LDC), as.logical(A.is.bm), 
    as.logical(B.is.bm), as.logical(C.is.bm), COFF, PACKAGE="bigalgebra")
}

# Vector addition and scaling
# Y := A * X  + Y
# A is a scalar double
# X is either a big.matrix or regular R matrix.
# Y is an optional matrix or vector of the same length as X.
# Returns a new matrix or big matrx with the same dimensions as X. If
# X is a dimension-less R vector, returns a column. Returned value type
# depends on the arguments and the value of the option
# options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]].
daxpy = function(A=1, X, Y)
{
  mixed = FALSE
  X.is.bm = check_matrix(X,classes=c('big.matrix','matrix','vector','numeric'))
# default to a column big matrix output
  M = length(X)
  L = M
  N = 1L
  D = dim(X)
  if(!is.null(D) && length(D)==2)
  {
    M = D[1]
    N = D[2]
  }
  Z = anon_matrix(M,N,val=0.0)
  if(!missing(Y))
  {
# Check conformity of Y and duplicate
    if(length(Y)!=length(X)) stop("Lengths of X and Y must match")
    mixed = (X.is.bm != check_matrix(Y,classes=c('big.matrix','matrix','vector','numeric')))
    Z[] = Y[]
  }
  ans = .Call("daxpy_wrapper", as.double(L), as.double(A), X, Z, X.is.bm,
              PACKAGE="bigalgebra")
  if(mixed) return(ans[])
  ans
}
