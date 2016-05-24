# SimplicialCubature - an R package to integrate over m-dimensional simplices in n-dimensions.
#   Is is required that n >= 2, while the simplices can be of dimension m, 1 <= m <= n.
#
# Written by John Nolan,  American University, Aug-Sep 2014. Modified in Spring 2015. 
#     e-mail jpnolan@american.edu
# 
# This research was supported by an agreement with Cornell University, Operations 
#    Research & Information Engineering, under contract W911NF-12-1-0385 from the Army 
#    Research Development and Engineering Command.
#
# There are two main functions in this package:
# 1. Function adaptIntegrateSimplex
#    This function adaptively integrates a function f(x) over a simplex in 
#    n-dimensional space, n >= 2.  
#    (a) If the the simplicies are n-dimensional, then the integration 
#    is performed by function adsimp(...), an R 
#    translation of Alan Genz's matlab function adsimp, which
#    adaptively integrates an arbitrary function over an n-dim. simplex.
#    See comments in function adsimp for references on the method.  The only
#    change to that method is that a new argument partitionInfo is added;
#    if that is TRUE, adsimp will return information about
#    the subdivision of the simplices.
#    (b) If the simplex is of dimension 1 <= m < n, then the problem is transformed
#    into an integral over m-dimensional space, and evaluated using adsimp(...)
#
#    In both cases, the adaptIntegrateSimplex function allows the integrand function
#    to have optional arguments by using "..." ; this is achieved by defining
#    a new integrand function locally that incorporates the extra arguments,
#    leaving adsimp unaware of that construct.
#
# 2. Function integrateSimplexPolynomial
#    R implementations of the Lasserre-Avrachenkov and Grundmann-Moller
#    methods for exactly integrating a polynomial over a simplex.
#    The polynomial p(x) must be defined by using function polyDefine.
#    It is specified as a list with fields:
#       p$coef[1:m] vector of doubles, the coefficients of the terms
#       p$k[1:m,1:n] matrix of integers
#       m and n are given implicitly in the sizes of these arrays
#    Precisely,  p(x[1],...,x[n])=sum (coef[i] * x[1]^k[i,1] * ... * x[n]^k[i,n]),
#       where the sum is over i=1,...,m.  
#
#    The functions polyDefine, polyEval, and polyPrint are used to define,
#    evaluate, and print a polynomial in coded form.  
#
#    The functions CanonicalSimplex, UnitSimplexV, SimplexVolume, and 
#    JacobianS2Canonical are utility routines to define the canonical
#    simplex, the unit simples, volume of an n-dim. simplex in R^n,
#    and the Jacobian of the transformation from an m-dim. simplex S2
#    to the m-dim. canonical simplex.
#
#    All other functions are internal and are used at your own risk.
#
################################################################################
adaptIntegrateSimplex <- function( f, S, fDim=1L, maxEvals=10000L, absError=0.0, 
    tol=1.0e-5, integRule=3L, partitionInfo=FALSE, ... ) {
# Adaptively integrate f(x) over one or more simplices in R^n, where 
# S[1:n,1:m,1:nS] may be a single simplex or a list of simplices. 
# The region of integration can be an m-dimensional simplex, 1 <= m <= n.
#
# When m=n, this is a wrapper function for Alan Genz's adsimp program to adaptively 
# integrate a function f(x) over an n-dimensional simplex S.  f(x) can be vector valued
# and S can be an array of simplices.  We change variable names to parallel
# the variable names in CRAN packages cubature and SphericalCubature.
# When m < n, the problem is transformed to integration in m-space and then
# adsimp is used to evaluate that m-dim. integral.
#
# Input variables:
#    f integrand function, returns a vector of length fDim
#    S either a single n x (m+1) matrix, with S[,i] being the i-th vertex 
#      of the simplex, or an n x (m+1) x nS array of simplices, with S[,,i] 
#      being the i-th simplex in the first form
#    fDim integer that gives the dimension of (the range of) function f
#    maxEvals integer that gives an upper limit for how many evaluations of
#      of f are allowed.  If this value is exceeded (or about to be exceeded),
#      the return variable will contain a non-zero value.
#    absError requested absolute error in the integral
#    tol requested relative error in the integral
#    integRule integer in the range 1:4 specifiying the integration rule used
#      (see function adsimp for more info)
#    partitionInfo TRUE/FALSE value indicating whether partition information
#      is be be returned.  This takes more memory (and time if m < n), but 
#      the results may be of use in some applications.  (see function adsimp for more info)
#    ...  optional arguments to f
#
# Returns a list with multiple fields:
#    integral estimated values of integrals, vector if fDim > 1
#    estAbsError estimated absolute error
#    functionEvaluations count of number of function evaluations
#    returnCode integer status: returnCode=0 is a successful return
#    message a text string describing what returnCode means
#    if partitionInfo=TRUE, other information about the final partition/subdivision 
#      of the region of integration is returned, see function adsimp for more details.

if( !is.function(f) ) { stop("first argument is not a function" ) }
if( is.matrix(S) )  { S <- array( S, dim=c(nrow(S),ncol(S),1)) }
if( !is.array(S) | (length(dim(S)) != 3) ) stop("S must be a single simplex or an array of simplices")
dimS <- dim(S); n <- dimS[1]; m <- dimS[2]-1; nS <- dimS[3]
if (n < 2) stop("dimension of the space is less than 2")
if( m > n ) stop("invalid simplices: m > n")

if( n == m ) { # region is a set of n-dim simplices inside R^n
   # define a new function to evaluate f using optional arguments in ...
   fNew <- function( x ) { f( x, ... ) }
   
   # call adsimp to do integration
   a <- adsimp( ND=n, VRTS=S, NF=fDim, F=fNew, MXFS=maxEvals, EA=absError, 
               ER=tol, KEY=integRule, partitionInfo )
} else {  
  # integration region is a set of m-dim. simplices inside R^n with m < n.
  # transform the problem to an m-dimensional one: map simplex S[,,i]
  # to a translation of the canonical simplex by (2*i,0,0,...,0).  
  # This translation allows the function transformedF to know 
  # which region of integration a point is in simply by looking at its first
  # coordinate.  This allows us to incorporate the change of variable
  # coefficient, JacobianCoef[i] for simplex S[,,i], into the transformed
  # integrand function.
  S0 <- CanonicalSimplex( m )
  transformedS <- array( 0.0, dim=c(m,m+1,nS) )
  JacobianCoef <- rep(0.0,nS)
  for (i in 1:nS) {
    JacobianCoef[i] <- JacobianS2Canonical( S[,,i] )
    S0[1, ] <- S0[1, ] + 2
    transformedS[,,i] <- S0
  }

  # integrate transformed f over transformed simplices
  if (m==1) {   # evaluate a line integral
    # define the transformed function which evaluates the original function f(x) at 
    # multiple points in R^n, returning multiple values.
    transformedF <- function( u ) { 
      # in 1-dim case, u is a vector of distinct points in R.
      # since the function integrate.vector.fn does one interval at at
      # time, all points in u must be in the same original parition interval.
      nu <- length(u)  
      y <- matrix(0.0,nrow=fDim,ncol=nu)
      for (i in 1:nu) {
        b <- original.coordinates( u[i], S )   
        y[,i] <- JacobianCoef[b$i] * f( b$x, ... )
      }
      return(y) }
            
    a <- integrate.vector.fn( transformedS, fDim, transformedF, maxEvals, absError, tol, partitionInfo )
    
  } else {  # evaluate m-dim. integral where 2 <= m < n
    # define the transformed function which evaluates the original function f(x) at a 
    # a single point in R^n, returning a single number
    transformedF <- function( u ) { 
      b <- original.coordinates( u, S )
      y <- JacobianCoef[b$i] * f( b$x, ... )
      return(y) }
            
    a <- adsimp( ND=m, VRTS=transformedS, NF=fDim, F=transformedF, MXFS=maxEvals, EA=absError, 
               ER=tol, KEY=integRule, partitionInfo )
  }
}

if (m==1) { # in univariate case, just return result of integrate.vector.fn 
  result <- a
} else { # rename results and return to caller
  if (partitionInfo) {
    result <-  list(integral=a$VL, estAbsError=a$AE, functionEvaluations=a$NV, 
            returnCode=a$FL, subsimplices=a$VRTS, subsimplicesIntegral=a$VLS, 
            subsimplicesAbsError=a$AES, subsimplicesVolume=a$VOL )
  } else {
    result <-  list(integral=a$VL, estAbsError=a$AE, functionEvaluations=a$NV, 
            returnCode=a$FL )
  }  
}

if(result$returnCode == 0) {
  if (partitionInfo & ( m < n )) {
    # the current partition info applies to transformed simplices,
    # translate those simplices back to the original coordinates. 
    n3 <- dim(result$subsimplices)[3]
    newSubsimplices <- array( 0.0, dim=c(n,m+1,n3) )
    for (i in 1:n3) {
      for (j in 1:(m+1)) {
        b <- original.coordinates( result$subsimplices[,j,i], S )
        newSubsimplices[,j,i] <- b$x
      }
      result$subsimplicesVolume[b$i] <- JacobianCoef[b$i]*result$subsimplicesVolume[b$i]
    }
    result$subsimplices <- newSubsimplices        
  }    
}

# convert returnCode to a message string
result$message <- adsimp.return.message( result$returnCode )  

return(result) }
################################################################################
adsimp.return.message <- function( rcode ) {
# return a text string explaining the return code recode from adsimp

if (rcode == 0) msg <- "OK"
if (rcode == 1) msg <- "error: maxEvals exceeded - too many function evaluations"
if (rcode == 2) msg <- "error: integRule < 0 or integRule > 4"
if (rcode == 3) msg <- "error: dimension n of the space < 2"
if (rcode == 4) msg <- "error: fDim = dimension of f < 1"
if (rcode == 5) msg <- "error: absError < 0 and tol < 0"
if ((rcode < 0) | (rcode > 5)) msg <- paste("error: unknown return code = ",rcode ) 

return(msg) }
################################################################################
integrate.vector.fn <- function( intervals, fDim, f, maxEvals, absError, tol, partitionInfo=FALSE ){
# Evaluate an integral of a vector valued function f(x)=(f[1](x),...,f[fDim](x)) of a real 
# variable x (not a vector) over a list of intervals.  The built-in R function integrate 
# uses quadpack function dqags when lower and upper are both finite, and that appears
# to evaluate f at 21 points each time. 
dqags.num.points <- 21L

if( !is.function(f) ) { stop("first argument is not a function" ) }

# define function to return a single coordinate of vector valued f(x)
# note that x will generally be a vector of points in this univariate example
newF <- function( x, which.coord ) { 
  if(length(x)!=dqags.num.points) { warning(paste("integrate.vector.fn: length(x) != ",dqags.num.points)) }  # debug
  f(x)[which.coord,] }

n.intervals <- dim(intervals)[3]
max.subdiv <- as.integer( maxEvals/dqags.num.points )
abs.tol <- absError/n.intervals
y <- rep(0.0,fDim)
if (partitionInfo) { 
  y.partition <- matrix( 0.0, nrow=fDim, ncol=n.intervals) 
  est.abs.error.partition <- y.partition
}
est.abs.error <- rep(0.0,fDim)
actual.subdiv <- rep( 0L, fDim )

# loop over coordinates of vector valued f
for (i in 1:fDim) { 
  # loop over intervals 
  for (j in 1:n.intervals) {
    tmp <- integrate( newF, lower=intervals[1,1,j], upper=intervals[1,2,j], rel.tol=tol, 
         abs.tol=abs.tol, subdivisions=max.subdiv, which.coord=i )
    if (tmp$message != "OK") { break }
    y[i] <- y[i] + tmp$value
    est.abs.error[i] <- est.abs.error[i] + tmp$abs.error
    actual.subdiv[i] <- actual.subdiv[i] + tmp$subdivisions
    if (partitionInfo) { 
      y.partition[i,j] <- y.partition[i,j] + tmp$value 
      est.abs.error.partition[i,j] <- est.abs.error.partition[i,j]+tmp$abs.error
    }    
  }
}

result <-  list(integral=y, estAbsError=est.abs.error, 
                functionEvaluations=dqags.num.points*sum(actual.subdiv), 
                returnCode=ifelse(tmp$message=="OK",0L,-1L), 
                message=tmp$message )
                
if(partitionInfo) {
  result$subsimplices <- intervals
  result$subsimplicesIntegral <- y.partition
  result$subsimplicesAbsError <- est.abs.error.partition
  result$subsimplicesVolume <- intervals[1,2,]-intervals[1,1,]
}
return(result)}
################################################################################
original.coordinates <- function( u, S ) {
# when m < n, this function is used to transform from the m-dimensional 
# new coordinates to original n-dimensional coordinates.  
# u should be a single point in R^m, e.g. a vector of length m.
# The basic idea is barycentric coordinates, after correcting for the 
# deliberate shift of the i-th simplex by (2*i,0,0,...,0).
#   x <- u[1]*S[,1,i] + u[2]*S[,2,i]+...+u[m]*S[,m,i]+ (1-sum(u))*S[,m+1,i]
#
# This function returns a list with two fields:
#   $i  for the index of the simplex
#   $x  for the point in R^n corresponding to u in R^m

i <- as.integer(u[1]/2.0) # figures out which simplex u is in
u[1] <- u[1]- 2*i   # subtract the shift to get barycentric coordinates
v <- c( u, 1-sum(u) )
return( list( i=i, x= S[,,i] %*% v ) ) }
################################################################################
integrateSimplexPolynomial <- function( p, S, method="GM" ) {
#
# Exactly integrate a polynomial p(x) over a simplex.
# 
# Input variables:
#   p a polynomial, which should be defined through function polyDefine.
#   S is either an n x (m+1) matrix with S[,j] the j-th vertex of the simplex.  
#      Alternatively, S can be an (n x (m+1) x k) array, with S[,,i] giving the 
#      i-th simplex as in the matrix case.
#   method is a string, one of:
#      "GM" for Grundmann-Moller method, works for 1 < m <= n.
#      "LA" for the Lasserre-Avrachenkov method, only works for m=n.
#
# The formulas used here are exact for polynomials, not approximate, but 
# use floating point arithmetic to evaluate the polynomial.  So there will be
# round-off due to this floating point evaluation.  If the coefficients of the
# polynomial are all rational numbers, and the vertices are all rational numbers,
# it is possible to evaluate the polynomials exactly as rational numbers, and
# the integral can be evaluated exactly.  We do not do this here, though it should
# be possible using CRAN packages gmp or Rmpfr.  Exact rational arithmetic is done 
# in the C program LattE, which is available from the website:
#     https://www.math.ucdavis.edu/~latte
#
# If method="GM", then the Grundmann-Moler algorithm of order ceiling( (degree(p)-1) / 2 ) 
# is used.  The GM method is exact for this order on a polynomial. See
# function grnmol( ) for details of the method and references.  
#
# If method="LA", this function splits the polynomial into a list of q-homogeneous 
# polynomials and calls function LasserreAvrachenkov( ) to exactly integrate a 
# q-homogeneous polynomial over each simplex S[,,i].  

# if S is a single simplex, make it an array of simplices of length 1
if (is.matrix(S)) { S <- array( S, dim=c(dim(S),1) ) }

# check input parameters
if( !is.array(S) ) stop("S must be a single simplex (a matrix) or a list of matrices (a 3-dim. array)")
# check input parameters
dimS <- dim(S)
if ( length(dimS) != 3) stop( "S must be a single simplex (a matrix) or a list of matrices (a 3-dim. array)" )
stopifnot( is.vector(p$coef), is.matrix(p$k), dim(p$k)[2]==dimS[1], method %in% c("GM","LA") )
n <- dimS[1]; m <- dimS[2]-1; k <- dimS[3]
if( m > n) stop("m must be <= n")
if( n < 1 ) stop("n=nrow(S) must be bigger >= 1")
if( m < 1 ) stop("m=ncol(S)-1 must be bigger than 1")

# compute the degrees of each of the terms in p(.)
deg <- rowSums(p$k)

# Grundmann-Moller algorithm
if (method=="GM") {
  gm.order <- ceiling( (max(deg)-1) / 2 ) # order of Grundmann-Moller method
  if (gm.order < 1) { gm.order <- 1 }
  integral <- 0; functionEvaluations <- 0
 
  if (m == n) {
    # polynomial interface to grnmol
    fnew <- function( x ) { evalPoly( x, p ) }

    for (i in 1:k) {   
      a <- grnmol( fnew, S[ , , i], gm.order ) 
      integral <- integral + a$Q[gm.order+1]
      functionEvaluations <- functionEvaluations + a$nv
    }
  } else {   # 1 <= m < n, map to an integral over m-dim. canonical simplex
  
    # integration region is a set of m-dim. simplices inside R^n with m < n.
    # transform the problem to an m-dimensional one: map simplex S[,,i]
    # to a translation of the canonical simplex by (2*i,0,0,...,0).  
    # This translation allows the function transformedF to know 
    # which region of integration a point is in simply by looking at its first
    # coordinate.  This allows us to incorporate the change of variable
    # coefficient, JacobianCoef[i] for simplex S[,,i], into the transformed
    # integrand function.
    S0 <- CanonicalSimplex( m )
    transformedS <- array( 0.0, dim=c(m,m+1,k) )
    JacobianCoef <- rep(0.0,k)
    for (i in 1:k) {
      JacobianCoef[i] <- JacobianS2Canonical( S[,,i] )
      S0[1, ] <- S0[1, ] + 2
      transformedS[,,i] <- S0
    }    

    # define the transformed function which evaluates the original function f(x) at a 
    # a single point in R^n, returning a single number
    transformedF <- function( u ) { 
      b <- original.coordinates( u, S )
      y <- JacobianCoef[b$i] * evalPoly( b$x, p )
      return(y) }
    
    for (i in 1:k) {   
      curS <- transformedS[ , , i]
      if (!is.matrix(curS)) { curS <- matrix( curS, nrow=1 ) }
      a <- grnmol( transformedF, curS, gm.order ) 
      integral <- integral + a$Q[gm.order+1]
      functionEvaluations <- functionEvaluations + a$nv
    }
  }
}

# Lasserre-Avrachenkov method
if (method=="LA") {
  # in principle, L-A method can be made to work when m < n, but it is
  # messy.  Since it requires that you separate the terms of the polynomial
  # into homogeneous terms, you have to express the transformed polynomial
  # in m-dimensions in polynomial form.  This requires some effort; we do
  # not see it worth the effort now.
  if (m != n) stop("Lasserre-Avranchenkov method only works with m=n")
  if (n == 1) stop("Lasserre-Avranchenkov method only works for n > 1")
  unique.deg <- unique(deg)

  # Loop through terms that are homog. of each degree q and all simplices.
  # Use a trick of setting some coef. to 0 when the corresponding term is 
  # not of degree q.
  integral <- 0.0; functionEvaluations <- 0
  for ( q in unique.deg ) {
    useTerm <- ( deg == q ) # useTerm[i] = 1 if term i is q-homog., 0 otherwise
    for (i in 1:k) { 
      tmp <- LasserreAvrachenkov( q, p, useTerm, S[ , , i] ) 
      integral <- integral + tmp$integral
      functionEvaluations <- functionEvaluations + tmp$functionEvaluations  
    }
  }
}

return( list(integral=integral,functionEvaluations=functionEvaluations )) }
################################################################################
definePoly <- function( coef, k ) {
# Define a polynomial in n variables with coefficients coef[1:m] and powers k[1:m,1:n]
# this is essentially an S3 method

# check sizes
stopifnot( is.vector(coef), is.numeric(coef), is.matrix(k), length(coef) == nrow(k) )

# as.integer drops the matrix dimensions
return( list( coef=coef, k=matrix( as.integer(k), nrow=nrow(k) ) ) ) }  
################################################################################
evalPoly <- function( x, p, useTerm=rep(TRUE,length(p$coef)) ) {
# Evaluate the polynomial y[j] = p(x[,j]) at j=1,...,nx
#
# The vector useTerm indicates which terms in polynomial p are
# used in this evaluation: if useTerm[i]=TRUE, include the i-th term
# in the function evaluation.  This is mainly for efficiency reasons:
# function integrateSimplexPolynomial splits the polynomial into terms of
# different degrees, and this vector of indicators allows us to use only
# the terms required.  This avoids repeated copying of terms, and will
# be simpler if these functions are coded in C.

if (is.vector(x)) { x <- matrix(x,ncol=1) }
stopifnot( nrow(x) == ncol(p$k) )
nx <- ncol(x);  m <- length(p$coef)

y <- rep(0.0,nx)
for (i in 1:m) {
  if (useTerm[i]) {   # skip terms when told to
    for (j in 1:nx) {
      y[j] <- y[j] + p$coef[i]*prod( x[,j]^p$k[i,] )
    }
  }
}
return(y) }
################################################################################
printPoly <- function( p, num.digits=5 ) {
# print the polynomial p in human readable form

# round coefficients for readability
coef <- round(p$coef,num.digits)

m <- length(coef)
k <- p$k
dimk <- dim(k)
if (m != dimk[1]) stop( "m != dimk[1]; not  valid polynomial" )
n <- dimk[2]
cat("Polynomial has",m,"terms and",n,"variables\n")
deg <- rowSums(k)
for (i in 1:m) {
  if (i > 1) { if (coef[i] > 0) cat ("\n  + ") else cat("\n  ") }
  cat( coef[i] )
  for (j in 1:n) {
    if (k[i,j] > 0 ) cat( " * x[",j,"]^", k[i,j], sep="" )
  }
  cat("       (degree=", deg[i],")" )  
}  
cat("\n")}
################################################################################
LasserreAvrachenkov <- function( q, p, useTerm, S ){
# Exactly integrate a homogeneous polynomial p over a single simplex S using the
# algorithm in Lasserre and Avrachenkov (2001), MAA Monthly, v 108, n 2, pg 151-154.
# We use the formula (18) of Baldoni, et. al. (2011), Math. of Comput., v 80, pg 297-325.
# We also incorporate the fact that you can eliminate half the terms in their inner sum and
# multiply the result by 2.  This is done below by restricting eps[1]=+1, or equivalently
# j[1]=0.
# The polynomial p(.) should be defined using function definePoly.
# Note that p may not be homogeneous, however argument useTerm is a boolean
# vector telling which terms of p to use.  It is assumed that all the terms
# in p with useTerm[i]=TRUE are homogeneous of order q.  Since this function
# is only called internally from function integrateSimplexPolynomial, we do 
# not check this here.

n <- nrow(S);   m <- ncol(S)-1
vol <- SimplexVolume(S) 

if (q == 0) { # constant terms are treated separately
  const <- vol
  total <- sum( p$coef[useTerm] )
  count <- 1
} else { # 
  const <- vol/ ( 2^(q-1) * factorial(q) * choose(q+n,n) )
  count <- 0
  total <- 0.0
  i <- rep( 1L, q )
  while ( !is.na( i[1] ) ) {
    j <- rep( 0L, q )
    while ( j[1]==0L ) {
      eps <- 1 - 2*j   # convert to +1s and -1s 
      tmp <- eps[1]*S[,i[1]]
      if( q > 1 ) {
        for (l in 2:q) { 
          tmp <- tmp + eps[l]*S[ ,i[l]]
        }
      }   
      val <- evalPoly( tmp, p, useTerm )
      total <- total + prod(eps) * val
      count <- count + 1
      j <- nextIntBaseB( j, b=2 )
    }
    i <- nextIndexLA( i, b=m+1 )
  }
} 
return( list(integral=const*total, functionEvaluations=count) ) }
################################################################################
SimplexVolume <- function( S ){
# compute the volume of the n-dim simplex S in R^n.  S is an n x (n+1) matrix, 
# with the columns of S giving the vertices of the simplex.  
# For S=CanonicalSimplex(n), the volume is 1/n!

stopifnot( is.matrix(S), nrow(S) >= 1, ncol(S) == nrow(S)+1 )
n <- nrow(S)
V <- S[ ,-1] - matrix( S[,1], nrow=n,ncol=n ) # n x n matrix
return( abs( det(V) )/factorial(n) ) }
################################################################################
SimplexSurfaceArea <- function( S3 ){
# compute the surface area of the (n-1)-dim simplex S3 in R^n.  S3 is an 
# n x n matrix, with the columns of S3 giving the vertices of the simplex.  
# For S3 = UnitSimplexV(n), 
#   SimplexSurfaceArea(UnitSimplexV(n))/SimplexVolume( CanonicalSimplex(n-1) ) 
#      = JacobianS2Canonical(SimplexSurfaceArea(UnitSimplexV(n)))= sqrt(n)

stopifnot( is.matrix(S3), nrow(S3) >= 1, ncol(S3) == nrow(S3) )
n <- nrow(S3)
return( JacobianS2Canonical(S3)*SimplexVolume(CanonicalSimplex(n-1)) ) }
################################################################################
JacobianS2Canonical <- function( S2 ){
# compute the Jacobian of transformation between the m-dimensional simplex S2 
# (sitting inside R^n) and the canonical simplex in R^m.
# S2 is an n x (m+1) matrix, with columns of S2 giving the vertices of the simplex.
# The method is to compute the Gram matrix G=t(V) %*% V, where V is as below.

stopifnot( is.matrix(S2), nrow(S2) > 1, ncol(S2) <= nrow(S2)+1 )
n <- nrow(S2); m <- ncol(S2)-1
V <- S2[ ,-1] - matrix( S2[,1], nrow=n,ncol=m ) # an n x m matrix

return( sqrt( abs( det( t(V) %*% V ) ) ) ) }
################################################################################
CanonicalSimplex <- function( n ) {
# return the vertices of the canonical simplex in R^n, the n-dim. simplex with
# vertices the origin (0,...,0) and the unit vectors e[1],...,e[n].
# It has n-dim. volume 1/n!

return( cbind( rep(0.0,n), diag(rep(1.0,n)) ) ) }
################################################################################
UnitSimplexV <- function( n ) {
# return the vertices (V-representaion) of the unit simplex in R^n, the (n-1) dim. simplex with
# vertices the unit vectors e[1],...,e[n].

return( diag(rep(1.0,n)) ) }
################################################################################
nextIntBaseB <- function( current.n, b ) {
# Compute the next integer in base b representation.  
# current.n is a vector of integers containing the base b representation of integer n: 
#     n=current.n[1]*b^k + current.n[2]*b^(k-1) + ... + current.n[2]*b + current.n[1].
# It is assumed that n >= 0.  
# Compute (n+1) in base b, by adding 1 and carrying.  If necessary, the 
# length of the integer vector is increased by 1 to add a new leading digit.

next.n <- current.n
k <- length(next.n)
repeat {
  next.n[k] <- next.n[k] + 1
  if (next.n[k] >= b) { next.n[k] <- 0;  carry <- TRUE
  } else { carry <- FALSE }
  k <- k-1
  if ( !carry | (k==0 ) ) break
}
if (carry) { next.n <- c(1L,next.n) }
return( next.n ) }
################################################################################
nextIndexLA <- function( current.n, b ) {
# Compute the next multi-index used in the outer sum in function LasserreAvrachenkov 
# Do this by incrementing current.n until we find one that meets the requirements.
# The number of values returned is choose( b+d, b ); after that
# a vector of NAs is returned to signal that there are no more indices.

d <- length( current.n ) 
ones <- rep( 1L, d )
next.n <- current.n - ones # subtract off the initial ones
repeat {
  next.n <- nextIntBaseB( next.n, b )
  if ( all( diff(next.n) >= 0) | (length(next.n) > d ) ) break
}
if (length(next.n) > d) { next.n <- rep( NA, d ) }
else { next.n <- next.n + ones }
return( next.n ) }
################################################################################

