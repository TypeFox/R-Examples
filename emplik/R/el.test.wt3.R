el.test.wt3 <- function(x, wt, mu, maxit=25, gradtol=1e-7, Hessian = FALSE, 
                        svdtol = 1e-9, itertrace=FALSE ){
x <- as.matrix(x)
n <- nrow(x)
p <- ncol(x)
mu <- as.vector(mu)
if( length(mu) !=p )
  stop("mu must have same dimension as observation vectors.")
if( n <= p )
  stop("Need more observations than length(mu) in el.test.wt2().")

z <- t( t(x) -mu )

wt <- as.vector(wt) 
if( length(wt) != n ) stop("length of wt must equal to n=nrow(x)")
if( any(wt < 0) ) stop("wt must be >= 0")

allw <- sum(wt)

#
#    Scale the problem, by a measure of the size of a 
# typical observation.  Add a tiny quantity to protect
# against dividing by zero in scaling.  Since z*lam is
# dimensionless, lam must be scaled inversely to z.
#
TINY <- sqrt( .Machine$double.xmin )
scale <- mean( abs(z) ) + TINY
z <- z/scale
##### if( !missing(lam) ){
##### lam <- as.vector(lam)
##### lam <- lam*scale
##### if( logelr(z,rep(0,p),lam)>0 )lam <- rep(0,p)
#####  }
##### if(  missing(lam)  )
  lam <- rep(0,p)
#
#     Take some precaution against users specifying
# tolerances too small.
#

if(svdtol < TINY ) svdtol <- TINY
if(gradtol < TINY) gradtol <- TINY

#
#    Preset the weights for combining Newton and gradient
# steps at each of 16 inner iterations, starting with
# the Newton step and progressing towards shorter vectors
# in the gradient direction.  Most commonly only the Newton
# step is actually taken, though occasional step reductions
# do occur.
#

nwts <- c( 3^-c(0:3), rep(0,12) )
gwts <- 2^( -c(0:(length(nwts)-1)))
gwts <- (gwts^2 - nwts^2)^.5
gwts[12:16] <- gwts[12:16] * 10^-c(1:5)

#
#    Iterate, finding the Newton and gradient steps, and
# choosing a step that reduces the objective if possible.
#

nits <- 0
gsize <- gradtol + 1
while(  nits<maxit && gsize > gradtol  ){

  arg  <- allw + z %*% lam
###  wts1 <- as.vector( llogp(arg, 1/n) )
  wts2 <- as.vector( -llogpp(arg, 1/n) )^.5
###  wtwts1 <- wt*wts1 
###  grad <- as.matrix( z*wtwts1 )
  #############grad <- as.vector( apply( grad, 2, sum ) )
###  grad <- as.vector(rowsum(grad, rep(1, nrow(grad)) ) )
  grad <- gradf(z, wt, lam) 
  gsize <- mean( abs(grad) )
  wtwts2 <- sqrt(wt)*wts2
  hess <- z*wtwts2
#                                   -1
#    The Newton step is -(hess'hess)    grad,
#  where the matrix hess is a sqrt of the Hessian.
#  Use svd on hess to get a stable solution.
#
## may try La.svd() in R (v. > 1.0) for better LAPACK.

##  svdh <- svd( hess )
####  svdh <- La.svd( hess )
##  if( min(svdh$d) < max(svdh$d)*svdtol )
##    svdh$d <- svdh$d + max(svdh$d)*svdtol
##  nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
#### nstep <- t(svdh$vt) %*% (t(svdh$u)/svdh$d)
##  nstep <- as.vector( nstep %*% matrix(wts1/wts2,n,1) )

svdh <- La.svd( hess )
if( min(svdh$d) < max(svdh$d)*svdtol )
    svdh$d <- svdh$d + max(svdh$d)*svdtol
nstep <- t(svdh$vt) %*% (svdh$vt/(svdh$d)^2)
nstep <- as.vector( nstep %*% grad )

  gstep <- grad
  if( sum(nstep^2) < sum(gstep^2) )
    gstep <- gstep*(sum(nstep^2)^.5/sum(gstep^2)^.5)

  ologelr <- logwelr( z, rep(0,p), wt, lam ) 

  ninner <- 0
  for(  i in 1:length(nwts) ){
   nlogelr <- logwelr(z,rep(0,p),wt, lam+nwts[i]*nstep+gwts[i]*gstep )
###    ngrad <- gradf(z,wt, lam+nwts[i]*nstep+gwts[i]*gstep )
###    ngsize <- mean( abs(ngrad) ) 
    if( nlogelr  < ologelr  ){
      lam <- lam+nwts[i]*nstep+gwts[i]*gstep
      ninner <- i
      break
    }
  }
  nits <- nits+1
  if( ninner==0 )nits <- maxit
  if( itertrace )
    print( c(lam, nlogelr, gsize, ninner) )
}

Hess <- NA 
if( Hessian ) Hess <- t(hess)%*%hess*scale^2

list( prob= as.vector(wt/as.vector(allw + z %*% lam)), 
       lambda = lam/scale, grad=grad*scale, hess=Hess, nits=nits )
}


logwelr <- function(x, mu, wt, lam){ 
x <- as.matrix(x)
n <- nrow(x)
p <- ncol(x)
if(  n <= p  )
  stop("Need more observations than variables in logelr.")
mu <- as.vector(mu)
if(  length(mu) != p  )
  stop("Length of mean doesn't match number of variables in logelr.")

allw <- sum(wt) 
z <- t( t(x) -mu )
arg <- allw/(allw + z %*% lam)
return( sum( wt* llog(arg,1/n) ) ) 
}


gradf <- function(z,wt,lam) {
allw <- sum(wt)
arg <- allw + z %*% lam
n <- length(wt) 
wts1 <- as.vector( llogp(arg, 1/n) )
wtwts1<- wt*wts1
grad <- as.matrix(z*wtwts1)
as.vector( rowsum( grad, rep(1, nrow(grad)) ) ) 
}

##########################################################
#    The function llog() is equal to the natural
#  logarithm on the interval from eps >0 to infinity.
#  Between -infinity and eps, llog() is a quadratic.
#  llogp() and llogpp() are the first two derivatives
#  of llog().  All three functions are continuous
#  across the "knot" at eps.
#
#    A variation with a second knot at a large value
#  M did not appear to work as well.
#
#    The cutoff point, eps, is usually 1/n, where n
#  is the number of observations.  Unless n is extraordinarily
#  large, dividing by eps is not expected to cause numerical
#  difficulty.
#
#### These functions have been defined inside el.test(). 
#### No need to repeat.
#
#llog <- function( z, eps ){
#
#ans <- z
#lo <- (z<eps)
#ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
#ans[ !lo ] <- log( z[!lo] )
#ans
#}
#
#llogp <- function( z, eps ){
#
#ans <- z
#lo <- (z<eps)
#ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
#ans[ !lo ] <- 1/z[!lo]
#ans
#}
#
#llogpp <- function( z, eps ){
#
#ans <- z
#lo <- (z<eps)
#ans[ lo  ] <- -1.0/eps^2
#ans[ !lo ] <- -1.0/z[!lo]^2
#ans
#}
##########################################################
