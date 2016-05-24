monhess <- function(x, Wfd, basislist, returnMatrix=FALSE)
{
#  MONHESS evaluates the second derivative of monotone fn. wrt coefficients
#  The function is of the form h[x] <- (D^{-1} exp Wfd)(x)
#  where  D^{-1} means taking the indefinite integral.
#  The interval over which the integration takes places is defined in
#       the basis object <- WFD.
#  The derivatives with respect to the coefficients in WFD up to order
#       NDERIV are also computed, max(NDERIV) <- 2.
#  Arguments:
#  X       argument values at which function and derivatives are evaluated
#             x[1] must be at lower limit, and x(n) at upper limit.
#  WFD     a functional data object
#  Returns:
#  D2H   values of D2 h wrt c
#  TVAL  Arguments used for trapezoidal approximation to integral
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Last modified 9 May 2012 by Jim Ramsay

#  set some constants

EPS    <- 1e-4
JMIN   <-  7
JMAX   <- 15

#  get coefficient matrix and check it

coef  <- Wfd$coefs
coefd <- dim(coef)
ndim  <- length(coefd)
if (ndim > 1 & coefd[2] != 1) stop("WFD is not a single function")

#  get the basis

basis    <- Wfd$basis
rangeval <- basis$rangeval
nbasis   <- basis$nbasis
nbaspr   <- nbasis*(nbasis+1)/2
onebaspr <- matrix(1,1,nbaspr)

#  set up first iteration

width <- rangeval[2] - rangeval[1]
JMAXP <- JMAX + 1
h     <- matrix(1,JMAXP,1)
h[2]  <- 0.25
#  matrix SMAT contains the history of discrete approximations to the
#    integral
smatD2h <- matrix(0,JMAXP,nbaspr)
#  array TVAL contains the argument values used <- the approximation
#  array FVAL contains the integral values at these argument values,
#     rows corresponding to argument values
#  the first iteration uses just the endpoints
iter  <- 1
xiter <- rangeval
tval  <- xiter
if (is.null(basislist[[iter]])) {
    bmat <- getbasismatrix(tval, basis, 0, returnMatrix)
    basislist[[iter]] <- bmat
} else {
    bmat <- basislist[[iter]]
}
fx   <- as.matrix(exp(bmat %*% coef))
D2fx <- matrix(0,2,nbaspr)
m <- 0
for (ib in 1:nbasis) {
   for (jb in 1:ib) {
      m <- m + 1
      D2fxij   <- as.matrix(fx*bmat[,ib]*bmat[,jb])
      D2fx[,m] <- D2fxij
   }
}
D2fval <- D2fx
smatD2h[1,] <- width*sum(D2fx)/2
tnm <- 0.5

#  now iterate to convergence

for (iter in 2:JMAX) {
   tnm   <- tnm*2
   del   <- width/tnm
   hdel  <- del/2
   xiter <- seq(rangeval[1]+del/2, rangeval[2]-del/2, del)
   tval  <- c(tval, xiter)
   if (is.null(basislist[[iter]])) {
      bmat <- getbasismatrix(xiter, basis, 0, returnMatrix)
      basislist[[iter]] <- bmat
   } else {
      bmat <- basislist[[iter]]
   }
   fx   <- as.matrix(exp(bmat%*%coef))
   D2fx <- matrix(0,length(xiter),nbaspr)
   m <- 0
   for (ib in 1:nbasis) {
      for (jb in 1:ib) {
         m <- m + 1
         D2fxij   <- as.matrix(fx*bmat[,ib]*bmat[,jb])
         D2fx[,m] <- D2fxij
      }
   }
   D2fval <- rbind(D2fval, D2fx)
   smatD2h[iter,] <- (smatD2h[iter-1,] + del*sum(D2fx))/2
   if (iter >= max(c(JMIN,5))) {
      ind <- (iter-4):iter
      result <- polintmat(h[ind],smatD2h[ind,],0)
      D2ss   <- result[[1]]
      D2dss  <- result[[2]]
      if (all(abs(D2dss) < EPS*max(abs(D2ss))) || iter >= JMAX) {
         # successful convergence
         # sort argument values and corresponding function values
         ordind  <- order(tval)
         tval    <- tval[ordind]
         D2fval  <- as.matrix(D2fval[ordind,])
         # set up partial integral values
         lval    <- outer(rep(1,length(tval)),D2fval[1,])
         del     <- tval[2] - tval[1]
         D2ifval <- del*(apply(D2fval,2,cumsum) - 0.5*(lval + D2fval))
         D2h     <- matrix(0,length(x),nbaspr)
         for (i in 1:nbaspr) D2h[,i] <- approx(tval, D2ifval[,i], x)$y
         return(D2h)
      }
    }
    h[iter+1] <- 0.25*h[iter]
  }
}

