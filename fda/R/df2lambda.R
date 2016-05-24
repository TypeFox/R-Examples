df2lambda <- function(argvals, basisobj, wtvec=rep(1,n), Lfdobj=0, df=nbasis)
{
#  Convert a degree of freedom DF for a smooth to the equivalent value
#    of the smoothing parameter lambda.

#  Last modified 26 October 2005

n <- length(argvals)
nbasis <- basisobj$nbasis
if (df >= nbasis) {
   lambda <- 0
   return(lambda)
}

TOL    <- 1e-3
GOLD   <- 1.0
GLIMIT <- 2.0
TINY   <- 1.0e-20

#  find machine precision
eps <- 1
tol1 <- 1 + eps
while (tol1 > 1) {
   eps  <- eps/2
   tol1 <- 1 + eps
}
eps <- sqrt(eps)

#  ------  initialization of lambda by finding bracketing values ------------
#             a < b < c such that  fb < fa  and  fb < fc
#  first use input value for lambda unless it is zero, in which case -1
bx <- -4.0
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lambda <- 10^(bx)
fb <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  now try bracketing the minimum by using a large value and a small
#  value.  If (this doesn't work, revert to the iterative method
#  at statement 5
if (bx >= -10 &&  bx <= 5) {
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   cx <- 5  #  the upper limit
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   lambda <- 10^(cx)
   fc <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ax <- -8  #  the lower limit
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   lambda <- 10^(ax)
   fa <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
}
#  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  check to see if minimum bracketed
#print(c(ax,bx,cx,fa,fb,fc, lambda))
if (fb >= fa || fb >= fc) {
  #  Failure to bracket minimum, proceed with iterative search for
  #    bracketing values.
  #  First, as an alternative value for ax, use the input value plus 0.1
  ax <- bx + 1
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  lambda <- 10^(ax)
  fa <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  now the bracketing process begins
  if (fb > fa) {
     #  exchange ax and bx
     dum <- ax
     ax  <- bx
     bx  <- dum
     dum <- fb
     fb  <- fa
     fa  <- dum
  }
  #  first guess at cx
  cx <- bx + GOLD*(bx - ax)
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  lambda <- 10^(cx)
  fc <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  check if (three values bracket minimum
  #print(c(ax,bx,cx,fa,fb,fc, lambda))
  while (fb >= fc) {
     r <- (bx - ax)*(fb - fc)
     q <- (bx - cx)*(fb - fa)
     u <- bx -
       ((bx - cx)*q - (bx - ax)*r)/(2.0*sign(max(c(abs(q-r),TINY)))*(q-r))
     ulim <- bx + GLIMIT*(cx - bx)
     if ((bx-u)*(u-cx) > 0.0) {
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        lambda <- 10^(u)
        fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (fu < fc) {
           #  success
           ax <- bx
           bx <- u
           fa <- fb
           fb <- fu
           break
        }
        if (fu > fb) {
           #  also success
           cx <- u
           fc <- fu
           break
        }
        #  failure:  fu >= fb
        u <- cx + GOLD*(cx - bx)
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        lambda <- 10^(u)
        fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     }
     if ((cx - u)*(u - ulim) > 0.0) {
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        lambda <- 10^(u)
        fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (fu < fc) {
           bx <- cx
           cx <-  u
           u  <- cx + GOLD*(cx - bx)
           fb <- fc
           fc <- fu
           #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           lambda <- 10^(u)
           fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
           #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }
     }
     if ((u-ulim)*(ulim-cx) >= 0.0) {
        u <- ulim
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        lambda <- 10^(u)
        fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     } else {
        u <- cx + GOLD*(cx - bx)
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        lambda <- 10^(u)
        fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     }
     ax <- bx
     bx <- cx
     cx <- u
     fa <- fb
     fb <- fc
     fc <- fu
     #print(c(ax,bx,cx,fa,fb,fc, lambda))
  }  #  end of while loop
}
#  ---------------------------------------------------------------------
#  --------------------  bracketing successful  ------------------------
#  ---------------------------------------------------------------------
a  <- min(c(ax,cx))
b  <- max(c(ax,cx))
v  <- bx
w  <- v
x  <- v
e  <- 0.0
fx <- fb
fv <- fx
fw <- fx
#  ---------------------------------------------------------------------
#  --------------------  main loop starts here -------------------------
#  ---------------------------------------------------------------------
xm   <- 0.5*(a + b)
tol1 <- eps*abs(x) + TOL/3
tol2 <- 2*tol1
crit <- abs(x - xm) - (tol2 - 0.5*(b - a))
#print(c(crit, lambda))
while (crit > 0) {
   #  is golden-section necessary?
   if (abs(e) > tol1) {
      #  fit parabola
      r <- (x - w)*(fx - fv)
      q <- (x - v)*(fx - fw)
      p <- (x - v)*q - (x - w)*r
      q <- 2.0*(q - r)
      if (q > 0.0)   p <- -p
      q <- abs(q)
      s <- e
      e <- d
      #  is parabola acceptable?
      if (abs(p) < abs(0.5*q*s) & p > q*(a - x) & p < q*(b - x)) {
         #  a parabolic interpolation step
         d <- p/q
         u <- x + d
         #  f must not be evaluated too close to a or b
         if ((u - a) < tol2 ||  b - u < tol2) {
            if (xm - x >= 0.0) d <- tol1 else d <- -tol1
         }
      } else {
         #  a golden-section step
         if (x >= xm) e <- a - x else e <- b - x
         d <- 0.382*e
      }
   } else {
      #  a golden-section step
      if (x >= xm) e <- a - x else e <- b - x
      d <- 0.382*e
   }
#  f must not be evaluated too close to x
   if (abs(d) >=  tol1) {
      u <- x + d
   } else {
      if (d >= 0.0) u <- x + tol1 else u <- x - tol1
  }
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  lambda <- 10^u
  fu <- (lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda) - df)^2
  #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  update  a, b, v, w, and x
  if (fu <= fx) {
     if (u  >= x) a <- x else b <- x
     v  <- w
     w  <- x
     x  <- u
     fv <- fw
     fw <- fx
     fx <- fu
  } else {
     if (u  < x) a <- u else b <- u
     if (fu <= fw || w == x) {
        v  <- w
        w  <- u
        fv <- fw
        fw <- fu
     }
     if (fu <= fv || v == x || v == w) {
        v  <- u
        fv <- fu
     }
  }
  xm   <- 0.5*(a + b)
  tol1 <- eps*abs(x) + TOL/3
  tol2 <- 2*tol1
  crit <- abs(x - xm) - (tol2 - 0.5*(b - a))
  #print(c(crit, lambda))
#  -------------------  } of main loop  ------------------------------
}
return(lambda)
}
