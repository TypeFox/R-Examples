# Function from package 'fda' (c) 2014
ppBspline <- function(t)
{
#PPBSPINE computes the coefficients of the polynomials
# of the piecewise polynomial B-spline of order norder=length(t)-1
# corresponding to the knot sequence T
# (i.e. B(k,T)(x) = (T(1+k)-T(1))[T(1),T(2),...,T(k+1)](.-x)_{+}^{k-1}),
# the coefficients of the polynomials each defines on
# distinct gapinvals of T.
# Say there are ngap distinct gapinvals in
# T [T(knots[1]),T(knots[2])), [T(knots[2]),T(knots[3])),...,
# T(knots[ngap]),T(knots[ngap+1])],
# then the coefficients are returned in the matrix COEFF as row vectors,
# the i-th row corresponding to the coefficients of the polynomial on the
# gapinval [T(knots[i]),T(knots[i+1])),
# such that, for x in [T(knots[i]),T(knots[i+1])),
# B(k,T(i))(x) =
#   COEFF(i,1)*(x-T(i))^{k-1} + ... + COEFF(i,k-1)*(x-T(i)) + COEFF(i,k).
# Note that we assume T(1) < T(k+1), i.e. T is not a
# sequence of the same knot.
# The function returns the ngap times k matrix COEFF and the vector
# INDEX indicating where in the sequence T we find the knots.
# This code uses part of the code from the function bspline.m.

#  Argument:  t ... a row vector of NORDER + 1 knot values,
#  where NORDER is the order of the spline.

# Last modified 26 January 2013

norder <- length(t) - 1
ncoef  <- 2*(norder-1)
if (norder > 1) {
   adds  <- rep(1,norder-1)
   tt    <- c(adds*t[1], t, adds*t[norder+1])
   gapin <- (1:(length(tt)-1))[diff(tt) > 0]
   ngap  <- length(gapin)
   iseq  <- (2-norder):(norder-1)
   ind   <- outer(rep(1,ngap),iseq) + outer(gapin,rep(1,ncoef))
   tx    <- matrix(tt[c(ind)],ngap,ncoef)
   ty    <- tx - outer(tt[gapin],rep(1,ncoef))
   b     <- outer(rep(1,ngap),(1-norder):0)  + outer(gapin,rep(1,norder))
   a     <- c(adds*0, 1, adds*0)
   d     <- matrix(a[c(b)],ngap,norder)
   for (j in 1:(norder-1)) {
      for (i in 1:(norder-j)) {
	      ind1 <- i + norder - 1;
	      ind2 <- i + j      - 1;
         d[,i] <-(ty[,ind1]*d[,i] - ty[,ind2]*d[,i+1])/
                 (ty[,ind1]       - ty[,ind2])
      }
   }
   Coeff <- d
   for (j in 2:norder) {
      factor <- (norder-j+1)/(j-1)
      ind <- seq(norder,j,-1)
      for (i in ind)
         Coeff[,i] <- factor*(Coeff[,i] - Coeff[,i-1])/ty[,i+norder-j]
   }
   ind   <- seq(norder,1,-1)
   if (ngap > 1) Coeff <- Coeff[,ind] else Coeff <- matrix(Coeff[,ind],1,norder)
   index <- gapin - (norder-1)
} else {
   Coeff <- as.matrix(1)
   index <- as.matrix(1)
}

list(Coeff,index)

}
