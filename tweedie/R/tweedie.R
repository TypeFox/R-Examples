#############################################################################
ptweedie.inversion <- function(q, mu, phi,  power, exact=FALSE ){ 
# Peter K Dunn 
# 08 Feb 2000--2005

y <- q
cdf <- array( dim=length(y) )
# Error checks
if ( power<1) stop("power must be greater than 1.")
if ( any(phi<= 0) ) stop("phi must be positive.")
if ( any(y<0) ) stop("y must be a non-negative vector.")
if ( any(mu<=0) ) stop("mu must be positive.")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
}
else {
   mu <- array( dim=length(y), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim=length(y), phi )
   # A vector of all phi's
}
  

its <- y
for (i in (1:length(y))) {

   # This has been added to avoid an issue with *very* small values of y
   # causing the FORTRAN to die when p>2; 
   # reported by Johann Cuenin 09 July 2013 (earlier, but that was the easiest email I could find about it :->)
   if ( ( power > 2 ) & (y[i] < 1.0e-300) ) {
   	### THE  e-300  IS ARBITRARY!!!!  ###
   	# Keep an eye on it; perhaps it needs changing
   
      # That is, y is very small: Use the limit as y->0 as the answer as FORTRAN has difficulty converging
      # I have kept the call to the FORTRAN (and for p>2, of course, the limiting value is 0).
      #  I could have done this differently
      # by redefining very small y as y=0.... but this is better methinks
      cdf[i] <- 0
   } 
   else
   {
		tmp <- .Fortran( "cdf",
				 as.double(power),
				 as.double(phi[i]),
				 as.double(y[i]),
				 as.double(mu[i]),
				 as.integer( exact ),
				 as.double(0), # funvalue
				 as.integer(0), # exitstatus
				 as.double(0), # relerr
				 as.integer(0), # its
				 PACKAGE="tweedie")

			 cdf[i] <- tmp[[6]]
	}
}

cdf

}
      


#############################################################################
dtweedie.dldphi.saddle <- function(phi, mu, power, y ){
# Calculates the derivative of log f wrt phi
# where the density is the saddlepoint density

# Peter Dunn
# 13 August 2002

dev <- tweedie.dev( power=power, y=y, mu=mu)
l <-  (-1)/(2*phi) + dev/(2*phi^2)

-2* sum(l)
}



 



#############################################################################
dtweedie.logl <- function(phi, y, mu, power) {
# Computes the log-likelihood for
# a Tweedie density.  

# Peter Dunn
# 26 April 2001

sum( log( dtweedie( y=y, mu=mu, phi=phi, power=power ) ) )

}


#############################################################################
dtweedie.logl.saddle <- function( phi, power, y, mu, eps=0){
# Calculates the log likelihood of Tweedie densities
# where the density is the saddlepoint density

# Peter Dunn
# 01 May 2001
sum( log( dtweedie.saddle(power=power, phi=phi, y=y, mu=mu, eps=eps) ) )

}


#############################################################################
dtweedie.logv.bigp <- dtweedie.logv.bigp <- function( y, phi, power){ 
# Peter K Dunn 
# 02 Feb 2000 
# 

#
# Error traps
#

if ( power<2) stop("power must be greater than 2.")
if ( any(phi<= 0) ) stop("phi must be positive.")
if ( any(y<=0) ) stop("y must be a strictly positive vector.")
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim = length(y), phi )
   # A vector of all phi's
}


#
# Set up
#

p <- power
a <- ( 2-p ) / ( 1-p )
a1 <- 1 - a

r <- -a1*log(phi) - log(p-2) - a*log(y) +
   a*log(p-1)
drop <- 37

#
# Now we find limits of summation, using Stirling's approximation
# to approximate the gamma terms, and find max as k.max.
#

logz <- max(r)
k.max <- max( y^(2-p) / ( phi * (p-2) ) )
k <- max( 1, k.max )

c <- logz + a1 + a*log(a)
vmax <- k.max * a1
estlogv <- vmax
#
# Now we search either side for when we can
# ignore terms as they are negligible.
# Remember we are using the envelope as the tests,
# not individual terms in the series.

# First:  the upper limit of k

while ( estlogv > (vmax - drop) ) {
	k <- k + 2
	estlogv <- k*( c - a1*log(k) )
}

hi.k <- ceiling(k)
#
# Now the lower limit of k
#
logz <- min(r)
k.max <- min( y^(2-p) / ( phi * (p-2) ) )
k <- max( 1, k.max )
c <- logz + a1 + a*log(a)
vmax <- k.max * a1
estlogv <- vmax

while ( (estlogv > (vmax-drop) ) && ( k>=2) ) {
	k <- max(1, k - 2)
	estlogv <- k*( c - a1*log(k) )
}

lo.k <- max(1, floor(k) )

#
# Now sum the series between established limits.
# We ensure it works for vector y.
#
k <- seq(lo.k, hi.k) 
o <- matrix( 1, nrow=length(y)) 

g <- matrix( lgamma( 1+a*k) - lgamma(1+k),
	nrow=1, ncol=length(k) )

og <- o %*% g
A <- outer(r, k) + og
C <- matrix( sin( -a*pi*k ) * (-1)^k,
	nrow=1, ncol=length(k) )

C <- o %*% C
m <- apply(A, 1, max)
ve <- exp(A - m)
sum.ve <- apply( ve*C, 1, sum )

# Now be careful!  Because of the +/- nature of the sin term,
# sum.ve can be very small but negative  due to subtractive
# cancellation.  Treat those carefully and separately.

neg.sum.ve <- (sum.ve<=0)
pos.sum.ve <- (sum.ve>0)

logv <- sum.ve
sum.ve[neg.sum.ve] <- 0
logv[neg.sum.ve] <- -Inf
logv[pos.sum.ve] <- log( sum.ve[pos.sum.ve] ) + m[pos.sum.ve]

list(lo=lo.k, hi=hi.k, logv=logv, k.max=k.max )

}


#############################################################################
dtweedie.logw.smallp <- function(y, phi, power){ 
#
# Peter K Dunn
# 02 Feb 2000
#

#
# Error traps
#

if ( power<1) stop("power must be between 1 and 2.")
if ( power>2) stop("power must be between 1 and 2.")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<=0) ) stop("y must be a strictly positive vector.")

#
# Set up
#
p <- power
a <- ( 2-p ) / ( 1-p )	 	# Note that a<0 for 1<p<2

a1 <- 1 - a
r <- -a*log(y) + a*log(p-1) - a1*log(phi) -
	log(2-p)		# All terms to power j

drop <- 37 			# Accuracy of terms: exp(-37)
#
# Find limits of summation using Stirling's approximation
# to approximate gamma terms, and find max as j.max
#
logz <- max(r) 			# To find largest  j  needed
j.max <- max( y^( 2-p ) / ( phi * (2-p) ) )
j <- max( 1, j.max )

cc <- logz + a1 + a*log(-a)	#
# This is all the terms to j-power, not needing any other j terms.
# The other terms are introduced when we know the values of j

wmax <- a1*j.max
estlogw <- wmax

# First, the upper limit of j
while(estlogw > (wmax-drop) ){
   j <- j + 2
   estlogw <- j*(cc - a1*log(j))
}

hi.j <- ceiling(j)

#
#
# Now the lower limit of j
#
#
logz <- min(r) 
j.max <- min( y^( 2-power ) / ( phi * (2-power) ) )

j <- max( 1, j.max)
wmax <- a1*j.max 
estlogw <- wmax 
 
# First, optimize to find the location of the maximum
while ( ( estlogw > (wmax-drop) ) && ( j>=2) ) {
   j <- max(1, j - 2)
	oldestlogw <- estlogw
   estlogw <- j*(cc-a1*log(j))
}

lo.j <- max(1, floor(j))
# 
# Now sum the series between established limits.
# We ensure it works for vector y.

j <- seq( lo.j, hi.j)      # sequence of j terms needed

o <- matrix( 1, nrow=length(y)) 	# matrix of ones
 
g <- matrix(lgamma( j+1 ) + lgamma( -a*j ),
	nrow=1, ncol=hi.j - lo.j + 1)

og <- o %*% g              # matrix of gamma terms
A <- outer(r, j) - og      # the series, almost ready to sum
m <- apply(A,1,max)        # avoid overflow; find maximum values
we <- exp( A - m )         # evaluate terms, less max.
sum.we <- apply( we,1,sum) # sum terms
logw <- log( sum.we ) + m  # now restore max.
 
list(lo=lo.j, hi=hi.j, logw=logw, j.max=j.max )

}


#############################################################################
dtweedie <- function(y, xi=power, mu, phi, power=NULL)
{
#
# This is a function for determining Tweedie densities.
# Two methods are employed:  cgf inversion (type=1)
# and series evaluation (type=2 if 1<p<2; type=3 if p>2)).
# This function uses bivariate interpolation to accurately
# approximate the inversion in conjunction with the series.
#
#
# FIRST, establish whether the series or the cgf
# # inversion is the better method.
#

#
# Here is a summary of what happens:
#
#   p                          |  Whats happpens
# -----------------------------+-------------------
#  1 -> p.lowest (=1.3)        |  Use series A
#  p.lowest (=1.3) -> p=2      |  Use series A if
#                              |        xix > xix.smallp (=0.8)
#                              |      inversion with rho if
#                              |        xix < xix.smallp (=0.8)
#  p=2 -> p=3                  |  Use series B if
#                              |        xix > xix.smallp (=0.8)
#                              |      inversion with rho if
#                              |        xix < xix.smallp (=0.8)
#  p=3 -> p.highest (=4)       |  Use series B if xix > xix.bigp (=0.9)
#                              |      inversion with 1/rho if
#                              |        xix < xix.bigp (=0.9)
#  p > 4                       |  Use series B if xix > 0.03
#                              | Saddlepoint approximation if xix<= 0.03 (NEEDS FIXING)
#

   # Sort out the xi/power notation
   if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
   xi.notation <- TRUE
   if ( is.null(power) ) power <- xi
   if ( is.null(xi) ) {
      xi <- power
      xi.notation <- FALSE
   }
   if ( xi != power ) {
      cat("Different values for xi and power given; the value of xi used.\n")
      power <- xi
   }
   index.par <- ifelse( xi.notation, "xi","p")


# Error checks
if ( power<1) stop("power must be greater than 1.\n")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<0) ) stop("y must be a non-negative vector.\n")
if ( any(mu<=0) ) stop("mu must be positive.\n")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.\n")
}
else {
   mu <- array( dim=length(y), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.\n")
}
else {
   phi <- array( dim=length(y), phi )
   # A vector of all phi's
}

density <- y

# Special Cases
if ( power==3 ){
	density <- dinvgauss(x=y, mean=mu, dispersion=phi)
	return(density)
}
if ( power==2 ) {
   density <- dgamma( rate=1/(phi*mu), shape=1/phi, x=y )
   return(density)
}
if ( power==0) {
   density <- dnorm( mean=mu, sd=sqrt(phi), x=y )
   return(density)
}
if ( (power==1) & (all(phi==1))) {
	# Poisson case
	density <- dpois(x=y/phi, lambda=mu/phi )
   return(density)
}

# Set up
id.type0 <- array( dim=length(y) )
id.series <- id.type0
id.interp <- id.type0

###   Now consider the cases 1<p<2   ###
id.type0 <- (y==0)
if (any(id.type0)) {
   if ( power > 2 ) {
      density[id.type0] <- 0
   }
   else {
      lambda <- mu[id.type0]^(2-power)/(phi[id.type0]*(2-power))
      density[id.type0] <- exp( -lambda )
   }
}


xi <- array( dim=length(y) )
xi[id.type0] <- 0
xi[!id.type0] <- phi[!id.type0] * y[!id.type0]^(power-2)
xix <- xi / ( 1 + xi )

if ( (power>1) && (power<=1.1) ) {
   id.series <- (!id.type0)
	if (any(id.series)){
      density[id.series] <- dtweedie.series(y=y[id.series],
                     mu=mu[id.series], phi=phi[id.series], power=power)
   }
   return(density=density)
}

if ( power==1 ) { # AND phi not equal to one here
	id.series <- rep(TRUE, length(id.series))
	id.interp <- rep(FALSE, length(id.series))
}
if ( (power>1.1) && (power<=1.2) ) {
   id.interp <- ( (xix>0) & (xix<0.1) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.1
      p.hi <- 1.2
      xix.lo <- 0
      xix.hi <- 0.1
      np <- 15
      nx <- 25
   }
}

if ( (power>1.2) && (power<=1.3) ) {
   id.interp <- ( (xix>0) & (xix<0.3) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.2
      p.hi <- 1.3
      xix.lo <- 0
      xix.hi <- 0.3
      np <- 15
      nx <- 25
   }
}
if ( (power>1.3) && (power<=1.4) ) {
   id.interp <- ( (xix>0) & (xix<0.5) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.3
      p.hi <- 1.4
      xix.lo <- 0
      xix.hi <- 0.5
      np <- 15
      nx <- 25
   }
}
if ( (power>1.4) && (power<=1.5) ) {
   id.interp <- ( (xix>0) & (xix<0.8) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.4
      p.hi <- 1.5
      xix.lo <- 0
      xix.hi <- 0.8
      np <- 15
      nx <- 25
   }
}
if ( (power>1.5) && (power<2) ) {
   id.interp <- ( (xix>0) & (xix<0.9) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 1.5
      p.hi <- 2
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
   }
}

### Cases p>2   ###
if ( (power>2) && (power<3) ) {
   id.interp <- ( (xix>0) & (xix<0.9) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 2
      p.hi <- 3
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
   }
}
if ( (power>=3) && (power<4) ) {
   id.interp <- ( (xix>0) & (xix<0.9) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 3
      p.hi <- 4
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
   }
}
if ( (power>=4) && (power<5) ) {
   id.interp <- ( (xix>0) & (xix<0.9) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 4
      p.hi <- 5
      xix.lo <- 0
      xix.hi <- 0.9
      np <- 15
      nx <- 25
   }
}
if ( (power>=5) && (power<7) ) {
   id.interp <- ( (xix>0) & (xix<0.5) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 5
      p.hi <- 7
      xix.lo <- 0
      xix.hi <- 0.5
      np <- 15
      nx <- 25
   }
}
if ( (power>=7) && (power<=10) ) {
   id.interp <- ( (xix>0) & (xix<0.3) )
   id.series <- (!(id.interp|id.type0))
   if ( any(id.interp)) {
      grid <- stored.grids(power)
      p.lo <- 7
      p.hi <- 10
      xix.lo <- 0
      xix.hi <- 0.3
      np <- 15
      nx <- 25
   }
}

if ( power>10) {
   id.series <- (y!=0)
   id.interp <- (!(id.series|id.type0))
}

if (any(id.series)) {
   density[id.series] <- dtweedie.series(y=y[id.series],
            mu=mu[id.series], phi=phi[id.series], power=power)
}

if (any(id.interp)) {
   dim( grid ) <- c( nx+1, np+1 )
   rho <- dtweedie.interp(grid, np=np, nx=nx,
               xix.lo=xix.lo, xix.hi=xix.hi,
               p.lo=p.lo, p.hi=p.hi,
               power=power, xix=xix[id.interp] )
   dev <- tweedie.dev(power = power, mu = mu[id.interp], y = y[id.interp])
	front <- rho/(y[id.interp] * sqrt(2 * pi * xi[id.interp]))
	density[id.interp] <- front * exp(-1/(2 * phi[id.interp]) * dev)
}


density

}


#############################################################################
dtweedie.saddle <- function(y, xi=power, mu, phi, eps=1/6, power=NULL) {
#
# Peter K Dunn
# 09 Jan 2002
#
#

   # Sort out the xi/power notation
   if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
   xi.notation <- TRUE
   if ( is.null(power) ) power <- xi
   if ( is.null(xi) ) {
      xi <- power
      xi.notation <- FALSE
   }
   if ( xi != power ) {
      cat("Different values for xi and power given; the value of xi used.\n")
      power <- xi
   }
   index.par <- ifelse( xi.notation, "xi","p")

# Error traps
#
	if( any(phi <= 0) )
		stop("phi must be positive.")
	if( (power >=1) & any(y < 0))
		stop("y must be a non-negative vector.")
	if( any(mu <= 0) )
		stop("mu must be positive.")
   if ( length(mu)==1 )
      mu <- array( dim=length(y), mu )
   if ( length(phi)==1 )
      phi <- array( dim=length(y), phi )
#
# END error checks
#
#   if ( any(y==0) ) y <- y + (1/6)
   y.eps <- y 
   
   if (power<2) y.eps <- y + eps
   # Nelder and Pregibon's idea for problems at y=0

   y0 <- (y == 0)
	density <- y
#	if(any(y == 0)) {
#		if(power >= 2) {
#         # This is a problem!
#			density[y0] <- 0 * y[y0]
#		}
#      else{

   dev <- tweedie.dev(y=y, mu=mu, power=power)
   density <- (2*pi*phi*y.eps^power)^(-1/2) * exp( -dev/(2*phi) )

#      }
#   }
   if ( any(y==0) ){
		if((power >= 1) && (power < 2)) {
			lambda <- mu[y0]^(2 - power)/(phi[y0] * (2 - power))
			density[y0] <- exp( -lambda)
		} else {
			density[y0] <- 0
		}
	}
	
#	if(any(y != 0)) {
#     dev <- tweedie.dev(y=y[yp], mu=mu[yp], power=power)
#      density[yp] <- (2*pi*phi[yp]*y[yp]^power)^(-1/2) * exp( -dev/(2*phi[yp]) )
#	}
	density
}


#############################################################################
dtweedie.series.bigp <- function(power, y, mu, phi){ 
 
# 
# Peter K Dunn 
# 02 Feb 2000 
# 

#
# Error traps
#

if ( power<2) stop("power must be greater than 2.")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<=0) ) stop("y must be a strictly positive vector.")
if ( any(mu<=0) ) stop("mu must be positive.")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
}
else {
   mu <- array( dim=length(y), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim=length(y), phi )
   # A vector of all phi's
}


result <- dtweedie.logv.bigp(power=power, y=y, phi=phi)
logv <- result$logv

theta <- mu^(1-power) / ( 1 - power )
kappa <- mu^(2-power) / ( 2 - power )

logfnew <- (y*theta-kappa)/phi - log( pi*y) + logv
f <- exp( logfnew )

list(density=f, logv=logv, lo=result$lo, hi=result$hi )

}


#############################################################################
dtweedie.series <- function(y, power, mu, phi){ 
#
# Peter K Dunn
# 09 Jan 2002
#

#
# Error traps
#

if ( power<1) stop("power must be between 1 and 2.")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<0) ) stop("y must be a non-negative vector.")
if ( any(mu<=0) ) stop("mu must be positive.")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
}
else {
   mu <- array( dim = length(y), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim = length(y), phi )
   # A vector of all phi's
}


y0 <- (y == 0 )
yp <- ( y!=0 )
density <- array( dim=length(y))

if ( (power == 2) | (power==1) ) { # Special cases
   if ( power == 2 ){
      density <- dgamma( y, shape=1/phi, rate=1/(phi * mu ) )
   }

   if ( (power == 1) ){
		density <- dpois(x=y/phi, lambda=mu/phi ) / phi  # Using identity: f(y; mu, phi) = c f(cy; c mu, c phi) for p=1.  Now set c = 1/phi
		if ( !all(phi==1)){
			warnings("The density computations when phi=1 may be subject to errors using floating-point arithmetic\n")
		}
   }
# CHANGED: 31 October, thanks to Dina Farkas' email 18 October 2012.  WAS:
#   if ( (power == 1) & (all(phi==1)) ){
#	density <- dpois(x=y/phi, lambda=mu/phi )
#   }
}
else{

   if ( any(y==0) ) {
      if ( power>2 ) {
         density[y0] <- 0*y[y0]
      }
      if ( ( power>1) && (power<2) ) {
         lambda <- mu[y0]^(2-power) / ( phi[y0] * (2-power) )
         density[y0] <- exp( -lambda )
      }
   }
   
   if ( any( y!=0 ) ) { 
      if ( power > 2 ) {
   	   density[yp] <- dtweedie.series.bigp(
   		   power=power,mu=mu[yp], y=y[yp], phi=phi[yp])$density
      }
   
      if ( ( power > 1 ) && ( power < 2 ) ) {
   	   density[yp] <- dtweedie.series.smallp(
   		   power=power,mu=mu[yp], y=y[yp], phi=phi[yp])$density
      }
   }
}

density

}
 



#############################################################################
dtweedie.series.smallp <- function(power, y, mu, phi){ 

#
# Peter K Dunn
# 02 Feb 2000
#

#
# Error traps
#

if ( power<1) stop("power must be between 1 and 2.")
if ( power>2) stop("power must be between 1 and 2.")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<=0) & (power>=2) ) stop("y must be a strictly positive vector.")
if ( any(y<0) & (power>1) & (power < 2) ) stop("y must be a non-negative vector.")
if ( any(mu<=0) ) stop("mu must be positive.")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.")
}
else {
   mu <- array( dim = length(y), mu )
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim=length(y), phi )
}

result <- dtweedie.logw.smallp( y=y, power=power, phi=phi)
logw <- result$logw

 
tau <- phi * (power-1) * mu^( power-1 )
lambda <- mu^( 2-power ) / ( phi * ( 2-power ) )
logf <- -y/tau - lambda - log(y) + logw
f <- exp( logf )

list(density=f, logw=logw, hi=result$hi, lo=result$lo)

}
 



#############################################################################
ptweedie <- function(q, xi=power, mu, phi, power=NULL) {
# Evaluates the cdf for
# Tweedie distributions.

# Peter Dunn
# 01 May 2001

if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
xi.notation <- TRUE
if ( is.null(power) ) power <- xi
if ( is.null(xi) ) {
   xi <- power
   xi.notation <- FALSE
}
if ( xi != power ) {
   cat("Different values for xi and power given; the value of xi used.\n")
   power <- xi
}
index.par <- ifelse( xi.notation, "xi","p")

y <- q

y.negative <- (y<0)
y.original <- y
y <- y[ !y.negative ]

# Error checks
if ( power<1) stop("power must be greater than 1.\n")
if ( any(phi<=0) ) stop("phi must be positive.")
# if ( any(y<0) ) stop("y must be a non-negative vector.\n")
if ( any(mu<=0) ) stop("mu must be positive.\n")
if ( length(mu)>1) {
   if ( length(mu)!=length(y.original) ) stop("mu must be scalar, or the same length as y.\n")
}else {
   mu <- array( dim=length(y.original), mu )
   # A vector of all mu's
}
mu.original <- mu
mu <- mu[ !y.negative ]

if ( length(phi)>1) {
   if ( length(phi)!=length(y.original) ) stop("phi must be scalar, or the same length as y.\n")
}else {
   phi <- array( dim=length(y.original), phi )
   # A vector of all phi's
}
phi.original <- phi
phi <- phi[ !y.negative ]



cdf.positives <- array( dim=length(y) )
cdf <- y.original

# Special Cases
if ( power==2 ) {
   f <- pgamma( rate=1/(phi*mu), shape=1/phi, q=y )
}
if ( power==0) {
   f <- pnorm( mean=mu, sd=sqrt(phi), q=y )
}
if ( power==1) {
   f <- ppois(q=y, lambda=mu/phi )
}

# Now, for p>2 the only option is the inversion
if ( power> 2 ) {

	# However, when y/q is very small, 
	# the function should return 0;
	# sometimes it fails (when *very* small...).
	# This adjustment is made in ptweedie.inversion()
   f <- ptweedie.inversion(power=power, mu=mu, q=y, phi=phi)
}

# Now for 1<p<2, the two options are the series or inversion.
# We avoid the series when p is near one, otherwise it is fine.
if ( (power>1) & (power<2) ) {
	if ( power <1.7 ) {
	   f <- ptweedie.series(power=power, q=y, mu=mu, phi=phi )
	}
	else{
	   f <- ptweedie.inversion( power=power, q=y, mu=mu, phi=phi)
	}
}

cdf[ !y.negative ] <- f
cdf[  y.negative ] <- 0

cdf[ is.infinite( cdf ) ] <- 1

return(cdf)

}



#############################################################################
ptweedie.series <- function(q, power, mu, phi) {
# 
# Peter K Dunn 
# 14 March 2000 
#
y <- q
# Error checks
if ( power<1) stop("power must be between 1 and 2.\n")
if ( power>2) stop("power must be between 1 and 2.\n")
if ( any(phi<= 0)) stop("phi must be positive.\n")
if ( any(y<0) ) stop("y must be a non-negative vector.\n")
if ( any(mu<=0) ) stop("mu must be positive.\n")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.\n")
}else {
   mu <- array( dim=length(y), mu )
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.\n")
}else {
   phi <- array( dim=length(y), phi )
}


# Set up
p <- power
lambda <- mu^(2-p) / ( phi * (2-p) )
tau <- phi * (p-1) * mu^( p-1 )
alpha <- (2-p) / ( 1-p )


# Now find the limits on N:
# First the lower limit on N

drop <- 39
lambda <- max(lambda )
logfmax <-  -log(lambda)/2
estlogf <- logfmax
N <- max( lambda )

while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
   estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
}

lo.N <- max(1, floor(N) )


# Now for the upper limit on N

lambda <- mu^(2-p) / ( phi * (2-p) )
lambda <- min( lambda )
logfmax <-  -log(lambda)/2
estlogf <- logfmax

N <- min( lambda )

while ( estlogf > (logfmax - drop) ) {
   N <- N + 1
   estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
}

hi.N <- max( ceiling(N) )

# Now evaluate between limits of N

cdf <- array( dim=length(y), 0 )

lambda <- mu^(2-p) / ( phi * (2-p) )
tau <- phi * (p-1) * mu^( p-1 )
alpha <- (2-p) / ( 1-p )


for (N in (lo.N:hi.N)) {

   # Poisson density
   pois.den <- dpois( N, lambda)

   # Incomplete gamma
   incgamma.den <- pchisq(2*y/tau, -2*alpha*N )

   # What we want
   cdf <- cdf + pois.den * incgamma.den

}

cdf <- cdf + exp( -lambda )
its <- hi.N - lo.N + 1


cdf

}



#############################################################################
stored.grids <- function(power){
# This S-Plus function contains all the stored interpolation
# grids for interpolating the Tweedie densities.

# Peter Dunn
# 18 April 2001

# Grids stored
#
##### PART 1:  1 < p < 2 ####
#
# 1.5 < p < 2 (A)
if ( (power>1.5) && ( power<2) ) {
grid <- c(0.619552027046892,
 -1.99653020187375,
 -8.41322974185725,
 -36.5410153849965,
 -155.674714429765,
 -616.623539139617,
 -2200.89580607229,
 -6980.89597754589,
 -19605.6031284706,
 -48857.4648924699,
 -108604.369946264,
 -216832.50696005,
 -391832.053251182,
 -646019.900543419,
 -979606.689484813,
 -1377088.40156071,
 -1808605.4504829,
 -2236021.71690435,
 -2621471.84422418,
 -2935291.72910167,
 -3160933.2895598,
 -3296019.98597692,
 -3350191.70672175,
 -3341188.71720722,
 -3290375.55904795,
 -3256251.04779701,
 0.55083671442337,
 2.79545515380255,
 6.43661784251532,
 0.0593276289419893,
 -88.7776267014812,
 -592.599945216246,
 -2715.1932676716,
 -9995.06110621652,
 -31004.2502694815,
 -82973.5646780374,
 -194644.329882664,
 -405395.126655972,
 -758085.608783886,
 -1285893.38788744,
 -1997446.83841278,
 -2866954.62503587,
 -3834679.24827227,
 -4818512.83125468,
 -5732403.6885235,
 -6504847.15089854,
 -7091578.60837028,
 -7479819.10447519,
 -7684747.8312771,
 -7739400.44941064,
 -7651256.22196861,
 -10780871.374931,
 -0.296443816521139,
 -0.740400436722078,
 13.5437223703132,
 132.570369303337,
 739.235236905108,
 3222.27476321793,
 11883.2578358337,
 37985.4069102717,
 106385.512302156,
 263186.759184411,
 579709.939191383,
 1146196.24715222,
 2051166.63217944,
 3349844.60761646,
 5033457.4421439,
 7014240.88296282,
 9135328.47624901,
 11203377.8720252,
 13031732.9028836,
 14478751.8662425,
 15470353.6878961,
 16004512.3987063,
 16147291.8371304,
 16087058.36173,
 17256553.6861257,
 -90196970.3516856,
 0.212508829604868,
 -1.25345600957563,
 -33.6344002424504,
 -229.413615894842,
 -961.2927083918,
 -3048.22512983433,
 -7865.94887016698,
 -16619.1472772416,
 -27426.4597246056,
 -29170.0602809743,
 6910.62618199269,
 133529.230915424,
 425031.289018802,
 962079.672230075,
 1805724.1265787,
 2970428.54627254,
 4407924.65252939,
 6009208.15039934,
 7623903.3639327,
 9090617.14760774,
 10273366.8474699,
 11112760.9034168,
 11783847.5998635,
 13967546.0605232,
 40580571.9968577,
 -1800585369.9502,
 -0.117096128441307,
 3.86814002139823,
 52.1297859521411,
 231.959187763769,
 297.55942409085,
 -2028.44588042558,
 -16035.7900988595,
 -69727.4207742442,
 -231001.098781009,
 -634711.96272568,
 -1500867.48368527,
 -3120952.565006,
 -5795192.65274283,
 -9728268.26803975,
 -14923011.3657198,
 -21127733.2735346,
 -27877621.6155107,
 -34629531.4516162,
 -40942366.5228353,
 -46618585.271072,
 -51693745.2168584,
 -56076194.3321918,
 -57846478.5394169,
 -39555239.6466308,
 244583291.9647,
 -19387464175.4177,
 -0.0251908119117217,
 -6.97888536964718,
 -57.2053622635109,
 -7.59149018702844,
 1940.74883722613,
 13776.3994355061,
 59688.9506688118,
 198967.267063524,
 554664.181209208,
 1342679.76760535,
 2880718.85057692,
 5556233.09746906,
 9746453.11803942,
 15705233.9839578,
 23444873.7106245,
 32630213.4195304,
 42473186.150293,
 51604301.4143779,
 57950838.5609574,
 58814629.2970318,
 51670532.360853,
 37073253.8293895,
 31112689.7825927,
 166639899.499491,
 2285908471.65836,
 -142300999941.473,
 0.239462332088407,
 9.7006247266445,
 24.6796376740648,
 -621.977746808668,
 -6061.14667092794,
 -28846.5231563469,
 -92974.9043842097,
 -229420.907671922,
 -455750.587752974,
 -720917.683643727,
 -798268.858897086,
 -142553.959799974,
 2227938.90487962,
 7722331.28668364,
 17902944.5216494,
 33799508.0633365,
 54727373.5621622,
 76699837.6804652,
 90928852.1062459,
 83791368.1450994,
 41369883.2049391,
 -33644574.3006822,
 -51466205.6969009,
 728822186.855893,
 12613725489.8081,
 -779054232341.425,
 -0.531163176360929,
 -9.62133388925249,
 83.5858879354598,
 1749.74435647472,
 10388.878276731,
 29739.3251394829,
 27228.8352302813,
 -137459.132798184,
 -797275.579291472,
 -2522656.85139472,
 -6002370.79311501,
 -11564879.4570243,
 -18293067.3260108,
 -23137574.7865943,
 -20926214.0959888,
 -6386212.74814566,
 20947334.9587757,
 49084025.5298454,
 45674483.7892477,
 -43357709.595878,
 -274284042.960864,
 -630264236.584573,
 -707988499.827136,
 2768071791.46789,
 55323323631.9224,
 -3346981844559.16,
 0.843274592042774,
 1.85634928776026,
 -303.233463115877,
 -2964.77206695662,
 -8035.18601516459,
 26623.9135533587,
 296203.748735318,
 1300209.48867761,
 4002931.81502221,
 10154091.5309971,
 23011763.2727463,
 48484290.5041288,
 96405481.1482104,
 181157702.817844,
 319621265.640707,
 523128884.608782,
 780656497.356576,
 1034030709.78307,
 1153235261.31943,
 934435310.723166,
 172736683.699636,
 -1059296923.94426,
 -1259513524.21048,
 11540211151.7843,
 202982012848.938,
 -11765811129127.3,
 -0.971763767126055,
 20.7278527481909,
 599.216068887268,
 2508.62108313741,
 -15799.3976275259,
 -190836.337882638,
 -873333.216664437,
 -2412364.04114457,
 -4372947.5863571,
 -4280053.82275206,
 4581124.92368874,
 37692021.4517723,
 125690565.165227,
 317577060.569135,
 672918414.800248,
 1235540175.83018,
 1981789079.14886,
 2741889556.9801,
 3118003228.13883,
 2478792089.71153,
 201910662.627241,
 -3440015764.71965,
 -3701941985.52277,
 36624484882.9033,
 631340208742.919,
 -35016783929182,
 0.465205915130579,
 -62.4878848976034,
 -713.029716677298,
 3463.89852940917,
 76620.1536606951,
 425050.797698706,
 1098749.49892239,
 830042.192958006,
 -3609253.07165929,
 -14019647.104919,
 -22008231.8565183,
 2879363.57573307,
 127763068.732814,
 471678257.276824,
 1199156009.13846,
 2449442490.02104,
 4177372240.464,
 5937970631.74612,
 6701001745.26596,
 4850964750.68393,
 -1252857624.6872,
 -11008101341.7195,
 -11613916592.1703,
 99208893729.6641,
 1723996766015.73,
 -90686535774445.6,
 1.35055309427014,
 109.632652005803,
 -15.551388984273,
 -19110.1792293505,
 -151248.805892165,
 -326520.136900253,
 1323004.57105885,
 11263794.9628682,
 40182893.7043995,
 99473459.937714,
 213007606.45782,
 462071289.943091,
 1036322096.46886,
 2233228139.73472,
 4392262214.22542,
 7778537815.0104,
 12326926829.6142,
 17068944624.2195,
 19342475640.096,
 14633223559.5287,
 -1238504576.84042,
 -24780929364.4666,
 -19774890981.5309,
 267278604099.047,
 4262049717379.08,
 -209109831966511,
 -4.90254965773785,
 -106.16559587471,
 2473.32253262609,
 39971.9971896279,
 105932.466006557,
 -940508.183476386,
 -7972770.02131952,
 -26249511.2854385,
 -38016600.0226778,
 28952783.4114192,
 296603175.897931,
 921440155.664983,
 2211825458.1676,
 4873052470.30472,
 9959226055.61093,
 18095018536.2409,
 28333607390.4192,
 37736361000.0061,
 42033152625.5253,
 35930405200.4038,
 11674689966.7752,
 -38568497842.5504,
 -74929319913.8016,
 440819854424.984,
 9367147142753.47,
 -436516234054193,
 9.21654006976202,
 -53.4970979572424,
 -6659.50032530198,
 -38474.8203617517,
 310488.103734294,
 3845463.04045484,
 14525744.8285833,
 14191442.4974447,
 -71954549.6112557,
 -299246254.252719,
 -416148031.634689,
 348026958.617362,
 2939406987.4056,
 7971636742.33878,
 16399104338.5733,
 31135394894.2163,
 55589206993.0977,
 84773664582.7313,
 93124980530.7634,
 37204300253.4301,
 -100809540619.219,
 -211523088540.365,
 124648961402.01,
 2454987217464.61,
 23136880987249.1,
 -833737986871973,
 -10.6044400547848,
 457.68066225078,
 9572.011521514,
 -37331.5997372033,
 -1178378.58119297,
 -5650741.83842416,
 1857995.74642602,
 109904372.796869,
 458568718.1367,
 915187897.2047,
 811019172.09732,
 92558199.5987161,
 1998081069.83546,
 13128498122.3811,
 36243363277.9486,
 61424177454.5951,
 73667085976.0338,
 86031421495.4807,
 166947042118.509,
 384532819727.471,
 591661667761.69,
 78943846328.2111,
 -2527580012466.73,
 -6916515802327.55,
 21469457257461.8,
 -1.4991232719638e+015,
 3.17664940253603,
 -1029.70926107305,
 -5007.37894552437,
 206353.644207492,
 1762192.60826409,
 -794281.059186644,
 -59885497.7612809,
 -267968880.144224,
 -316522642.887633,
 1411438830.94685,
 7247583561.26104,
 15913831315.8407,
 19845848486.0143,
 17971388397.4505,
 42579150295.9815,
 155992918698.532,
 356278537036.23,
 414110925668.561,
 -161408422221.427,
 -1685255831372.5,
 -3094542673159.69,
 -53087024607.4896,
 16845038275179.6,
 63212141053606.3,
 202776604073488,
 -2.43199634134677e+015 )
}

# 1.4 < p < 1.5 (B)
if ( (power>1.4) && (power<=1.5) ) {
grid <- c(0.52272502143588,
 -2.37970276267028,
 -4.50105382434541,
 2.01559335438892,
 26.7080425402173,
 53.1211365756307,
 56.5864772151065,
 83.6551659699119,
 218.128756577377,
 414.067827117006,
 586.302276594556,
 925.084114296478,
 1629.94052297329,
 2248.04335981282,
 2325.84860904575,
 3056.70827313074,
 5936.34902374365,
 7512.68000065931,
 -494.046409795432,
 -14632.889451627,
 2555.32119471521,
 103793.057384242,
 239732.321730146,
 57045.261742077,
 -1151133.66663602,
 -3657091.74594141,
 0.775913134802564,
 4.5833536253197,
 -6.51027140852265,
 -110.820592889565,
 -302.27161508231,
 -170.893952600532,
 597.013908647114,
 1036.46479866453,
 362.328670850611,
 1143.67350470972,
 5427.58842419152,
 7481.1938127866,
 3916.80746756334,
 10919.4206700089,
 37992.348893648,
 35578.2283061684,
 -48344.4290328334,
 -70796.935934355,
 307989.153472163,
 965852.855353596,
 280146.00264422,
 -4462973.64355295,
 -12247413.3846698,
 -7748551.19364459,
 47938914.6705982,
 196377667.248434,
 -1.15101738661085,
 -3.82183996501946,
 79.2652077448879,
 447.712841790259,
 183.895349671323,
 -4007.2875514913,
 -11554.3202406694,
 -11119.3558179111,
 -3321.49044452746,
 -24694.4788687522,
 -79648.9320167266,
 -73209.1251877769,
 -2021.74771903288,
 -154892.733325762,
 -580785.160917161,
 -291546.650641339,
 1473778.38160541,
 1519601.02525224,
 -7163122.30868759,
 -21848349.4729583,
 -6440961.41109336,
 106830967.379257,
 319671808.797567,
 306212673.299803,
 -944079002.967604,
 -4580578568.46053,
 2.02189579962441,
 -10.0772337050698,
 -303.024532759887,
 -704.351532760201,
 4973.76150820855,
 26142.4723427423,
 31405.9126438149,
 -53806.275525627,
 -149359.118457136,
 16092.9360420975,
 206115.379917112,
 -509055.061292639,
 -1505746.99603438,
 440398.507094803,
 4047882.74286314,
 -3759763.21945456,
 -28362810.7587706,
 -20497981.7777362,
 115049844.871681,
 326041126.501133,
 57527117.4877412,
 -1804448247.28013,
 -5607533850.24561,
 -6742827026.7323,
 11967660233.3111,
 137394861109.439,
 -2.64873193579361,
 68.3499218097747,
 718.463817413271,
 -1571.25771413107,
 -28554.9807730767,
 -66684.6392549349,
 112682.126175196,
 683976.718235879,
 693804.971402528,
 -977072.064073228,
 -991329.886033511,
 5636906.19208592,
 7342984.59011112,
 -20840829.3577062,
 -49193222.8209866,
 42836025.7323296,
 238149829.126778,
 -55474045.455026,
 -1800027569.31187,
 -4174218574.22397,
 -583545105.082449,
 22566852342.6238,
 72170080841.219,
 99488552940.42,
 -95607027136.1,
 2344789183218.79,
 0.265173429404147,
 -224.175784222991,
 -799.783571636822,
 14465.1992899084,
 83422.4433648738,
 -33747.4954822623,
 -1182034.20503755,
 -2562724.00673901,
 1327150.31866679,
 9096932.76169983,
 -8147841.7628163,
 -67997474.423206,
 -84381101.1698542,
 1167786.45920554,
 -290262963.143079,
 -1966559598.21217,
 -4495910335.48392,
 -2778049317.5256,
 9570853857.51408,
 20045770280.4768,
 -39053969152.1498,
 -304481715916.404,
 -861160659351.957,
 -1248601137158.62,
 1149200212746.68,
 152670892302636,
 12.5290983345227,
 499.005760177699,
 -1884.50478980683,
 -52292.5606105439,
 -102007.716028595,
 974234.695542423,
 4265255.61960463,
 -518423.706368479,
 -37847638.3377804,
 -94994579.0086004,
 -141430225.033227,
 -554077045.707268,
 -2621090414.90387,
 -7585723619.93884,
 -14489129113.1742,
 -23451598610.8142,
 -52316888100.9028,
 -154472187133.545,
 -386170628207.441,
 -654732126120.445,
 -458996711714.355,
 1207352849976.42,
 5003368112913.01,
 8162669646404.62,
 10014475965674.9,
 4.30099151054367e+015,
 -49.3184821796504,
 -659.190055418151,
 13760.6854640598,
 117281.576980704,
 -286688.694056081,
 -4514692.65697881,
 -11833557.5599173,
 -13306480.8841683,
 -177569989.290577,
 -1501492416.77122,
 -6595461136.24562,
 -20356400466.8674,
 -54254546277.6263,
 -140855615592.576,
 -350477425656.988,
 -772978728940.34,
 -1445553693725.57,
 -2349378847823.61,
 -3825995343158.07,
 -7818306940077.94,
 -19731827478910.9,
 -48657002077699.6,
 -98504762092397.2,
 -121989778390269,
 454260903660862,
 1.09114573020388e+017,
 123.410260081819,
 -197.403953092516,
 -40101.3176876088,
 -52951.2659905612,
 2144698.8673286,
 1739654.2351533,
 -141996380.772718,
 -1289342757.59316,
 -7156566175.76133,
 -32325427550.7015,
 -126169094996.878,
 -419281859332.1,
 -1174950784066.64,
 -2822039940568,
 -6037668177832.46,
 -12096555834995.5,
 -23612374044250.6,
 -45167080717259.5,
 -82416904456952.1,
 -138051430229770,
 -204193982413479,
 -251800695384389,
 -208815124278779,
 288829117097310,
 9.54395313223954e+015,
 2.24604455857544e+018,
 -210.739560523616,
 5859.20505890453,
 140328.667144469,
 278664.828817554,
 -10534350.9218635,
 -206287796.083656,
 -2616524539.48725,
 -22630598386.6276,
 -138968736141.657,
 -644711350455.876,
 -2387358710732.05,
 -7357036386363.44,
 -19420572655537.7,
 -44732883747785.5,
 -91030332840422.6,
 -165757883253373,
 -276320429460859,
 -441737780262524,
 -728934745124743,
 -1.31361122221926e+015,
 -2.49168560525055e+015,
 -4.31067214964099e+015,
 -4.70964682712407e+015,
 9.73624002713054e+015,
 2.04847011281059e+017,
 4.09466705815479e+019,
 278.504340481087,
 1227.95202954615,
 609240.949141676,
 8398669.53224352,
 -122964349.27965,
 -4120134389.77069,
 -51789735957.0207,
 -413475908382.113,
 -2400754425461.51,
 -10787324023927.9,
 -38888467331936.8,
 -115207923139742,
 -285196579464986,
 -596235892662438,
 -1.05638430488185e+015,
 -1.57370592252881e+015,
 -1.90815996658209e+015,
 -1.68347182504384e+015,
 -470951786697564,
 2.27515760918722e+015,
 8.65004127903261e+015,
 2.98774925717651e+016,
 1.1577555109681e+017,
 4.88806957320561e+017,
 3.5377171188424e+018,
 6.72365027054536e+020,
 1302.48019924681,
 259560.540734449,
 8128828.4793897,
 25490162.9304377,
 -3094217693.64666,
 -72611077275.6044,
 -860362188931.341,
 -6697568428379.51,
 -37943648370441.8,
 -165469263635562,
 -574639366352108,
 -1.61975362474037e+015,
 -3.72208462865178e+015,
 -6.84085667350701e+015,
 -9.24907873336246e+015,
 -5.63230537269553e+015,
 1.49215241150529e+016,
 7.00536184155271e+016,
 1.84186220493898e+017,
 3.92386129591378e+017,
 7.72461936211387e+017,
 1.59366458882732e+018,
 3.85548874292868e+018,
 1.14456878085735e+019,
 5.75359458832251e+019,
 1.02900591814312e+022,
 13286.5118009426,
 2546171.88500447,
 80098454.5079186,
 -575521939.430233,
 -60032299026.9678,
 -1202606432540.75,
 -13452077960416.8,
 -101362184753825,
 -559725106993583,
 -2.37495806056525e+015,
 -7.95202306660676e+015,
 -2.11888175159205e+016,
 -4.40669841594096e+016,
 -6.47957726367261e+016,
 -3.24083222588877e+016,
 1.83291719462516e+017,
 8.40260894845939e+017,
 2.37093537889737e+018,
 5.4317336077623e+018,
 1.10375221140593e+019,
 2.11436475962051e+019,
 4.09291799720987e+019,
 8.65527998606027e+019,
 2.14451316905219e+020,
 8.54313682198285e+020,
 1.49594792908599e+023,
 203942.274027855,
 26468105.5883892,
 462332761.175384,
 -25284736506.3552,
 -1094062620100.06,
 -19258710333321.7,
 -203938699172942,
 -1.48348319911623e+015,
 -7.94195851466556e+015,
 -3.25428940850727e+016,
 -1.03863076023785e+017,
 -2.5561068913892e+017,
 -4.48713805903153e+017,
 -3.40534549726198e+017,
 1.18853947957601e+018,
 6.7824313161356e+018,
 2.1555870304696e+019,
 5.41093701368525e+019,
 1.17698643723088e+020,
 2.32697280763342e+020,
 4.35104831467264e+020,
 8.06863397336643e+020,
 1.57397090527128e+021,
 3.47588746977727e+021,
 1.24037541169807e+022,
 2.082483901283e+024,
 2425722.75874371,
 168829967.462896,
 -6124967796.72205,
 -684794159351.057,
 -20099159890195.9,
 -312088468164664,
 -3.08510744950999e+015,
 -2.12874367646947e+016,
 -1.08016089267394e+017,
 -4.13278775037366e+017,
 -1.18371276048634e+018,
 -2.32798990878473e+018,
 -1.62965158901462e+018,
 9.87177193402913e+018,
 5.68548659646942e+019,
 1.93636952907764e+020,
 5.21964303390126e+020,
 1.21047675802231e+021,
 2.51556116094745e+021,
 4.81371375246211e+021,
 8.68998842506269e+021,
 1.52222181674142e+022,
 2.68233759254245e+022,
 5.09349061051762e+022,
 1.81405407714602e+023,
 2.73961792646156e+025,
 18727182.1067274,
 -3856365344.26257,
 -418983873891.394,
 -17333732195464.3,
 -378965106591800,
 -5.04766648385378e+015,
 -4.45572615875926e+016,
 -2.73473979347584e+017,
 -1.17835431873962e+018,
 -3.3155488196652e+018,
 -3.27497940087677e+018,
 2.43638130489573e+019,
 1.80581815169857e+020,
 7.62727068811908e+020,
 2.48326825663032e+021,
 6.79744241369944e+021,
 1.63174880753305e+022,
 3.52624193788357e+022,
 6.98893205641466e+022,
 1.28897681504313e+023,
 2.23908254454823e+023,
 3.69904872503628e+023,
 5.82169304070961e+023,
 8.85953055589028e+023,
 3.25024671656323e+024,
 3.52294080805527e+026 )
}
# 1.3 < p < 1.4 (C)

if ( (power>1.3) && (power<=1.4) ) {
grid <- c(0.897713960567798,
 -0.47040810564112,
 -1.43992915137813,
 -5.23984219701719,
 -15.0352102570745,
 -13.4481420654518,
 87.0680569141198,
 319.051289202574,
 -18.4506632683117,
 -2319.89441748794,
 -3118.4617838911,
 13713.8934655348,
 46257.8044960381,
 -39135.1988030987,
 -498172.85786169,
 -811816.300243115,
 2568821.88452726,
 15957456.1105785,
 30538856.3230172,
 -38265076.8403202,
 -445019432.131326,
 -1535324610.183,
 -2685571226.57687,
 1844994500.64329,
 33799853247.9255,
 148548951511.487,
 -0.0411924170235428,
 -0.280923312875219,
 1.85631524787589,
 35.4548526987637,
 182.842029466341,
 186.535484241998,
 -1885.62536085722,
 -7050.5577782269,
 3914.68802807428,
 72775.307234475,
 81966.2319220996,
 -546556.025881189,
 -1736496.33116182,
 2018406.08511184,
 22431395.8530421,
 38371009.900033,
 -121432141.656335,
 -823775133.379304,
 -1829474775.72627,
 1064685431.79,
 23026018217.409,
 92967101827.417,
 209244156825.792,
 108176144299.254,
 -1494179320052.26,
 -8820798547246.88,
 0.335385258677322,
 4.50345806518915,
 10.608217004451,
 -198.77147800601,
 -1523.07566142071,
 -1500.09456259484,
 23139.2326521245,
 83055.2779182394,
 -113519.71875397,
 -1189391.83133768,
 -880016.414378636,
 11594447.2485996,
 32944023.5477499,
 -57906962.3154716,
 -530418894.202135,
 -898586763.98809,
 3160587598.72435,
 22164670922.5832,
 54438271034.5952,
 -7405444239.37665,
 -615716148031.174,
 -2843896714598.56,
 -7637460249121.44,
 -9706633868213.88,
 27251052973973.8,
 220225910667289,
 -0.792263177969423,
 -21.2045420582639,
 -91.4401569487452,
 1101.7598394385,
 10118.5419251243,
 5726.65994681374,
 -215851.229167288,
 -669893.217282033,
 1905585.55556376,
 13591395.7207961,
 2340171.0149139,
 -169616173.653862,
 -405049039.820392,
 1178374485.90503,
 8651034672.11598,
 13041124217.199,
 -63654520224.3913,
 -427730047807.293,
 -1124157151846.8,
 -274204341784.276,
 11179066056046.7,
 58371933813335.8,
 179156520740047,
 326715205946158,
 -110966263637577,
 -7.30089254118998e+015,
 1.98088010280869,
 88.7205095644789,
 451.25451837431,
 -6238.38424401488,
 -56974.0505836968,
 23643.1529279025,
 1692687.00100628,
 4116879.50794305,
 -20945171.4628138,
 -107415978.096988,
 132557253.129546,
 2070091722.93035,
 3657706483.5209,
 -19429026052.7314,
 -118401385300.158,
 -176359135255.114,
 835686139242.973,
 5703797966574.08,
 14884329900045.1,
 -552681886101.364,
 -188818869070189,
 -1.00878535432539e+015,
 -3.36911029646898e+015,
 -7.50185785257662e+015,
 -5.12493009299165e+015,
 -1.75430627966896e+017,
 -6.01808947590792,
 -358.631844263425,
 -1646.82372824411,
 35116.9756113627,
 280756.645438638,
 -507546.635114422,
 -10020156.6805911,
 -4945705.75502759,
 272336260.80756,
 1104711765.97144,
 -568290834.123902,
 -14088675880.1234,
 -12950574108.0922,
 210192537963.341,
 859065438623.435,
 -747462235105.65,
 -20708041247987.7,
 -104755203973610,
 -298469477349205,
 -363165453972938,
 1.27630293527266e+015,
 1.00317977358032e+016,
 3.80411831558389e+016,
 9.62493626170414e+016,
 1.15786608399486e+017,
 -1.03118547220491e+019,
 22.952229421138,
 1405.69142020857,
 3967.15419654569,
 -181648.014611196,
 -997642.253786536,
 8248812.38538573,
 91248769.6906494,
 269940828.639906,
 457556637.301566,
 7822037651.92945,
 72574078872.7822,
 285799266923.8,
 186491437806.251,
 -3589451347123.35,
 -20096636648266.9,
 -55011152917522.8,
 -75141293052079.6,
 -68199466557531.5,
 -1.03828417550274e+015,
 -1.01420393701492e+016,
 -5.89973268530095e+016,
 -2.54452415456125e+017,
 -8.83842320399224e+017,
 -2.50090967000384e+018,
 -5.75293047799529e+018,
 -3.0288264126223e+020,
 -99.9878158346858,
 -5210.49423696155,
 4761.55706178008,
 1014983.38061518,
 6678222.58484158,
 18138060.5674288,
 523215144.725566,
 8388725586.42611,
 65308743699.6125,
 317410830648.742,
 1134002934158.94,
 3498250050884.25,
 9275186197008.01,
 5168202521004.28,
 -177786455579982,
 -1.541966007589e+015,
 -8.17352077554882e+015,
 -3.26417185677986e+016,
 -1.04838259404326e+017,
 -2.78718323122462e+017,
 -6.27055158391557e+017,
 -1.27711947908717e+018,
 -2.98493427948817e+018,
 -1.08753143842026e+019,
 -6.81489155059389e+019,
 -7.09103794715832e+021,
 445.7016632619,
 18118.8922546317,
 -59575.916393475,
 -1546683.20227499,
 66577416.4833025,
 1581334778.6952,
 18270427567.616,
 159620726230.676,
 1179819333919.81,
 7136701581206.33,
 33141064579502.4,
 108828204724357,
 187282794130913,
 -416027319900596,
 -5.25364928574569e+015,
 -2.85113559834757e+016,
 -1.22522620595716e+017,
 -4.79778248465207e+017,
 -1.80908798010203e+018,
 -6.6584955339892e+018,
 -2.38169794720652e+019,
 -8.23438092498963e+019,
 -2.75142145405642e+020,
 -8.98926942417278e+020,
 -3.30676128118007e+021,
 -1.39506678075182e+023,
 -1872.21347826633,
 -46283.567888389,
 1588637.59435705,
 58147592.3245715,
 1251197774.35208,
 23222717220.2576,
 321130227736.287,
 3181320366971.15,
 23454746025596.1,
 134242541498815,
 608809039493410,
 2.13579959983314e+015,
 4.94319739507647e+015,
 -642441041478663,
 -8.45372610438658e+016,
 -5.78902547364894e+017,
 -2.75526324472539e+018,
 -1.0840919422109e+019,
 -3.76905193553513e+019,
 -1.20991644215191e+020,
 -3.73703358932369e+020,
 -1.16027081029794e+021,
 -3.77617115231973e+021,
 -1.34025753565226e+022,
 -5.99737264585084e+022,
 -2.43233695351887e+024,
 7556.65517982232,
 215106.469973897,
 8930450.66728664,
 660119294.92026,
 21206807863.5722,
 398339343315.102,
 5190034087210.52,
 50568961091221.7,
 380791925478234,
 2.24722396598682e+015,
 1.03784283558222e+016,
 3.62656117204867e+016,
 8.1707753418827e+016,
 -1.42136634601931e+016,
 -1.3434681243398e+018,
 -8.98526482243229e+018,
 -4.27270839102888e+019,
 -1.72138580695355e+020,
 -6.27800609966133e+020,
 -2.15417480180286e+021,
 -7.16024766698911e+021,
 -2.36639068137229e+022,
 -7.98151888085422e+022,
 -2.8511432263888e+023,
 -1.2325305554076e+024,
 -4.09104391632781e+025,
 -22430.9879586296,
 1389052.54923454,
 207349675.587708,
 10195666599.2626,
 296796553336.983,
 5643928776333.53,
 74892895567209.3,
 732432041468099,
 5.49414382342695e+015,
 3.23706890283122e+016,
 1.50367841248151e+017,
 5.32587551175168e+017,
 1.22108658524644e+018,
 -2.40776349600474e+017,
 -2.10569831857027e+019,
 -1.42079098762757e+020,
 -6.7576612032549e+020,
 -2.69444367649145e+021,
 -9.63469651591202e+021,
 -3.22501888704947e+022,
 -1.04839963336164e+023,
 -3.43153835335537e+023,
 -1.17225158559623e+024,
 -4.36196591662948e+024,
 -1.98653397314499e+025,
 -6.88448478890627e+026,
 113277.793293176,
 18517911.460342,
 2450455304.46694,
 133732757996.437,
 3939467103342.68,
 73996025955050,
 973922420426797,
 9.51229943083948e+015,
 7.14628383239514e+016,
 4.211386770595e+017,
 1.94758654233723e+018,
 6.79901662283943e+018,
 1.47957074222858e+019,
 -9.95844999508158e+018,
 -3.10086406715516e+020,
 -2.02215387815306e+021,
 -9.57663684541226e+021,
 -3.84223024829774e+022,
 -1.39073983190476e+023,
 -4.72891356250685e+023,
 -1.5634530908607e+024,
 -5.19829965024042e+024,
 -1.79911666610008e+025,
 -6.76349804228038e+025,
 -3.07798865328979e+026,
 -1.20960393431303e+028,
 262418.415973756,
 250096978.583001,
 31172920464.8184,
 1642871008200.88,
 48188034311242.9,
 904074332252990,
 1.18607968360778e+016,
 1.15420861185206e+017,
 8.6517873709592e+017,
 5.09510325906077e+018,
 2.3530063632105e+019,
 8.13726112964088e+019,
 1.67367130237733e+020,
 -2.20436191917162e+020,
 -4.39432362856901e+021,
 -2.78056375326784e+022,
 -1.30625919423068e+023,
 -5.22563427007622e+023,
 -1.88958355853632e+024,
 -6.42962247956913e+024,
 -2.13314085393945e+025,
 -7.14689344739656e+025,
 -2.50319784231518e+026,
 -9.53595733424881e+026,
 -4.39966873221631e+027,
 -2.25020366975876e+029,
 5301456.42169776,
 2874720507.00564,
 367880928787.512,
 19643816240801.3,
 579437137354404,
 1.09310772186823e+016,
 1.44543936097009e+017,
 1.42253933395633e+018,
 1.08266855728833e+019,
 6.5042335654303e+019,
 3.08534291740846e+020,
 1.11258933788151e+021,
 2.54876590864897e+021,
 -1.16381734289056e+021,
 -5.15203133618314e+022,
 -3.49775856982996e+023,
 -1.70121626254362e+024,
 -6.97428046783589e+024,
 -2.57225393421212e+025,
 -8.89902171064717e+025,
 -2.99264095778419e+026,
 -1.01251337109722e+027,
 -3.56015066420899e+027,
 -1.34824053722105e+028,
 -6.3098335548906e+028,
 -4.46843594631534e+030,
 51969153.4318949,
 37970527243.0619,
 5016911719944.11,
 280591474635532,
 8.66499271311873e+015,
 1.71435134742076e+017,
 2.38968675963669e+018,
 2.4977218053678e+019,
 2.03817113129374e+020,
 1.33075091271939e+021,
 7.03087446839742e+021,
 2.99159383585797e+022,
 9.87178121143328e+022,
 2.15294869859065e+023,
 -1.70137308294605e+022,
 -3.20794310177486e+024,
 -2.12464634368172e+025,
 -9.96697641227305e+025,
 -3.96959854025141e+026,
 -1.44146335908902e+027,
 -5.00266337924097e+027,
 -1.73044335238055e+028,
 -6.20780217337989e+028,
 -2.41346445464175e+029,
 -1.2385425771224e+030,
 -9.58840383704782e+031 )
}

# 1.2 < p < 1.3 (D)
if ( (power>1.2) && (power<=1.3) ) {
grid <- c(0.959582735112296,
 -0.196662499936635,
 -0.301224766165899,
 -0.502046137683381,
 -1.5895374772965,
 -7.89601110410117,
 11.3307142810595,
 762.363953933895,
 7332.01631877668,
 29138.4582145098,
 -76339.1153989754,
 -1796970.2346991,
 -11726396.1851003,
 -33450515.2165733,
 98652180.9916008,
 1769595410.03626,
 11618980670.3263,
 49132193375.8717,
 121727250506.162,
 -108084821334.712,
 -3366026267162.29,
 -24251373433568,
 -125561039968382,
 -532052675965940,
 -1.87778088137625e+015,
 -5.20253340562286e+015,
 -0.0117920451301422,
 -0.0865135588789834,
 -0.305001740311124,
 0.805339553985313,
 45.1127229820372,
 543.393289708699,
 2052.1244600063,
 -23479.8112741772,
 -378633.827624817,
 -2160148.0589237,
 -71997.7261576385,
 94765697.6597504,
 790575791.875542,
 3071423100.48738,
 -1263803153.75136,
 -106897577346.909,
 -858303914821.876,
 -4302053589524.5,
 -14378154556128.4,
 -18792439949861.5,
 158250693996112,
 1.65848181037028e+015,
 1.01561192369251e+016,
 4.93870814001288e+016,
 2.03707862081918e+017,
 7.32155754220541e+017,
 0.0276487092578896,
 0.273727106195384,
 3.10493645078303,
 6.74883014953071,
 -791.517414710189,
 -14659.1601086036,
 -97513.5598522334,
 253826.352485024,
 9659750.39303433,
 72426555.4027922,
 108077398.518517,
 -2603530107.50454,
 -27142206424.4266,
 -129010093038.665,
 -115592228710.678,
 3373975177466.36,
 33069166217378.1,
 190970279120807,
 775088407219876,
 1.95531526988697e+015,
 -1.08844022433567e+015,
 -5.03041904941278e+016,
 -3.88737032793942e+017,
 -2.17996177307102e+018,
 -1.03733702308722e+019,
 -4.56687817243336e+019,
 0.00460342184603284,
 -1.4191835457503,
 -61.0896531126379,
 -722.851866770932,
 5719.49208464532,
 231804.146910697,
 2166106.69539304,
 160219.206657996,
 -174724415.264779,
 -1623587498.26411,
 -4008185392.86052,
 51909707843.9719,
 660864914267.222,
 3692989887945.52,
 7286367447296.07,
 -66279529836476.3,
 -829978185197398,
 -5.3925210853509e+015,
 -2.46096213563743e+016,
 -7.65535837790057e+016,
 -7.14762806859001e+016,
 1.1170095179102e+018,
 1.09817546139889e+019,
 6.89213311167972e+019,
 3.45878357931231e+020,
 9.50469739120499e+019,
 -0.135631224023216,
 9.82565018790751,
 752.94430372433,
 13626.2467481671,
 2250.45356087176,
 -2801646.63778468,
 -34431818.5075219,
 -60859775.5058763,
 2414520669.31045,
 26210017396.4224,
 71434279662.7282,
 -950461823661.22,
 -12870038897298.3,
 -75651275138190.2,
 -153087109622533,
 1.57838782399087e+015,
 2.08458508751872e+016,
 1.47352057035772e+017,
 7.65309425661461e+017,
 3.07056896036942e+018,
 8.78357943721663e+018,
 7.63366506587975e+018,
 -1.20231415827737e+020,
 -1.19381796423493e+021,
 -8.4234840567679e+021,
 -2.34505645250979e+023,
 2.86398306366288,
 4.76124711367818,
 -6457.28547283359,
 -164758.27150296,
 -589766.226217846,
 28067458.5191206,
 401139336.264519,
 416911989.44211,
 -40222720768.9054,
 -459374207097.421,
 -1682302721723.89,
 12458076906678.7,
 216651931418543,
 1.51961566257658e+015,
 5.4584519020853e+015,
 -5.19685396766197e+015,
 -2.27022332259378e+017,
 -1.85664800967564e+018,
 -9.86237368369173e+018,
 -3.61931853333836e+019,
 -6.07700089132706e+019,
 3.81186434434807e+020,
 4.88459347203027e+021,
 3.36478940057831e+022,
 1.56888688507016e+023,
 -1.17704898582546e+025,
 -30.8405895817803,
 -700.185732152437,
 47091.7507004153,
 1596708.4034919,
 7159288.38315759,
 -321689179.440809,
 -5431369518.38533,
 -21502235882.8679,
 277109738167.276,
 3653526301164.72,
 5678876796467.48,
 -241265821098598,
 -2.84217859105037e+015,
 -1.41758365590388e+016,
 1.445542591378e+016,
 9.32012909683014e+017,
 9.89863985847271e+018,
 7.14065483941308e+019,
 4.10045064062402e+020,
 1.97534191394671e+021,
 8.07769891709407e+021,
 2.70685067720764e+022,
 6.02927395033448e+022,
 -9.22337544984228e+022,
 -3.1648353300078e+024,
 -5.61603302303069e+026,
 250.323709544488,
 8545.64176293525,
 -338614.039444803,
 -14362564.9062655,
 -92763398.7027217,
 2174220601.06076,
 29854756860.9811,
 -212387816186.495,
 -8028244379539.38,
 -76503721502058.2,
 -140612055233030,
 4.75984421187865e+015,
 6.87179074752808e+016,
 5.46755810778416e+017,
 3.00178419969484e+018,
 1.15075325799871e+019,
 2.56621820643217e+019,
 -9.99967497065107e+018,
 -2.40812092499658e+020,
 1.000269417608e+021,
 2.80929218333346e+022,
 2.76237750806896e+023,
 1.92887184219875e+024,
 1.03524743837727e+025,
 2.54120088308415e+025,
 -1.67374835060249e+028,
 -1788.89737939434,
 -72709.9323496133,
 2484233.26866509,
 105824573.878487,
 46058548.2567914,
 -41008709091.9079,
 -661049991528.956,
 -2869882451384.02,
 34629600071672.9,
 584592196626057,
 3.85506398038373e+015,
 1.85280581626924e+016,
 2.40594389635909e+017,
 4.39015768159735e+018,
 5.68500042585068e+019,
 5.46101970218862e+020,
 4.19859227844407e+021,
 2.71529223079585e+022,
 1.52448364339692e+023,
 7.57284365298019e+023,
 3.350130760701e+024,
 1.30047897758996e+025,
 4.0904969542102e+025,
 5.75746965876781e+025,
 -8.76067257210248e+026,
 -4.18601579421029e+029,
 12278.7581956566,
 496755.477226686,
 -22174016.5110608,
 -982354213.175663,
 -7651083547.2257,
 74418384151.9456,
 -158221447812.226,
 -43635118674338.1,
 -493258537684304,
 3.84848286761028e+015,
 1.75887871994533e+017,
 2.68669123829776e+018,
 2.71502811657687e+019,
 2.09549316384593e+020,
 1.32297232667709e+021,
 7.18370559114638e+021,
 3.53542513917885e+022,
 1.68168650207952e+023,
 8.25241080182299e+023,
 4.30466100873181e+024,
 2.32738116938991e+025,
 1.23477038775956e+026,
 6.04371816144202e+026,
 2.416150272357e+027,
 6.92811011890208e+026,
 -8.30613806320748e+030,
 -79244.2377537492,
 -2985145.62680809,
 144884969.0038,
 3426964337.87191,
 -120514004452.734,
 -4568263757378.07,
 -41894084093912.4,
 490584024082356,
 1.9241818444463e+016,
 3.11732137423038e+017,
 3.73577725478084e+018,
 3.90759600611774e+019,
 3.77869191059759e+020,
 3.36829721948857e+021,
 2.71562326419553e+022,
 1.95921942858408e+023,
 1.26492353680313e+024,
 7.3517449416279e+024,
 3.87469051660304e+025,
 1.86123499753071e+026,
 8.13297579530344e+026,
 3.16817378492111e+027,
 1.01489843927767e+028,
 1.5807400240241e+028,
 -2.021313503304e+029,
 -1.45957971542625e+032,
 521788.722135468,
 9472176.27503808,
 -1834442451.42211,
 -67207743391.7279,
 -811715015878.492,
 -2843970947585.28,
 171797242776594,
 1.07249545636155e+016,
 3.34114929712701e+017,
 6.59909611038221e+018,
 9.29872040038554e+019,
 1.00743611185641e+021,
 8.83945953164831e+021,
 6.53568861391816e+022,
 4.21089667672736e+023,
 2.43823913494771e+024,
 1.30691447936152e+025,
 6.66339503265359e+025,
 3.30039394627076e+026,
 1.60447868150269e+027,
 7.63096240611393e+027,
 3.48620982415311e+028,
 1.46850577313574e+029,
 5.0140288773257e+029,
 -3.94001569880189e+029,
 -2.26766193173512e+033,
 -2825375.34945871,
 -61496352.0334524,
 2493173658.3079,
 -333237181754.782,
 -18996407047060.9,
 -185565155812803,
 9.2705029708914e+015,
 3.99039220138775e+017,
 8.61232487663729e+018,
 1.32051416047768e+020,
 1.60969350838808e+021,
 1.6490151337151e+022,
 1.46466548357974e+023,
 1.14951997029669e+024,
 8.07500862310061e+024,
 5.1281421080201e+025,
 2.96949619046387e+026,
 1.57946761642313e+027,
 7.76011998009655e+027,
 3.52921742493438e+028,
 1.47769012971021e+029,
 5.55807112615323e+029,
 1.72069255611094e+030,
 2.45634325498066e+030,
 -4.33879732862619e+031,
 -3.40784734034525e+034,
 17658315.6693612,
 -1121888367.87085,
 -181270807554.95,
 -7242038986871.72,
 -120082281227313,
 2.31291589803502e+015,
 2.20848277428698e+017,
 7.48543822457755e+018,
 1.61954010343843e+020,
 2.55218480608558e+021,
 3.13237500947909e+022,
 3.13138431695314e+023,
 2.6379745983487e+024,
 1.92552803357783e+025,
 1.24720050890185e+026,
 7.32106996964179e+026,
 3.9678324994807e+027,
 2.01741972525514e+028,
 9.74388679831022e+028,
 4.5059355340519e+029,
 1.99801024333608e+030,
 8.4241281541234e+030,
 3.28635082828553e+031,
 1.06272567840542e+032,
 -2.16825023052698e+032,
 -5.26834661824363e+035,
 -48488127.6322903,
 -9941368958.95411,
 -1701495949545.85,
 -115729090727230,
 -2.46271860575519e+015,
 5.70594335232465e+016,
 4.91770401202199e+018,
 1.54976753508251e+020,
 3.17533901633655e+021,
 4.84563283115304e+022,
 5.87562438473976e+023,
 5.8931457323464e+024,
 5.03145371691109e+025,
 3.73914198714188e+026,
 2.46316680755782e+027,
 1.46055114567583e+028,
 7.89759555011506e+028,
 3.9367543253666e+029,
 1.82414003560443e+030,
 7.89299615934093e+030,
 3.18195160616357e+031,
 1.17447804036354e+032,
 3.74155783420197e+032,
 7.41309337505642e+032,
 -9.98125051400098e+033,
 -9.74447097008623e+036,
 575794042.269378,
 -359696847349.006,
 -46582766175582.5,
 -2.5172826170398e+015,
 -4.85788166089857e+016,
 1.35929935720063e+018,
 1.12759088306295e+020,
 3.61755309466896e+021,
 7.52761720304016e+022,
 1.15671836837518e+024,
 1.39982433553636e+025,
 1.39034166332529e+026,
 1.16802699889137e+027,
 8.49867342754942e+027,
 5.46129812043187e+028,
 3.15144234524267e+029,
 1.65678656217278e+030,
 8.03580116507344e+030,
 3.63383187207197e+031,
 1.54337891859975e+032,
 6.16359372138549e+032,
 2.28426586046228e+033,
 7.47381246679803e+033,
 1.68250380067942e+034,
 -1.71346139855961e+035,
 -2.23719793136107e+038 )
}

# 1.1 < p < 1.2 (E)
if ( (power>=1.1) && (power<=1.2) ) {
grid <- c(0.989973425448979,
 -0.111802566436515,
 -0.128138198419178,
 -0.144020752096856,
 -0.156865568330961,
 -0.139657194482135,
 3.536695146237,
 303.956020562169,
 16971.561119878,
 711511.636211499,
 23703218.8754648,
 650925619.121534,
 15122827515.9695,
 302907419807.887,
 5300237420118.51,
 81627637427388.5,
 1.10626268382253e+015,
 1.30088329870981e+016,
 1.26350384665605e+017,
 8.34748022753782e+017,
 -1.48193266830101e+018,
 -1.83141490787536e+020,
 -4.36156639238075e+021,
 -8.28126332215202e+022,
 -1.08728181988638e+024,
 -1.54431621438905e+027,
 -0.00313462642593353,
 -0.0388702138075973,
 -0.0903770763700018,
 -0.178015741979943,
 -0.442864053335794,
 -13.4654296356606,
 -879.397544291817,
 -44326.6781688309,
 -1719945.84474206,
 -52165024.459267,
 -1245203281.74841,
 -23434690273.4911,
 -352258956554.21,
 -4722829656798.21,
 -83507513328313.6,
 -2.39585027904263e+015,
 -7.50927529719736e+016,
 -2.07270273612729e+018,
 -4.97162068932377e+019,
 -1.06766044761913e+021,
 -2.12429857352596e+022,
 -4.05944617085845e+023,
 -7.80448328338183e+024,
 -1.62547597498508e+026,
 -3.41127169961574e+027,
 5.72214427762114e+029,
 0.00507283741509471,
 0.0611283656924362,
 0.118275683282882,
 0.0638053314994744,
 -22.662058007833,
 -2973.06775411472,
 -243880.688361423,
 -13804096.0539272,
 -566905197.509389,
 -17310061654.4572,
 -391792916235.597,
 -6189711302728.5,
 -49123570530885.4,
 601194433055370,
 2.98825837615e+016,
 5.47657960626698e+017,
 3.5405751519937e+018,
 -1.26138966768788e+020,
 -5.9951048018011e+021,
 -1.65556052103626e+023,
 -3.74652939321509e+024,
 -7.7582256629277e+025,
 -1.58711907280473e+027,
 -3.48648576340569e+028,
 -7.77485808458045e+029,
 1.03410922602955e+031,
 0.000415849616352672,
 0.0106967682516999,
 0.158203224877345,
 4.49743806433715,
 -1145.2512510092,
 -257762.508706388,
 -25867475.247582,
 -1637247291.63787,
 -72907884189.7384,
 -2409044555520.72,
 -60655326006420.5,
 -1.17024127884226e+015,
 -1.69507534186056e+016,
 -1.73353934726338e+017,
 -1.19768287172967e+018,
 -1.6925948580369e+019,
 -8.40082861320748e+020,
 -3.19703011126093e+022,
 -9.01680307946117e+023,
 -2.10464678914191e+025,
 -4.37986494585875e+026,
 -8.62023394900792e+027,
 -1.70595448526874e+029,
 -3.65540979450112e+030,
 -7.98948970235937e+031,
 -8.72229109997141e+033,
 -0.000277296401191326,
 0.00785581341291909,
 7.70945326818073,
 1047.85637415002,
 12617.6081671164,
 -11318694.4184715,
 -1521632559.16444,
 -109651752055.073,
 -5365129636121.41,
 -195049584142469,
 -5.55257586362403e+015,
 -1.28875628635706e+017,
 -2.54416168432133e+018,
 -4.52586741722226e+019,
 -7.85882741468167e+020,
 -1.43675527082477e+022,
 -2.81191719043273e+023,
 -5.63847556272856e+024,
 -1.10380771755939e+026,
 -2.06538042850477e+027,
 -3.71079218491912e+028,
 -6.55949719267708e+029,
 -1.19216339970128e+031,
 -2.37651568510983e+032,
 -4.91119938960528e+033,
 -9.19298686871738e+035,
 -4.0453365352258e-005,
 0.808120615148473,
 462.33674259018,
 75430.0290598068,
 4331914.30073851,
 -205449452.816475,
 -53008737555.3822,
 -4488085605111.43,
 -241208396247575,
 -9.50346382864154e+015,
 -2.94900174959011e+017,
 -7.58627631848485e+018,
 -1.69435784683881e+020,
 -3.4411087434676e+021,
 -6.63519114738198e+022,
 -1.2514488828429e+024,
 -2.32881811283864e+025,
 -4.24297462063675e+026,
 -7.48476399731269e+027,
 -1.27258719106238e+029,
 -2.09917093736239e+030,
 -3.42920519941086e+031,
 -5.77122270489562e+032,
 -1.06549974420896e+034,
 -2.07816556556317e+035,
 -4.87536422889045e+037,
 0.00675158161047212,
 30.709935211898,
 17928.1999465246,
 3251282.07045102,
 259150408.193801,
 4065970794.89839,
 -1081487969698.15,
 -120396596572924,
 -7.24129443659267e+015,
 -3.08345729970917e+017,
 -1.02588642970565e+019,
 -2.83199105144724e+020,
 -6.79734935044741e+021,
 -1.47579853708826e+023,
 -2.99031990839464e+024,
 -5.7671579517374e+025,
 -1.06688395883241e+027,
 -1.89309297237781e+028,
 -3.21672741997238e+029,
 -5.24303649293225e+030,
 -8.26753906534395e+031,
 -1.28566458271899e+033,
 -2.0446711266049e+034,
 -3.53884600906098e+035,
 -6.61800844109463e+036,
 -1.75572778206776e+039,
 0.258679354840137,
 840.360258770242,
 507444.100015697,
 99442112.1901258,
 9368084176.4115,
 380375494385.308,
 -9258609923520.86,
 -2.23106090141392e+015,
 -1.58174211439346e+017,
 -7.3607104165464e+018,
 -2.61924232181302e+020,
 -7.67289085556648e+021,
 -1.94419837765708e+023,
 -4.42140613397016e+024,
 -9.26766241853947e+025,
 -1.82037390050501e+027,
 -3.37830777373629e+028,
 -5.94348917058151e+029,
 -9.9356636154327e+030,
 -1.58518493151072e+032,
 -2.43604478698234e+033,
 -3.67029806957723e+034,
 -5.60458388706875e+035,
 -9.22449768201003e+036,
 -1.6904092353277e+038,
 -4.89465564411877e+040,
 6.41800151101403,
 18182.0558161604,
 11380508.3514703,
 2361269590.73996,
 245103026658.95,
 12999718384259.7,
 144508294149080,
 -3.1075671302527e+016,
 -2.81566498473887e+018,
 -1.44690887425048e+020,
 -5.49152466538323e+021,
 -1.692387361606e+023,
 -4.47164792383022e+024,
 -1.05093110630748e+026,
 -2.25304900366832e+027,
 -4.47763141879538e+028,
 -8.32719512434137e+029,
 -1.45736744448634e+031,
 -2.41132989810688e+032,
 -3.79399183745648e+033,
 -5.72975750077454e+034,
 -8.44163520665995e+035,
 -1.2506421655368e+037,
 -1.98161951684487e+038,
 -3.63563587421314e+039,
 -1.14232297834824e+042,
 130.4805948356,
 319626.889155575,
 206944151.13871,
 44870875123.7641,
 4947555746666.1,
 294990039461152,
 6.86125098068619e+015,
 -3.63786637655628e+017,
 -4.49266130472303e+019,
 -2.52923657396585e+021,
 -1.01259704639731e+023,
 -3.24472934374372e+024,
 -8.83324503429516e+025,
 -2.1213437136394e+027,
 -4.60888186519505e+028,
 -9.21064839386007e+029,
 -1.71143638740235e+031,
 -2.97850907929798e+032,
 -4.88496458063687e+033,
 -7.6009042549272e+034,
 -1.13248736207305e+036,
 -1.64007876899463e+037,
 -2.37417054431763e+038,
 -3.65864990028412e+039,
 -6.84510916139291e+040,
 -2.35049544393324e+043,
 963.245715385098,
 4945125.74645911,
 3444934228.87347,
 766212865526.459,
 85513510174417.3,
 5.0967420057197e+015,
 1.1088205871047e+017,
 -7.45318428114027e+018,
 -8.70993014114392e+020,
 -4.89483272136854e+022,
 -1.96863338240187e+024,
 -6.33450769555172e+025,
 -1.72614090113824e+027,
 -4.13046184681991e+028,
 -8.89861663528251e+029,
 -1.7563102367483e+031,
 -3.21437889043872e+032,
 -5.50329625743619e+033,
 -8.87758069741474e+034,
 -1.35900361585279e+036,
 -1.99194876368027e+037,
 -2.83400436317181e+038,
 -4.01653931243524e+039,
 -6.05394644066378e+040,
 -1.16801065759799e+042,
 -4.42572447936007e+044,
 51233.824363294,
 71261780.1129107,
 41401713028.4753,
 8934258413489.44,
 988188857001118,
 5.68035904215963e+016,
 8.47175769846526e+017,
 -1.32231153172396e+020,
 -1.32772648545315e+022,
 -7.30083365623776e+023,
 -2.93724557054933e+025,
 -9.51930409068657e+026,
 -2.6176588474417e+028,
 -6.31768365981176e+029,
 -1.37023407131264e+031,
 -2.71586983800496e+032,
 -4.97904709617111e+033,
 -8.52052408526674e+034,
 -1.37148020128378e+036,
 -2.09209859073969e+037,
 -3.05157815956964e+038,
 -4.31198694357118e+039,
 -6.0516153334399e+040,
 -9.04746158745842e+041,
 -1.83089018411669e+043,
 -7.85792264499214e+045,
 -1497353.00381781,
 723041183.21408,
 798570121442.665,
 179799692875761,
 1.75400882707701e+016,
 5.90805036793554e+017,
 -4.49390336205661e+019,
 -7.12757601895464e+021,
 -4.99886165418825e+023,
 -2.40376807666057e+025,
 -8.92056369816489e+026,
 -2.71249252736027e+028,
 -7.02691152031126e+029,
 -1.59610221811831e+031,
 -3.25079147118034e+032,
 -6.04362056723916e+033,
 -1.04033863298669e+035,
 -1.67721342078464e+036,
 -2.55641247592603e+037,
 -3.71487059097827e+038,
 -5.19156012641781e+039,
 -7.05767849053118e+040,
 -9.54079997480186e+041,
 -1.37682173010483e+043,
 -2.81944608373897e+044,
 -1.35568431120779e+047,
 50890062.0717611,
 14367608874.8046,
 -2441428268667.74,
 -1.1038760824906e+015,
 -1.39748101423647e+017,
 -9.39850732312024e+018,
 -3.91099001869365e+020,
 -1.10065830730325e+022,
 -2.72302043182617e+023,
 -1.10817009473656e+025,
 -6.4123402563569e+026,
 -3.20104135366767e+028,
 -1.26609617693228e+030,
 -4.08925485679342e+031,
 -1.11710339860563e+033,
 -2.65351306408463e+034,
 -5.59921841894524e+035,
 -1.06786709570634e+037,
 -1.86801008380294e+038,
 -3.03739730383572e+039,
 -4.65152896853984e+040,
 -6.8118589248791e+041,
 -9.79544921447061e+042,
 -1.51514184913643e+044,
 -3.61763707621482e+045,
 -2.33635753945665e+048,
 -1632017496.10562,
 110317600316.58,
 330180395951750,
 6.7231021628206e+016,
 3.97828600685863e+018,
 -2.95784859710759e+020,
 -6.98592705045258e+022,
 -6.25659244735647e+024,
 -3.64861751928678e+026,
 -1.5833753360599e+028,
 -5.45880499049355e+029,
 -1.55850222438659e+031,
 -3.79783004095228e+032,
 -8.08967352013915e+033,
 -1.53636362179642e+035,
 -2.64622016202475e+036,
 -4.19608871034953e+037,
 -6.20828443819136e+038,
 -8.67464415143325e+039,
 -1.15730260817286e+041,
 -1.48916975227492e+042,
 -1.86658981818677e+043,
 -2.31349724473074e+044,
 -3.02888770555798e+045,
 -5.99835061959248e+046,
 -4.23733375004766e+049,
 38417214435.0839,
 -726039946748.135,
 -7.31431915372476e+015,
 -1.56745545878788e+018,
 -1.02395774227146e+020,
 4.67055145073008e+021,
 1.3399291038475e+024,
 1.20997806353346e+026,
 6.95795172708743e+027,
 2.95024225282919e+029,
 9.86887837596031e+030,
 2.71349946947451e+032,
 6.30942685965434e+033,
 1.26704450603725e+035,
 2.23282870445361e+036,
 3.49345248902422e+037,
 4.88842278953211e+038,
 6.12750579137622e+039,
 6.82428986341838e+040,
 6.55589778514714e+041,
 4.93603193443648e+042,
 1.70522518557545e+043,
 -3.25119471661806e+044,
 -1.19081432094687e+046,
 -5.38145156298819e+047,
 -9.5337915116242e+050 )
}


# 2 < p < 3 (F)

if ( (power>2) && (power<3) ) {
grid <- c(
0.999320859444135,
-0.00339863996520666,
-0.017405834743935,
-0.094415885413374,
-0.474943122220438,
-2.12217526114901,
-8.28657088560929,
-28.1983730239296,
-83.9866444485264,
-220.714016257597,
-516.840853257854,
-1089.80906886319,
-2091.04653592782,
-3687.84899560698,
-6035.14688398691,
-9245.16287671588,
-13364.8056799522,
-18368.6312623177,
-24170.5284208701,
-30652.0675548239,
-37701.6909476229,
-45257.9298753625,
-53351.6680939871,
-62148.0795745321,
-71944.5580351856,
-92584.6307265222,
0.283341460167403,
1.4198595014972,
7.26323206149768,
39.3470798317175,
197.737872455163,
882.775423162057,
3442.34491552515,
11680.2981292535,
34586.2312661717,
89941.2772399677,
207071.437085722,
425804.06636173,
789124.09455765,
1329886.46550453,
2055949.50037776,
2940437.73707706,
3922432.53173063,
4918772.14066394,
5842643.29870822,
6622148.2558007,
7212994.5161919,
7602717.43947442,
7807280.48778773,
7863250.68212219,
7817746.83340438,
7725050.76243782,
-0.0600764380680103,
-0.392050079124554,
-1.58937633747626,
-6.08664376399068,
-21.2367358176583,
-63.2532442004164,
-151.966790509734,
-262.308768568861,
-171.730813184702,
848.887604105088,
4456.38515028148,
13551.4810243792,
32150.4572347722,
64571.5085421313,
114045.803910683,
181267.784390237,
263555.257000779,
355057.226722376,
447959.727920055,
534210.69005017,
607115.455356772,
662364.32312859,
698095.516091012,
715453.61053991,
709404.020526808,
746938.05963745,
0.0245245775014973,
0.155539966613345,
0.314103723867304,
-0.974507538399082,
-14.3303629251534,
-93.9696323047437,
-452.198879088332,
-1754.67665143518,
-5704.65115146005,
-15896.7164855003,
-38614.4509020089,
-82918.5184841515,
-159371.089358326,
-277309.878496798,
-441509.504984657,
-649717.382851132,
-892380.682977137,
-1154918.2173388,
-1421750.28004424,
-1680648.33849184,
-1926151.40276689,
-2161369.37925071,
-2399069.23604879,
-2658615.89551391,
-2996584.32616847,
-3244309.2211475,
-0.0114761862276732,
-0.0592619852303128,
0.149317524117004,
3.1187422858429,
23.1527784014401,
124.86200951222,
542.858009682394,
1973.11003584607,
6116.61892273702,
16406.537857693,
38528.1580303481,
80007.829774633,
148150.992289522,
246211.713959957,
368704.211125567,
497587.798271825,
601013.469645754,
635274.165995275,
549071.310137826,
288010.646988658,
-203166.790795375,
-982564.60054751,
-2119680.41798567,
-3704741.92691329,
-5931154.48660554,
-8809844.1184798,
0.00575921327240797,
0.0165592415712064,
-0.310297398420674,
-3.50279461877254,
-22.6809711897805,
-113.076858198932,
-465.191104437209,
-1624.40007498999,
-4906.66394460468,
-13042.9156711335,
-31056.3757863278,
-67485.3530473821,
-136357.774862024,
-260555.840068835,
-476923.168849248,
-842508.117841961,
-1440283.41961548,
-2382568.94608435,
-3811813.17259831,
-5900950.25877304,
-8858145.80320564,
-12942579.2179964,
-18500106.3718362,
-26024539.6191355,
-36343654.5689473,
-51760152.8969368,
-0.00300526348269945,
0.00332582644494308,
0.352539219455832,
3.25105798234768,
19.0591425935349,
87.5918992144986,
333.900689176084,
1071.0953680497,
2874.83712283637,
6250.68546381489,
9860.08999845749,
5357.73872252884,
-33652.9231630469,
-170144.676973302,
-528295.164157734,
-1318787.55146601,
-2857009.40475253,
-5568572.63242854,
-9983746.86426784,
-16731892.6976634,
-26555077.5194214,
-40365482.5555169,
-59377234.8688406,
-85360882.1668992,
-121173410.369109,
-177831829.550814,
0.00156864913600111,
-0.0131184236007616,
-0.348599078434079,
-2.77780309250851,
-14.4361149519266,
-58.8320138727439,
-203.569099773387,
-645.03374745931,
-2061.89411085477,
-7067.36215375172,
-25074.4673862501,
-85193.394389473,
-263487.697425328,
-729232.04752068,
-1807739.26844015,
-4046984.43990115,
-8268652.76656352,
-15588975.4539116,
-27412600.48323,
-45427572.5997714,
-71652140.9796545,
-108601600.144071,
-159667231.52998,
-229874415.837102,
-327427911.37296,
-486267279.031203,
-0.000751098308492496,
0.0183607429610277,
0.32560155770289,
2.20414845756873,
9.30362608159624,
27.198466156795,
41.1328701889614,
-147.601617296233,
-1778.48998961713,
-10615.1337836503,
-48121.4524064771,
-180117.812148503,
-576102.516202511,
-1608834.11717014,
-3989140.08790138,
-8909582.40404006,
-18156969.4268232,
-34163633.7942921,
-60005088.1200331,
-99404892.9765117,
-156858144.123099,
-238028512.439587,
-350642272.707905,
-506310271.081071,
-724411047.556449,
-1085756135.90953,
0.00023268753144256,
-0.0214573124116083,
-0.290521160221242,
-1.52935064098677,
-3.59013983096331,
3.92703520664888,
60.9297355136416,
133.554643751598,
-1021.17664549104,
-11590.6683029842,
-66685.879961975,
-281254.304289336,
-960925.358517155,
-2788713.9792639,
-7077864.08967289,
-16041484.6298937,
-33004923.3507173,
-62506118.7968919,
-110299815.326543,
-183381109.539535,
-290240024.52127,
-441647425.450367,
-652423811.084398,
-945103927.138418,
-1358058723.14185,
-2050647607.71471,
0.000141239642755454,
0.0232583731015363,
0.239064459228162,
0.680275724339044,
-3.14242038852734,
-42.3056828515275,
-274.070438394209,
-1492.11348836561,
-7656.84591539787,
-36265.7134048327,
-152004.565173214,
-553517.300322874,
-1753199.12822194,
-4878532.09154446,
-12079077.4876383,
-26959794.0499863,
-54922340.2612081,
-103328174.078604,
-181512150.271794,
-300836965.808012,
-475136824.235409,
-722051548.703595,
-1066029150.12698,
-1544638613.28437,
-2222964632.64711,
-3372045968.09847,
-0.000444860575100238,
-0.0236847905695715,
-0.161058411288511,
0.407895924404698,
10.7548689475349,
73.0538463066281,
281.046556675514,
383.530300181306,
-3421.68161502129,
-33452.5094386606,
-180847.896775903,
-735888.039250597,
-2456210.90539502,
-7010955.43700536,
-17571067.8749036,
-39423967.3829222,
-80439386.6066522,
-151267929.917988,
-265326130.629948,
-438855452.428628,
-691566344.59827,
-1048598947.11475,
-1544960914.68021,
-2234995984.89575,
-3214527134.98847,
-4893933770.89838,
0.000703732261226084,
0.0219231227951077,
0.0425149877267029,
-1.78013554662686,
-18.8024427670158,
-106.878166315772,
-491.309303563418,
-2410.2059871961,
-13085.0171426766,
-67405.3817964102,
-299027.541071522,
-1121165.69781058,
-3591096.47322529,
-9998617.39241346,
-24622492.5304239,
-54479632.7305074,
-109838493.819578,
-204359729.041694,
-354963747.711344,
-581822254.35559,
-909137795.541055,
-1367637645.50103,
-2000281379.29399,
-2874669331.61444,
-4113040594.89626,
-6286557382.99022,
-0.000959178477233191,
-0.0172549225992488,
0.119521088004283,
3.25076886753573,
24.469396043305,
77.1613526002622,
-413.359872941446,
-8197.32770301892,
-68568.1716842124,
-394782.555626955,
-1746273.13505017,
-6261284.87873029,
-18848727.5557636,
-48909488.5602613,
-111724871.937695,
-228668100.175633,
-425755681.178738,
-730883280.355038,
-1170935158.54954,
-1770335830.80924,
-2552404036.30653,
-3544443902.10624,
-4787694880.17175,
-6355795903.3852,
-8395778630.02711,
-11697767237.6242,
-0.000133829436888704,
-0.00402825444711644,
-0.39467858862294,
-4.78902212824088,
-23.5283465754233,
-495.238264820242,
-11644.3071017881,
-143948.075333492,
-1133949.74994485,
-6423704.38453159,
-28106118.540934,
-99563713.8931762,
-295413197.299012,
-753554043.615336,
-1687633184.19881,
-3377063524.57688,
-6129762082.98725,
-10226991209.9182,
-15871393620.5447,
-23160885564.1831,
-32101906559.4075,
-42662039575.6347,
-54852876075.286,
-68836142772.9544,
-85067835018.7585,
-105221455769.378,
-0.0159347723438438,
-0.139677220958585,
-0.501567031055687,
-1.15800906623169,
-0.28145406205746,
-5860.29004848185,
-152251.522259703,
-1873731.64284049,
-14679841.3179456,
-82892307.4814804,
-361969279.219145,
-1280474002.74372,
-3795015931.20803,
-9670673302.28459,
-21636448308.9111,
-43251166651.4325,
-78419385616.2198,
-130679713013.754,
-202537118746.243,
-295128647709.38,
-408391269373.725,
-541723400932.74,
-695004000950.63,
-869822131123.467,
-1070898637903.45,
-1311173239757.92)
}

# 3<=p<4 (G)
if ( (power>=3) && (power<4) ) {
grid <- c(
0.807045029092288,
-0.742817271116966,
-3.33628159647649,
-17.1892563957248,
-83.606924206003,
-365.071192253248,
-1400.71893678865,
-4694.58855569616,
-13766.4944192918,
-35518.0066568666,
-81236.2491747527,
-166111.053093699,
-306340.130694687,
-514013.401873071,
-791488.479479209,
-1127820.45893735,
-1499215.00485969,
-1873687.95856514,
-2218215.22065063,
-2505773.30210073,
-2720001.08609826,
-2856885.35011282,
-2921562.04807504,
-2942860.2782209,
-2724727.65113592,
-6393599.83992389,
-0.136539772699055,
-0.3712773738108,
-1.34060216897536,
-6.31281195417687,
-28.7743296376622,
-120.000150401644,
-444.643043043044,
-1450.18256577969,
-4160.32884893442,
-10542.0837079396,
-23750.546686822,
-47945.825282172,
-87450.1861981138,
-145332.756352579,
-221913.780561476,
-313883.835136119,
-414535.3380698,
-515112.581103503,
-606768.028120295,
-682457.909013515,
-738016.337582899,
-773197.016326358,
-785451.127656309,
-832107.719453664,
-65576.1860675511,
-16681644.1359413,
0.0390190571539757,
0.217823357411596,
1.05307617320898,
5.46628587502473,
26.7079674195229,
116.777871575447,
448.243317291106,
1502.23344046383,
4404.21226928853,
11360.1463323601,
25977.6589655352,
53114.6692363956,
97962.628202229,
164425.342777476,
253336.682488931,
361326.799670235,
480960.413361712,
602200.055474944,
714666.469016104,
809833.902052876,
882644.532433897,
931195.500052511,
963253.25246436,
909728.615872677,
2118955.81118029,
-28223901.1772564,
-0.0135344333535702,
-0.105123408597516,
-0.597134560351953,
-3.30138792730478,
-16.7678509371076,
-75.2601777869174,
-294.380322822538,
-1000.62051658493,
-2965.26791173967,
-7710.19873496618,
-17729.923018914,
-36365.1109880616,
-67103.1049888601,
-112342.914312332,
-172033.611258013,
-242814.982348664,
-318156.29116042,
-389547.350581509,
-448309.499169329,
-487358.607042042,
-502288.280654725,
-491896.568009773,
-454167.367329454,
-436702.595676197,
592554.019288957,
-28175511.1369623,
0.00516979878298409,
0.049186369552098,
0.313361035619778,
1.82380477698879,
9.55014857006874,
43.7632364583752,
173.85444871189,
599.000056509861,
1801.41110174687,
4775.62836370548,
11292.324476147,
24116.159116646,
47097.8558813758,
85097.7219846563,
143732.570485074,
228950.891134358,
346463.780768116,
501132.836280051,
696519.783935997,
934868.31070775,
1217737.35915428,
1547476.14417311,
1929941.85114461,
2371249.26369004,
3259938.37856758,
-10481780.2845814,
-0.00210569527115553,
-0.0231405024519351,
-0.161238923976278,
-0.980926414519322,
-5.27476113570811,
-24.596921426023,
-98.4679251076239,
-335.875641752254,
-966.009103975749,
-2295.01927050472,
-4285.53743080301,
-5305.54520561608,
298.828932844618,
25527.1544820149,
93698.0628522327,
239219.913855273,
504807.80240039,
935434.395727312,
1571191.70516859,
2442224.38428192,
3568464.03604572,
4965686.14710397,
6658844.98464334,
8715913.57478651,
11333531.2576485,
12135648.0516836,
0.000904992444600065,
0.0111168052401889,
0.0833664125811974,
0.527134997928565,
2.89808864704583,
13.7770488849601,
57.6536611318659,
222.31978415136,
830.973533709237,
3052.29525500565,
10683.2000738177,
34231.1489571498,
98080.965980282,
249787.613056241,
567853.696586215,
1162730.51664402,
2167379.38902465,
3718634.07876582,
5936186.98691647,
8908850.23901353,
12696377.3865465,
17351472.7681057,
22966767.528756,
29776590.0504127,
38478019.7505486,
50867749.5850808,
-0.000409107208596102,
-0.00550111018040516,
-0.0439376208280782,
-0.288505996551621,
-1.63212078318116,
-7.75520903054617,
-27.7954092957704,
-45.8656633398937,
259.464351955424,
2854.74539011967,
15714.5555070392,
63230.6123114516,
204776.271508092,
558927.74579377,
1323826.85797549,
2780007.90439764,
5266252.25523076,
9131789.73343168,
14682560.3322589,
22146472.0731672,
31679055.6650164,
43421701.98937,
57625070.1381395,
74910568.9517981,
97163770.0643605,
127754334.17767,
0.000194307002727792,
0.00281756674493309,
0.0237136332062326,
0.159405381495387,
0.898230867797415,
4.60345984485858,
27.6400397896721,
205.661534317615,
1448.56631969149,
8291.67719559354,
37903.115045638,
141316.649850398,
441087.624260472,
1180252.50521153,
2763322.5415243,
5761704.22812121,
10864619.9821838,
18782495.180654,
30138958.1231858,
45402276.5830579,
64900090.6565401,
88943149.521878,
118085299.880759,
153669984.92336,
199638879.192562,
263053537.257288,
-9.6813970562699e-005,
-0.00149818930643216,
-0.0133055656398246,
-0.0945115298691158,
-0.579908914141087,
-2.47603853062905,
4.13987484493922,
197.011666853647,
1993.36630427103,
12842.763799835,
61521.7609859012,
234369.520358223,
739791.084258406,
1992499.07342665,
4684730.93725107,
9797208.48563922,
18517157.6071782,
32074859.2918227,
51560378.2551523,
77807929.4600213,
111424466.116615,
153009280.153805,
203615359.468966,
265712174.566699,
346298256.693629,
458843429.853769,
5.0510737242787e-005,
0.000825454278033001,
0.00754236535471132,
0.0518083089767653,
0.274621994126378,
1.98114382300995,
29.5576814444073,
370.652270203056,
3159.13901640118,
19407.1216982271,
91464.1099352158,
346424.449467708,
1091742.85694851,
2941228.74414017,
6923668.63289443,
14504109.9212339,
27468203.8817623,
47684788.4333904,
76836997.6005713,
116252578.587634,
166950731.275237,
229981255.740226,
307139843.974733,
402454023.216751,
526912234.233937,
702818878.857313,
-2.7418070102492e-005,
-0.000471538370173029,
-0.0046165376568393,
-0.036587699748026,
-0.259306297425945,
-0.513545608507554,
24.5167655002064,
422.079080450406,
3864.47708330576,
24364.9020805635,
116380.447641915,
444742.438017592,
1411324.94709018,
3824627.71731913,
9050677.77134831,
19052334.7104253,
36248377.3563104,
63208294.4047604,
102300026.230609,
155465995.050065,
224290351.769315,
310469434.35741,
416810360.956061,
549298562.77482,
723618453.467978,
973334765.794742,
1.5473116821586e-005,
0.000276604798913003,
0.00264969418261052,
0.0174900279365758,
0.0766314928362542,
1.35761435509042,
35.91128539334,
507.476924676718,
4502.64745318921,
28323.9262461144,
135970.357806147,
523202.778300938,
1672239.7082978,
4563033.2399971,
10868160.2019504,
23016985.4182579,
44040524.526026,
77210463.32594,
125613056.798018,
191876531.064486,
278262231.602734,
387270853.524674,
522928872.391834,
693428189.404835,
919453049.166152,
1247727868.46718,
-8.8935422652274e-006,
-0.000166030604502527,
-0.00180320520193969,
-0.0158064250790231,
-0.112690635415925,
0.605792100613844,
35.8920609726225,
549.294618680751,
5007.87803105007,
32039.3397103949,
155789.404879026,
605600.459856101,
1951459.22132761,
5359784.8952449,
12832094.7037116,
27287162.7391152,
52377387.4480437,
92054947.2665312,
150060697.345423,
229605705.993322,
333504156.464628,
464936857.308676,
629016062.372982,
835873005.440438,
1110473643.66999,
1511045682.99279,
5.2836771853153e-006,
0.000101503163727336,
0.000993429510955962,
0.00768536258973886,
0.0762004564807366,
2.5289608455215,
61.1254255902454,
833.768195052616,
7354.03786958471,
46348.1408509856,
223057.664768771,
859057.342688984,
2741866.57069041,
7454083.6797878,
17650308.942397,
37089910.5528786,
70295592.6389546,
121895207.707387,
195915436.6333,
295398279.000812,
422647796.67442,
580244039.797793,
772834327.8546,
1010119774.0892,
1316484399.07765,
1741035939.30733,
-3.064377771159e-006,
-5.980922030261e-005,
-0.000636343983276073,
-0.00078615973385692,
0.104702137517939,
4.49267742988801,
99.9926447089914,
1292.20378864764,
10963.7639617,
66756.6908405897,
310381.696918354,
1152600.74257852,
3536665.37191601,
9209232.17803997,
20794869.4119386,
41458235.608168,
74105805.4180075,
120360968.097685,
179752435.234823,
249520137.896054,
325147964.097433,
401275805.233735,
472066391.64127,
529481285.471072,
561635505.992595,
419706664.414131)
}

# 4 <= p < 5 (H)
if ( (power>=4) && (power<5) ) {
grid <- c(
0.699925346022531,
-0.967077543569014,
-4.01132613691835,
-20.1399473550733,
-96.2422051602767,
-415.365980035549,
-1580.21518068211,
-5262.33523148009,
-15354.0312963513,
-39454.0262740431,
-89938.4982463926,
-183391.698932919,
-337402.996179853,
-564968.193375439,
-868383.335714499,
-1235426.26432415,
-1639939.91489061,
-2046986.48980844,
-2420647.49346844,
-2731722.04945314,
-2962469.78638038,
-3110166.85555429,
-3167407.14273507,
-3347866.5224723,
220066.29635677,
-101623416.405294,
-0.0856827146015816,
-0.129194389632916,
-0.275741081827042,
-0.995185808644477,
-3.44304072656882,
-11.1738124760052,
-32.3161811060365,
-82.0820212947264,
-181.186657587385,
-345.1190334895,
-560.487576581505,
-755.308564356227,
-780.02789704195,
-420.782646552896,
549.203902376783,
2292.75919352824,
4844.65304486468,
8078.19818210918,
11752.1076384033,
15489.0074612448,
19295.1080639715,
20576.2965127781,
43458.6905027808,
-272065.299693591,
7025098.66546232,
-233107902.179985,
0.0168969338431571,
0.0658060047775105,
0.252941889440879,
1.1930373881614,
5.46017018978347,
22.7950867197324,
84.5102486368491,
275.66150811094,
790.803379055727,
2003.61163553412,
4513.25595683332,
9109.342888956,
16611.6941580461,
27601.4176540827,
42137.2784788037,
59588.6192203015,
78680.9563361091,
97750.8199011234,
115126.506416694,
129453.386054272,
140101.469157844,
145829.366024648,
160224.556969174,
-43787.4381679082,
5529231.86708544,
-215530714.274696,
-0.0040335782993408,
-0.0234817794259036,
-0.110587472703672,
-0.561343113249565,
-2.70232831818259,
-11.6848904632598,
-44.472343214546,
-148.047891896767,
-431.69906130895,
-1108.54940295543,
-2525.46917566759,
-5147.26241339181,
-9467.85961656983,
-15855.3598573824,
-24383.2818096163,
-34724.30268989,
-46163.6433390437,
-57734.3395387668,
-68421.6708128944,
-77354.6697154979,
-83913.2169220258,
-87816.2036212983,
-87020.6389051429,
-128645.623351391,
1605434.61145123,
-93593764.193113,
0.00106368144379509,
0.0078728133185123,
0.0423052301671125,
0.226374767743507,
1.12595759435215,
4.97840665743307,
19.2538811192632,
64.8572506879996,
190.743003851787,
492.585067098186,
1125.23986845133,
2291.87964977926,
4195.71676979692,
6958.1794888724,
10533.1890357286,
14663.2366154889,
18914.8130296304,
22801.0859452267,
25968.0905843622,
28413.4033257206,
30721.5445983317,
34434.2563191554,
41999.1412856565,
63751.8453401499,
142903.993721121,
-14773309.479337,
-0.0002988416478197,
-0.0026128755213848,
-0.0155064767814191,
-0.0866816879592185,
-0.442392050536327,
-1.99055510244982,
-7.80955850461462,
-26.7322098215555,
-80.4281213186107,
-215.042835786294,
-517.452034952092,
-1134.82818218957,
-2293.86917654107,
-4310.82691914404,
-7575.43567664562,
-12486.4405956532,
-19320.884708294,
-28032.1935354864,
-37976.4690140615,
-47537.7814603878,
-53537.2289629646,
-50176.3090563924,
-27111.3184862328,
35023.1543546029,
97772.386910872,
1945675.4370755,
8.8111852916931e-005,
0.000874298083941498,
0.00562299485434772,
0.0327030242649852,
0.171726872656,
0.790684116809639,
3.1130199548177,
10.1710449772459,
26.346202329904,
48.3548936424336,
31.8771254839422,
-186.395662514681,
-1016.84516919046,
-3258.27728769355,
-8162.68516490687,
-17312.2380768225,
-32239.1706818026,
-53802.4297241544,
-81382.1575234866,
-111877.895287325,
-138254.603738189,
-146999.311848418,
-113518.232821567,
3726.73086187838,
273817.624118751,
757199.553621758,
-2.7085879774208e-005,
-0.000297498761276123,
-0.00204245827241299,
-0.0121381151934212,
-0.061733063855247,
-0.26699858234589,
-1.1122126687644,
-5.31569083447484,
-27.8731682359897,
-134.874202708192,
-553.327324835573,
-1903.59891170908,
-5566.04749679175,
-14086.3410351293,
-31399.2007566238,
-62614.6499198601,
-113234.321754671,
-187884.675320979,
-288785.120075029,
-414032.58145669,
-555329.364059289,
-694159.031925558,
-795349.360834586,
-800767.324268166,
-634936.677499812,
-1390026.09670781,
8.6547130615233e-006,
0.000103571257343756,
0.000763282155189074,
0.00503495689305218,
0.0339503192962971,
0.201182583372745,
0.679200600696824,
-1.40781922364333,
-34.7607658472872,
-251.167735291424,
-1214.87163835565,
-4539.92881673297,
-13946.9496782576,
-36547.4380325321,
-83834.1134982485,
-171726.961316446,
-319298.861981019,
-546376.82435059,
-870690.377723771,
-1305228.69338149,
-1855917.38613467,
-2519159.50584102,
-3280158.49857154,
-4124134.63935406,
-5101301.89955304,
-8040658.58499969,
-2.863420992543e-006,
-3.6705711047603e-005,
-0.000267736627702962,
-0.00118987631537025,
0.00500374982239699,
0.0906349021665064,
-0.03964765353084,
-8.86995948576699,
-96.4179103903694,
-622.102159599423,
-2925.97709670939,
-10897.5978947285,
-33690.7023380288,
-89258.4433448129,
-207545.201776322,
-431710.061827054,
-816252.874623963,
-1422300.99639821,
-2311804.80050625,
-3543091.25767127,
-5169984.23754927,
-7246834.8839518,
-9846927.88218808,
-13128306.876303,
-17545564.4897858,
-25760072.0632545,
9.862455038939e-007,
1.3687405837123e-005,
0.000130756634353942,
0.00161205842445782,
0.0260978720622342,
0.205692870824142,
-0.217742919739228,
-19.0670830309275,
-201.334161636633,
-1297.1205613785,
-6144.28137107098,
-23122.4216546073,
-72314.544985341,
-193853.872583445,
-455974.877979017,
-958991.306803016,
-1832262.21329924,
-3224343.55883332,
-5290401.34780242,
-8183538.60025336,
-12057410.3744032,
-17089089.6467557,
-23542978.7525315,
-31949315.2967005,
-43603277.3073714,
-61461357.0483208,
-3.3898916523125e-007,
-4.553631043329e-006,
-6.3146792385749e-006,
0.000996685778798012,
0.0283638961892905,
0.205028928564685,
-1.5209878776305,
-40.9193475021547,
-399.642312169516,
-2535.0167059544,
-12007.7742204405,
-45413.275846636,
-142947.623599072,
-385712.615552065,
-912728.368992133,
-1929683.59832148,
-3702976.31241841,
-6539242.01547221,
-10759257.0350438,
-16681448.9101165,
-24632750.0692673,
-35008236.3642549,
-48422810.8521226,
-66089231.9656645,
-90790259.3017384,
-124782634.831774,
1.3670317399056e-007,
2.8922795451404e-006,
8.0311153292009e-005,
0.00209218748382464,
0.042288144829039,
0.224071380381849,
-3.88595477558196,
-80.3828018669491,
-756.261098584473,
-4760.28166814639,
-22556.6461988553,
-85542.2228609317,
-270094.51932018,
-730737.212456536,
-1732553.94919498,
-3667198.87154186,
-7040018.87161964,
-12429346.6875643,
-20436837.3840823,
-31660964.6132857,
-46731285.933596,
-66447851.8147418,
-92106261.376575,
-126229910.642334,
-174276152.679859,
-236992787.643919,
-2.9515402119269e-008,
1.4230937296803e-006,
0.000121964247166344,
0.00353296955768749,
0.065444187405576,
0.228124670269768,
-8.75709483474552,
-158.883665823477,
-1458.29696267044,
-9116.76392174209,
-43122.5389213505,
-163475.6881624,
-516064.788911388,
-1395522.20150944,
-3305660.01902088,
-6987365.36000612,
-13391031.6773161,
-23598208.7876355,
-38732820.725386,
-59928676.182718,
-88431093.3806788,
-125922659.671926,
-175220715.536157,
-241690986.941136,
-336090539.887824,
-458661737.622016,
4.2195953578419e-008,
4.7262596571073e-006,
0.000268409542694463,
0.00711771507088045,
0.117709590133548,
0.185527243324104,
-20.5126719620795,
-341.724384735611,
-3067.30039938492,
-18999.4745727416,
-89434.205638909,
-338004.516140821,
-1064613.61087423,
-2873576.28161018,
-6796130.72493253,
-14346593.3785094,
-27467963.4091078,
-48381210.7295654,
-79426953.5719077,
-123043357.109127,
-182047880.504625,
-260412755.530795,
-364810110.446681,
-507417982.521783,
-710389534.723263,
-981948479.261186,
3.2858256896933e-008,
8.6262465133508e-006,
0.000619619047250753,
0.0177352317421679,
0.286522133653506,
0.930762013518072,
-34.6477998067267,
-621.723655825581,
-5696.25104328862,
-35724.7343454217,
-170023.838772441,
-649865.670726276,
-2071812.51834388,
-5666452.90654193,
-13596350.0025411,
-29159538.9365487,
-56804496.6608088,
-101967267.133786,
-170892032.312482,
-270727256.793951,
-410267304.386682,
-601763640.206147,
-864300376.261109,
-1229059105.40151,
-1743431919.71737,
-2419331563.0187)
}

# 5 <= p < 7 (I)
if ( (power>=5) && (power<7) ) {
grid <- c(
0.786434222068313,
-0.346012079815851,
-0.111403879165884,
-0.37684834413327,
-0.227942819345689,
-0.748816199446288,
-0.429657425382675,
-1.74160208112774,
-0.380029940723467,
-4.75431964095465,
3.04669074354556,
-18.7254127758988,
42.2595203965277,
-127.075610501907,
597.354383371277,
-1783.48701130305,
11845.8231537032,
-56908.6196171327,
443014.347945159,
-3727125.91186479,
41915265.5772462,
-613145856.515179,
12314442848.2845,
-353565282768.735,
15056060852441.6,
-935512014119462,
-0.0372367474176292,
-0.0300330290598776,
0.0194320567856015,
-0.0286632698352744,
0.0445003394042309,
-0.060948400702475,
0.132495375767898,
-0.213544390269728,
0.423844494115997,
-1.61814802662431,
-1.11118707905738,
-26.5227940118242,
-67.7293316745046,
-425.458598662828,
-1001.46373152249,
-5317.32161731875,
-6139.36538655163,
-71044.3658618384,
194667.687528656,
-2788059.62676948,
30424957.0834719,
-488372280.977764,
10689386756.2126,
-340771077349.598,
16472908240969.7,
-1.18753639215006e+015,
0.00272199211678462,
0.00600457935696839,
0.000771571810584139,
0.00448332513732787,
0.00189163733536341,
0.00589886871984612,
0.00239247793636055,
-0.0286338123481378,
-0.315494675308072,
-2.47958187250564,
-14.5308793302688,
-72.977265238482,
-303.738213023504,
-1105.5766522905,
-3387.71681429456,
-9110.18620747135,
-18243.3848321049,
-29784.0538288672,
75227.02534185,
-48264.5575731943,
8222854.19107454,
-103843212.592982,
2899046709.27306,
-108457074214.473,
6391125630551.74,
-572642997373694,
-0.000235725585525655,
-0.000939712997116092,
-0.000628940476056143,
-0.000732868263909586,
-0.00123615938011299,
0.000322710464360014,
0.0114320282205539,
0.137507245882712,
1.1142446434786,
8.13970036772955,
52.4668442088235,
302.39429825616,
1551.75362597238,
7136.34444891635,
29595.0973156719,
111934.511057314,
390212.830424631,
1272109.66953064,
3921732.69419177,
11664373.480828,
33562473.6192614,
96275902.7551133,
386971758.390603,
-7779194558.35237,
869116980003.18,
-120470494442105,
2.2131289395219e-005,
0.000143553411002132,
0.000170793457061453,
0.000168846275141174,
0.000879434613056662,
0.00980484450010389,
0.12801774812558,
1.33913284198253,
11.5110566068305,
82.8906614784825,
510.178293855477,
2728.20581284727,
12861.2337022356,
54151.1413496475,
206133.384228595,
717759.486483854,
2312738.60846531,
6977878.20433018,
19960158.1406384,
54872807.7633121,
147027968.000616,
392467947.201468,
1007010930.32944,
4277401271.27259,
-32175521793.3006,
-5944774644892.04,
-2.1158714077244e-006,
-2.2245962793571e-005,
-3.7048515147495e-005,
2.6733273608329e-005,
0.00159538838855198,
0.0313650954321681,
0.434832189395928,
4.6653031612727,
40.2540676573929,
287.502036146911,
1739.44171649193,
9089.28709664162,
41709.0038241202,
170553.998156786,
629707.080488873,
2124898.97405362,
6629823.09908051,
19347946.9909511,
53449119.9055938,
141603511.259636,
365139375.94298,
931639294.908808,
2399554030.90057,
6365379770.8051,
6629077555.34454,
1325257347799.35,
1.8706218898433e-007,
3.5595867805602e-006,
1.0461142077657e-005,
0.000125822395739863,
0.00343944847191744,
0.0663161997681272,
0.932504916856121,
10.0628718510742,
86.680128866696,
614.451310218334,
3674.06045960846,
18916.695687826,
85354.4120703696,
342720.990999314,
1241358.93947498,
4106484.2634202,
12551665.0365457,
35852478.1688061,
96827356.785738,
250430599.497147,
629495711.021131,
1564844953.08431,
3922964834.69878,
10060076758.7741,
28625631416.721,
59232587864.0314,
-1.0736419206779e-008,
-5.5682318067157e-007,
1.521132088195e-006,
0.000169816881450556,
0.00513471228094518,
0.101321240932149,
1.43787147451134,
15.5501270999697,
133.573873042078,
940.8573571264,
5575.65420221155,
28399.6411897312,
126608.618970022,
501858.833477652,
1793488.15346414,
5851297.04271422,
17631543.2949301,
49624288.1787248,
131965198.165583,
335793721.545266,
829834387.493952,
2027957061.72163,
5001954036.02077,
12663935117.7047,
32689069667.4124,
18114334235.1189,
-1.2346951657342e-009,
1.2206227284162e-007,
3.9297895016142e-006,
0.000197419028183901,
0.00602851603167602,
0.120332460066855,
1.71628779074478,
18.5634802061232,
158.929000041028,
1113.04191204279,
6547.27067245127,
33064.4777866304,
146043.954634495,
573317.997343912,
2028724.86709794,
6553279.09218201,
19550737.044832,
54474555.29211,
143383996.095597,
361033832.339612,
882831234.550636,
2136130396.56957,
5224797761.40674,
13104715569.8387,
32786630095.1341,
37594565903.4031,
7.6808927594202e-010,
3.8884985498993e-009,
3.1693128416737e-006,
0.000185283673344334,
0.00577829446251914,
0.116358132374561,
1.6638196322789,
17.9690172702516,
153.197902136379,
1066.5313667415,
6229.28866525285,
31215.0838976393,
136765.982020265,
532563.383330351,
1869646.06197352,
5993641.65402298,
17752031.0930607,
49122405.022279,
128443983.218799,
321382731.452622,
781396537.897085,
1882476272.22521,
4594300256.1123,
11490195036.7978,
27999004292.5526,
22293146654.9453,
-2.3198818654552e-010,
1.7768432128195e-008,
2.5170398599975e-006,
0.000147600250383574,
0.0046726297592754,
0.0946608920777015,
1.35378607562751,
14.5687692163386,
123.471744448172,
853.19624444289,
4941.88962460754,
24548.7134403604,
106623.582662544,
411715.309842279,
1434104.29957124,
4564806.01096974,
13434953.6203005,
36971805.6020912,
96215692.5085437,
239814563.118944,
581599219.830851,
1400939948.6891,
3430310154.31951,
8605849617.10048,
20412486855.5244,
8789001626.70916,
7.1680465084153e-011,
8.2425411536541e-009,
1.643130414875e-006,
0.000101010537904274,
0.00324768919113633,
0.0660426805148569,
0.941790746358288,
10.06409355469,
84.4700029565598,
577.065123627681,
3301.27606364533,
16190.5323305796,
69436.106953296,
264898.991604471,
912482.279496566,
2875705.99281888,
8390894.88493373,
22923038.4921442,
59301250.3489504,
147168407.593844,
356301731.33071,
860860680.309354,
2129050981.40544,
5406333032.3096,
12477283740.0382,
-626655436.574466,
-1.6592791250243e-011,
4.0243341048155e-009,
9.144573335573e-007,
5.9220160549e-005,
0.00193782694015436,
0.039445147224364,
0.55777780092105,
5.874590816796,
48.3891511020367,
323.416138439675,
1806.1076174669,
8634.35391245784,
36073.4833351782,
134086.227392197,
450406.616132974,
1386158.30589757,
3956483.9824662,
10592449.579616,
26908399.3711585,
65779957.3387473,
157925265.23842,
383776594.467057,
976109761.141651,
2581703121.52396,
5773630017.85034,
-4997652480.24311,
3.0254155247433e-012,
6.480434669686e-010,
3.8358388175966e-007,
2.805883798994e-005,
0.000946081848661616,
0.0191810731660941,
0.265048170865662,
2.69054536709911,
21.1078472872048,
132.823226564838,
690.014628342345,
3029.32274273717,
11460.2837105301,
37977.6328286557,
111736.598497042,
294820.966835864,
701306.659651741,
1501566.94543448,
2865072.81069675,
4832518.62805497,
7836886.63058089,
19150010.6873593,
89185710.2082667,
421423170.310072,
903955723.061903,
-4583125093.25688,
-4.6611488044825e-012,
-1.238013643329e-009,
6.0868733097722e-008,
8.742211953457e-006,
0.000328731139452505,
0.00659863314758232,
0.0844309956857507,
0.742711442123046,
4.57966556822816,
18.4867246619517,
26.504351580449,
-278.749024654208,
-3001.65994744203,
-18452.141557694,
-87884.1942798782,
-354357.650582793,
-1264077.22326902,
-4102494.78750277,
-12352555.6345389,
-34956221.4649815,
-93464824.2667665,
-234738633.43955,
-541385094.991335,
-1108582688.25269,
-2320229645.40791,
-923285695.803153,
-3.6747431215211e-012,
-1.4003818727403e-009,
3.2110392985838e-008,
8.0692472963221e-006,
0.000344137055565625,
0.0076331695212749,
0.108809443680536,
1.09803763922457,
8.2754406180277,
47.8674513539372,
212.505496482315,
684.095100622171,
1139.77502759503,
-3574.08105308545,
-43796.4476172793,
-252629.656994052,
-1124806.21665939,
-4304377.68043652,
-14822703.8520193,
-47098602.6891216,
-140045686.50728,
-391792106.570437,
-1027557566.96914,
-2504246941.83045,
-5877038428.29484,
-2114725574.7047)
}

# 7 <= p <= 10 (J)
if ( (power>=7) && (power<=10) ) {
grid <- c(
0.77536534618895,
-0.430563489144848,
0.248933090844648,
-0.776691071939754,
1.42435320848802,
-3.96378215199569,
10.6082763732752,
-32.6646389815812,
107.594392018696,
-385.875807619758,
1554.08106263416,
-6562.13654153062,
33311.2180271377,
-175799.39356045,
1100307.32618603,
-8076936.20134624,
64818006.8415247,
-672272627.675492,
7951846880.99446,
-121489926465.507,
2346719230615.01,
-60597733349803.8,
2.17252324837994e+015,
-1.13969296410995e+017,
9.19134050303596e+018,
-1.13643544485344e+021,
-0.0228287393493532,
-0.0194681344424241,
0.0389473545390899,
-0.084197161744077,
0.22422473111649,
-0.619022478530716,
1.92391840210546,
-5.54851462125942,
31.3124540235083,
32.9263463638783,
1526.62316779162,
8457.00783887994,
82602.3748106012,
435383.429491403,
3109338.73107016,
12375814.1074866,
91668768.6318141,
102755151.774348,
4154874222.20171,
-39332097778.8331,
956298556688.97,
-25984409388026.1,
1.02539880527099e+015,
-6.00621364640884e+016,
5.55452727721829e+018,
-8.10234314733144e+020,
0.000990449241892617,
0.00245658725606492,
-0.00243739393253416,
0.00484724481040902,
-0.00926456300347046,
0.0281306359896192,
0.156239298430668,
4.17858979355244,
62.7520656235757,
780.015382056436,
7975.4456209287,
68984.4386566011,
515762.532144664,
3384991.53977911,
19858782.3438055,
105265045.358707,
514394379.656016,
2317182902.4062,
10175938926.1775,
37656980976.6547,
252620644930.558,
-2174477649470.25,
136934481128099,
-9.59053217215135e+015,
1.11722695207039e+018,
-2.10701957364907e+020,
-4.739402348122e-005,
-0.000231852992426791,
7.5710640928626e-005,
-0.000134980675522712,
0.000354522998021474,
0.0192668478637093,
0.520510022906938,
10.7916294399501,
169.819567472824,
2141.74949227175,
22251.8592082094,
195421.11468649,
1479886.25999519,
9836509.12146506,
58265320.1095556,
312023644.738493,
1530894452.40315,
6980944196.86386,
29986172284.6235,
123833899892.611,
492927274255.622,
2078899780231.5,
6850572520142.16,
-127562576321897,
5.83782301712251e+016,
-2.12174087566134e+019,
2.1161776328178e-006,
2.0906300316289e-005,
7.2258740877989e-006,
-5.5308991127345e-006,
0.000518091195594951,
0.0211849717747649,
0.651125893562394,
14.0232486435453,
229.888540942629,
2996.41818368473,
32119.513737,
290589.577507762,
2266622.81206172,
15520120.3234281,
94786645.0963208,
523940347.205171,
2658194740.25114,
12552187877.3901,
55962360223.2959,
239035772367.959,
991951859205.356,
4061174998315.7,
15708445458792.5,
118144223622463,
-5.07332493452949e+015,
6.4479142073582e+016,
-5.2809354056961e-008,
-1.8733028739588e-006,
-2.173112075266e-006,
-1.4999243577998e-006,
-1.6130806107812e-005,
0.0062857185977922,
0.291729587926022,
7.72340499458176,
144.463644174106,
2082.77967708294,
24311.1198183042,
237511.034262544,
1991173.97650715,
14616262.9410656,
95580256.036951,
565523735.27937,
3071700431.36731,
15529390455.3399,
74043372155.9458,
336948400501.14,
1477696012185.8,
6275079442094.89,
25790455939029,
96328050938412.8,
190447249714597,
1.00992445696896e+017,
-6.18972021646e-009,
1.6716212953349e-007,
1.819435975603e-007,
-1.4425153664549e-005,
-0.000610612063278709,
-0.0167634347715261,
-0.311332599080969,
-3.94138325688247,
-29.7238495417058,
6.2555652775239,
4042.6061664409,
72350.5669137632,
851581.13806843,
7901256.5282858,
61774104.4395478,
422241976.017542,
2587726006.34361,
14498413734.2669,
75443525203.5096,
369254749083.404,
1715435173133.93,
7593311590379.09,
31824142191324.9,
121720773218198,
396847554120315,
-4.20744529251827e+015,
1.6263773828611e-009,
-1.5358949677752e-008,
-2.4571305966315e-007,
-1.8428794010876e-005,
-0.000937415911513019,
-0.0306800211195246,
-0.708070117922537,
-12.1299638962403,
-159.03668317182,
-1618.62979619491,
-12693.3036944218,
-72284.2761884275,
-217166.845678661,
1044589.95460791,
23105099.1246563,
228679026.770271,
1721827459.30268,
11024395805.3666,
62955703260.573,
329358475506.641,
1604477935025.79,
7336606897800.17,
31379245002770.4,
121522052593940,
347347465989858,
370932937610640,
-2.6808728840779e-010,
7.9020135767862e-010,
-1.5514883197396e-007,
-1.6276332336879e-005,
-0.000876562698279347,
-0.0303398354210813,
-0.738225701613061,
-13.3408979363142,
-185.677983898719,
-2038.39425331329,
-17876.6499694802,
-124905.284696409,
-670434.749401926,
-2347535.07095901,
618376.316323485,
94377618.194901,
987830337.872776,
7299222860.15746,
45155167569.2116,
248286662532.449,
1249349399406.29,
5836370519195.45,
25333886840723.4,
99430753191700.9,
299121530322441,
625102006387724,
4.0255310414341e-011,
-4.0075628102411e-010,
-1.0747626938796e-007,
-1.1032436842021e-005,
-0.000614076074677421,
-0.0218656078360402,
-0.54523287925633,
-10.0746614750899,
-143.319283331831,
-1611.54363682135,
-14561.2179931157,
-106231.861503575,
-615725.868813294,
-2629632.88870542,
-5405277.88267556,
35284319.2277274,
545097935.998791,
4473496821.19931,
29049471850.145,
164111610940.043,
839519077635.141,
3964583399582.53,
17362504452757.4,
69068653702018.9,
218680025952334,
491443538772321,
-5.5715190324839e-012,
-1.6407435736789e-010,
-5.7888379466468e-008,
-6.1468686512493e-006,
-0.000349575517522456,
-0.0126681684563627,
-0.320310340444328,
-5.984181303379,
-85.8933670684766,
-973.163367650677,
-8853.94643946026,
-65051.8256679697,
-380371.836571965,
-1650092.93269853,
-3630335.63566798,
19513096.3150577,
321805563.17688,
2676104393.70507,
17469092428.2355,
98920174459.9976,
506594361070.731,
2395294913189.68,
10527983468557,
42415349208286.2,
142117239446790,
424957076182122,
9.0943122447293e-013,
-7.6362405646957e-011,
-2.713727877671e-008,
-2.928550895988e-006,
-0.000169272201914165,
-0.00620868754605586,
-0.158241844155348,
-2.96853021751442,
-42.6197249239258,
-480.920925303812,
-4332.65386909201,
-31222.9264278609,
-175486.089577641,
-686688.578405898,
-700984.246915129,
17606789.7093272,
210224332.113958,
1618919804.83086,
10205343393.9837,
56638977400.7164,
286412959136.382,
1344441166878.75,
5905112905298.04,
24102001115608.5,
86044212091524.3,
330694952809613,
-7.2206827308685e-014,
-2.9539633014077e-011,
-1.1047488592785e-008,
-1.2213639725909e-006,
-7.167213125224e-005,
-0.00265286298445461,
-0.0678318219234301,
-1.2685871768583,
-18.0238255068467,
-199.279184183516,
-1731.66032795093,
-11672.9085946704,
-56607.411359395,
-125551.570376765,
1014540.31262727,
16556288.6301798,
145428947.733874,
1003290094.11754,
5972683400.78656,
32007315379.0508,
158215899080.749,
732613701543.087,
3206486334222.12,
13279647471500.9,
50798414623688.9,
237038150364079,
4.3609819166476e-014,
-8.2019208583295e-012,
-3.8681869833068e-009,
-4.4736540762345e-007,
-2.681092767123e-005,
-0.00100096482939058,
-0.02552849391404,
-0.470191072698517,
-6.46906795194434,
-67.4507414346289,
-524.848782738057,
-2753.53939768566,
-4110.14385245369,
104460.727972183,
1538900.68620935,
14024098.498795,
101919977.04296,
638471700.403706,
3583252863.16384,
18458609364.3726,
88808087125.7936,
404428978510.69,
1762294434820.44,
7420587100023.73,
30447350490820.6,
161257587722473,
1.102927814593e-014,
-2.9833499013762e-013,
-1.0237201386522e-009,
-1.357500505407e-007,
-8.5194681661709e-006,
-0.000321263444569207,
-0.00803242863440226,
-0.139760568118627,
-1.70845294566464,
-13.7392013834384,
-43.3312232463909,
672.210515304308,
14566.0279110099,
170927.78078521,
1525821.93755474,
11366312.5854366,
73868184.8847336,
430610301.867199,
2297726058.33628,
11406906364.7078,
53429633348.726,
239204245061.543,
1037569229691.63,
4440670923491.02,
19385786140661,
108390210481645,
2.4429367439691e-015,
1.4097646088969e-012,
-8.1610167897946e-011,
-2.2790490204496e-008,
-1.5951868755541e-006,
-5.4515698129209e-005,
-0.000902748062402921,
0.000476641583177264,
0.403718060688101,
11.3250139182509,
195.75060974824,
2526.10978407101,
26230.959211931,
228725.678544486,
1724448.556825,
11494565.1027429,
68980892.2243655,
378540210.929366,
1926018147.60076,
9203415540.67664,
41820439362.7481,
183069601908.131,
783979860688.641,
3359830680474.02,
15066614162007.1,
80531744425420.1)
}

grid
}






#############################################################################
tweedie.dev <- function(y, mu, power)
{
# 
# Peter K Dunn 
# 29 Oct 2000 
# 

	p <- power
	if(p == 1) {
      dev <- array( dim=length(y) )
      mu <- array( dim=length(y), mu )
		dev[y!=0] <- y[y!=0] * log(y[y!=0]/mu[y!=0]) - (y[y!=0] - mu[y!=0])
      dev[y==0] <- mu[y==0]
	}
   else{
   	if(p == 2) {
	   	dev <- log(mu/y) + (y/mu) - 1
   	}
      else{
      	if (p== 0) {
   		   dev <- (y-mu)^2
      		dev <- dev/2
         }
         else{
      		dev <- (y^(2 - p))/((1 - p) * (2 - p)) - 
		      	(y * (mu^(1 - p)))/(1 - p) + 
      			(mu^(2 - p))/(2 - p)
         }
	  	}
   }
	dev * 2
}


#############################################################################
dtweedie.dldphi <- function(phi, mu, power, y ){
# Calculates the log-likelihood
# function, wrt phi, for p>2.  In particular, it returns
#    sum{ d(log f)/d(phi) } = d( log-likelihood )/d(phi).
# The mle of phi can be found, therefore, by setting this to zero.
# y  is generally a vector of observed values.

#
# Peter Dunn
# 31 Jan 2001
#

if ( (power != 2 ) & ( power != 1 ) ) {

   k <- phi^(1/(power-2))
# cat("k=",k,"\n")
# cat("phi=",phi,"\n")
# cat("power=",power,"\n")
   if ( k < 1  & k > 0 ) {
      # Use the transform f(y; mu, phi) = c f(c*y; c*mu, c^(2-p)*phi)
      # and differentiate with c=phi^(1/(p-2)):
      #    d log f / d phi = c^(2-p) * {df(cy; c*mu, 1)/dphi} / f(cy; c*mu, 1)
      f <- dtweedie( y=k*y, power=power, mu=k*mu, phi=1 )
      d <- dtweedie.dlogfdphi( y=k*y, power=power, mu=k*mu, phi=1 )
         # Note:  We need dlogf/dphi = dlogf.dphi * f
      top <- d * f
      d <- -2* sum( top / f * k^(2-power) )

   }
   else{
    # Compute directly
      d <- -2*sum( dtweedie.dlogfdphi(y=y, power=power, mu=mu, phi=phi) )
   }
}
else{
   # Cases p==1 and  p=2
   d <- -2*sum( dtweedie.dlogfdphi(y=y, power=power, mu=mu, phi=phi) )
}
d
}



#############################################################################
dtweedie.dlogfdphi <- function(y, mu, phi, power)
{
#
# Calculates d(log f)/d(phi) for the Tweedie
# densities.
# It is used, for example, in mle fitting of phi.  We would then
# sum over  y  and set this function to 0.
#
#
# Peter Dunn
# 31 Jan 2001
#

	p <- power
	a <- (2 - p)/(1 - p)
	if(length(phi) == 1) {
		phi <- array(dim = length(y), phi)
	}
	if(length(mu) == 1) {
		mu <- array(dim = length(y), mu)
	}
	A <- (y * mu^(1 - p))/(phi^2 * (p - 1))
	B <- mu^(2 - p)/(phi^2 * (2 - p))

	if(power > 2) {
		f <- array(dim = c(length(y)))
      # Here, we evaluate logv and everything as normal.
      # If logv has infinite values then we resort to other tactics.

		kv <- dtweedie.kv.bigp(power = power, phi = phi, y = y)$kv
		dv.dphi <- (kv * (a - 1))/phi
		out.logv <- dtweedie.logv.bigp(power = power, phi = phi, y = y)

      # Now see if this causes problems.
      logv <- out.logv$logv

      # Now detect problem computing  logv  and remedy them
      probs <- (is.infinite(logv)) | (is.nan(logv)) | (y<1)

 		if(any(probs)) {

         # OK then:  Troubles computing log(V).
         # best we can do is use definition I think.
         delta <- 1.0e-5
         a1 <- dtweedie(power=power, phi=phi[probs], mu=mu[probs], y=y[probs])
         a2 <- dtweedie(power=power, phi=phi[probs]+delta, mu=mu[probs], y=y[probs])
         f[probs] <- (log(a2) - log(a1) ) / delta

		}

		f[!probs] <- A[!probs] + B[!probs] + dv.dphi[!probs]/exp(logv[!probs])
	}
#  END p>2

	if(power == 2) {
      f <-   -log(y)  + ( y/mu ) + digamma(1/phi) - 1 + log( mu*phi )
      f <- f / (phi^2)
	}

	if(power == 1) {

      f <- mu - y - y*log(mu/phi) + y*digamma(1+(y/phi))
      f <- f / (phi^2)
	}

	if((power > 1) && (power < 2)) {
      # We need to treat y==0 separately.
      # We fill  f  with the Y=0 case
   	f <- array(  dim = length(y), mu^(2-power) / ( phi^2 * (2-power) )  )
		jw <- dtweedie.jw.smallp(power=power, phi=phi[y > 0], y=y[y > 0])$jw
		dw.dphi <- (jw * (a - 1)) / phi[y > 0]
		logw <- dtweedie.logw.smallp(power=power, phi=phi[y > 0], y=y[y > 0])$logw
		f[y>0] <- A[y > 0] + B[y > 0] + dw.dphi/exp(logw)
	}
	f
}


#############################################################################
dtweedie.interp <- function(grid, nx, np, xix.lo, xix.hi,
                            p.lo, p.hi, power, xix) {
# Does the interpolation calculation

# Peter Dunn
# 17 April 2001

# Set-up
# xi-dimension
jt <- seq(0, nx, 1)
ts <- cos((2 * jt + 1)/(2 * nx + 2) * pi)
ts <- ((xix.hi + xix.lo) + ts * (xix.hi - xix.lo))/2	#

# p-dimension
jp <- seq(0, np, 1)
ps <- cos((2 * jp + 1)/(2 * np + 2) * pi)
ps <- ((p.hi + p.lo) + ps * (p.hi - p.lo))/2	#

rho <- array(dim = nx + 1)
dd1 <- array(dim = nx + 1)	#

# interpolate to find divided differences
for(i in 1:(nx + 1)) {
   dd1[i] <- grid[i, np + 1]
   for(k in seq(np, 1, -1)) {
      dd1[i] <- dd1[i] * (power - ps[k]) + grid[i, k]
	}
}

# Now use divided difference to interpolate rho or 1/rho
rho <- dd1[nx + 1]
for(k in seq(nx, 1, -1)) {
   rho <- rho * (xix - ts[k]) + dd1[k]
}

if ( power >= 3) {
   rho <- 1 / rho
}
rho

}



#############################################################################
dtweedie.inversion <- function(y, power, mu, phi, exact=TRUE, method=3){ 
# 
# Peter K Dunn 
# 06 Aug 2002
#

# Error checks
if ( power<1) stop("power must be greater than 1.")
if ( any(phi <= 0)) stop("phi must be positive.")
#if ( power>2) if ( any(y <= 0) ) stop("y must be a positive vector.")
#if ( (power>1 & (power<2) ) if ( any(y < 0) ) stop("y must be a non-negative vector.")
if ( any(mu <= 0) ) stop("mu must be positive.")
if ( length(mu)>1) {
   if ( length(mu)!=length(y) ) {
      stop("mu must be scalar, or the same length as y.")
   }
}
else {
   mu <- array( dim=length(y), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim=length(y), phi )
   # A vector of all phi's
}
save.method <- method
if ( !is.null(method)){
   if ( length(method)>1 ) {
      method <- save.method <- method[1]
   }
   if ( !(method %in% c(1,2,3)) ) stop("method must be 1, 2 or 3 (or left empty).")
}

y.len <- length(y)
density <- y
its <- y

if ( is.null(method)){
   method <- array( dim=length(y))
} else {
   method <- array( method, dim=length(y))
}

for (i in (1:y.len)) {

   # There are three approaches, each a product of a simple bit
   # and a complicated bit computed in FORTRAN
   # We choose Method 3 if no other is requested.
   #
   # The methods are documented in Dunn and Smyth (2008):
   # Method 1: Evaluate a: compute a(y, phi) = f(y; 1, phi)
   # Method 2: Rescale the mean to 1
   # Method 3: Rescale y to 1 and evaluate b.
   
   if ( y[i] <= 0 ) {
   	if ( (power>1) & (power<2) ) {
			if ( y[i]==0 ) {
		      density[i] <- exp( -mu[i] ^ (2-power) / ( phi[i] * (2-power) ) )
   		} else {
				density[i] <- 0
   		}
	   } else {
   		density[i] <- 0
   	}
   } else {
      # Here, y > 0
      m2 <- 1/mu[i]
      
      theta <- ( mu[i]^(1-power) - 1 ) / ( 1 - power )
      if ( ( abs(power - 2 ) ) < 1.0e-07 ){
         kappa <- log(mu[i]) + (2 - power) * ( log(mu[i])^2 ) / 2
      } else {
         kappa <- ( mu[i]^(2-power) - 1 ) / ( 2 - power )
      }
      m1 <- exp( (y[i]*theta - kappa )/phi[i] )
   
      dev <- tweedie.dev(y=y[i], mu=mu[i], power=power )
      m3 <- exp( -dev/(2*phi[i]) ) / y[i]
   
      min.method <- min( m1, m2, m3 )
   
      # Now if no method requested, find the notional "optimal"
      if ( is.null(method[i]) ) {
         if ( min.method==m1 ){
            use.method <- 1
         } else {
            if ( min.method==m2 ) {
               use.method <- 2
            } else {
               use.method <- 3
            }
         }
      } else {
         use.method <- method[i]
      }
   
      # Now use the method
      # NOTE: FOR ALL  METHODS, WE HAVE mu=1
      verbose <- FALSE
      
      if ( use.method==2 ) {
         tmp <- .Fortran( "pdf",
            as.double(power),
            as.double(phi[i] / (mu[i]^(2-power)) ), # phi
            as.double(y[i]/mu[i]), # y
            as.double(1), # mu
            as.integer( exact ), #exact as an integer
            as.integer( verbose ), #verbose as an integer
            as.double(0), # funvalue
            as.integer(0), # exitstatus
            as.double(0), # relerr
            as.integer(0), # its
            PACKAGE="tweedie")
   
         den <- tmp[[7]]
         density[i] <- den * m2
   
      } else {
         if ( use.method==1 ) {
            tmp <- .Fortran( "pdf",
               as.double(power),
               as.double(phi[i]), # phi
               as.double(y[i]), # y
               as.double(1), # mu
               as.integer( exact ), #exact as an integer
               as.integer( verbose ), #verbose as an integer
               as.double(0), # funvalue
               as.integer(0), # exitstatus
               as.double(0), # relerr
               as.integer(0), # its
               PACKAGE="tweedie")
   
            den <- tmp[[7]]
            density[i] <- den * m1
   
         } else { # use.method==3
            tmp <- .Fortran( "pdf",
               as.double(power),
               as.double(phi[i]/(y[i]^(2-power))), # phi
               as.double(1), # y
               as.double(1), # mu
               as.integer( exact ), #exact as an integer
               as.integer( verbose ), #verbose as an integer
               as.double(0), # funvalue
               as.integer(0), # exitstatus
               as.double(0), # relerr
               as.integer(0), # its
               PACKAGE="tweedie")
   
            den <- tmp[[7]]
            density[i] <- den * m3
         }
      }
   }
}

density

}
      

#############################################################################
dtweedie.jw.smallp <- function(y, phi, power ){ 
#
# Peter K Dunn
# 18 Jun 2002
#

#
# Error traps
#

if ( power<1) stop("power must be between 1 and 2.")
if ( power>2) stop("power must be between 1 and 2.")
if ( any(phi<=0) ) stop("phi must be strictly positive.")
if ( any(y<=0) ) stop("y must be a strictly positive vector.")

#
# Set up
#
p <- power
a <- ( 2-p ) / ( 1-p )	 	# Note that a<0 for 1<p<2

a1 <- 1 - a
r <- -a*log(y) + a*log(p-1) - a1*log(phi) -
	log(2-p)		# All terms to power j

drop <- 37 			# Accuracy of terms: exp(-37)
#
# Find limits of summation using Stirling's approximation
# to approximate gamma terms, and find max as j.max
#
logz <- max(r) 			# To find largest  j  needed
j.max <- max( y^( 2-p ) / ( phi * (2-p) ) )
j <- max( 1, j.max )

c <- logz + a1 + a*log(-a)
wmax <- a1*j.max
estlogw <- wmax

# First, the upper limit of j

while(estlogw > (wmax-drop) ){
   j <- j + 2
   estlogw <- j*(c - a1*log(j))
}

hi.j <- ceiling(j)

# Now the lower limit of j
logz <- min(r) 
j.max <- min( y^( 2-power ) / ( phi * (2-power) ) )

j <- max( 1, j.max)
wmax <- a1*j.max 
estlogw <- wmax 
 
while ( ( estlogw > (wmax-drop) ) && ( j>=2) ) {
   j <- max(1, j-2)
   estlogw <- j*(c-a1*log(j))
}

lo.j <- max(1, floor(j))
# 
# Now sum the series between established limits.
# We ensure it works for vector y.

j <- seq( lo.j, hi.j) 				# sequence of j terms needed
o <- matrix( 1, nrow=length(y)) 	# matrix of ones
 
g <- matrix(lgamma( j+1 ) + lgamma( -a*j ),
	nrow=1, ncol=hi.j - lo.j + 1) 
logj <- matrix(log(j),
	nrow=1, ncol=hi.j - lo.j + 1) 
og <- o %*% g 							# matrix of gamma terms
ologj <- o %*% logj 
A <- outer(r,j) - og + ologj		# the series, almost ready to sum
m <- apply(A,1,max)					# avoid overflow; find maximum values
we <- exp( A - m )					# evaluate terms, less max.
sum.we <- apply( we,1,sum)			# sum terms
jw <- sum.we * exp( m )				# now restore max.
# Since derivs may be negative, can't use log-scale

list(lo=lo.j, hi=hi.j, jw=jw, j.max=j.max )

}
 


#############################################################################
dtweedie.kv.bigp <- function(y, phi, power){ 
# 
# Peter K Dunn 
# 18 Jun 2002
# 

#
# Error traps
#

if ( power<2) stop("power must be greater than 2.")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( any(y<=0) ) stop("y must be a strictly positive vector.")
if ( length(phi)>1) {
   if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.")
}
else {
   phi <- array( dim = length(y), phi )
   # A vector of all phi's
}



#
# Set up
#

p <- power
a <- ( 2 - p ) / ( 1 - p )
a1 <- 1 - a

r <- -a1*log(phi) - log(p-2) - a*log(y) +
    a*log(p-1)
drop <- 37

#
# Now we find limits of summation, using Stirling's approximation
# to approximate the gamma terms, and find max as k.max.
#

logz <- max(r)
k.max <- max( y^(2-p) / ( phi * (p-2) ) )
k <- max( 1, k.max )

c <- logz + a1 + a*log(a)
vmax <- k.max * a1
estlogv <- vmax
#
# Now we search either side for when we can
# ignore terms as they are negligible.
# Remember we are using the envelope as the tests,
# not individual terms in the series.

# First:  the upper limit of k

while ( estlogv > (vmax - drop) ) {
	k <- k+2
	estlogv <- k*( c -a1*log(k) )
}

hi.k <- ceiling(k)
#
# Now the lower limit of k
#
logz <- min(r)
k.max <- min( y^(2-p) / ( phi * (p-2) ) ) 
k <- max( 1, k.max )
c <- logz + a1 + a*log(a)
vmax <- k.max * a1
estlogv <- vmax

while ( (estlogv > (vmax-drop) ) && ( k>=2) ) {
	k <- max(1, k-2)
	estlogv <- k*( c - a1*log(k) )
}

lo.k <- max(1, floor(k) )

#
# Now sum the series  between established limits.
# We ensure it works for vector y.
#
k <- seq(lo.k, hi.k) 
o <- matrix( 1, nrow=length(y)) 

g <- matrix( lgamma( 1+a*k) - lgamma(1+k),
	nrow=1, ncol=length(k) )
logk <- matrix( log(k),
	nrow=1, ncol=length(k) )

og <- o %*% g
ologk <- o %*% logk

A <- outer(r,k) + og + ologk

C <- matrix( sin( -a*pi*k ) * (-1)^k,
	nrow=1, ncol=length(k) )
C <- o %*% C

m <- apply(A, 1, max)
ve <- exp(A - m)
sum.ve <- apply( ve*C, 1, sum )
kv <- sum.ve * exp( m )
# Since derivs may be negative, can't use log-scale

list(lo=lo.k, hi=hi.k, kv=kv, k.max=k.max )

}


#############################################################################
qtweedie <- function(p, xi=power, mu, phi, power=NULL){

   if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
   xi.notation <- TRUE
   if ( is.null(power) ) power <- xi
   if ( is.null(xi) ) {
      xi <- power
      xi.notation <- FALSE
   }
   if ( xi != power ) {
      cat("Different values for xi and power given; the value of xi used.\n")
      power <- xi
   }
   index.par <- ifelse( xi.notation, "xi","p")

# Error checks
if ( any(power<1) ) stop("power must be greater than 1.\n")
if ( any(phi <= 0) ) stop("phi must be positive.")
if ( any(p<0) ) stop("p must be between zero and one.\n")
if ( any(p>1) ) stop("p must be between zero and one.\n")
if ( any(mu <= 0) ) stop("mu must be positive.\n")
if ( length(mu)>1) {
   if ( length(mu)!=length(p) ) stop("mu must be scalar, or the same length as p.\n")
}
else {
   mu <- array( dim=length(p), mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=length(p) ) stop("phi must be scalar, or the same length as p.\n")
}
else {
   phi <- array( dim=length(p), phi )
   # A vector of all phi's
}
# Now, if mu is a vector, make power correct length
if ( length(power) == 1 ) {
   if ( length(mu) > 1 ) {
      power <- rep( power, length(mu) ) 
   }
}

#
# End error checks
#

len <- length(p) 

# Some monkeying around to explicitly account for the cases p=1 and p=0
ans <- ans2 <- rep( NA, length=len )
if ( any(p==1) ) ans2[p==1] <- Inf
if ( any(p==0) ) ans2[p==0] <- 0

ans     <-  ans[ ( (p>0) & (p<1) ) ]
mu.vec  <-  mu[ ( (p>0) & (p<1) ) ]
phi.vec <-  phi[ ( (p>0) & (p<1) ) ]
p.vec   <- p[ ( (p>0) & (p<1) ) ]

for (i in (1:length(ans)) ) {

   mu.1 <- mu.vec[i]
   phi.1 <- phi.vec[i]
   p.1 <- p.vec[i]
   pwr <- power[i]
   
   prob <- p.1 # Rename p to avoid confusion with  pwr
   
   if ( pwr<2 ) {
      qp <- qpois(prob, lambda=mu.1/phi.1)
      if ( pwr == 1 ) ans[i] <- qp   
   }
   
   qg <- qgamma(prob,  rate=1/(phi.1*mu.1), shape=1/phi.1 )
   if ( pwr==2 ) ans[i] <- qg
   
   # Starting values
   # - for 1<pwr<2, linearly interpolate between Poisson and gamma
   if ( (pwr>1) & ( pwr<2) ) {
      start <- (qg - qp)*pwr + (2*qp - qg)
   }
   
   # - for pwr>2, start with gamma
   if ( pwr>2 ) start <- qg
   
   # Solve!
   if ( ( pwr>1) & (pwr<2) ) { # This gives a *lower* bound on the value of the answer (if y>0)
      step <- dtweedie(y=0, mu=mu.1, phi=phi.1, power=pwr)
      # This is P(Y=0), the discrete "step"
      if ( prob <= step ) {
         ans[i] <- 0
      }
   }
   
   if ( is.na(ans[i]) ) { # All cases except Y=0 when 1 < pwr < 2
         
      pt2 <- function( q, mu, phi, pwr, p.given=prob ){ 

            ptweedie(q=q, mu=mu, phi=phi, power=pwr ) - p.given
      }
      pt <- pt2( q=start, mu=mu.1, phi=phi.1, pwr=pwr, p.given=prob)
      
      if ( pt == 0 ) ans2[i] <- start
      
      if ( pt > 0 ) { 
		    loop <- TRUE
				start.2 <- start
	        while ( loop ) {
               # Try harder
               start.2 <- 0.5*start.2
               if (pt2( q=start.2, mu.1, phi.1, pwr, p.given=prob )<0 ) loop=FALSE
					# RECALL:  We are only is this part of the loop if  pt>0
            }
      }
      
		if ( pt < 0) {
         loop <- TRUE
         start.2 <- start
         
         while ( loop ) {
            # Try harder
            start.2 <- 1.5*(start.2 + 2)
            if (pt2( q=start.2, mu.1, phi.1, pwr, p.given=prob )>0 ) loop=FALSE
				# RECALL:  We are only is this part of the loop if  pt<0
         }
      }

       ans[i] <- uniroot(pt2, c(start, start.2), mu=mu.1, phi=phi.1, p=pwr, 
                        p.given=prob )$root
   }
    
}

ans2[ is.na(ans2) ] <-  ans
ans2
}



#############################################################################
rtweedie <- function(n, xi=power, mu, phi, power=NULL){

   if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
   xi.notation <- TRUE
   if ( is.null(power) ) power <- xi
   if ( is.null(xi) ) {
      xi <- power
      xi.notation <- FALSE
   }
   if ( xi != power ) {
      cat("Different values for xi and power given; the value of xi used.\n")
      power <- xi
   }
   index.par <- ifelse( xi.notation, "xi","p")

# Error checks
if ( power<1) stop("power must be greater than or equal to 1.\n")
if ( any(phi<=0) ) stop("phi must be positive.")
if ( n<1 ) stop("n must be a positive integer.\n")
if ( any(mu<=0) ) stop("mu must be positive.\n")
if ( length(mu)>1) {
   if ( length(mu)!=n ) stop("mu must be scalar, or of length n.\n")
}
else {
   mu <- array( dim=n, mu )
   # A vector of all mu's
}
if ( length(phi)>1) {
   if ( length(phi)!=n ) stop("phi must be scalar, or of length n.\n")
}
else {
   phi <- array( dim=n, phi )
   # A vector of all phi's
}
    if (power==1) {
        rt <- phi * rpois(n, lambda=mu / ( phi  ) )
    }
    
    if (power==2) {
        alpha <- (2-power)/(1-power)
        gam <- phi * (power-1) * mu ^(power-1)
        rt <- rgamma( n, shape=1/phi, scale=gam )
    }
    
    if ( power>2) {
        rt <- qtweedie( runif(n), mu=mu ,phi=phi , power=power)
    }
    
    if ( (power>1) & (power<2) ) {
        # Two options:  As above or directly.
        # Directly is faster
        rt <- array( dim=n, NA)
    
        lambda <- mu^(2-power) / ( phi * (2-power) )
        alpha <- (2-power)/(1-power)
        gam <- phi * (power-1) * mu ^(power-1)
    
        N <- rpois(n, lambda=lambda)
        for (i in (1:n) ){
           rt[i] <- sum( rgamma(N[i], shape = -alpha, scale=gam[i]) )
        }
    }
    as.vector(rt)
}


#############################################################################
tweedie.profile <- function(formula, p.vec=NULL, xi.vec=NULL, link.power = 0, 
		data, weights, offset, fit.glm=FALSE, 
      do.smooth=TRUE, do.plot=FALSE, 
      do.ci=do.smooth, eps=1/6,
      control=list( epsilon=1e-09, maxit=glm.control()$maxit, trace=glm.control()$trace ),
      do.points=do.plot, method="inversion", conf.level=0.95, 
      phi.method=ifelse(method=="saddlepoint","saddlepoint","mle"), verbose=FALSE, add0=FALSE) {
# verbose gives feedback on screen:
#    0 : minimal (FALSE)
#    1 : small amount (TRUE)
#    2 : lots
# Determine the value of  p  to use to fit a Tweedie distribution
# A profile likelihood is used (can can be plotted with do.plot=TRUE)

# The plot can be (somewhat) recovered using (where  out  is the
# returned object, with smoothing):
#     plot( out$x, out$y)
# The actual points can then be added:
#     points( out$p, out$L)
# (Note that out$L and out$y are the same if smoothing is  not requested.  
# Likewise for out$p and out$x.)

# Peter Dunn
# 07 Dec 2004


if ( is.logical( verbose ) ) {
   verbose <- as.numeric(verbose)
}

if (verbose >= 1 ) {
	cat("---\n This function may take some time to complete;\n")
	cat(" Please be patient.  If it fails, try using  method=\"series\"\n")
	cat(" rather than the default  method=\"inversion\"\n")
	cat(" Another possible reason for failure is the range of p:\n")
	cat(" Try a different input for  p.vec\n---\n")
}


cl <- match.call()
mf <- match.call()
m <- match(c("formula", "data", "weights","offset"), names(mf), 0L)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
Y <- model.response(mf, "numeric")
X <- if (!is.empty.model(mt))
	model.matrix(mt, mf, contrasts)
else matrix(, NROW(Y), 0)
weights <- as.vector(model.weights(mf))
if (!is.null(weights) && !is.numeric(weights))
	stop("'weights' must be a numeric vector")
offset <- as.vector(model.offset(mf))
if (!is.null(weights) && any(weights < 0))
	stop("negative weights not allowed")
if (!is.null(offset)) {
	if (length(offset) == 1)
		offset <- rep(offset, NROW(Y))
	else if (length(offset) != NROW(Y))
		stop(gettextf("number of offsets is %d should equal %d (number of observations)",
				length(offset), NROW(Y)), domain = NA)
}



### NOW some notation stuff
xi.notation <- TRUE
if ( is.null(xi.vec) & !is.null(p.vec) ){ # If p.vec given, but not xi.vec
	xi.vec <- p.vec
	xi.notation <- FALSE
}

if ( is.null(p.vec) & !is.null(xi.vec) ){ # If xi.vec given, but not p.vec
	p.vec <- xi.vec
	xi.notation <- TRUE
}

if ( is.null( p.vec ) & is.null(xi.vec)) { # Neither xi.vec or p.vec are given
	if ( any(Y==0) ){
		p.vec <- seq(1.2, 1.8, by=0.1)
	} else {
		p.vec <- seq(1.5, 5, by=0.5)
	}
	xi.vec <- p.vec
	xi.notation <- TRUE
}

### AT THIS POINT, we have both xi.vec and p.vec declared, and both are the same
### but we stick with using xi.vec hereafter


# Determine notation to use in output (consistent with what was supplied by the user)
index.par <- ifelse( xi.notation, "xi","p")

# Find if values of p/xi have been specified that are not appropriate
xi.fix <- any( xi.vec<=1 & xi.vec!=0, na.rm=TRUE)

if ( xi.fix ) {
	xi.fix.these <- (xi.vec<=1 & xi.vec!=0)
	xi.vec[ xi.fix.these ] <- NA 
	if ( length(xi.vec) == 0 ) {
		  stop(paste("No values of",index.par,"between 0 and 1, or below 0, are possible.\n"))
	} else {
		cat("Values of",index.par,"between 0 and 1 and less than zero have been removed: such values are not possible.\n")
	}
}

# Warnings
if ( any( Y == 0 ) & any( (xi.vec >= 2) | (xi.vec <=0) ) ) {
	xi.fix.these <- (xi.vec>=2 | xi.vec<=0)
	xi.vec[ xi.fix.these ] <- NA 
	if ( length(xi.vec) == 0 ) {
		  stop(paste("When the response variale contains exact zeros, all values of",index.par,"must be between 1 and 2.\n"))
	} else {
		cat("When the response variale contains exact zeros, all values of",index.par,"must be between 1 and 2; other values have been removed.\n")
	}
}


# Remove any NAs in the p/xi values
xi.vec <- xi.vec[ !is.na(xi.vec) ] 

if ( do.smooth & (length(xi.vec) < 5) ) {
   warning(paste("Smoothing needs at least five values of",index.par,".") )
   do.smooth <- FALSE
   do.ci <- FALSE
}
if ( (conf.level >= 1) | (conf.level <=0)  ){
   stop("Confidence level must be between 0 and 1.")
}

if ( !do.smooth & do.ci ) {
   do.ci <- FALSE
	warning("Confidence intervals only computed if  do.smooth=TRUE\n")
}

if ( add0 ) {
	p.vec  <- c( 0, p.vec )
	xi.vec <- c( 0, xi.vec )
}


cat(xi.vec,"\n")




# Some renaming
ydata <- Y
model.x <- X



# Now, fit models!  Need the Tweedie class

# First define the function to *minimize*;
# since we want a *maximum* likeihood estimator,
# define the negative.
dtweedie.nlogl <- function(phi, y, mu, power) {
    ans <- -2 * sum( log( dtweedie( y=y, mu=mu, phi=phi, power=power ) ) )
    if ( is.infinite( ans ) ) {
      # If infinite, replace with saddlepoint estimate?
        ans <- sum( tweedie.dev(y=y, mu=mu, power=power) )/length( y )
     }
    attr(ans, "gradient") <- dtweedie.dldphi(y=y, mu=mu, phi=phi, power=power)
    
    ans
}


# Set up some parameters
xi.len <- length(xi.vec)
phi <- NaN

L <- array( dim = xi.len )
phi.vec <- L

# Now most of these are for debugging:
b.vec <- L
c.vec <- L
mu.vec <- L
b.mat <- array( dim=c(xi.len, length(ydata) ) )

for (i in (1:xi.len)) {

   if ( verbose>0) {
      cat( paste(index.par," = ",xi.vec[i],"\n", sep="") )
   } else {
		cat(".")
	}


   # We find the mle of phi.
   # We try to bound it as well as we can, 
   # then use  uniroot to solve for it.

   # Set things up
   p <- xi.vec[i]

   phi.pos <- 1.0e2
   bnd.neg <- -Inf
   bnd.pos <- Inf

    if (verbose==2) cat("* Fitting initial model:")

   # Fit the model with given p
   catch.possible.error <- try(
      fit.model <- glm.fit( x=model.x, y=ydata, weights=weights, offset=offset,
      							 control=control,
                            family=tweedie(var.power=p, link.power=link.power)),
      silent = TRUE
   )


	skip.obs <- FALSE
   if ( class( catch.possible.error )=="try-error" ) {
      skip.obs <- TRUE 
   }
   
   if( skip.obs ) {
      warning(paste("  Problem near ",index.par," = ",p,"; this error reported:\n     ",
                     catch.possible.error,
                    " Examine the data and function inputs carefully.") )
   
      # NOTE:  We need epsilon to be very small to make a
      #        smooth likelihood plot
      mu <- rep(NA, length(ydata) )
      
    } else {
      mu <- fitted( fit.model )
    }
    if (verbose==2) cat(" Done\n")

    ### -- Start:  Estimate phi
    if (verbose>=1) cat("* Phi estimation, method: ", phi.method)
    
    if( skip.obs ) {
      if (verbose>=1) cat("; but skipped for this obs\n")
      phi.est <- phi <- phi.vec[i] <- NA
    } else {
      if ( phi.method=="mle"){
   
         if (verbose>=1) cat(" (using optimize): ")
      
            # Saddlepoint approx of phi:
            phi.saddle <- sum( tweedie.dev(y=ydata, mu=mu, power=p) )/length( ydata )
      
            if ( is.nan(phi) ) {
      
               # NOTE:  This is only used for first p.
               #        The remainder use the previous estimate of phi as a starting point
      
               # Starting point for phi:  use saddlepoint approx/EQL.
            phi.est <- phi.saddle
      
            } else {
               phi.est <- phi
            }
            low.limit <- min( 0.001, phi.saddle/2)
            
      #    ans <- nlm(p=phi.est, f=dtweedie.nlogl, 
      #                hessian=FALSE,
      #                power=p, mu=mu, y=data)
			if ( p!= 0 ) {
				ans <- optimize(f=dtweedie.nlogl, maximum=FALSE, interval=c(low.limit, 10*phi.est),
							power=p, mu=mu, y=ydata )
				phi <- phi.vec[i] <- ans$minimum
			} else {
				phi <- phi.vec[i] <- sum( (ydata-mu)^2 ) / length(ydata)
			}
		if (verbose>=1) cat(" Done (phi =",phi.vec[i],")\n")
   
      } else{ # phi.method=="saddlepoint")
   
         if (verbose>=1) cat(" (using mean deviance/saddlepoint): ")
         phi <- phi.est <- phi.vec[i] <- sum( tweedie.dev(y=ydata, mu=mu, power=p) )/length( ydata )
         if (verbose>=1) cat(" Done (phi =",phi,")\n")
      }
    }
    
    ### -- Done:  estimate phi

    # Now determine log-likelihood at this p
    # Best to use the same density evaluation for all:  series or inversion
    
    if (verbose>=1) cat("* Computing the log-likelihood ")

   # Now compute the log-likelihood
    if (verbose>=1) cat("(method =",method,"):")

    if ( skip.obs ) {
      if (verbose>=1) cat(" but skipped for this obs\n")
      L[i] <- NA
    } else {   
      if ( method=="saddlepoint") {
         L[i] <- dtweedie.logl.saddle(y=ydata, mu=mu, power=p, phi=phi, eps=eps)
      } else {
      if (p==2) {
            L[i] <- sum( log( dgamma( rate=1/(phi*mu), shape=1/phi, x=ydata ) ) )
         } else {
            if ( p == 1 ) {
            	if ( phi==1 ){
						# The phi==1 part added 30 Sep 2010.
	               L[i] <- sum( log( dpois(x=ydata/phi, lambda=mu/phi ) ) )
	            } else { # p=1, but phi is not equal to zero
						# As best as we can, we determine if the values
						# of  ydata  are multiples of  phi
						y.on.phi <- ydata/phi
						close.enough <- array( dim=length(y.on.phi))
						for (i in (1:length(y.on.phi))){
							if (isTRUE(all.equal(y.on.phi, as.integer(y.on.phi)))){
		          	     L[i] <- sum( log( dpois(x=round(y/phi), lambda=mu/phi ) ) )
							} else {
								L[i] <- 0
							}
						}
	            }
            } else {
				if ( p == 0 ) {
					L[i] <- sum( dnorm(x=ydata, mean=mu, sd=sqrt(phi), log=TRUE) )
				} else {
					L[i] <- switch(
					pmatch(method, c("interpolation","series", "inversion"), 
                                 nomatch=2),
					"1" = dtweedie.logl( mu=mu, power=p, phi=phi, y=ydata),
					"2" = sum( log( dtweedie.series( y=ydata, mu=mu, power=p, phi=phi) ) ),
					"3" = sum( log( dtweedie.inversion( y=ydata, mu=mu, power=p, phi=phi) ) ) )
				}
            }
         }
      }
   }

    if (verbose>=1) {
       cat(" Done: L =",L[i],"\n")
    } 
}

if ( verbose == 0 ) cat("Done.\n")
### Smoothing
   # If there are infs etc in the log-likelihood,
	# the smooth can't happen.  Perhaps this can be worked around.
	# But at present, we note when it happens to ensure no future errors
	# (eg in computing confidence interval for p).

# y and x are the smoothed plot produced; here we set them
# as NA so they are defined when the function is returned to
# the workspace
y <- NA
x <- NA

if ( do.smooth ) {
    L.fix <- L
    xi.vec.fix <- xi.vec
	 phi.vec.fix <- phi.vec
	 if ( any( is.nan(L) ) | any( is.infinite(L) ) | any( is.na(L) ) ) {
       retain.these <- !( ( is.nan(L) ) | ( is.infinite(L) ) | ( is.na(L) ) )
	    xi.vec.fix <- xi.vec.fix[ retain.these ]
		 phi.vec.fix <- phi.vec.fix[ retain.these ]
	    L.fix <- L.fix[ retain.these ]
		 
# 	    p.vec.fix <- p.vec.fix[ !is.infinite(L.fix) ]
# 	    phi.vec.fix <- phi.vec.fix[ !is.infinite(L.fix) ]
# 	    L.fix <- L.fix[ !is.infinite(L.fix) ]
	 
        if (verbose>=1) cat("Smooth perhaps inaccurate--log-likelihood contains  Inf  or  NA.\n")
    }
	 #else {
    
    if ( length( L.fix ) > 0 ) {
	   if (verbose>=1) cat(".")
      if (verbose>=1) cat(" --- \n")
      if (verbose>=1) cat("* Smoothing: ")
    # Smooth the points
       # - get smoothing spline
    ss <- splinefun( xi.vec.fix, L.fix )

      # Plot smoothed data
      xi.smooth <- seq(min(xi.vec.fix), max(xi.vec.fix), length=50 )
      L.smooth <- ss(xi.smooth )
      
      if ( do.plot) {
         keep.these <- is.finite(L.smooth) & !is.na(L.smooth)
         L.smooth <- L.smooth[ keep.these ] 
         xi.smooth <- xi.smooth[ keep.these ] 
         if ( verbose>=1 & any( !keep.these ) ) {
            cat(" (Some values of L are infinite or NA for the smooth; these are ignored)\n")
         }
        
         yrange <- range( L.smooth, na.rm=TRUE )
         
         plot( yrange ~ range(xi.vec),
            type="n",
				las=1,
            xlab=ifelse(xi.notation, expression(paste( xi," index")), expression(paste( italic(p)," index")) ),
            ylab=expression(italic(L)))
         lines( xi.smooth, L.smooth,
            lwd=2)
         rug( xi.vec )
			if (do.points) {
			   points( L ~ xi.vec, pch=19)
			}
		 if (add0) lines(xi.smooth[xi.smooth<1], L.smooth[xi.smooth<1], col="gray", lwd=2)
      }
      x <- xi.smooth
      y <- L.smooth
      
   } else {
      cat("  No valid values of the likelihood computed: smooth aborted\n",
          "  Consider trying another value for the input  method.\n")
   }
}
else {
   if ( do.plot) {
         
      keep.these <- is.finite(L) & !is.na(L)
      xi.vec <- xi.vec[ keep.these ] 
      L <- L[ keep.these ] 
      phi.vec <- phi.vec[ keep.these ]
      if ( verbose>=1 & any( keep.these ) ) {
         cat(" Some values of L are infinite or NA, and hence ignored\n")
      }

      yrange <- range( L, na.rm=TRUE )
      # Plot the data we have
      plot( yrange ~ range(xi.vec),
         type="n",     
         las=1,
         xlab=ifelse( xi.notation, expression(paste(xi," index")), expression(paste(italic(p)," index")) ),
         ylab=expression(italic(L)))
      lines( L ~ xi.vec, lwd=2)
	
      rug( xi.vec )
		if (do.points) {
		   points( L ~ xi.vec, pch=19)
		}

   }
   x <- xi.vec
   y <- L
}
if (verbose>=2) cat(" Done\n")


### Maximum likelihood estimates
if ( do.smooth ){
   # WARNING:  This spline does not necessarily pass
   #           through the given points.  This should
   #           be seen as a deficiency in the method.

   if (verbose>=2) cat(" Estimating phi:  ")
   
   # Find maximum from profile

   L.max <- max(y, na.rm=TRUE)
   xi.max <- x[ y==L.max ]
   
	# In some unusual cases, when add0=TRUE, the maximum of p/xi
	# may be between 0 and 1.
	# In this situation, we... set the mle to either 0 or 1,
	# whichever gives the larger values of the log-likelihood
	if ( (xi.max > 0)  & ( xi.max < 1 ) ) {
		L.max <- max( c( y[xi.vec==0], y[xi.vec==1]) )
		xi.max <- xi.vec[ L.max == y ]
		cat("MLE of",index.par,"is between 0 and 1, which is impossible.",
		    "Instead, the MLE of",index.par,"has been set to",xi.max,
		    ".  Please check your data and the call to tweedie.profile().")
	}

   # Now need to find mle of  phi  at this very value of  p
   # - Find bounds
   
   phi.1 <-   2 * max( phi.vec.fix, na.rm=TRUE )
   phi.2 <- 0.5 * min( phi.vec.fix, na.rm=TRUE )
   
	if ( phi.1 > phi.2 ) {
	   phi.hi <- phi.1
		phi.lo <- phi.2 
	} else {
	   phi.hi <- phi.2
		phi.lo <- phi.1
	}
   
   if (verbose>=2) cat(" Done\n")
	
   # Now, if the maximum happens to be at an endpoint,
   # we have to do things differently:
	if ( (xi.max==xi.vec[1]) | (xi.max==xi.vec[length(xi.vec)]) ) {
   
   	if ( xi.max==xi.vec[1]) phi.max <- phi.vec[1]
   	if ( xi.max==xi.vec[length(xi.vec)]) phi.max <- phi.vec[length(xi.vec)]
#      phi.max <- phi.lo # and is same as phi.max
      warning("True maximum possibly not detected.")

   } else {

      # - solve

      if ( phi.method=="saddlepoint"){
         
         mu <- fitted( glm.fit( y=ydata, x=model.x, weights=weights, offset=offset,
                                family=tweedie(xi.max, link.power=link.power)))
         phi.max <- sum( tweedie.dev(y=ydata, mu=mu, power=xi.max) )/length( ydata )
         
      } else {
#        phi.max <- nlm(p=phi.est, f=dtweedie.nlogl, 
#                    hessian=FALSE,
#                    power=p, mu=mu, y=data, 
#                    gradient=dtweedie.dldphi)$estimate
        mu <- fitted( glm.fit( y=ydata, x=model.x, weights=weights, offset=offset,
                               family=tweedie(xi.max, link.power=link.power)))
        phi.max <- optimize( f=dtweedie.nlogl, maximum=FALSE, interval=c(phi.lo, phi.hi ), 
            # set lower limit phi.lo???
                             power=xi.max, mu=mu, y=ydata)$minimum
       
#       Note that using the Hessian often produces computational problems
#       for little (no?) advantage.
      
      }
   }
} else {
   # Find maximum from computed data
   if (verbose>=2) cat(" Finding maximum likelihood estimates:  ")
  
   L.max <- max(L)

   xi.max   <- xi.vec  [ L == L.max ]
   phi.max <- phi.vec[ L == L.max ]

   if (verbose>=2) cat(" Done\n")
}

# Now report
if ( verbose >= 2 ) {
   cat( "ML Estimates:  ",index.par,"=",xi.max," with phi=",phi.max," giving L=",L.max,"\n")
   cat(" ---\n")
}

# Now find approximate, nominal 95% CI info
ht <- L.max - ( qchisq(conf.level, 1) / 2 )
ci <- array( dim=2, NA )


if ( do.ci ) {
    if (verbose==2) cat("* Finding confidence interval:")
   if ( !do.smooth ) {
      warning("Confidence interval may be very inaccurate without smoothing.\n")
		y <- L
		x <- xi.vec
   }

   if ( do.plot ) {
      abline(h=ht, lty=2)
      title( sub=paste("(",100*conf.level,"% confidence interval)", sep="") )
   }

   # Now find limits on p:
   # - fit a smoothing spline the other way:  based on L predicting p
   #   so we can ensure that the limits on  p  are found OK

   # --- Left side ---
   cond.left <- (y < ht ) & (x < xi.max )
   if ( all(cond.left==FALSE) ) {
      warning("Confidence interval cannot be found: insufficient data to find left CI.\n")
   }else{
		# Now find it!
      approx.left <- max( x[cond.left] )
      index.left <- seq(1, length(cond.left) )
      index.left <- index.left[x==approx.left]
      left.left <- max( 1, index.left - 5)
      left.right <- min( length(cond.left), index.left + 5)

      # Left-side spline
      ss.left <- splinefun( y[left.left:left.right], x[left.left: left.right] )
      ci.new.left <- ss.left( ht )

      ci[1] <- ci.new.left
   }

   # --- Right side ---
   cond.right <- (y < ht ) & (x > xi.max )
   if ( all( cond.right==FALSE ) ) {
      warning("Confidence interval cannot be found: insufficient data to find right CI.\n")
   }else{
     approx.right <- min( x[cond.right] )
      index.right <- seq(1, length(cond.right) )
      index.right <- index.right[x==approx.right]
      right.left <- max( 1, index.right - 5)
      right.right <- min( length(cond.left), index.right + 5)

      # Right-side spline
      ss.right <- splinefun( y[right.left:right.right], x[right.left: right.right] )
      ci.new.right <- ss.right(ht )

      ci[2] <- ci.new.right
   }
    if (verbose==2) cat(" Done\n")
}

if ( fit.glm ) {
   out.glm <- glm.fit( x=model.x, y=ydata, weights=weights, offset=offset,
         family=tweedie(var.power=xi.max, link.power=link.power) )

	if ( xi.notation){
		out <- list( y=y, x=x, ht=ht, L=L, xi=xi.vec, xi.max=xi.max, L.max=L.max, 
           phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method,
           glm.obj = out.glm)
	} else {
		out <- list( y=y, x=x, ht=ht, L=L, p=p.vec, p.max=xi.max, L.max=L.max, 
					phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method,
					glm.obj = out.glm)
	}
} else {
	if (xi.notation ){
		out <- list( y=y, x=x, ht=ht, L=L, xi=xi.vec, xi.max=xi.max, L.max=L.max, 
					phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method)
	
	} else {
		out <- list( y=y, x=x, ht=ht, L=L, p=p.vec, p.max=xi.max, L.max=L.max, 
					phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method)
	}

}

if ( verbose ) cat("\n")

invisible( out )

}



#############################################################################
dtweedie.stable <- function(y, power, mu, phi)
{
         # Error checks
         if ( power<1) stop("power must be greater than 2.\n")
         if ( any(phi<=0) ) stop("phi must be positive.")
         if ( any(y<0) ) stop("y must be a non-negative vector.\n")
         if ( any(mu<=0) ) stop("mu must be positive.\n")
         if ( length(mu)>1) {
            if ( length(mu)!=length(y) ) stop("mu must be scalar, or the same length as y.\n")
         }
         else {
            mu <- array( dim=length(y), mu )
            # A vector of all mu's
         }
         if ( length(phi)>1) {
            if ( length(phi)!=length(y) ) stop("phi must be scalar, or the same length as y.\n")
         }
         else {
            phi <- array( dim=length(y), phi )
            # A vector of all phi's
         }
         
        density <- y
	alpha <- (2-power)/(1-power)
	beta <- 1
	k <- 1  # The parameterization used
	delta <- 0
	gamma <- phi * (power-1) * 
                    ( 1/(phi*(power-2)) * cos( alpha * pi / 2 ) ) ^ (1/alpha)
        
        require(stabledist) # Needed for  dstable
        ds <- dstable(y, alpha=alpha, beta=beta, gamma=gamma, delta=delta, pm=k)
        density <- exp((y*mu^(1-power)/(1-power)-mu^(2-power)/(2-power))/phi)*ds

        
        density
}


#############################################################################
tweedie.plot <- function(y, xi=power, mu, phi, type="pdf", power=NULL, 
                         add=FALSE, ...) {

   if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
   xi.notation <- TRUE
   if ( is.null(power) ) power <- xi
   if ( is.null(xi) ) {
      xi <- power
      xi.notation <- FALSE
   }
   if ( xi != power ) {
      cat("Different values for xi and power given; the value of xi used.\n")
      power <- xi
   }
   index.par <- ifelse( xi.notation, "xi","p")

   if ( ( power < 0 ) | ( ( power > 0 ) & ( power < 1 ) ) ) {
      stop("Plots cannot be produced for power =",power,"\n")
   }   
   
   is.pg <- ( power > 1 ) & ( power < 2 )
   
   if ( type=="pdf") {
      fy <- dtweedie( y=y, power=power, mu=mu, phi=phi)
   } else {
      fy <- ptweedie( q=y, power=power, mu=mu, phi=phi)
   }
   
   if ( !add ) {
      if ( is.pg ) {
         plot( range(fy) ~ range( y ), 
            type="n", ...)
         if ( any( y==0 ) ) { # The exact zero
            points( fy[y==0] ~ y[y==0], pch=19, ... )
         }
         if ( any( y>0 ) ) { # The exact zero
            lines( fy[y>0]   ~ y[y>0], pch=19, ... )
         }
      } else {  # Not a Poison-gamma dist
         plot( range(fy) ~ range( y ), 
            type="n", ...)
         lines( fy ~ y, pch=19, ... )
      }
   } else {# Add; no new plot
      if ( is.pg ) {
         if ( any( y==0 ) ) { # The exact zero
            points( fy[y==0] ~ y[y==0], pch=19, ... )
         }
         if ( any( y>0 ) ) { # The exact zero
            lines( fy[y>0]   ~ y[y>0], pch=19, ... )
         }
      } else {  # Not a Poison-gamma dist
         lines( fy ~ y, pch=19, ... )
      }
   
   }
   return(invisible(list(y=fy, x=y) ))
   
}
#############################################################################
AICtweedie <- function( glm.obj, k=2){
	
	wt <- glm.obj$prior.weights
	n <- length(glm.obj$residuals)
	edf <- n - glm.obj$df.residual

	dev <- deviance(glm.obj)
	disp <- dev/n

	mu <- fitted( glm.obj )
	y  <- glm.obj$y
	p <- get("p", envir = environment(glm.obj$family$variance))

	den <- dtweedie( y=y, mu=mu, phi=disp, power=p)
	AIC <- -2*sum( log(den) * wt) 

	return( AIC + k*(edf+1) )

}



