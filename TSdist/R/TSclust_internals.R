
#Internal functions taken directly from TSclust for pairwise 
#calculations of some functions.


#Internal function needed for pairwise calculation of the pairwise.interp.specDistance
#function.

interp.SPEC.LS <- function(x, n)  {
  pgram <- spec.pgram(x, plot=FALSE)
  Yk <- pgram$spec
  lambdas <- pgram$freq
  h <- dpill(lambdas, Yk)
  interplambdas <- seq(min(lambdas), max(lambdas), length.out=n)
  ys <- leastsquares.spec(Yk, lambdas, h, interplambdas)
  ys[ys<0.0001] <- 0.0001 #no zeroes allowed
  list( x = interplambdas, y = ys)
}

leastsquares.spec <- function( Yk, lambdas, h, lambdaeval) {
  d <- data.frame(Yk)
  d$lambdas <- lambdas
  lp <- locpol(Yk~lambdas, d, bw=h,kernel=gaussK, xeval=lambdaeval )
  ord <- order(lambdaeval) #trick to get the original order of the lambas, locpol sorts the input vector
  ord2 <- order(ord)   #second part of the trick
  lp$lpFit$Yk[ord2]
}

interp.W.LK <- function(x, n) {
  interps <- interp.SPEC.LOGLIKELIHOOD(x, n)
  interps$y <- exp(interps$y)
  interps
}

interp.SPEC.LOGLIKELIHOOD <- function(x, n) {
  pgram <- spec.pgram(x, plot=FALSE)
  Yks <- log(pgram$spec)
  lambdas <- pgram$freq
  interplambdas <- seq(min(lambdas), max(lambdas), length.out=n)
  hX <- 0.93*dpill(lambdas, Yks)
  approx( interplambdas, .vectorized.lk.optim( interplambdas, Yks, lambdas, hX ) )
}

interp.SPEC.GLK <- function(x, n) {
  pgram <- spec.pgram(x, plot=FALSE)
  Yk <- log(pgram$spec)
  lambdas <- pgram$freq
  h <- 0.93*dpill(lambdas, Yk)
  ys <- .vectorized.lk.optim( lambdas, Yk, lambdas, h )
  list( x = lambdas, y = list(mu = ys, Z = Yk) )
}

#####Integration functions

#Integration function for GLK
integrate.GLK <- function( base, x, y) {
  Z <- x$Z - y$Z
  mu <- x$mu - y$mu
  sum(Z - mu - 2*log(1 + exp(Z - mu))) - sum( Z - 2*log(1 + exp(Z)))
}

#Integration function for ISD
integrate.ISD <- function(base, x, y) {
  trapez(base, (x - y)^2)
}

#trapezoind integration
trapez <- function(x,y) {
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

#Trapezoid integration of W divergence.
integrate.divergenceW <- function(base, x, y, alpha ) {
  val <- simetricDivergenceW( x / y, alpha)
  trapez(base, val)
}

#Divergence functions

#W divergence
divergenceW <- function(x, alpha) {
  if ((alpha > 0) & (alpha< 1)) {
    log(alpha*x + 1 - alpha) - alpha*log(x)
  }
  else {
    stop("condition 0 < alpha < 1 not fulfilled")
  }
}

#Symmetric divergence function.
simetricDivergenceW <- function(x,alpha) {
  divergenceW(x,alpha) + divergenceW(1/x,alpha)  
}


#Maximum likelihood functions functions needed for calculation of distances based
#on spectra. 

Spectral.AB <- function (  ABvec, lambda, Yks, lambdas, h) {
  acum <- 0
  a <- ABvec[1]
  b <- ABvec[2]
  -sum ( ( -exp(Yks - a -b*(lambdas - lambda) ) + Yks - a -b*(lambdas - lambda) ) * funcionKh( lambdas - lambda, h) )
}

likelihood.optim <- function(  lambda, Yks, lambdas, h) {
  startA = lambdas[1]
  startB = 0
  optim(c(startA,startB), Spectral.AB, lambda=lambda, Yks = Yks, lambdas=lambdas, method="L-BFGS-B",lower = c(min(Yks) , -101), upper= c(max(Yks),101), h =h)$par[1]
}

.vectorized.lk.optim <- Vectorize(likelihood.optim,"lambda")  #needed for likelihood.optim to accept a vector, required for integrate


#Gaussian kernel function
funcionKh <- function( value, h) {
  value <- value/h
  dnorm( -(value**2) ) / h
}


#Internal functions neccessary to calculate predDistance of databases. 

#Integrated L1 distance
L1dist <- function( density.x, density.y, bw.x, bw.y) {
  r <- range(density.x[,1], density.y[,1])
  a <- r[1] - 0.025*(r[2]-r[1]) ; b <- r[2] + 0.025*(r[2]-r[1])
  #L1 function between forecasts.
  integrand_L1 <- function(u) {
    abs(estimator.density( density.x, u, bw.x ) - estimator.density( density.y, u, bw.y ) )
  }
  DL1 <- 2
  tryCatch( {
    DL1 <- integrate( integrand_L1, lower=a, upper=b)$value }, error = function(e) {
      plot.default( density.x, type="l", col="red", lwd=2, main="Error, problem intergrating the L1 distance of these densities, approximating...", xlab="", ylab="Density",
                    xlim=c(min(density.x[,1], density.y[,1]), max(density.x[,1], density.y[,1])), ylim=c(0, max(density.x[,2], density.y[,2]) ) )
      lines( density.y, col="blue", lwd=2 )
      legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
      DL1 <- integrate( integrand_L1, lower=a, upper=b, stop.on.error=FALSE)$value
    })
  DL1
}

# From a sample "data" and a sequence of bandwidths "h.seq", this procedure 
# computes the kernel density function evaluated on s set of points "x".
estimator.density <- function( data, x, h.seq, kind.of.kernel="dnorm"){
  
  result <- approx(data[,1], data[,2], x)$y
  result[ is.na(result)] <- 0
  result
  
}