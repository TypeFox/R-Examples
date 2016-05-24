
# A Poisson Subordinated Distribution By Lihn (LIHNPSD)

library(moments)
library(BB) # spg
library(Bolstad2)
library(optimx)
library(Rmpfr) # Multiple Precision Float, erf, erfc are defined here !
library(sn) # skew normal

.N <- function(x,p=128) mpfr(x, precBits = p)


# library(codetools)
# findGlobals(calcqq)
# ------------------------------------------
# LIHNPSD_UnitTest(): for unit testing

prepare <- function(d) UseMethod("prepare", d)
rawmean <- function(d) UseMethod("rawmean", d)
rawmu1  <- function(d) UseMethod("rawmu1",  d)
rawmu2  <- function(d) UseMethod("rawmu2",  d)
rawmu3  <- function(d) UseMethod("rawmu3",  d)
rawmu4  <- function(d) UseMethod("rawmu4",  d)

rawsn       <- function(d, type, x, k) UseMethod("rawsn", d)
density     <- function(d,x) UseMethod("density",   d)
rawdensity  <- function(d,x) UseMethod("rawdensity", d)
rawdensity_kth   <- function(d,x,k) UseMethod("rawdensity_kth", d)
rawdensityslope  <- function(d,x) UseMethod("rawdensityslope", d)
rawdensity0 <- function(d)   UseMethod("rawdensity0", d)
rawcdf      <- function(d,x) UseMethod("rawcdf", d)
rawcdfinv <- function(d,c,xinit) UseMethod("rawcdfinv", d)

tailindex   <- function(d,x) UseMethod("tailindex", d)
tailindex_plot <- function(d,xmin,xmax,ymax=0.4) UseMethod("tailindex_plot", d)

psdkernel    <- function(d,k) UseMethod("psdkernel", d)
psdmagnitude <- function(d,r,step=0.1) UseMethod("psdmagnitude", d)
poisson_sum  <- function(d,fn) UseMethod("poisson_sum", d)
poisson_sum_kth <- function(d,fn,k) UseMethod("poisson_sum_kth", d)


psdunittest  <- function(d) UseMethod("psdunittest", d)
generatepdf <- function(d,NS,NT,raw=1) UseMethod("generatepdf", d)
calcqq <- function(d,hq,step=5, debug=0) UseMethod("calcqq",  d)

mu1_analytic <- function(d) UseMethod("mu1_analytic", d)
mu2_analytic <- function(d) UseMethod("mu2_analytic", d)
mu3_analytic <- function(d) UseMethod("mu3_analytic", d)
mu4_analytic <- function(d) UseMethod("mu4_analytic", d)
mu_n_core    <- function(d,n) UseMethod("mu_n_core", d)
psdvariance <- function(d) UseMethod("psdvariance", d)
psdskewness <- function(d) UseMethod("psdskewness", d)
psdkurtosis <- function(d) UseMethod("psdkurtosis", d)

standardfit <- function(d, r, hist, trace, iter, plotqq, weights, merge_tails) UseMethod("standardfit", d)

# when this is in a package
#formals(sd.default) <- c(formals(sd.default), alist(... = ))

# ------------------------------------------
# simple PSD constructor
SPSD <- function(sigma, alpha, gamma, beta=0, mpfr=0) {
  d <- list( sigma= sigma, alpha= alpha, gamma= gamma, beta= beta, mpfr= mpfr )
  class(d) <- "LIHNPSD"
  rawmean(d)
}

# the input can be numeric or character
# notice that 0.05 in double is not exactly that in mpfr
prepare.LIHNPSD <- function(d) {
  # you have to provide sigma, alpha and gamma
  stopifnot( ! is.null(d$sigma) )
  stopifnot( ! is.null(d$alpha) )
  stopifnot( ! is.null(d$gamma) )

  if (is.null(d$mpfr)) d$mpfr <- 0
  D <- function(x) {
    if ( d$mpfr > 0 ) {
      if ( class(x) != "mpfr" ) mpfr(x, precBits = floor(d$mpfr))
      else x # already mpfr
    }
    else as.numeric(x)
  }

  if (is.null(d$location)) d$location <- 0.0
  if (is.null(d$beta)) d$beta <- 0.0
  if (is.null(d$lambda)) d$lambda <- 1.0
  if (is.null(d$epsilon)) d$epsilon <- 1e-10
  d$pi <- if ( d$mpfr > 0 ) Const("pi", prec= d$mpfr) else d$pi <- pi

  d$sigma <- D(d$sigma)
  d$alpha <- D(d$alpha)
  d$gamma <- D(d$gamma)
  d$location <- D(d$location)
  d$beta  <- D(d$beta)
  d$lambda <- D(d$lambda)
  d$epsilon <- D(d$epsilon)

  d
}

# ------------------------------------------

# PoissonTerm <- function(lambda,k) dpois(k,lambda) # doesn't work in MPFR 
PoissonTerm <- function(lambda,k) lambda^k*exp(-lambda)/factorial(k)  

LihnFunctionKth <- function(alpha,lambda,k) {
  PoissonTerm(lambda,k)*(k+1)^alpha
}
LihnFunctionSum <- function(alpha,x,epsilon = 1e-10) {
  # s: sum, k: poisson index, v: value of k-th term
  s <- k <- v <- 0.0 * alpha # multiply by alpha to inherit its type

  while ( k < 20 || v > epsilon*s || v == 0 ) {
    v <- LihnFunctionKth(alpha,x,k)
    s <- s + v
    k <- k + 1
  }
  s
}

LihnFunctionAnalytic <- function(alpha,x) {
  if ( x <= 0 ) stopifnot( paste("Unsupported x:",as.numeric(x)) == "STOP" )
  if ( alpha == 1 ) (x+1)
  else if ( alpha == 2 ) (x^2+3*x+1)
  else if ( alpha == 3 ) (x^3+6*x^2+7*x+1)
  else if ( alpha == 4 ) (x^4+10*x^3+25*x^2+15*x+1)
  else if ( alpha == 0 ) 1
  else if ( alpha == -1 && x != 0 ) (1-exp(-x))/x
  else stopifnot( paste("Unsupported alpha:",as.numeric(alpha)) == "STOP" )
}

LihnBetaPoly <- function(N,b) {
  if ( N == 1 ) b
  else if ( N == 2 ) 1
  else if ( N == 4 ) 3
  else if ( N == 6 ) 15
  else {
    ppi <- if ( class(b) != "mpfr" ) pi 
           else Const("pi", prec= getPrec(b))
    if ( N == 3 ) (-ppi/2*b^3+3*b)
    else if ( N == 5 ) (3*ppi^2/4*b^5-5*ppi*b^3+15*b)
    else stopifnot( paste("Unsupported N:",as.numeric(N)) == "STOP" )
  }
}

LihnTildeFunction <- function(x,alpha,p,epsilon = 1e-10) {
  # s: sum, k: poisson index, v: value of k-th term
  s <- k <- v <- 0.0 * alpha # multiply by alpha to inherit its type
  while ( k < 20 || v > epsilon*s || v == 0 ) {
    sigma <- (k+1)^alpha
    v <- (k+1)^p/factorial(k) * exp(-x^2/(2*sigma^2))
    s <- s + v
    k <- k + 1
  }
  C <- if ( class(alpha) == "mpfr" ) 
          exp(-1+0*alpha)/sqrt(2* Const("pi", prec= getPrec(alpha))) 
       else exp(-1)/sqrt(2*pi)
  s * C
}

ParetoDiff2 <- function(x,alpha) {
  L5 <- LihnTildeFunction(x,alpha,-5*alpha)
  L3 <- LihnTildeFunction(x,alpha,-3*alpha)
  A <- (x^2*L5- L3)
  L1 <- LihnTildeFunction(x,alpha,-alpha)
  L3 <- LihnTildeFunction(x,alpha,-3*alpha);
  A/L1 - x^2 * L3^2/L1^2
}

LihnFunctionValidate <- function() {
  # we will validate Lihn function at alpha =-1,0,1,2,3,4
  err <- 1e-6
  for (a in c(-1,0,1,2,3,4)) {
    for( x in c(1,2,3,4,5) ) {
      # print(paste("LihnFunction: validating for a=",
      #   as.numeric(a),"x=",as.numeric(x)))
      stopifnot( abs(LihnFunctionSum(a,x) - LihnFunctionAnalytic(a,x)) < err )  
    }
  }
  print("LihnFunction is validated")

  stopifnot( abs( LihnTildeFunction(0,0,0) - 1/sqrt(2*pi) ) < err ) 
  stopifnot( abs( LihnTildeFunction(1,0,0) - 0.24197072451914 ) < err ) 
  stopifnot( abs( LihnTildeFunction(1,0.5,-0.5) - 0.21950993766309 ) < err ) 
  stopifnot( abs( ParetoDiff2(0,0)   + 1 ) < err ) 
  stopifnot( abs( ParetoDiff2(1,0)   + 1 ) < err ) 
  stopifnot( abs( ParetoDiff2(1,0.5) + 0.57402540354604 ) < err ) 
  stopifnot( abs( ParetoDiff2(10,0.62) + 0.018150532857661 ) < err ) 

  print("LihnTildeFunction is validated");
}


# ------------------------------------------

psdkernel.LIHNPSD <- function(d,k) ((k+1)^d$alpha) * ((1+d$gamma)^k)

# return the k magnitude of the return series
psdmagnitude.LIHNPSD <- function(d, r, step=0.1) {
  ks <- c()
  for ( i in 1:length(r) ) {
    j <- 0
    repeat {
      if ( abs(r[i]) <= psdkernel(d,j) * d$sigma ) {
        ks[i] <- j
        break
      }
      j <- j + step
    }
  }
  ks
}

# This is the utility to build moments
poisson_sum_kth.LIHNPSD <- function(d, fn, k) {
  PoissonTerm(d$lambda, k) * fn(d,k)
}
poisson_sum.LIHNPSD <- function(d, fn) {
  # sm: sum, m: k-th term of sum
  sm <- k <- m <- 0.0*d$sigma # multiply by sigma to inherit its type
  kmax <- if ( class(d$sigma) == "mpfr" ) 100000 else 160
  cond <- TRUE
  while ( cond && k <= kmax ) {
    m  <- poisson_sum_kth(d, fn, k)
    sm <- sm + m
    cond <- ( k < 20 
            || ( abs(m) >= abs(d$epsilon*sm) && sm != 0.0 )
            || ( sm == 0.0 && k < 20 ) )
    # print(paste(cond,"k=",as.numeric(k),"sm=",as.numeric(sm),
    #  "m=",as.numeric(m),"err=",as.numeric(d$epsilon)))
    k  <- k + 1
  }
  sm
}

# ------------------------------------------
rawmu1.LIHNPSD <- function(d) {
  mu1core <- function(d,k) ( psdkernel(d, k) * d$sigma )
  poisson_sum(d, mu1core) * d$beta
}
rawmean.LIHNPSD <- function(d) {
  d <- prepare(d)
  d$rawmean <- rawmu1(d)
  d
}
rawmu2.LIHNPSD <- function(d) {
  mu2core <- function(d,k) {
    sigma2 <- psdkernel(d, k) * d$sigma 
    sigma2^2
  }
  poisson_sum(d, mu2core)
}
rawmu3.LIHNPSD <- function(d) {
  mu3core <- function(d,k) {
    sigma2 <- psdkernel(d, k) * d$sigma 
    sigma2^3
  }
  poisson_sum(d, mu3core) * LihnBetaPoly(3,d$beta)
}
rawmu4.LIHNPSD <- function(d) {
  mu4core <- function(d,k) {
    sigma2 <- psdkernel(d, k) * d$sigma 
    sigma2^4
  }
  poisson_sum(d, mu4core) * LihnBetaPoly(4,d$beta)
}
# ------------------------------------------
# also raw moments, but expressed in Lihn functions
mu_n_core.LIHNPSD <- function(d, n) {
  L <- LihnFunctionSum(n*d$alpha, d$lambda*(1+d$gamma)^n, d$epsilon)
  d$sigma^n * exp(d$lambda*((1+d$gamma)^n-1.0)) * L # muN core component
}

mu1_analytic.LIHNPSD <- function(d) {
  d <- prepare(d)
  L <- LihnFunctionSum(d$alpha, d$lambda*(1+d$gamma), d$epsilon)
  mu1 <- d$beta * d$sigma * exp(d$lambda*d$gamma) * L
  mu1
}
mu2_analytic.LIHNPSD <- function(d) {
  d <- prepare(d)
  L <- LihnFunctionSum(2*d$alpha, d$lambda*(1+d$gamma)^2, d$epsilon)
  mu2 <- d$sigma^2 * exp(d$lambda*d$gamma*(d$gamma+2)) * L
  mu2
}
mu3_analytic.LIHNPSD <- function(d) {
  d <- prepare(d)
  mu3 <- LihnBetaPoly(3,d$beta) * mu_n_core(d,3)
  mu3
}
mu4_analytic.LIHNPSD <- function(d) {
  d <- prepare(d)
  mu4 <- LihnBetaPoly(4,d$beta) * mu_n_core(d,4)
  mu4
}

# ------------------------------------------
psdvariance.LIHNPSD <- function(d) {
  rawmu2(d) - d$rawmean^2
}
# http://en.wikipedia.org/wiki/Moment_about_the_mean#Relation_to_moments_about_the_origin
psdskewness.LIHNPSD <- function(d) {
  mu1 <- d$rawmean
  (rawmu3(d) - 3*mu1*rawmu2(d) + 2*mu1^3 ) / psdvariance(d)^(3/2)
}
psdkurtosis.LIHNPSD <- function(d) {
  mu1 <- d$rawmean
  (rawmu4(d) - 4*mu1*rawmu3(d) + 6*mu1^2*rawmu2(d) - 3*mu1^4 ) / psdvariance(d)^2 - 3.0
}

# ------------------------------------------
# PDF offset by mu1+d$location
density.LIHNPSD <- function(d, x) {
  d2 <- d
  d2$location <- 0
  rawdensity(d2, x - d$rawmean - d$location)
}
# PDF without any offset (not even d$location)
# k >= 0 (int) is used to study the kth term contribution
rawdensity_kth.LIHNPSD <- function(d, x, k) {
  gauss <- function(d,k2) {
    sigma2 <- psdkernel(d, k2) * d$sigma 
    mn <- d$beta * sigma2
    gs <- exp(-(x-mn)^2/(2*sigma2^2))/sqrt(2*d$pi)/sigma2
    gs
  }
  poisson_sum_kth(d, gauss, k)
}
# rawsn is an internal function
# type 1: for pdf, 2: for 2nd term of dP/dx
rawsn.LIHNPSD <- function(d, type, x, k) {
  sigma2 <- psdkernel(d, k) * d$sigma 
  delta <- d$beta * sqrt(d$pi/2);
  if ( abs(delta) >= 1 ) stopifnot( paste("d$beta is out of range:",as.numeric(d$beta)) == "STOP" ) 

  if ( type == 1 ) {
    sqrt2 <- sqrt( d$lambda / d$lambda * 2 ) # accomodate mpfr
    sn_a  <- delta/sqrt(1-delta^2);
    gs <- exp(-x^2/(2*sigma2^2))/sqrt(2*d$pi)/sigma2
    gcdf <- (1+erf(sn_a * x/sqrt2/sigma2))/2
    gs * 2*gcdf
  }
  else if ( type == 2 ) {
    x2 <- x * sqrt( 2/d$pi ) / sqrt( 2/d$pi - d$beta^2 )
    gs <- exp(-x2^2/(2*sigma2^2))/sqrt(2*d$pi)/sigma2
    gcdf <- (1+erf(0))/2
    gs * 2*gcdf
  }
  else stopifnot( paste("Unsupported type:",as.numeric(type)) == "STOP" )
}
rawdensity.LIHNPSD <- function(d, x) {
  sqrt2 <- sqrt( d$lambda / d$lambda * 2 ) # accomodate mpfr
  delta <- d$beta * sqrt(d$pi/2);
  if ( abs(delta) >= 1 ) stopifnot( paste("d$beta is out of range:",as.numeric(d$beta)) == "STOP" ) 
  sn_a  <- delta/sqrt(1-delta^2);
  gauss <- function(d,k) rawsn(d, 1, x, k)
  poisson_sum(d, gauss)  
}
rawdensityslope.LIHNPSD <- function(d, x) { 
  dpk_dx <- function(d,k) {
    sigma2 <- psdkernel(d, k) * d$sigma 
    t1 <- -rawsn(d, 1, x, k) * ( x / sigma2^2 ) 
    b2 <- d$beta / sigma2 * sqrt( 2/d$pi ) / sqrt( 2/d$pi - d$beta^2 )
    t2 <- b2 * rawsn(d, 2, x, k)
    t1 + t2
  }
  poisson_sum(d, dpk_dx)
}

rawcdf.LIHNPSD <- function(d, x) {
  # erf <- function(x) 2*pnorm(x*sqrt(2))-1
  # we are using erf from rmpfr package
  sqrt2 <- sqrt( d$lambda / d$lambda * 2 ) # accomodate mpfr
  delta <- d$beta * sqrt(d$pi/2);
  if ( abs(delta) >= 1 ) stopifnot( paste("d$beta is out of range:",as.numeric(d$beta)) == "STOP" ) 
  sn_a  <- delta/sqrt(1-delta^2);
  cdf <- function(d,k) {
    sigma2 <- psdkernel(d, k) * d$sigma 
    gcdf <- (1+erf(x/sqrt2/sigma2))/2
    gcdf - 2*T.Owen(x/sigma2, sn_a)
  }
  poisson_sum(d, cdf)
}
# xinit: if you don't know what to supply, use d$rawmean
rawcdfinv.LIHNPSD <- function(d, c, xinit) {
  err <- d$epsilon
  x  <- xinit
  x1 <- d$rawmean # a fake number
  repeat {
    pdf <- rawdensity(d,x)
    cdf <- rawcdf(d,x)
    x1 <- x - (cdf-c)/pdf
    if ( abs(x1-x) < err ) break
    x <- x1
  }
  x1
}

rawdensity0.LIHNPSD <- function(d) {
  gs <- 1/sqrt(2*d$pi)/d$sigma 
  f  <- exp(-d$lambda*d$gamma/(d$gamma+1))
  L <- LihnFunctionSum( -d$alpha, d$lambda/(d$gamma+1) )
  gs * f * L
}

tailindex.LIHNPSD <- function(d,x) {
  cdf <- rawcdf(d,x)
  pdf <- rawdensity(d,x)
  dff <- rawdensityslope(d,x)
  -1-(1-cdf)*dff/pdf^2
}

tailindex_plot.LIHNPSD <- function(d,xmin,xmax,ymax=0.4) {
  N <- 1000
  x <- xmin+(0:(N-1))*(xmax-xmin)/N
  idx <- tailindex(d,x)
  ymax <- min(abs(ymax), max(abs(as.numeric(idx))*2))
  plot( x, idx, pch = ".", col="red",
    ylim=c(-ymax,+ymax),
    ylab="Tail Index(x)", xlab="x", main="Tail Index")
  abline(h=0)
  # text to describe d
}
# -----------------------------------------
generatepdf.LIHNPSD <- function(d, NS, NT, raw=1) {
  d <- prepare(d)
  # NS: number of sigma, NT: number of tick sample per unit of sigma
  XM <- d$sigma*NS # max of x
  N  <- NS*NT # number of samples in x-axis
  x  <- seq(-XM, XM, length=N+1)
  dx <- (2*XM)/N
  pb <- if ( raw==1 ) rawdensity(d,x) else density(d,x)

  list( NS = NS, NT = NT, XM = XM,
    N = N, x  = x, dx = dx, pb = pb
  )
}

# hq = list( qhx = h$mids, qhy = h$counts )
calcqq.LIHNPSD <- function(d, hq, step=5, debug=0 ) {

  tm0 <- unclass(Sys.time())

  hq$qhyc <- cumsum(hq$qhy) / sum(hq$qhy)
  qqp <- list ( x = c(), xq = c(), y = c(), yq = c() )
  xmin <- min(hq$qhx)
  xmax <- max(hq$qhx)
  xabs <- max(abs(xmin),abs(xmax))
  tx <- d$location # starting point of theoretical x
  dx <- ( hq$qhx[2] - hq$qhx[1] )/step
  while ( tx >= hq$qhx[1] ) tx <- tx-dx
  stopifnot( hq$qhyc[1] != 0 )

  tm1 <- unclass(Sys.time()) - tm0

  while( rawcdf(d, tx-d$rawmean-d$location) >= hq$qhyc[1] 
         && abs(tx) < 5*xabs ) tx <- tx-dx
  if ( debug > 0 ) print(paste("tx=",as.numeric(tx)))

  tm2 <- unclass(Sys.time()) - tm0

  tcdf <- 0
  cdfcnt <- 0
  for ( i in 1:(length(hq$qhy)-1) ) {
    qqp$x[i] <- hq$qhx[i] # data's x
    qqp$xq[i] <- qval <- hq$qhyc[i] # data's cdf(x)
    # the goal is to find the min(x) of fit whose cdf >= data's cdf
    repeat {
      tcdf <- rawcdf(d, tx-d$rawmean-d$location)
      cdfcnt <- cdfcnt + 1
      if ( tcdf >= qval || abs(tx) > 5*xabs ) break
      tx <- tx+dx
    }
    qqp$y[i]  <- tx # fit's x
    qqp$yq[i] <- tcdf # fit's cdf(x)
    if ( debug == 2 ) print(paste("qx=",as.numeric(qqp$x[i]),
      "qval=",as.numeric(qval),"tcdf=",as.numeric(tcdf),"tx=",as.numeric(tx)))
  }
  if ( debug > 0 ) print(paste("qqp=","tcdf=",as.numeric(tcdf),"tx=",as.numeric(tx)))
  tm3 <- unclass(Sys.time()) - tm0
  qqp$tm1 <- tm1
  qqp$tm2 <- tm2
  qqp$tm3 <- tm3
  print(paste("qqp tm1=",sprintf("%.2f",tm1),
    "tm2=",sprintf("%.2f",tm2),
    "tm3=",sprintf("%.2f",tm3),
    "cdfcnt=",cdfcnt))

  qqp
}

# ------------------------------------------
# Unit Test Functions
# ------------------------------------------
psdunittest.LIHNPSD <- function(d) {
  d <- prepare(d)
  d <- rawmean(d)
  err <- 1e-6
  errn <- 1e-6 # for the extreme case
  dgs <- 7 # round
  Q <- function(x) round(as.numeric(x),dgs) # helper
  TM <- function() format(Sys.time(), "%X %x") # helper

  print(paste("--Unit Test on: sigma=",as.numeric(d$sigma),
    "alpha=",Q(d$alpha),
    "gamma=",Q(d$gamma),
    "beta=",Q(d$beta),
    "lambda=",Q(d$lambda),
    "mpfr=",Q(d$mpfr),"--"))

  # when symmetric, test a few obvious facts 
  if ( d$beta == 0 ) {
    # cdf at x=0 is 1/2
    cdf0 <- rawcdf(d,0)
    print(paste("cdf_at_0 =",Q(cdf0)))
    stopifnot( abs(cdf0 - 0.5) < 1e-4 )

    # d pdf/dx at x=0 is 0
    stopifnot( rawdensityslope(d,0) == 0 )
  }
  # test location (only applicable to density, not rawdensity)
  d$location <- d$sigma
  stopifnot( abs( density(d,d$sigma) - rawdensity(d,d$sigma-d$rawmean-d$location) ) < err )
  d$location <- 0
  stopifnot( abs( density(d,d$sigma) - rawdensity(d,d$sigma-d$rawmean) ) < err )

  NS <- 20  # number of sigma
  repeat {
    pdf <- rawdensity(d,d$sigma*NS)
    if ( pdf <= 1e-20 || NS >= 400 ) break
    NS <- NS+20
    if ( NS %% 80 == 0 ) print(paste("Pdf test: # sigma=",NS,"pdf=",as.numeric(pdf),"time=",TM() ))
  }
  NT <- floor(40000/NS) # number of tick sample per unit of sigma
  print(paste("Integral: #sigma=",NS,"#ticks=",NT,"pdf=",as.numeric(pdf), "time=",TM() ))

  rs <- generatepdf(d, NS, NT)
  N  <- rs$N
  x  <- rs$x 
  pb <- rs$pb
  print(paste("PDF generated at ",TM() ))

  mu0_num <- sintegral( x, pb, N/2 )$int 
  mu1_num <- sintegral( x, x * pb, N/2 )$int
  mu2_num <- sintegral( x, x^2 * pb, N/2 )$int
  mu3_num <- sintegral( x, x^3 * pb, N/2 )$int
  mu4_num <- sintegral( x, x^4 * pb, N/2 )$int

  # total of probably density should be one
  if ( !( abs(1.0 - mu0_num) < err ) ) {
    print(paste("mu0_num=",Q(mu0_num)))
    stopifnot( abs(1.0 - mu0_num) < err )
  }
  # validate cdf, seem a bit off that I need to use err=1e-4
  L <- length(x)
  for( i in floor(L* c(0.45, 0.5-NT/L/2, 0.5, 0.5+NT/L*2, 0.55, 1.0)) ) { 
    xm <- x[i]
    cdf_raw <- rawcdf(d,xm)
    cdf_num <- sintegral( x[1:i], pb[1:i], i/2 )$int 
    print(paste("cdf[",i,"](",Q(xm),")_raw=",
      Q(cdf_raw),"_num=",Q(cdf_num)
    ))
    stopifnot( abs(cdf_raw - cdf_num) < 1e-4 )

    # pdf slope testing
    slope_num <- (rawdensity(d,x[i])-rawdensity(d,x[i-1]))/(x[i]-x[i-1])
    slope_raw <- rawdensityslope(d,x[i])
    print(paste("slope[",i,"](",Q(xm),")_raw=",
      Q(slope_raw),"_num=",Q(slope_num)
    ))
    stopifnot( abs(slope_raw - slope_num) < 5e-3 ) # very loose

  }
  pb0_num <- rawdensity(d,0)
  pb0_raw <- rawdensity0(d)

  print(paste("pb0_raw=",Q(pb0_raw),Q(pb0_num)))
  stopifnot( abs(pb0_raw - pb0_num)/pb0_raw < err )


  mu2_raw <- rawmu2(d)
  mu3_raw <- rawmu3(d)
  mu4_raw <- rawmu4(d)

  var_raw <- psdvariance(d)
  sk_raw <- psdskewness(d)
  kt_raw <- psdkurtosis(d)

  var_num <- mu2_num - mu1_num^2
  sk_num <- (mu3_num - 3*mu1_num*mu2_num + 2*mu1_num^3 ) / var_num^(3/2)
  kt_num <- (rawmu4(d) - 4*d$rawmean*rawmu3(d) + 6*d$rawmean^2*rawmu2(d) - 3*d$rawmean^4 ) / var_num^2 - 3.0

  print(paste("var_raw=",Q(var_raw),Q(var_num),
    "err=", as.numeric(abs(var_raw-var_num)/var_raw)
  ))
  print(paste("sk_raw=",Q(sk_raw),Q(sk_num),
    "err=", as.numeric(abs(sk_raw-sk_num)/max(1,as.numeric(abs(sk_raw))))
  ))
  print(paste("kt_raw=",Q(kt_raw),Q(kt_num),
    "err=", as.numeric(abs(kt_raw-kt_num)/max(1,as.numeric(abs(kt_raw))))
  ))

  
  mu1_a <- mu1_analytic(d) 
  mu2_a <- mu2_analytic(d)
  mu3_a <- mu3_analytic(d)
  mu4_a <- mu4_analytic(d)

  print(paste(
    "rawmean=",Q(d$rawmean)
    ,"mu1_a=",Q(mu1_a)
    ,"mu1_n=",Q(mu1_num)
  ))
  if ( kt_raw > 10 ) errn <- 1e-4
  stopifnot( abs(d$rawmean - mu1_a) < err )
  stopifnot( abs(d$rawmean - mu1_num) < errn )

  print(paste(
    "mu2_raw=",Q(mu2_raw)
    ,"mu2_a=",Q(mu2_a)
    ,"mu2_n=",Q(mu2_num)
  ))
  if ( kt_raw > 10 ) errn <- 1e-3
  stopifnot( abs(mu2_raw - mu2_a)/mu2_raw < err )
  stopifnot( abs(mu2_raw - mu2_num)/mu2_raw < errn )
  stopifnot( abs(var_raw - var_num)/var_raw < errn )

  print(paste(
    "mu3_raw=",Q(mu3_raw)
    ,"mu3_a=",Q(mu3_a)
    ,"mu3_n=",Q(mu3_num)
  ))
  if ( kt_raw > 10 ) errn <- 1e-2
  dn <- max(abs(mu3_raw),mu2_raw)
  stopifnot( abs(mu3_raw - mu3_a)/dn < err )
  stopifnot( abs(mu3_raw - mu3_num)/dn < errn )
  stopifnot( abs(sk_raw - sk_num) < errn )

  print(paste(
    "mu4_raw=",Q(mu4_raw)
    ,"mu4_a=",Q(mu4_a)
    ,"mu4_n=",Q(mu4_num)
  ))
  if ( kt_raw > 10 ) errn <- 1/20
  stopifnot( abs(mu4_raw - mu4_a)/mu4_raw < err )
  stopifnot( abs(mu4_raw - mu4_num)/mu4_raw < errn )
  stopifnot( abs(kt_raw - kt_num)/max(0.1,as.numeric(kt_raw)) < errn )

  # CDF^-1 test
  CDFCUT <- 0.01 # 1 percent
  xmax <- 0
  repeat {
    if ( rawcdf(d,xmax-1) < 2*CDFCUT ) break
    xmax <- xmax-1
  }
  xmin <- xmax-1
  repeat {
    if ( rawcdf(d,xmin-1) < 0.3*CDFCUT ) break
    xmin <- xmin-1
  }
  xi <- seq(xmin,xmax,by=0.001)
  ci <- rawcdf(d,xi)
  x_num <- xi[max(which(ci<CDFCUT))] # solved by brute force
  x_ntn <- rawcdfinv(d, CDFCUT, d$rawmean) # solved by Newton's method
  print(paste(
    "CDF=",CDFCUT
    ,"x_num=",Q(x_num)
    ,"x_ntn=",Q(x_ntn)
  ))
  stopifnot( abs(x_num - x_ntn)/abs(x_ntn) < 1e-3 )
  
  print(paste("----","----",TM(),"----","----"))

}

LIHNPSD_UnitTest <- function(mpfr=0) {
  tm <- format(Sys.time(), "%X %x")
  print(paste("----","----",tm,"----","----"))
  LihnFunctionValidate()
  
  # In mpfr, fractional number needs to be quoted to be precise
  # simple case
  dist <- list( sigma= 1.0, alpha= 0.0, gamma= 0, beta= "-0.1", mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)

  dist <- list( sigma= 1.0, alpha= 0.5, gamma= "0.2", beta= 0, mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)

  dist <- list( sigma= 1.0, alpha= 0.5, gamma= 0, beta= "-0.1", mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)

  # symmetric case
  dist <- list( sigma= 1.0, alpha= 1.0, gamma= 0, beta= 0, mpfr= mpfr)
  class(dist) <- "LIHNPSD"
  psdunittest(dist)

  dist <- list( sigma= 1.0, alpha= 1.0, gamma= 0, beta= -0.25, mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)

  dist <- list( sigma= 1.0, alpha= 1.0, gamma= "0.3", beta= 0.0, mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)
  
  # Approximately DJI 
  # We need to study why this is not converging? 
  # var err ~ 8e-6, can't get much better no matter NT, NS, epsilon
  dist <- list( sigma= 1.0, alpha= 1.0, gamma= "0.3", beta= "-0.2", mpfr= mpfr )
  class(dist) <- "LIHNPSD"
  psdunittest(dist)
  # Maxima: rawmu1(0.002,2,0.3,-0.1,1,300),numer; -0.00177911390838

  print("--LIHNPSD_UnitTest is completed successfully--")
}
# ------------------------------------------
# Fitting Functions
# ------------------------------------------
LIHNPSD_standardfit_test <- function(d, r, hist, plotqq= 1, 
    weights=list(), merge_tails=c(0,0) ) {
  data_stats <- c( mean(r), sqrt(var(r)), skewness(r), kurtosis(r)-3.0 )
  psd <- c( d$location/data_stats[2], d$sigma/data_stats[2], d$alpha, d$gamma, d$beta )
  LIHNPSD_standardfit_fn(psd, data_stats, hist, weights=weights, plotqq= plotqq, merge_tails= merge_tails)
}

LIHNPSD_standardfit_fn <- function( psd, data_stats, hist, plotqq=1, 
    weights=list(), merge_tails=c(0,0), debug=0 ) {

  tm0 <- unclass(Sys.time())

  d <- list( location= psd[1]*data_stats[2],
    sigma= psd[2]*data_stats[2], alpha= psd[3], gamma= psd[4], beta= psd[5] )
  class(d) <- "LIHNPSD"
  d <- rawmean(d)

  if ( abs(d$beta) >= sqrt(2/pi) ) {
    delta <- if ( debug == 1 ) list( delta= abs(d$beta)*1000 ) else abs(d$beta)*1000
    return(delta)
  } 

  if (is.null(weights$m2)) weights$m2 <- 1
  if (is.null(weights$m3)) weights$m3 <- 1
  if (is.null(weights$m4)) weights$m4 <- 1
  if (is.null(weights$mmx)) weights$mmx <- 1
  if (is.null(weights$pdf_df)) weights$pdf_df <- 1
  if (is.null(weights$qq_df))  weights$qq_df  <- 0.1
  if (is.null(weights$a_bump)) weights$a_bump <- 1
   
  st <- c( 0, sqrt(psdvariance(d)), psdskewness(d), psdkurtosis(d) )

  tm1 <- unclass(Sys.time()) - tm0

  NS <- 4
  xs <- seq( d$rawmean-NS*d$sigma, d$rawmean+NS*d$sigma, by=2*NS*d$sigma/200 )
  ds <-rawdensity(d,xs)
  max_density <- max(hist$density)
  mmx <- abs( max(ds)-max_density ) / max_density

  tm2 <- unclass(Sys.time()) - tm0

  NS <- 4
  pdf_df <- 0.0
  xpdf <- c() 
  tpdf <- c()
  j <- 1
  for ( i in 1:length(hist$mids) ) {
    x <- hist$mids[i]
    if ( abs(x-d$rawmean-d$location) <= NS*d$sigma ) {
      xpdf[j] <- x
      tpdf[j] <- rawdensity(d, x-d$rawmean-d$location)
      df <- (log(hist$density[i]) - log(tpdf[j]))^2
      pdf_df <- pdf_df + if ( abs(df) <= 1000.0 ) df else 0
      j <- j + 1
    }
  }
  pdf_df <- sqrt(abs(pdf_df))
  pdf_df <- if ( pdf_df == Inf ) 10000.0 else pdf_df

  tm3 <- unclass(Sys.time()) - tm0

  hq <- MergeTailHistogram( list(qhx=hist$mids, qhy=hist$counts), merge_tails)
  qqp <- calcqq(d, hq, debug=debug)
  qq_df <- sqrt(sum(abs(qqp$y-qqp$x)^2))/max(hist$mids)

  tm4 <- unclass(Sys.time()) - tm0

  alpha_bump <- if ( d$alpha > 1.5 ) abs(d$alpha-1.5)/10 else 0
  m2 <- abs( st[2]-data_stats[2] ) / data_stats[2] # var pct
  m3 <- abs( st[3]-data_stats[3] )
  m4 <- abs( st[4]-data_stats[4] ) / data_stats[4] # kurt pct

  delta <- 0
  delta <- delta + (m2 * weights$m2)^2
  delta <- delta + (m3 * weights$m3)^2
  delta <- delta + (m4 * weights$m4)^2
  delta <- delta + (alpha_bump * weights$a_bump)^2
  delta <- delta + (mmx * weights$mmx)^2
  delta <- delta + (qq_df * weights$qq_df)^2
  delta <- delta + (pdf_df * weights$pdf_df)^2
  delta <- sqrt(delta)


  if ( plotqq == 1 ) {
    par(mfcol=c(2,2)); # 1 plot
    par(oma = c(1, 1, 4, 1))

    par(mar=c(4,4,2,2)) 
    xmin <- min(xs+d$rawmean+d$location)
    xmax <- max(xs+d$rawmean+d$location)
    plot(hist$mids, hist$density, type="o", pch=22, lty=2, col="red", 
      main="PDF", xlab="r", ylab="PDF",
      xlim=c(xmin,xmax) )
    lines(xpdf, tpdf, pch = ".", col="black")
    lines(xs+d$rawmean+d$location, ds, pch = ".", col="green")

    par(mar=c(4,4,2,2)) 
    plot(hist$mids, log(hist$density), type="o", pch=22, lty=2, col="red", 
      main="Log PDF", xlab="log(r)", ylab="Log PDF")
    lines(xpdf, log(tpdf), pch = ".", col="black")

    xmin <- min(hist$mids)
    xpos <- xmin*0.5
    ypos <- max(log(hist$density))
    sc <- 0.8
    text(xmin*0.5,ypos,labels=c("data"),cex=sc)
    text(xmin*0.5,ypos-1,labels=paste("mean", sprintf("%.6f",d$rawmean)),cex=sc)
    text(xmin*0.5,ypos-2,labels=paste("std",  sprintf("%.6f",sqrt(psdvariance(d)))),cex=sc)
    text(xmin*0.5,ypos-3,labels=paste("skew", sprintf("%.6f",psdskewness(d))),cex=sc)
    text(xmin*0.5,ypos-4,labels=paste("kurt", sprintf("%.6f",psdkurtosis(d))),cex=sc)


    par(mar=c(4,4,2,2))
    plot(qqp$x, qqp$y, pch = "o", col="red", 
      main=paste("QQ Plot",Sys.time()), 
      ylim=c( min(qqp$y), max(qqp$y)*1.5 ),
      xlab="Observed Quantile", ylab="Theoretical Quantile")
    abline(h=0,col="yellow")
    abline(v=0,col="yellow")
    lines(qqp$x, abs(qqp$y-qqp$x), pch = ".", col="green")
    abline(0,1,col="blue")

    tm9 <- unclass(Sys.time()) - tm0

    xmin <- min(qqp$x)
    xpos <- xmin*0.5
    ypos <- max(qqp$y)*1.5
    ydf <- (max(qqp$y)-min(qqp$y))/10
    sc <- 0.8

    text(xpos,ypos-1*ydf,labels=c("timing"),cex=sc)
    text(xpos,ypos-2*ydf,labels=paste("tm1",   sprintf("%.2f",tm1)),cex=sc)
    text(xpos,ypos-3*ydf,labels=paste("tm2",   sprintf("%.2f",tm2)),cex=sc)
    text(xpos,ypos-4*ydf,labels=paste("tm3",   sprintf("%.2f",tm3)),cex=sc)
    text(xpos,ypos-5*ydf,labels=paste("tm4",   sprintf("%.2f",tm4)),cex=sc)
    text(xpos,ypos-6*ydf,labels=paste("tm9",   sprintf("%.2f",tm9)),cex=sc)

    par(mar=c(4,4,2,2)) 
    plot(qqp$x, qqp$y, pch = "o", col="green", 
      main=paste("QQ Data",Sys.time()), 
      ylim=c( min(qqp$y), max(qqp$y)*1.5 ),
      xlab="Observed Quantile", ylab="Theoretical Quantile")



    text(xpos,ypos-1*ydf,labels=c("fit"),cex=sc)
    text(xpos,ypos-2*ydf,labels=paste("delta",   sprintf("%.8f",delta)),cex=sc)
    text(xpos,ypos-4*ydf,labels=paste("alpha", sprintf("%.8f",d$alpha)),cex=sc)
    text(xpos,ypos-5*ydf,labels=paste("gamma", sprintf("%.8f",d$gamma)),cex=sc)
    text(xpos,ypos-6*ydf,labels=paste("beta",  sprintf("%.8f",d$beta)),cex=sc)
    text(xpos,ypos-7*ydf,labels=paste("sigma",  sprintf("%.6f %.8f",d$sigma,psd[2])),cex=sc)
    text(xpos,ypos-8*ydf,labels=paste("loc",  sprintf("%.6f %.8f",d$location,psd[1])),cex=sc)
    xmax <- max(qqp$x)
    xpos <- xmax*0.4
    text(xpos,ypos-1*ydf,labels=c("delta parts"),cex=sc)
    text(xpos,ypos-2*ydf,labels=paste("max_pdf", sprintf("%.6f",max_density)),cex=sc)
    text(xpos,ypos-3*ydf,labels=paste("mmx", sprintf("%.6f %.2f",mmx, weights$mmx)),cex=sc)
    text(xpos,ypos-4*ydf,labels=paste("m2",  sprintf("%.6f %.2f",m2, weights$m2)),cex=sc)
    text(xpos,ypos-5*ydf,labels=paste("m3",  sprintf("%.6f %.2f",m3, weights$m3)),cex=sc)
    text(xpos,ypos-6*ydf,labels=paste("m4",  sprintf("%.6f %.2f",m4, weights$m4)),cex=sc)
    text(xpos,ypos-7*ydf,labels=paste("a_bump",  sprintf("%.6f %.2f",alpha_bump, weights$a_bump)),cex=sc)
    text(xpos,ypos-8*ydf,labels=paste("pdf_df",  sprintf("%.6f %.2f",pdf_df, weights$pdf_df)),cex=sc)
    text(xpos,ypos-9*ydf,labels=paste("qq_df", sprintf("%.8f %.2f",qq_df, weights$qq_df)),cex=sc)

  }
  if ( debug == 1 ) list( delta= delta,
    mmx= mmx, m2=m2, m3=m3, m4=m4, alpha_bump=alpha_bump,
    max_density=c(max(ds),max_density), qq_df= qq_df ) 
  else delta
}

# optimx: http://cran.r-project.org/web/packages/optimx/optimx.pdf
# spg: http://cran.r-project.org/web/packages/BB/BB.pdf

standardfit.LIHNPSD <- function(d, r, hist, trace= 0, iter=5000, 
  plotqq= 1, weights=list(), merge_tails=c(0,0) ) {

  print(paste("standardfit:","iter=",iter,"plotqq=",plotqq,
    "weights=(",paste(weights),")"))

  data_stats <- c( mean(r), sqrt(var(r)), skewness(r), kurtosis(r)-3.0 )
  psd <- c( d$location/data_stats[2], d$sigma/data_stats[2], d$alpha, d$gamma, d$beta )

  psdout <- optimx( psd, LIHNPSD_standardfit_fn, method=c("spg"), 
    itnmax=iter, control=list(trace=trace,eps=1e-5), 
    data_stats= data_stats, hist= hist, plotqq= plotqq, 
    weights=weights, merge_tails=merge_tails )
  par <- psdout$par$par
  d2 <- list( location=par[1]*data_stats[2],
    sigma=par[2]*sqrt(var(r)), alpha=par[3], gamma=par[4], beta=par[5] )
  class(d2) <- "LIHNPSD"
  d2 <- rawmean(d2)
  list( dist= d2, psdout= psdout)
}

# ------------------------------------------
# Misc Helper Functions
# ------------------------------------------

# converts the time series to log returns
# to calculates the log-return distribution

TimeSeriesLogReturn <- function( pr, days ) {
  n <- length(pr)
  s <- seq(1, n, by=days)
  pr2 <- pr[s]
  r <- c(0.0)
  for (i in 2:length(pr2)) {
    if ( pr[i-1] > 0.0 && pr[i] ) { 
      r[i-1] <- log(pr[i]/pr[i-1])
    }
  }
  r
}
ConvertDistToSample <- function( x, hx, N, ngc = 200 ) {
  NR <- N / sum(hx)
  hx_r <- c()
  for( i in 1:length(hx) ) {
    xi <- x[i]
    xs <- seq(xi,xi,length=hx[i]*NR)
    hx_r <- c(hx_r, xs)
    if ( i %% ngc ) gc()
  }
  gc()
  hx_r
}

# q = list( qhx = h$mids, qhy = h$counts )
MergeTailHistogramOneSide <- function( q, allowed_merge ) {
  qhx <- q$qhx
  qhy <- q$qhy

  owed <- 0
  while ( allowed_merge > 0 ) {
    if ( qhy[1] > 0 ) {
      qhy <- qhy[-1]
      qhx <- qhx[-1]
      allowed_merge <- allowed_merge - 1
      owed <- owed + 1
      while ( qhy[1] == 0 ) {
        qhy <- qhy[-1]
        qhx <- qhx[-1]
      }
    }
    else break
  }
  qhy[1] <- qhy[1] + owed
  list( qhx= qhx, qhy= qhy )
}

MergeTailHistogram <- function ( q, merge_tails ) {
  allowed_negative_merge <- merge_tails[1]
  allowed_positive_merge <- merge_tails[2]
  q2 <- MergeTailHistogramOneSide ( q, allowed_negative_merge )
  q2 <- list( qhx= rev(q2$qhx), qhy= rev(q2$qhy) )
  q3 <- MergeTailHistogramOneSide ( q2, allowed_positive_merge )
  list( qhx= rev(q3$qhx), qhy= rev(q3$qhy) )
}

# ------------------------------------------------------------------
# Standard plot facility
# ------------------------------------------------------------------
# h: histogram
# tx: theoretical x, tpdf/tcdf: theoretical PSD
LIHNPSD_plotpdf <- function (dist, h, tx, tpdf, xlab="log(r)", main="PSD PDF" ) {
  par(mar=c(4,4,2,1)) 
  s <- sqrt(psdvariance(dist))
  i <- 3
  repeat {
     if ( rawcdf(dist, dist$rawmean+i*s) >= 0.995 || i > 20 ) break
     i <- i + 1
  }
  xmin <- -s*i+dist$rawmean
  ymax <- max(max(tpdf),max(h$density))

  plot( h$mids, h$density, pch = "o", col="red",
    xlim=c(-s*i, s*i)+dist$rawmean,
    ylim=c(0,ymax),
    ylab="PDF", xlab=xlab, main=main)
  lines(tx, tpdf, pch = ".", col="blue")

  sc=0.8
  yseg <- 0.20*max(h$density)
  segments( dist$location+dist$rawmean-s, yseg, 
            dist$location+dist$rawmean+s, yseg, lwd=3 ) 
  text( dist$location, yseg*0.5, "std dev", cex=sc )
  # legends
  legend(xmin, ymax, c("data","fit"), cex=sc, 
   col=c("red","blue"), pch=c("o",NA), lty=c(NA,1));
}
LIHNPSD_plotlogpdf <- function (dist, h, tx, tpdf, xlab="log(r)", main="PSD Log PDF" ) {
  par(mar=c(4,4,2,1)) 

  fnt <- log(h$density) [ which(is.finite(log(h$density))) ]

  xmin <- min(h$mids)
  xmax <- max(h$mids)
  plot( h$mids, log(h$density), pch = "o", col="red", 
    ylim=c( min(max(fnt)-8,min(fnt)), max(fnt)+0.5 ),
    xlim=c( xmin, xmax ),
    ylab="log(PDF)", xlab=xlab, main=main)
  lines(tx, log(tpdf), pch = ".", col="blue")

  # legends

  # parameters
  xpos <- 0.7 * if ( abs(xmin) > abs(xmax) ) xmin else xmax 
  ypos <- max(log(h$density))+0.5
  sc <- 0.8
  sp <- 0.8
  text(xpos,ypos-1*sp,labels=c("psd fit"),cex=sc)
  text(xpos,ypos-2*sp,labels=paste("location", sprintf("%.6f",dist$location)),cex=sc)
  text(xpos,ypos-3*sp,labels=paste("sigma",    sprintf("%.6f",dist$sigma)),cex=sc)
  text(xpos,ypos-4*sp,labels=paste("alpha", sprintf("%.4f",dist$alpha)),cex=sc)
  text(xpos,ypos-5*sp,labels=paste("gamma", sprintf("%.4f",dist$gamma)),cex=sc)
  text(xpos,ypos-6*sp,labels=paste("beta",  sprintf("%.4f",dist$beta)),cex=sc)
}

LIHNPSD_plotcdf <- function (dist, h, data_st, tx, tcdf, xlab="log(r)", main="PSD CDF" ) {
  par(mar=c(4,4,2,1)) 
  s <- sqrt(psdvariance(dist))
  i <- 3
  repeat {
     if ( rawcdf(dist, dist$rawmean+i*s) >= 0.995 || i > 20 ) break
     i <- i + 1
  }

  plot( h$mids, cumsum(h$counts) / sum(h$counts), lwd=2, type="s", col="red", 
    xlim=c(-s*i, s*i)+dist$rawmean,
    ylab="CDF", xlab=xlab, main=main)
  dx <- tx[2]-tx[1]
  # lines( h$mids, cumsum(h$counts) / sum(h$counts), col="red") 
  lines(tx, tcdf, pch = ".", col="blue")
  xmin <- -s*i+dist$rawmean
  xmax <- s*i+dist$rawmean
  sc <- 0.8
  text(xmin*0.5,0.90,labels=c("data"),cex=sc)
  text(xmin*0.5,0.75,labels=paste("mean", sprintf("%.6f",data_st[1])),cex=sc)
  text(xmin*0.5,0.60,labels=paste("std",  sprintf("%.4f",data_st[2])),cex=sc)
  text(xmin*0.5,0.45,labels=paste("skew", sprintf("%.4f",data_st[3])),cex=sc)
  text(xmin*0.5,0.30,labels=paste("kurt", sprintf("%.4f",data_st[4])),cex=sc)

  text(xmax*0.5,0.75,labels=c("psd fit"),cex=sc)
  text(xmax*0.5,0.60,labels=paste("mean", sprintf("%.6f",dist$rawmean+dist$location)),cex=sc)
  text(xmax*0.5,0.45,labels=paste("std",  sprintf("%.4f",sqrt(psdvariance(dist)))),cex=sc)
  text(xmax*0.5,0.30,labels=paste("skew", sprintf("%.4f",psdskewness(dist))),cex=sc)
  text(xmax*0.5,0.15,labels=paste("kurt", sprintf("%.4f",psdkurtosis(dist))),cex=sc)

}
LIHNPSD_plotqq <- function (dist, qqp, merge_tails, main="PSD QQ-Plot" ) {
  par(mar=c(4,4,2,1))
 
  dx <- qqp$x[2] - qqp$x[1]
  plot(qqp$x+dx/2, qqp$y, pch="o", col="red", main=main, 
    xlab="Observed Quantile", ylab="Theoretical Quantile")
  abline(0,1,col="blue")
  abline(h=0,col="yellow")
  abline(v=0,col="yellow")
  lines(qqp$x, abs(qqp$y-qqp$x), pch = ".", col="green")

  # legends
  sc <- 0.8
  xmin <- min(qqp$x)
  ymax <- max(qqp$y)
  legend(xmin, ymax, c("qq data","45 degree","error"), cex=sc, 
   col=c("red","black","green"), pch=c("o",NA,NA), lty=c(NA,1,1));

  xmax <- max(qqp$x)
  ymin <- min(qqp$y)
  text(xmax*0.5, ymin*0.4, "tail dropped:",cex=sc);
  text(xmax*0.5, ymin*0.6, paste("left",  merge_tails[1]),cex=sc);
  text(xmax*0.5, ymin*0.8, paste("right", merge_tails[2]),cex=sc);
}
LIHNPSD_plot_std4gr <- function ( th, dt, EPS=FALSE, file=NA ) {
  if ( EPS ) {
    postscript(file=file,
      paper="special", width=8, height=8, horizontal=FALSE)
  }

  par(mfcol=c(2,2)); # 4 plots
  par(oma = c(1, 1, 4, 1))

  LIHNPSD_plotpdf    (th$dist, dt$h, th$tx, th$tpdf)
  LIHNPSD_plotlogpdf (th$dist, dt$h, th$tx, th$tpdf)
  LIHNPSD_plotcdf    (th$dist, dt$h, dt$stats, th$tx, th$tcdf)
  LIHNPSD_plotqq     (th$dist, th$qqp, th$merge_tails )

  if ( EPS ) dev.off()

}

LIHNPSD_prepare_data <- function( logr, breaks, merge_tails ) {
  dt <- list( logr= logr, N= length(logr), breaks=breaks, merge_tails=merge_tails )
  dt$stats <- c( mean(logr), sqrt(var(logr)), skewness(logr), kurtosis(logr)-3.0 )
  h <- hist(logr, breaks = breaks, plot = FALSE)
  h$cuml <- cumsum(h$counts) / sum(h$counts)
  dt$hq <- MergeTailHistogram( list(qhx=h$mids, qhy=h$counts), merge_tails )
  dt$h  <- h

  dt
}

LIHNPSD_theoretical_result <- function ( dist, dt, N=5000 ) {
  th <- list( dist=dist, N= N, qqp= calcqq(dist, dt$hq), merge_tails=dt$merge_tails )
  tx  <- seq(min(dt$logr), max(dt$logr), length=N+1)
  th$dx   <- tx[2]-tx[1]
  th$tpdf <- density(dist,tx)
  th$tcdf <- rawcdf(dist, tx-dist$rawmean-dist$location)
  th$tx <- tx

  #k  <- seq(0,10,by=1)
  #kn <- psdkernel(dist, k)
  #kp <- PoissonTerm( dist$lambda, k)

  th
}

