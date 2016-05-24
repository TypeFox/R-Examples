# Checks for lower <= Dx <= upper, where
# mean (d x 1), sigma (d x d), D (r x d), x (d x 1), lower (r x 1), upper (r x 1)
# Uses partly checks as in mvtnorm:::checkmvArgs!
# 
checkTmvArgs2 <- function(mean, sigma, lower, upper, D)
{
	if (is.null(lower) || any(is.na(lower))) 
		stop(sQuote("lower"), " not specified or contains NA")
	if (is.null(upper) || any(is.na(upper))) 
		stop(sQuote("upper"), " not specified or contains NA")
	if (!is.numeric(mean) || !is.vector(mean)) 
		stop(sQuote("mean"), " is not a numeric vector")
	if (is.null(sigma) || any(is.na(sigma))) 
		stop(sQuote("sigma"), " not specified or contains NA")
	if (is.null(D) || any(is.na(D))) 
		stop(sQuote("D"), " not specified or contains NA")
	
	if (!is.matrix(sigma)) {
		sigma <- as.matrix(sigma)
	}
	
	if (!is.matrix(D)) {
		D <- as.matrix(D)
	}
	
	if (NCOL(lower) != NCOL(upper)) {
		stop("lower and upper have non-conforming size")
	}
	
	checkSymmetricPositiveDefinite(sigma)
	
	d <- length(mean)
	r <- length(lower)
	
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
	}
	
	if (length(lower) != NROW(D) || length(upper) != NROW(D)) {
		stop("D (r x d), lower (r x 1) and upper (r x 1) have non-conforming size")
	}
	
	if (length(mean) != NCOL(D)) {
		stop("D (r x d) and mean (d x 1) have non-conforming size")
	}
	
	if (any(lower>=upper)) {
		stop("lower must be smaller than or equal to upper (lower<=upper)")
	}
	
	# checked arguments
	cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper, D=D)
	return(cargs)
}
# Gibbs sampling with general linear constraints a <= Dx <= b
# with x (d x 1), D (r x d), a,b (r x 1) requested by Xiaojin Xu [xiaojinxu.fdu@gmail.com]
# which allows for (r > d) constraints!

# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (d x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (d x d) der Normalverteilung
# @param lower unterer Trunkierungsvektor (d x 1) mit lower <= Dx <= upper
# @param upper oberer Trunkierungsvektor (d x 1) mit lower <= Dx <= upper
# @param D Matrix for linear constraints, defaults to (d x d) diagonal matrix
# @param H Precision matrix (d x d) if given
# @param algorithm c("rejection", "gibbs", "gibbsR")
rtmvnorm2 <- function(n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)),
    D = diag(length(mean)),
    algorithm=c("gibbs", "gibbsR", "rejection"), ...) {
  algorithm <- match.arg(algorithm)
  
  # check of standard tmvtnorm arguments
  # Have to change check procedure to handle r > d case
  cargs <- checkTmvArgs2(mean, sigma, lower, upper, D)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  D <- cargs$D
  
  # check of additional arguments
  if (n < 1 || !is.numeric(n) || n != as.integer(n) || length(n) > 1) {
	  stop("n must be a integer scalar > 0")
  }
  
  if (!identical(D,diag(length(mean)))) {
    # D <> I : general linear constraints
    if (algorithm == "gibbs") {
      # precision matrix case H vs. covariance matrix case sigma will be handled inside method 
      retval <- rtmvnorm.gibbs2.Fortran(n, mean=mean, sigma=sigma, D=D, lower=lower, upper=upper, ...)
    } else if (algorithm == "gibbsR") {
      # covariance matrix case sigma
      retval <- rtmvnorm.gibbs2(n, mean=mean, sigma=sigma, D=D, lower=lower, upper=upper, ...)  
    } else if (algorithm == "rejection") {
      retval <- rtmvnorm.rejection(n, mean=mean, sigma=sigma, D=D, lower=lower, upper=upper, ...)
    }
    return(retval)
  } else {
    # for D = I (d x d) forward to normal rtmvnorm() method
    retval <- rtmvnorm(n, mean=mean, sigma=sigma, lower=lower, upper=upper, D=D, ...)
    return(retval)
  }
  return(retval)
}


# Gibbs sampler implementation in R for general linear constraints
# lower <= Dx <= upper where D (r x d), x (d x 1), lower, upper (r x 1)
# which can handle the case r > d.
#
# @param n
# @param mean
# @param sigma
# @param D
# @param lower
# @param upper
# @param burn.in.samples
# @param start.value
# @param thinning
rtmvnorm.gibbs2 <- function (n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)), 
    D = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), 
    upper = rep(Inf, length = length(mean)), 
    burn.in.samples = 0, 
    start.value = NULL, 
    thinning = 1) 
{
    if (thinning < 1 || !is.numeric(thinning) || length(thinning) > 1) {
        stop("thinning must be a integer scalar > 0")
    }
    d <- length(mean)
    S <- burn.in.samples
    if (!is.null(S)) {
        if (S < 0) 
            stop("number of burn-in samples must be non-negative")
    }
    if (!is.null(start.value)) {
        if (length(mean) != length(start.value)) 
            stop("mean and start value have non-conforming size")
        if (any(D %*% start.value < lower || D %*% start.value > upper)) 
            stop("start value does not suffice linear constraints lower <= Dx <= upper")
        x0 <- start.value
    }
    else {
        x0 <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
    }
    if (d == 1) {
      X <- rtnorm.gibbs(n, mu = mean[1], sigma = sigma[1, 1], 
            a = lower[1], b = upper[1])
      return(X)
    }
    # number of linear constraints lower/a <= Dx <= upper/b, D (r x n), a,b (r x 1), x (n x 1)
    r <- nrow(D)
    
    X <- matrix(NA, n, d)
    U <- runif((S + n * thinning) * d)
    l <- 1
    sd <- list(d)
    P <- list(d)
    
    # [ Sigma_11  Sigma_12 ] = [ sigma_{i,i}   sigma_{i,-i}  ]
    # [ Sigma_21  Sigma_22 ]   [ sigma_{-i,i}  sigma_{-i,-i} ]
    for (i in 1:d) {
        Sigma_11 <- sigma[i, i]     # (1 x 1)
        Sigma_12 <- sigma[i, -i]    # (1 x (d - 1))
        Sigma_22 <- sigma[-i, -i]   # ((d - 1) x (d - 1))
        P[[i]] <- t(Sigma_12) %*% solve(Sigma_22)
        sd[[i]] <- sqrt(Sigma_11 - P[[i]] %*% Sigma_12)
    }
    x <- x0
    # for all draws
    for (j in (1 - S):(n * thinning)) {
      # for all x[i]
      for (i in 1:d) {
        lower_i <- -Inf
        upper_i <- +Inf
            
        # for all linear constraints k relevant for variable x[i]. 
        # If D[k,i]=0 then constraint is irrelevant for x[i]
        for (k in 1:r) {
          if (D[k,i] == 0) next
          bound1 <- lower[k]/D[k, i] - D[k,-i] %*% x[-i] /D[k, i]
          bound2 <- upper[k]/D[k, i] - D[k,-i] %*% x[-i] /D[k, i]

          if (D[k, i] > 0) {
            lower_i <- pmax(lower_i, bound1)
            upper_i <- pmin(upper_i, bound2)
          } else {
            lower_i <- pmax(lower_i, bound2)
            upper_i <- pmin(upper_i, bound1)
          }
        }
            
        mu_i <- mean[i] + P[[i]] %*% (x[-i] - mean[-i])
        F.tmp <- pnorm(c(lower_i, upper_i), mu_i, sd[[i]])
        Fa <- F.tmp[1]
        Fb <- F.tmp[2]
        x[i] <- mu_i + sd[[i]] * qnorm(U[l] * (Fb - Fa) + Fa)
        l <- l + 1
      }
      if (j > 0) {
        if (thinning == 1) {
          X[j, ] <- x
        }
        else if (j%%thinning == 0) {
          X[j%/%thinning, ] <- x
        }
      }
    }
    return(X)
}

rtmvnorm.gibbs2.Fortran <- function(n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)),
    D     = diag(length(mean)), 
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)), 
		burn.in.samples = 0, 
    start.value = NULL, 
    thinning = 1)
{
  # No checks of input arguments, checks are done in rtmvnorm()
  
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	  if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
    if (NCOL(D) != length(start.value) || NROW(D) != length(lower) || NROW(D) != length(upper))  stop("D, start.value, lower, upper have non-conforming size")
	  if (any(D %*% start.value < lower || D %*% start.value > upper)) stop("start value must lie in simplex defined by lower <= Dx <= upper") 
	  x0 <- start.value 
  } else {
    stop("Must give start.value with lower <= D start.value <= upper")
  }
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1) {
    X <- rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    return(X)
  }
      
  # Ergebnismatrix (n x d)
  X <- matrix(0, n, d)
  
  # number of linear constraints lower/a <= Dx <= upper/b, D (r x n), a,b (r x 1), x (n x 1)
  r <- nrow(D)
  
  # Call to Fortran subroutine
  # TODO: Aufpassen, ob Matrix D zeilen- oder spaltenweise an Fortran übergeben wird!
  # Bei sigma ist das wegen Symmetrie egal.
  ret <- .Fortran("rtmvnormgibbscov2",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              r     = as.integer(r),
                              mean  = as.double(mean),
                              sigma = as.double(sigma),
                              C     = as.double(D),
                              a     = as.double(lower), 
                              b     = as.double(upper),
                              x0    = as.double(x0),
							                burnin   = as.integer(burn.in.samples),
							                thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
  X <- matrix(ret$X, ncol=d, byrow=TRUE)
  return(X)
}

if (FALSE) {
# dimension d=2
# number of linear constraints r=3 > d 
# linear restrictions a <= Dx <= b with x (d x 1); D (r x d); a,b (r x 1)
D <- matrix(
      c(  1,  1,
          1, -1,
        0.5, -1), 3, 2, byrow=TRUE)
a <- c(0, 0, 0)
b <- c(1, 1, 1)

# mark linear constraints as lines
plot(NA, xlim=c(-0.5, 1.5), ylim=c(-1,1))
for (i in 1:3) {
  abline(a=a[i]/D[i, 2], b=-D[i,1]/D[i, 2], col="red")
  abline(a=b[i]/D[i, 2], b=-D[i,1]/D[i, 2], col="red")
}

# Gibbs sampling:
# determine lower and upper bounds for each index i given the remaining variables: x[i] | x[-i]

### Gibbs sampling for general linear constraints a <= Dx <= b
x0 <- c(0.5, 0.2)
sigma <- matrix(c(1, 0.2, 
                  0.2, 1), 2, 2)
X <- rtmvnorm.gibbs2(n=1000, mean=c(0, 0), sigma, D, lower=a, upper=b, start.value=x0)
points(X, pch=20, col="black")

X2 <- rtmvnorm.gibbs2(n=1000, mean=c(0, 0), sigma, D, lower=a, upper=b, start.value=x0)
points(X2, pch=20, col="green")

# Rejection sampling (rtmvnorm.rejection) funktioniert bereits mit beliebigen Restriktionen (r > d)
X3 <- rtmvnorm.rejection(n=1000, mean=c(0, 0), sigma, D, lower=a, upper=b)
points(X3, pch=20, col="red")


rtmvnorm.gibbs2(n=1000, mean=c(0, 0), sigma, D, lower=a, upper=b, start.value=c(-1, -1))


colMeans(X)
colMeans(X2)
}

