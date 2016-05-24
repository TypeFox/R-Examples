# estimation methods for the parameters of the truncated multivariate normal distribution
#
# Literatur:
#
# Amemiya (1974)   : Instrumental Variables estimator
# Lee (1979)
# Lee (1983)
# Griffiths (2002) : 
# "Gibbs Sampler for the parameters of the truncated multivariate normal distribution"
#
# Stefan Wilhelm, wilhelm@financial.com
#library(tmvtnorm)
library(stats4)

# Hilfsfunktion : VECH() Operator
vech=function (x)
{
  # PURPOSE: creates a column vector by stacking columns of x
  #          on and below the diagonal
  #----------------------------------------------------------
  # USAGE:  v = vech(x)
  # where:  x = an input matrix
  #---------------------------------------------------------
  # RETURNS:
  #         v = output vector containing stacked columns of x
  #----------------------------------------------------------

  # Written by Mike Cliff, UNC Finance  mcliff@unc.edu
  # CREATED: 12/08/98
  
  #if(!is.matrix(x))
  #{
  #   
  #}

  rows    = nrow(x)
  columns = ncol(x);
  v = c();
  for (i in 1:columns)
  {
   v = c(v, x[i:rows,i]);
  }
 v
}

# Hilfsfunktion : Operator für Namensgebung sigma_i.j (i <= j), d.h. wie vech(), nur Zeilenweise
vech2 <- function (x)
{
  # PURPOSE: creates a column vector by stacking columns of x
  #          on and below the diagonal
  #----------------------------------------------------------
  # USAGE:  v = vech2(x)
  # where:  x = an input matrix
  #---------------------------------------------------------
  # RETURNS:
  #         v = output vector containing stacked columns of x
  #----------------------------------------------------------

  # Written by Mike Cliff, UNC Finance  mcliff@unc.edu
  # CREATED: 12/08/98
  
  rows    = nrow(x)
  columns = ncol(x);
  v = c();
  for (i in 1:rows)
  {
   v = c(v, x[i,i:columns]);
  }
 v
}

# Hilfsfunktion : Inverser VECH() Operator
inv_vech=function(v)
{
  #----------------------------------------------------------
  # USAGE:  x = inv_vech(v)
  # where:  v = a vector
  #---------------------------------------------------------
  # RETURNS:
  #         x = a symmetric (m x m) matrix containing de-vectorized elements of v
  #----------------------------------------------------------
  
  # Anzahl der Zeilen
  m = -0.5+sqrt(0.5^2+2*length(v))
  x = matrix(0,nrow=m,ncol=m)
  if (length(v) != m*(m+1)/2)
  {
    # error
    stop("v must have m*(m+1)/2 elements")
  }
  
  for (i in 1:m)
  {
    #cat("r=",i:m," c=",i,"\n")
    x[ i:m, i]   = v[((i-1)*(m-(i-2)*0.5)+1) : (i*(m-(i-1)*0.5))]
    x[ i,   i:m] = v[((i-1)*(m-(i-2)*0.5)+1) : (i*(m-(i-1)*0.5))]
  }
  x  
}

# 1. Maximum-Likelihood-Estimation of mu and sigma when truncation points are known
#
# TODO/Idee: Cholesky-Zerlegung der Kovarianzmatrix als Parametrisierung
#  
# @param X data matrix (T x n)
# @param lower, upper truncation points
# @param start list of start values for mu and sigma
# @param fixed a list of fixed parameters
# @param method
# @param cholesky flag, if TRUE, we use the Cholesky decomposition of sigma as parametrization
# @param lower.bounds lower bounds for method "L-BFGS-B"
# @param upper.bounds upper bounds for method "L-BFGS-B"
mle.tmvnorm <- function(X, 
 lower=rep(-Inf, length = ncol(X)), 
 upper=rep(+Inf, length = ncol(X)), 
 start=list(mu=rep(0,ncol(X)), sigma=diag(ncol(X))), 
 fixed=list(), method="BFGS", 
 cholesky=FALSE, 
 lower.bounds=-Inf,
 upper.bounds=+Inf,
 ...)
{
  # check of standard tmvtnorm arguments
  cargs <- checkTmvArgs(start$mu, start$sigma, lower, upper)
  start$mu    <- cargs$mean
  start$sigma <- cargs$sigma
  lower       <- cargs$lower
  upper       <- cargs$upper
  
  # check if we have at least one sample
  if (!is.matrix(X) || nrow(X) == 0) {
    stop("Data matrix X with at least one row required.")
  }
  
  # verify dimensions of x and lower/upper match
  n <- length(lower)
  if (NCOL(X) != n) {
		stop("data matrix X has a non-conforming size. Must have ",length(lower)," columns.")
	}
	
	# check if lower <= X <= upper for all rows
	ind <- logical(nrow(X))
  for (i in 1:nrow(X))
  {
    ind[i] = all(X[i,] >= lower & X[i,] <= upper)
  }
  if (!all(ind)) {
    stop("some of the data points are not in the region lower <= X <= upper")
  }
  
  if ((length(lower.bounds) > 1L || length(upper.bounds) > 1L || lower.bounds[1L] != 
        -Inf || upper.bounds[1L] != Inf) && method != "L-BFGS-B") {
        warning("bounds can only be used with method L-BFGS-B")
        method <- "L-BFGS-B"
  }
  
  # parameter vector theta = mu_1,...,mu_n,vech(sigma)
  if (cholesky) {
    # if cholesky == TRUE use Cholesky decomposition of sigma
    # t(chol(sigma)) returns a lower triangular matrix which can be vectorized using vech() 
    theta <- c(start$mu, vech2(t(chol(start$sigma))))
  } else {
    theta <- c(start$mu, vech2(start$sigma))
  }
  # names for mean vector elements : mu_i
  nmmu     <- paste("mu_",1:n,sep="")
  # names for sigma elements : sigma_ij
  nmsigma  <- paste("sigma_",vech2(outer(1:n,1:n, paste, sep=".")),sep="")
  names(theta) <- c(nmmu, nmsigma)
  
  # negative log-likelihood-Funktion dynamisch definiert mit den formals(), 
  # damit mle() damit arbeiten kann
  #
  # Eigentlich wollen wir eine Funktion negloglik(theta) mit einem einzigen Parametersvektor theta.
  # Die Methode mle() braucht aber eine "named list" der Parameter (z.B. mu_1=0, mu_2=0, sigma_1=2,...) und entsprechend eine 
  # Funktion negloglik(mu1, mu2, sigma1,...)
  # Da wir nicht vorher wissen, wie viele Parameter zu schätzen sind, definieren wir die formals()
  # dynamisch um
  # 
  # @param x dummy/placeholder argument, will be overwritten by formals() with list of skalar parameters
  negloglik <- function(x)
  {
    nf <- names(formals())
    
    # recover parameter vector from named arguments (mu1=...,mu2=...,sigma11,sigma12 etc).
    # stack all named arguments to parameter vector theta
    theta <- sapply(nf, function(x) {eval(parse(text=x))})
    
    # mean vector herholen
    mean <- theta[1:n]
        
    # Matrix für sigma bauen
    if (cholesky) {
      L <- inv_vech(theta[-(1:n)])
      L[lower.tri(L, diag=FALSE)] <- 0  # L entspricht jetzt chol(sigma), obere Dreiecksmatrix
      sigma <- t(L) %*% L
    } else {
      sigma <- inv_vech(theta[-(1:n)])
    }
    
    # if sigma is not positive definite, return MAXVALUE
    if (det(sigma) <= 0 || any(diag(sigma) < 0)) {
      return(.Machine$integer.max)
    }
    
    # Log-Likelihood
    # Wieso hier nur dmvnorm() : Wegen Dichte = Conditional density
    f <- -(sum(dmvnorm(X, mean, sigma, log=TRUE)) - nrow(X) * log(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)))
    
    if (is.infinite(f) || is.na(f)) {
      # cat("negloglik=",f," for parameter vector ",theta,"\n")
      # "L-BFGS-B" requires a finite function value, other methods can handle infinte values like +Inf
      # return a high finite value, e.g. integer.max, so optimize knows this is the wrong place to be
      # TODO: check whether to return +Inf or .Machine$integer.max, certain algorithms may prefer +Inf, others a finite value
      #return(+Inf)
      return(.Machine$integer.max)
    }
    
    f
  }
  formals(negloglik) <- theta
  
  # for method "L-BFGS-B" pass bounds parameter "lower.bounds" and "upper.bounds"
  # under names "lower" and "upper"
  if ((length(lower.bounds) > 1L || length(upper.bounds) > 1L || lower.bounds[1L] != 
        -Inf || upper.bounds[1L] != Inf) && method == "L-BFGS-B") {
    mle.fit <- eval.parent(substitute(mle(negloglik, start=as.list(theta), fixed=fixed, method = method, lower=lower.bounds, upper=upper.bounds, ...)))
    
    #mle.call <- substitute(mle(negloglik, start=as.list(theta), fixed=fixed, method = method, lower=lower.bounds, upper=upper.bounds, ...)) 
    #mle.fit <- mle(negloglik, start=as.list(theta), fixed=fixed, method = method, lower=lower.bounds, upper=upper.bounds, ...)
    #mle.fit@call <- mle.call
    return (mle.fit)
  } else {
    # we need evaluated arguments in the call for profile(mle.fit)
    mle.fit <- eval.parent(substitute(mle(negloglik, start=as.list(theta), fixed=fixed, method = method, ...))) 
    
    #mle.call <- substitute(mle(negloglik, start=as.list(theta), fixed=fixed, method = method, ...)) 
    #mle.fit <- mle(negloglik, start=as.list(theta), fixed=fixed, method = method, ...)
    #mle.fit@call <- mle.call
    return (mle.fit)
  }
}

# Beispiel:
if (FALSE) {

lower=c(-1,-1)
upper=c(1, 2)
mu   =c(0, 0)
sigma=matrix(c(1, 0.7,
               0.7, 2), 2, 2)
               
# generate random samples               
X <- rtmvnorm(n=500, mu, sigma, lower, upper)
method <- "BFGS"

# estimate mu and sigma from random samples
# Standard-Startwerte
mle.fit1 <- mle.tmvnorm(X, lower=lower, upper=upper)
mle.fit1a <- mle.tmvnorm(X, lower=lower, upper=upper, cholesky=TRUE)

mle.fit1b <- mle.tmvnorm(X, lower=lower, upper=upper, method="L-BFGS-B", 
  lower.bounds=c(-1, -1, 0.001, -Inf, 0.001), upper.bounds=c(2, 2, 2, 2, 3))


Rprof("mle.profile1.out") 
mle.profile1 <- profile(mle.fit1, X, method="BFGS", trace=TRUE)
Rprof(NULL)
summaryRprof("mle.profile1.out")
confint(mle.profile1)
par(mfrow=c(2,2))
plot(mle.profile1)
summary(mle.fit1)
logLik(mle.fit1)
vcov(mle.fit1)
#TODO: confint(mle.fit1)
#profile(mle.fit1)

# andere Startwerte, näher am wahren Ergebnis
mle.fit2 <- mle.tmvnorm(x=X, lower=lower, upper=upper, start=list(mu=c(0.1, 0.1), 
  sigma=matrix(c(1, 0.4, 0.4, 1.8),2,2)))
# --> funktioniert jetzt besser...
summary(mle.fit2)

# andere Startwerte, nimm mean und Kovarianz aus den Daten (stimmt zwar nicht, ist aber sicher
# ein besserer Startwert als 0 und diag(n).
mle.fit3 <- mle.tmvnorm(x=X, lower=lower, upper=upper, start=list(mu=colMeans(X), 
  sigma=cov(X)))
summary(mle.fit3)  
}