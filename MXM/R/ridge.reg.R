####################
#### Multivariate (and univariate) ridge regression
####################

### usage: ridge.reg(target, dataset, lambda, B = 1, newdata = NULL) 


ridge.reg <- function(target, dataset, lambda, B = 1, newdata = NULL) {
  ## target is the dependent variable and can be a matrix as well
  ## However we only use it with a univariate target variable
  ## dataset contains the independent, CONTINUOUS ONLY, variables
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## newdata is the new independent variables values 
  ## whose values of y you want to estimate
  ## by default newdata is NULL
  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  n <- length(target)  ## sample size
  p <- ncol(dataset)  ## dimensionality of dataset
  my <- mean(target)
  yy <- target - my  ## center the dependent variables
  xx <- scale(dataset)[1:n, ]  ## standardize the independent variables
  W <- solve( crossprod(xx) + lambda * diag(p) )
  beta <- W %*% crossprod(xx, yy)
  est <- xx %*% beta + my
  va <- var(target - est) * (n - 1) / (n - p - 1)
  vab <- kronecker(va, W %*% crossprod(xx) %*% W  ) 
  seb <- matrix( sqrt( diag(vab) ), nrow = p )
  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
       id <- sample(1:n, n, replace = T)
       yb <- yy[id]  ;  xb <- xx[id, ]
       W <- solve( crossprod(xb) + lambda * diag(p) )
       be[i, ] <- W %*% crossprod(xb, yb)
    }
    seb <- matrix( apply(be, 2, sd), nrow = p ) ## bootstrap standard errors of betas
  } 
  ## seb contains the standard errors of the coefficients
  if ( is.null( colnames(dataset) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(beta) <- paste("X", 1:p, sep = "")
  } else  names(seb) <- names(beta) <- colnames(dataset)
  if ( !is.null(newdata) ) {
    mx <- colMeans(dataset)
    s <- apply(dataset, 2, sd)
    s <- diag(1/s)
    newdata <- as.matrix(newdata)
    newdata <- matrix(newdata, ncol = p)
    nu <- nrow(newdata)
    newdata <- ( newdata - rep( mx, rep(nu, p) ) ) %*% s ## standardize the newdata values 
    est <- newdata %*% beta + my 
  } else est <- dataset %*% beta + my 
  est <- as.vector(est) 
  list(beta = beta, seb = seb, est = est)
}



### References
### Hoerl, A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems". Technometrics, 12(1):55?67
### Brown, P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

