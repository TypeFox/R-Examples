


#' The loss function to perform a least mean square regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
lmsRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  
  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- 0.5*(f-y)^2
    grad <- f-y
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a least absolute deviation regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
ladRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- abs(f-y)
    grad <- sign(f-y)
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a logistic regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
logisticRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- log(1+exp(-y*f))
    grad <- -y/(1+exp(-y*f))
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a quantile regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param q a numeric value in the range [0-1] defining quantile value to consider
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
quantileRegressionLoss <- function(x,y,q=0.5) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')    
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (length(q)!=1 || q<0 || q>1) stop('q must be a length one numeric in the range [0-1]')

  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- pmax(q*(f-y),(1-q)*(y-f))
    grad <- ifelse(f>y,q,q-1)
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a epsilon-insensitive regression (Vapnik et al. 1997)
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param epsilon a numeric value setting tolerance of the epsilon-regression
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
epsilonInsensitiveRegressionLoss <- function(x,y,epsilon) {
  
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- pmax(0,abs(f-y)-epsilon)
    grad <- ifelse(abs(f-y)<epsilon,0,sign(f-y))
    val <- sum(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}





