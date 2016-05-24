#' CV for ISS
#' 
#' Cross-validation method to tuning the parameter t for ISS.
#' 
#' K-fold cross-validation method is used to tuning the parameter $t$ for ISS.
#' Mean square error is used as prediction error.
#'
#' @param X An n-by-p matrix of predictors
#' @param y Response Variable
#' @param K Folds number for CV. Default is 5.
#' @param t A vector of predecided tuning parameter.
#' @param intercept If TRUE, an intercept is included in the model (and not 
#' penalized), otherwise no intercept is included. Default is TRUE.
#' @param normalize if TRUE, each variable is scaled to have L2 norm 
#' square-root n. Default is TRUE.
#' @param plot.it Plot it? Default is TRUE
#' @param se Include standard error bands? Default is TRUE
#' @param \dots Additonal arguments passing to lb
#' @return A list is returned. The list contains a vector of parameter t, 
#' crossvalidation error cv.error, and the estimated standard deviation for it cv.sd
#' @references Ohser, Ruan, Xiong, Yao and Yin, Sparse Recovery via Differential
#'  Inclusions, \url{http://arxiv.org/abs/1406.7728}
#' @author Feng Ruan, Jiechao Xiong and Yuan Yao
#' @keywords Cross-validation
#' @examples
#' #Examples in the reference paper
#' library(MASS)
#' n = 200;p = 100;k = 30;sigma = 1
#' Sigma = 1/(3*p)*matrix(rep(1,p^2),p,p)
#' diag(Sigma) = 1
#' A = mvrnorm(n, rep(0, p), Sigma)
#' u_ref = rep(0,p)
#' supp_ref = 1:k
#' u_ref[supp_ref] = rnorm(k)
#' u_ref[supp_ref] = u_ref[supp_ref]+sign(u_ref[supp_ref])
#' b = as.vector(A%*%u_ref + sigma*rnorm(n))
#' cv.iss(A,b,intercept = FALSE,normalize = FALSE)
#' 

cv.iss <-function(X, y, K = 5, t, intercept = TRUE,
          normalize = TRUE,plot.it = TRUE, se = TRUE,...)
{
  if (!is.matrix(X) || !is.vector(y)) stop("X must be a matrix and y must be a vector!")
  if (nrow(X) != length(y)) stop("Number of rows of X must equal the length of y!")
  n <- dim(X)[1]
  folds <- split(sample(1:n), rep(1:K, length = n))
  if(missing(t)) {
      one <- rep(1, n)
      if(intercept){
          meanx <- drop(one %*% X)/n
          X <- scale(X, meanx, FALSE)
      }
      if(normalize){
          normx <- sqrt(drop(one %*% (X^2))/n)
          X <- scale(X, FALSE, normx)
      }
      t <- seq(100)/max(abs(y%*%X))*n
  }
  
  residmat <- matrix(0, length(t), K)
  for(i in seq(K)) {
      omit <- folds[[i]]
      fit <- iss(X[ - omit,,drop=FALSE  ], y[ - omit],intercept = intercept,
                 normalize=normalize ,...)
      fit <- predict(fit, X[omit,  ,drop=FALSE], t)$fit
      if(length(omit)==1)fit<-matrix(fit,nrow=1)
      residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
  }
  cv.error <- apply(residmat, 1, mean)
  cv.sd <- sqrt(apply(residmat, 1, var)/K)
  object<-list(t = t, cv.error = cv.error, cv.sd = cv.sd)
  if (plot.it){
      plot(t, cv.error,pch="*",type = "b", ylim = range(cv.error, cv.error + cv.sd, 
                cv.error - cv.sd), xlab=("t"),ylab="Cross-Validated MSE")
      if(se)
          segments(t, cv.error + cv.sd, t, cv.error - cv.sd, col="red")
  }
  object
}

