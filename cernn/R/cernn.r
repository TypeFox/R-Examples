#' Compute the regularization path for Covariance Estimate Regularized by Nuclear Norms (CERNN)
#' 
#' \code{cernn} performs stable covariance estimation over a grid of regularization parameters.
#' 
#' @param X The data matrix whose rows are observations and columns are covariates.
#' @param lambda vector of regularization parameters controling amount of shrinkage towards the target.
#' @param alpha Parameter that controls mixture between the trace and inverse trace penalties.
#' @seealso \code{get_alpha}, \code{shrink_eigen}, \code{select_lambda}
#' @references Eric C. Chi and Kenneth Lange, Stable estimation of a covariance matrix guided by nuclear norm penalties,
#' Computational Statistics and Data Analysis, 80:117-128, 2014.
#' @export
#' @examples
#' n <- 10
#' p <- 5
#' set.seed(12345)
#' X <- matrix(rnorm(n*p),n,p)
#' alpha <- get_alpha(X)
#' lambda <- 10**(seq(-1,4,length.out=100))
#' sol_path <- cernn(X,lambda,alpha)
#' df <- t(sol_path$e)
#' 
#' ## Plot regularization paths of eigenvalues
#' matplot(x=log10(lambda),y=df,type='l',ylab='shrunken eigenvalue')
#' grand_mean <- (norm(scale(X,center=TRUE,scale=FALSE),'f')**2)/(n*p)
#' abline(h=grand_mean)
cernn <- function(X,lambda,alpha) {
  n <- nrow(X)
  p <- ncol(X)
  nLambda <- length(lambda)
  X <- scale(X,center=TRUE,scale=FALSE)
  svd_sol <- svd((1/sqrt(n))*X,nv=p)
  U <- svd_sol$v
  d <- svd_sol$d**2
  if (length(d) < p) {
    d <- c(d,double(p-length(d)))
  }
  e <- shrink_eigen(d,lambda,alpha,n)
  return(list(e=e,U=U))
}

#' Selection of penalty parameter based on cross-validation
#' 
#' \code{select_lambda} selects the best regularization parameter from a grid of values based
#' on minimal predictive negative log-likelihood.
#' 
#' @param X n-by-p data matrix
#' @param lambda vector of penalties for cross-validation
#' @param fold number of folds for cross-validation
#' @export
#' @examples
#' n <- 30
#' p <- 30
#' set.seed(12345)
#' X <- matrix(rnorm(n*p),n,p)
#' alpha <- get_alpha(X)
#' lambda_max <- get_lambda_max(svd(X)$d**2,alpha,n)
#' lambda <- 10**(seq(-1,log10(lambda_max),length.out=100))
#' sol_path <- cernn(X,lambda,alpha)
#' df <- t(sol_path$e)
#' 
#' ## Plot regularization paths of eigenvalues
#' matplot(x=log10(lambda),y=df,type='l',ylab='shrunken eigenvalue')
#' grand_mean <- (norm(scale(X,center=TRUE,scale=FALSE),'f')**2)/(n*p)
#' abline(h=grand_mean)
#' 
#' ## Plot selected lambda
#' abline(v=log10(select_lambda(X,lambda)$lambda))
select_lambda <- function(X, lambda, fold=min(nrow(X),10)) {
  
  n <- nrow(X)
  p <- ncol(X)
  g <- length(lambda)
  
  sections <- cut(1:n,breaks=fold,labels=c(1:fold))
  
  negloglikelihood <- matrix(0, fold, g);
  
  for (i in c(1:fold)){
    tsindx <- which(sections==i)
    trindx <- which(sections!=i)
    
    Xtr <- scale(X[trindx,,drop=FALSE],center=TRUE,scale=FALSE)
    Xts <- scale(X[tsindx,,drop=FALSE],center=TRUE,scale=FALSE)
    
    ntr <- nrow(Xtr)
    nts <- nrow(Xts)
    
    alpha <- get_alpha(Xts)
    
    soln <- cernn( Xtr, lambda, alpha)
    
    z <- matrix(apply(Xts %*% soln$U,2,FUN=function(x){sum(x**2)}),nrow=1)
    negloglikelihood[i,] <- c((1/nts)*z%*%(1/soln$e) + apply(log(soln$e),2,sum))
  }
  
  nL <- colSums(negloglikelihood)
  minind <- floor(median(which(nL==min(nL))))
  lambda_max <- lambda[minind]
  
  alpha <- get_alpha(X)
  mean_eigenvalue <- (norm(scale(X),'f')**2)/(n*p)
  alpha <- 1 + mean_eigenvalue**2
  alpha <- 1/alpha
  sol <- cernn(X,lambda[minind],alpha)
  e <- sol$e
  U <- sol$U
  S <- U%*%diag(c(e))%*%t(U)
  return(list(S=S,e=e,U=U,lambda=lambda[minind],negL=nL,alpha=alpha))
}

#' Nonlinear shrinkage of sample eigenvalues
#' 
#' \code{shrink_eigen} shrinks the sample eigenvalues.
#' 
#' @param d Vector of sample eigenvalues to shrink. These must be nonnegative.
#' @param lambda Regularization parameter controling amount of shrinkage towards the target.
#' @param alpha Parameter that controls mixture between the trace and inverse trace penalties.
#' @param n The number of observations.
#' @export
#' @return Vector of shrunken eigenvalues.
#' @examples
#' set.seed(12345)
#' nLambda <- 100
#' lambda <- 10**seq(-2,2,length.out=nLambda)
#' alpha <- 0.5
#' n <- 10
#' p <- 5
#' d <- sort(2*runif(p))
#' e <- shrink_eigen(d,lambda,alpha,n)
#' 
#' ## Plot regularization paths of eigenvalues
#' matplot(x=log10(lambda),y=t(e),type='l',ylab='shrunken eigenvalue')
shrink_eigen <- function(d, lambda, alpha, n) {
  if (any(d < 0))
    stop("Sample eigenvalues must be nonnegative.")
  p <- length(d)
  nLambda <- length(lambda)
  e <- matrix(0,p,nLambda)
  for (i in 1:nLambda) {
    e[,i] <- (-n + sqrt(n^2 + 4*lambda[i]*alpha*(n*d + lambda[i]*(1-alpha))))/(2*lambda[i]*alpha)
  }
  return(e)
}

#' Compute alpha parameter for covariance regularization.
#' 
#' \code{get_alpha} computes the alpha parameter that shrinks eigenvalues of the sample covariance to their grand mean.
#' 
#' @param X The data matrix whose rows are observations and columns are covariates.
#' @export
#' @examples
#' n <- 10
#' p <- 5
#' set.seed(12345)
#' X <- matrix(rnorm(n*p),n,p)
#' get_alpha(X)
get_alpha = function(X) {
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center=TRUE, scale=FALSE)
  mean_eigenvalue <- (norm(X,'f')**2)/(n*p)
  alpha <- 1 + mean_eigenvalue**2
  alpha <- 1/alpha
  return(alpha)
}

#' Compute lambda_max parameter for covariance regularization.
#' 
#' \code{get_lambda_max} computes a maximum lambda value that will shrink eigenvalues nearly to the grand mean.
#' 
#' @param d Vector of sample eigenvalues to shrink. These must be nonnegative.
#' @param alpha Parameter that controls mixture between the trace and inverse trace penalties.
#' @param n The number of observations.
#' @param eps tolerance
#' @export
#' @examples
#' n <- 10
#' p <- 5
#' set.seed(12345)
#' X <- matrix(rnorm(n*p),n,p)
#' d <- svd(X)$d**2
#' alpha <- get_alpha(X)
#' get_lambda_max(d,alpha,n)
get_lambda_max = function(d,alpha,n,eps=1e-2) {
  beta <- sqrt((1 - alpha)/alpha)
  return(max(abs(0.5*(beta*n*d/(1-alpha) - n/alpha)))/(eps*beta))
}

#' Quadratic Loss
#' 
#' \code{loss_quadratic} computes the quadratic loss.
#' 
#' @param S Covariance Estimate
#' @param Sinv Reference Precision Matrix
#' @export
#' @examples
#' set.seed(12345)
#' p <- 20
#' d <- sort(abs(rcauchy(p)),decreasing=TRUE)
#' sigma <- diag(d)
#' n <- 20
#' X <- scale(matrix(rnorm(n*p),n,p),center=FALSE,scale=1/sqrt(d))
#' alpha <- get_alpha(X)
#' lambda <- 10**(seq(-2,2,length.out=100))
#' sol_cv <- select_lambda(X,lambda)
#' loss_quadratic(sol_cv$S,solve(sigma))
loss_quadratic <- function(S,Sinv) {
  eye <- diag(1,ncol(S))
  return(norm(S%*%Sinv - eye,'f')**2)
}

#' Entropy Loss
#' 
#' \code{loss_entropy} computes the entropy loss, which is also known as Stein's loss.
#' 
#' @param S Covariance Estimate
#' @param Sinv Reference Precision Matrix
#' @export
#' @examples
#' set.seed(12345)
#' p <- 20
#' d <- sort(abs(rcauchy(p)),decreasing=TRUE)
#' sigma <- diag(d)
#' n <- 20
#' X <- scale(matrix(rnorm(n*p),n,p),center=FALSE,scale=1/sqrt(d))
#' alpha <- get_alpha(X)
#' lambda <- 10**(seq(-2,2,length.out=100))
#' sol_cv <- select_lambda(X,lambda)
#' loss_entropy(sol_cv$S,solve(sigma))
loss_entropy <- function(S,Sinv) {
  p <- nrow(S)
  SiS <- Sinv %*% S
  return(sum(diag(SiS)) - sum(log(svd(SiS)$d)) - p)
}
