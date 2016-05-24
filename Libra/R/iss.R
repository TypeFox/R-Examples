#' ISS solver for linear model with lasso penalty
#' 
#' Solver for the entire solution path of coefficients for ISS.  
#' 
#' The ISS solver computes the whole regularization path for 
#' lasso-penalty for linear model. It gives the piecewise constant
#'  solution path for Bregman Inverse Scale Space Differential 
#'  Inclusion. It is the asymptotic limit of LB method with kaapa 
#'  goes to infinity and alpha goes to zero.
#' 
#' @param X An n-by-p matrix of predictors
#' @param y Response Variable
#' @param intercept if TRUE, an intercept is included in the model (and 
#'  not penalized), otherwise no intercept is included. Default is TRUE.
#' @param normalize if normalize, each variable is scaled to have L2 norm
#'  square-root n. Default is TRUE.
#' @param nvar Maximal number of variables allowed in the model.
#' @return An "LB" class object is returned. The list contains the call, 
#'  the family, the path, the intercept term a0 and value for alpha, kappa,
#'  iter, and meanvalue, scale factor of X, meanx and normx.
#' @references Ohser, Ruan, Xiong, Yao and Yin, Sparse Recovery via Differential
#'  Inclusions, \url{http://arxiv.org/abs/1406.7728}
#' @author Feng Ruan, Jiechao Xiong and Yuan Yao
#' @keywords regression
#' @examples
#' #Examples in the reference paper
#' library(MASS)
#' library(lars)
#' library(MASS)
#' library(lars)
#' n = 80;p = 100;k = 30;sigma = 1
#' Sigma = 1/(3*p)*matrix(rep(1,p^2),p,p)
#' diag(Sigma) = 1
#' A = mvrnorm(n, rep(0, p), Sigma)
#' u_ref = rep(0,p)
#' supp_ref = 1:k
#' u_ref[supp_ref] = rnorm(k)
#' u_ref[supp_ref] = u_ref[supp_ref]+sign(u_ref[supp_ref])
#' b = as.vector(A%*%u_ref + sigma*rnorm(n))
#' lasso = lars(A,b,normalize=FALSE,intercept=FALSE,max.steps=100)
#' par(mfrow=c(3,2))
#' matplot(n/lasso$lambda, lasso$beta[1:100,], xlab = bquote(n/lambda), 
#'         ylab = "Coefficients", xlim=c(0,3),ylim=c(range(lasso$beta)),type='l', main="Lasso")
#' object = iss(A,b,intercept=FALSE,normalize=FALSE)
#' plot(object,xlim=c(0,3),main=bquote("ISS"))
#' kappa_list = c(4,16,64,256)
#' alpha_list = 1/10/kappa_list
#' for (i in 1:4){
#'   object <- lb(A,b,kappa_list[i],alpha_list[i],family="gaussian",group=FALSE,
#'                trate=20,intercept=FALSE,normalize=FALSE)
#'   plot(object,xlim=c(0,3),main=bquote(paste("LB ",kappa,"=",.(kappa_list[i]))))
#' }
#'

iss <- function(X, y,intercept = TRUE, normalize = TRUE,nvar = min(dim(X)))
{
	call <- match.call()
	if (!is.matrix(X) || !is.vector(y)) stop("X must be a matrix and y must be a vector!")
	if (nrow(X) != length(y)) stop("Number of rows of X must equal the length of y!")
	np <- dim(X)
  	n <- np[1]
  	p <- np[2]
  	one <- rep(1, n)
	if(intercept){
    	meanx <- drop(one %*% X)/n
    	X <- scale(X, meanx, FALSE)
    	mu <- mean(y)
    	y <- drop(y - mu)
    }
  	else {
    	meanx <- rep(0,p)
    	mu <- 0
    	y <- drop(y)
  	}
	if(normalize){
	    normx <- sqrt(drop(one %*% (X^2))/n)
    	X <- scale(X, FALSE, normx)
  	}else normx <- rep(1,p)
	
	maxitr <- 20*min(n,p)
	res <- y
	rho <- rep(0,p)
	active <- rep(FALSE,p)
	t <- 0
	hist_t <- rep(0,maxitr+1)
	hist_path <- matrix(0,maxitr+1,p)
	hist_rho <- matrix(0,maxitr+1,p)
	for (i in 1:maxitr){
		ng <- t(X)%*%res
		delta <- pmax((1-rho)/ng,(-1-rho)/ng)
		delta[active] <- Inf
		delta_t <- min(delta)
		add <- which(delta == delta_t)
		
		t <- t + delta_t
		rho[!active] <- rho[!active] + delta_t*ng[!active]
		rho[add] <- round(rho[add])
		active[add] <- TRUE
		hist_rho[i+1,] <- rho
		hist_t[i+1] <- t

		obj_nnls <- nnnpls(X[,active,drop=FALSE],y,rho[active])
		hist_path[i+1,] <- hist_path[i,]
		hist_path[i+1,active] <- obj_nnls$x		
		res <- obj_nnls$residuals
		active[active] <- (obj_nnls$x!=0)
		if (sum(active)>=min(nvar,n,p))
			break
	}
	hist_path <- hist_path[seq(i+1), ,drop=FALSE]
	#re-scale
	hist_path <- t(scale(hist_path, FALSE, normx))
	a0 <- mu - meanx%*%hist_path
	hist_rho <- t(hist_rho[seq(i+1), ,drop=FALSE])
	hist_t <- hist_t[seq(i+1)]*n
	object <- list(call = call, family="gaussian",group=FALSE, kappa = Inf, alpha = 0, path = hist_path, t = hist_t,iter = Inf, normx = normx,meanx = meanx,a0 = a0)
	class(object) <- "lb"
	return(object)
}