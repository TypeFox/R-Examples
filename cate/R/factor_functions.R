## This file contains all the factor analysis methods.
## R comes with a factor analysis function factanal()
## There is also fa() in the package psych.
## Here we implement some recent algorithms in the case of n < p and build a high level wrapper of them.
## For p >> n, empirical results show that using all three methods give similar results

#' Factor analysis
#'
#' The main function for factor analysis with potentially high dimensional variables.
#' Here we implement some recent algorithms that is optimized for the high dimensional problem where
#' the number of samples n is less than the number of variables p.
#' 
#' @param Y data matrix, a n*p matrix
#' @param r number of factors
#' @param method algorithm to be used
#'
#' @details The three methods are quasi-maximum likelihood (ml),
#' principal component analysis (pc),
#' and factor analysis using an early stopping criterion (esa).
#'
#' The ml is iteratively solved the Expectation-Maximization algorithm
#' using the PCA solution as the initial value.
#' See Bai and Li (2012) and for more details. For the esa method, see
#' Owen and Wang (2015) for more details.
#' 
#' @return a list of objects
#' \describe{
#' \item{Gamma}{estimated factor loadings}
#' \item{Z}{estimated latent factors}
#' \item{Sigma}{estimated noise variance matrix}
#' }
#'
#' @references {
#' Bai, J. and Li, K. (2012). Statistical analysis of factor models of high dimension. \emph{The Annals of Statistics 40}, 436-465.
#' Owen, A. B. and Wang, J. (2015). Bi-cross-validation for factor analysis. \emph{arXiv:1503.03515}.
#' }
#' 
#' @examples
#' ## a factor model
#' n <- 100
#' p <- 1000
#' r <- 5
#' Z <- matrix(rnorm(n * r), n, r)
#' Gamma <- matrix(rnorm(p * r), p, r)
#' Y <- Z %*% t(Gamma) + rnorm(n * p)
#' 
#' ## to check the results, verify the true factors are in the linear span of the estimated factors.
#' pc.results <- factor.analysis(Y, r = 5, "pc") 
#' sapply(summary(lm(Z ~ pc.results$Z)), function(x) x$r.squared)
#' 
#' ml.results <- factor.analysis(Y, r = 5, "ml") 
#' sapply(summary(lm(Z ~ ml.results$Z)), function(x) x$r.squared)
#' 
#' esa.results <- factor.analysis(Y, r = 5, "esa")
#' sapply(summary(lm(Z ~ esa.results$Z)), function(x) x$r.squared)
#' 
#' @seealso \code{\link{fa.pc}}, \code{\link{fa.em}}, \code{\link[esaBcv]{ESA}}
#'
#' @import esaBcv
#' @export
#'
factor.analysis <- function(Y,
                            r,
                            method = c("ml", "pc", "esa")) {

                                        # match the arguments
    if (r == 0) {
    	return(list(Gamma = NULL, 
    				Z = NULL,
    				Sigma = apply(Y, 2, function(v) mean(v^2))))
    }
    method <- match.arg(method, c("ml", "pc", "esa"))

    if (method == "pc") {
        fa.pc(Y, r)
    } else if (method == "ml") {
        fa.em(Y, r)
    } else if(method == "esa") {
    	result <- ESA(Y, r)
    	return(list(Gamma = result$estV %*% diag(result$estD, r,r),
    				Z = sqrt(nrow(Y)) * result$estU, 
    				Sigma = result$estSigma))
    }
}

#' Factor analysis via principal components
#'
#' @inheritParams factor.analysis
#' 
#' @seealso \code{\link{factor.analysis}} for the main function.
#' 
#' @export
#'
fa.pc <- function(Y, r) {

    svd.Y <- svd(Y)
    
    Gamma <- svd.Y$v[, 1:r] %*% diag(svd.Y$d[1:r], r, r) / sqrt(nrow(Y))
    Z <- sqrt(nrow(Y)) * svd.Y$u[, 1:r]
    
    Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) mean(x^2))
    
    return(list(Gamma = Gamma,
    			 Z = Z,
                Sigma = Sigma))

}


#' Factor analysis via EM algorithm to maximize likelihood
#'
#' @inheritParams factor.analysis
#' @param tol a tolerance scale of change of log-likelihood for convergence in the EM iterations
#' @param maxiter maximum iterations
#' 
#' @seealso \code{\link{factor.analysis}} for the main function.
#' 
#' @references {
#' Bai, J. and Li, K. (2012). Statistical analysis of factor models of high dimension. \emph{The Annals of Statistics 40}, 436-465.
#' See \url{http://www.mathworks.com/matlabcentral/fileexchange/28906-factor-analysis/content/fa.m} for a MATLAB implementation of the EM algorithm.
#' }
#' 
#' @export
#'
fa.em <- function(Y, r, tol = 1e-6, maxiter = 1000) {

    ## The EM algorithm:
    ##
    ## Y = Z Gamma' + E Sigma^{1/2}                        (n * p)
    ## mle to estimate Gamma (p * r) and Sigma             (p * p)

    ## E step:
    ## EZ = Y (Gamma Gamma' + Sigma)^{-1} Gamma            (n * r)
    ## VarZ = I - Gamma' (Gamma Gamma' + Sigma)^{-1} Gamma (r * r)
    ## EZ'Z = n VarZ + EZ' * EZ                            (r * r)

    ## M step:
    ## update Gamma: Gamma = (Y' EZ)(EZ'Z)^{-1}
    ## update Sigma:
    ## Sigma = 1/n diag(Y'Y - Gamma EZ' Y - Y' EZ' Gamma' + Gamma EZ'Z Gamma')

    ## The log-likelihood (simplified)
    ## llh <- -log det(Gamma Gamma' + Sigma) - tr((Gamma Gamma' + Sigma)^{-1} S)
    ##

    ## For details, see http://cs229.stanford.edu/notes/cs229-notes9.pdf

    p <- ncol(Y)
    n <- nrow(Y)

    ## initialize parameters
	init <- fa.pc(Y, r)
    Gamma <- init$Gamma
	#Gamma <- matrix(runif(p * r), nrow = p)
    #invSigma <- 1/colMeans(Y^2)
    invSigma <- 1/init$Sigma
	#invSigma <- 1/colMeans((Y - sqrt(n) * start.svd$u %*% t(Gamma))^2)
    llh <- -Inf # log-likelihood

    ## precompute quantitites
    I <- diag(rep(1, r))
    ## diagonal of the sample Covariance S
    #sample.var <- apply(Y, 2, var)
	sample.var <- colMeans(Y^2)

    ## compute quantities needed
    tilde.Gamma <- sqrt(invSigma) * Gamma
    M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, symmetric = TRUE)
    YSG <- Y %*% (invSigma * Gamma)

	logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
    B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
    logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
    llh <- -logdetY - logtrY

	converged <- FALSE
    for (iter in 1:maxiter) {

        ## E step:
        ## Using Woodbury matrix identity:
        ## tilde.Gamma = Sigma^{-1/2} Gamma
        ## VarZ = (I + tilde.Gamma' tilde.Gamma)^{-1}
        ## EZ = Y Sigma^{-1} Gamma VarZ
        varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
        EZ <- YSG %*% varZ
        EZZ <- n * varZ + t(EZ) %*% EZ

        ## M step:
        eigenEZZ <- eigen(EZZ, symmetric = TRUE)
        YEZ <- t(Y) %*% EZ
        ## EZ'Z = G'G
        G <- sqrt(eigenEZZ$values) * t(eigenEZZ$vectors)
        ## updating invSigma
        invSigma <- 1/(sample.var - 2/n * rowSums(YEZ * Gamma) +
                           1/n * rowSums((Gamma %*% t(G))^2))
        ## updating Gamma
        Gamma <- YEZ %*% eigenEZZ$vectors %*%
            (1/eigenEZZ$values * t(eigenEZZ$vectors))

        ## compute quantities needed
        tilde.Gamma <- sqrt(invSigma) * Gamma
        M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
        eigenM <- eigen(M, T)
        YSG <- Y %*% (invSigma * Gamma)

        ## compute likelihood and check for convergence
        old.llh <- llh
        ## Ussing Woodbury matrix identity
        ## log det(Gamma Gamma' + Sigma) =
        ## log [det(invSigma^{-1})det(M)]
        ## tr((Gamma Gamma' + Sigma)^{-1} S) =
        ## tr(invSigma S) - tr(B'B)/n
        ## where B = (M')^{-1/2} YSG'
        logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
        B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
        logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
        llh <- -logdetY - logtrY

        if (abs(llh - old.llh) < tol * abs(llh)) {
            converged <- TRUE
            break
        }
    }
    
    # GLS to estimate factor loadings
    svd.H <- svd(t(Gamma) %*% (invSigma * Gamma))
    Z <- Y %*% (invSigma * Gamma) %*% (svd.H$u %*% (1/svd.H$d * t(svd.H$v))) 

    return(list(Gamma = Gamma,
                Sigma = 1/invSigma,
                Z = Z,
                niter = iter,
                converged = converged))

}

#' @describeIn est.confounder.num Estimate the number of factors
#' 
#' @export
#' 
#' @examples
#' ## example for est.factor.num
#' n <- 50
#' p <- 100
#' r <- 5
#' Z <- matrix(rnorm(n * r), n, r)
#' Gamma <- matrix(rnorm(p * r), p, r)
#' Y <- Z %*% t(Gamma) + rnorm(n * p)
#' 
#' est.factor.num(Y, method = "ed")
#' est.factor.num(Y, method = "bcv")
#' 
est.factor.num <- function(Y,
                           method = c("bcv", "ed"),
                           rmax = 20,                           
                           nRepeat = 12,
                           bcv.plot = TRUE, log = "") {
  
  method <- match.arg(method, c("bcv", "ed"))
  
  n <- nrow(Y)
  p <- ncol(Y)
    
  if (identical(method, "bcv")) {
    result <- EsaBcv(Y, nRepeat = nRepeat, r.limit = rmax)
    r <- result$best.r
    avg.sigma2 <- mean(result$estSigma)
    errors <- colMeans(result$result.list)/avg.sigma2
    if (bcv.plot) {
      plot(0:(length(errors)-1), errors, log = log, xlab = "r",
           ylab = "bcv MSE relative to the noise", type = "o")
    }
    return(list(r = r, errors = errors))
  } else if (identical(method, "ed")) {
 	  return(EigenDiff(Y, rmax))
  }

}

#' the ed method to estimate the number of factors proposed by Onatski(2010)
#' @inheritParams est.confounder.num
#' 
#' @keywords internal
EigenDiff <- function(Y, rmax = 20, niter = 10) {
	n <- nrow(Y)
	p <- ncol(Y)
	ev <- svd(Y)$d^2 / n 
	n <- length(ev)
	if (is.null(rmax)) 
		rmax <- 3 * sqrt(n)
	j <- rmax + 1
	
	diffs <- ev - c(ev[-1], 0)
	
	for (i in 1:niter) {
		y <- ev[j:(j+4)]
		x <- ((j-1):(j+3))^(2/3)
		lm.coef <- lm(y ~ x)
		delta <- 2 * abs(lm.coef$coef[2])
		idx <- which(diffs[1:rmax] > delta)
		if (length(idx) == 0) 
			hatr <- 0
		else hatr <- max(idx)
		
		newj = hatr + 1
		if (newj == j) break
		j = newj
	}
	
	return(hatr)
	
}
