#' Estimate Latent Factor Matrix
#'
#' Find out the best number of factors using Bi-Cross-Validation (BCV) with 
#' Early-Stopping-Alternation (ESA) and then estimate the factor matrix.
#'
#' The model is 
#' \deqn{Y = 1 \mu' + X \beta + n^{1/2}U D V' + E \Sigma^{1/2}}
#' where \eqn{D} and \eqn{\Sigma} are diagonal matrices, \eqn{U} and \eqn{V} 
#' are orthogonal and \eqn{mu'}
#' and \eqn{V'} represent _mu transposed_ and _V transposed_ respectively.
#' The entries of \eqn{E} are assumed to be i.i.d. standard Gaussian.
#' The model assumes heteroscedastic noises and especially works well for
#' high-dimensional data. The method is based on Owen and Wang (2015).  Notice that
#' when nonnull \code{X} is given or centering the data is required (which is essentially
#' adding a known covariate with all \eqn{1}), for identifiability, it's required that
#' \eqn{<X, U> = 0} or \eqn{<1, U> = 0} respectively. Then the method will first make a rotation
#' of the data matrix to remove the known predictors or centers, and then use 
#' the latter \code{n - k} (or \code{n - k - 1} if centering is required) samples to 
#' estimate the latent factors. The rotation idea first appears in Sun et.al. (2012).
#'
#' @param Y observed data matrix. p is the number of variables and 
#' n is the sample size. Dimension is \code{c(n, p)}
#' @param X the known predictors of size \code{c(n, k)} if any. Default is NULL (no known predictors). 
#' \code{k} is the number of known covariates.
#' @param r.limit the maximum number of factor to try. Default is 20. 
#' Can be set to Inf.
#' @param niter the number of iterations for ESA. Default is 3.
#' @param nRepeat number of repeats of BCV. In other words, the random partition of \eqn{Y}
#' will be repeated for \code{nRepeat} times. Default is 12.
#' @param only.r whether only to estimate and return the number of factors.
#' @param svd.method either "fast", "propack" or "standard".
#' "fast" is using the \code{\link[corpcor]{fast.svd}} function in package corpcor to compute SVD, "propack" is using the \code{\link[svd]{propack.svd}} to compute SVD and "standard"
#' is using the \code{\link[base]{svd}} function in the base package. Because of PROPACK 
#' issues, "propack" fails for some matrices, and when that happens, 
#' the function will use "fast" to
#' compute the SVD of that matrix instead. Default method is "fast".
#' @param center logical, whether to add an intercept term in the model.
#' Default is False. 
#' @return \code{EsaBcv} returns an obejct of \code{\link[base]{class}} "esabcv" 
#' The function \code{plot} plots the cross-validation results and points out the  
#' number of factors estimated
#' An object of class "esabcv" is a list containing the following components:
#'   \item{best.r}{the best number of factor estimated}
#'   \item{estSigma}{the diagonal entries of estimated \eqn{\Sigma} 
#'   which is a vector of length \code{p}} 
#'   \item{estU}{the estimated \eqn{U}. Dimension is \code{c(n, r)}}
#'   \item{estD}{the estimated diagonal entries of \eqn{D} 
#'   which is a vector of length \code{r}}
#'   \item{estV}{the estimated \eqn{V}. Dimension is \code{c(p, r)}}
#'   \item{beta}{the estimated \eqn{\beta} which is a matrix of size \code{c(k, p)}. 
#'   Return NULL if the argument \code{X} is NULL.}
#'   \item{estS}{the estimated signal(factor) matrix \eqn{S} where
#'              \deqn{S = 1 \mu' + X \beta + n^{1/2}U D V'}}
#'   \item{mu}{the sample centers of each variable which is a vector of length
#'   \code{p}. It's an estimate of \eqn{\mu}. Return 
#'   NULL if the argument \code{center} is False.}
#'   \item{max.r}{the actual maximum number of factors used. For the details of how this is decided,
#'   please refer to Owen and Wang (2015)}
#'   \item{result.list}{a matrix with dimension \code{c(nRepeat, (max.r + 1))} storing 
#'   the detailed BCV entrywise MSE of each repeat for r from 0 to \code{max.r}}
#' @examples
#' Y <- matrix(rnorm(100), nrow = 10)
#' EsaBcv(Y)
#' @references Art B. Owen and Jingshu Wang(2015), Bi-cross-validation for factor analysis,
#' \url{http://arxiv.org/abs/1503.03515}
#' @references Yunting Sun, Nancy R. Zhang and Art B. Owen, Multiple hypothesis testing adjusted for latent variables, with an application to the AGEMAP gene expression data. The Annuals of Applied Statistics, 6(4): 1664-1688, 2012
#' @seealso \code{\link{ESA}}, \code{\link{plot.esabcv}}
#' @importFrom corpcor fast.svd
#' @importFrom svd propack.svd
#' @export
EsaBcv <- function (Y, X = NULL, r.limit = 20, niter = 3, nRepeat = 12,
                    only.r = F, svd.method = "fast", center = F){ 
	Y <- as.matrix(Y);
    p <- ncol(Y);
    n <- nrow(Y);
    X.ori <- X;
    if (center) 
    	X <- cbind(rep(1, n), X);
    qr.X <- NULL;
    Y.ori <- Y;
    if(!is.null(X)) {
		X <- as.matrix(X);
		k <- ncol(X);
		if (k >= n)
			stop("Too many predictors!! The number of predictors 
			     is expected to be much less than the sample size!")
		qr.X <- qr(X);
		Y <- qr.qty(qr.X, Y)[-(1:k), ];
		n <- n-k;
    } 
    ## decide the held-in size
    gamma <- p/n;
    bar.gamma <- ((sqrt(gamma) + sqrt(1/gamma))/2)^2;
    sqrt.rho <- sqrt(2)/(sqrt(bar.gamma) + sqrt(bar.gamma + 3));
    held.in.size.small <- min(round(sqrt.rho*sqrt(p * n)), p - 1, n - 1);
    held.in.size.large <- round(sqrt.rho^2 * p * n/held.in.size.small);
    if (n < p) {
        n1 <- held.in.size.small;
        p1 <- held.in.size.large;
    } else {
        n1 <- held.in.size.large;
        p1 <- held.in.size.small;
    }
    max.r <- min(n1, p1);
    if(max.r > r.limit)
        max.r <- r.limit;
    if (is.null(nRepeat)) {
        nRepeat <- max(round(p/p1), round(n/n1), round(p/(p-p1)), round(n/(n-n1)));
    } 

    result.list <- NULL;
    for (i in 1:nRepeat) {
       	Y.resample <- Y[sample(1:n, n), sample(1:p, p)]; 
       	
        result.list <- rbind(result.list,
                             sapply(0:max.r, function(r)
                                    BcvGivenRank(Y.resample, r, niter, n1, p1, 
                                    svd.method)));
    }
    not.na.result <- which(!is.na(colMeans(result.list)));
    max.r <- sum(not.na.result == 1:length(not.na.result)) - 1;
    result.list <- result.list[, 1:(max.r + 1)];
	colnames(result.list) <- 0:max.r;

    best.r <- as.numeric(which.min(colMeans(result.list)) - 1);

    if (only.r)
        return(list(result.list = result.list, best.r = best.r));
    est <- ESA(Y.ori, best.r, X.ori, center, niter, svd.method);
    result <- list(best.r = best.r, estSigma = est$estSigma, estU = est$estU, 
				estD = est$estD, estV = est$estV, beta = est$beta,
				estS = est$estS, mu = est$mu, max.r = max.r, 
				result.list = result.list);
	class(result) <- "esabcv"
    return(result);
}


## return a vector of BCV MSE from the four folds for a given rank
BcvGivenRank <- function (Y, r, niter, n1, p1, svd.method) {
    p <- ncol(Y);
    n <- nrow(Y);
    hp <- p - p1;
    hn <- n - n1;
    bcv.result <- BcvPerFoldGivenRank(Y[1:hn, 1:hp, drop = F], Y[1:hn, -(1:hp), drop = F],
                                        Y[-(1:hn), 1:hp, drop = F], Y[-(1:hn), -(1:hp), drop = F],
                                        r, niter, svd.method);

    return(bcv.result);
}

## calculate BCV MSE for one fold
## inputs:
# Y00, Y01, Y10, Y11: the four folds of data and Y00 is the held-out fold
# r: given number of factor to try
# niter: number of iteration steps for ESA.
## outputs:
# the average entrywiseestimation error of the held-out block. 
# If the estimate of Sigma is unreasonable, return NA
BcvPerFoldGivenRank <- function(Y00, Y01, Y10, Y11, r, niter, svd.method,
                                  tol.Sigma.scale = 6){

    if (r ==0)
        return(mean(Y00^2));
    result <- try(ESA(Y11, r, niter = niter, svd.method = svd.method));
    if (class(result)== "try-error") {
        save(Y11, r, file = "problematic.RData")
        print("Encounter problematic data!")
    }
    Sigma1 <- result$estSigma;

    if ((mean(-log10(Sigma1)) > tol.Sigma.scale - log10(max(diag(Sigma1)))) |
        (max(diag(Sigma1)) < .Machine$double.eps)) {
        return(NA);
    }

    Y11.tilde <- t(t(Y11)/ sqrt(Sigma1));
    Y01.tilde <- t(t(Y01) / sqrt(Sigma1));
    pseudo.inv.Y11.tilde <- PseudoInv(Y11.tilde, r, svd.method);
    Y00.est <- Y01.tilde %*% pseudo.inv.Y11.tilde %*% Y10;
    held.out.res <- Y00 - Y00.est;
    est.error <-  mean((held.out.res)^2);

    return (est.error);
}

## Moore-Penrose Pseudo Inverse of A Matrix
#
# Calculate the Moore-Penrose pseudo inverse of matrix \code{Y} up to a given rank
#
# This function does similar things as \code{\link[MASS]{ginv}} in package MASS while controls 
# the rank.
#
# Y the given matrix
# k the given rank
# A pseudo-inverse matrix of rank \code{k}
PseudoInv <- function(Y, k, svd.method = "fast") {
    svd.trunc <- ComputeSVD(Y, svd.method, k);
      tol <- sqrt(.Machine$double.eps);
      pseudo.inv.d <- sapply(svd.trunc$d, function(x) 
                             ifelse(x > svd.trunc$d[1] * tol, x^(-1), 0));
    pseudo.inv.Y <- svd.trunc$v %*% diag(pseudo.inv.d[1:k], k, k)    %*% 
                    t(svd.trunc$u);
      return(pseudo.inv.Y)
}

#' Estimate Latent Factor Matrix With Known Number of Factors
#'
#' Estimate the latent factor matrix and noise variance using early stopping 
#' alternation (ESA) given the number of factors.
#'
#' The model used is 
#' \deqn{Y = 1 \mu' + X \beta + n^{1/2}U D V' + E \Sigma^{1/2}}
#' where \eqn{D} and \eqn{\Sigma} are diagonal matrices, \eqn{U} and \eqn{V} 
#' are orthogonal and \eqn{\mu'}
#' and \eqn{V'} mean _mu transposed_ and _V transposed_ respectively.
#' The entries of \eqn{E} are assumed to be i.i.d. standard Gaussian.
#' The model assumes heteroscedastic noises and especially works well for
#' high-dimensional data. The method is based on Owen and Wang (2015). Notice that
#' when nonnull \code{X} is given or centering the data is required (which is essentially
#' adding a known covariate with all \eqn{1}), for identifiability, it's required that
#' \eqn{<X, U> = 0} or \eqn{<1, U> = 0} respectively. Then the method will first make a rotation
#' of the data matrix to remove the known predictors or centers, and then use 
#' the latter \code{n - k} (or \code{n - k - 1} if centering is required) samples to 
#' estimate the latent factors.
#'
#' @inheritParams EsaBcv
#' @param r The number of factors to use
#' @return The returned value is a list with components
#'   \item{estSigma}{the diagonal entries of estimated \eqn{\Sigma} 
#'   which is a vector of length \code{p}} 
#'   \item{estU}{the estimated \eqn{U}. Dimension \code{c(n, r)}}
#'   \item{estD}{the estimated diagonal entries of \eqn{D} 
#'   which is a vector of length \code{r}}
#'   \item{estV}{the estimated \eqn{V}. Dimension is \code{c(p, r)}}
#'   \item{beta}{the estimated \eqn{beta} which is a matrix of size \code{c(k, p)}. 
#'   Return NULL if the argument \code{X} is NULL.}
#'   \item{estS}{the estimated signal (factor) matrix \eqn{S} where 
#'               \deqn{S = 1 \mu' + X \beta + n^{1/2}U D V'}}
#'   \item{mu}{the sample centers of each variable which is a vector of length
#'   \code{p}. It's an estimate of \eqn{\mu}. Return 
#'   NULL if the argument \code{center} is False.}
#' @examples
#' Y <- matrix(rnorm(100), nrow = 10) + 3 * rnorm(10) %*% t(rep(1, 10))
#' ESA(Y, 1)
#' @references Art B. Owen and Jingshu Wang(2015), Bi-cross-validation for factor analysis,
#' \url{http://arxiv.org/abs/1503.03515}

#' @export
ESA <- function(Y, r, X = NULL, center = F, niter = 3, svd.method = "fast"){

    Y <- as.matrix(Y);
    p <- ncol(Y);
    n <- nrow(Y);
    Y.ori <- Y;
    if (center) 
    	X <- cbind(rep(1, n), X);
    qr.X <- NULL;
    if(!is.null(X)) {
		X <- as.matrix(X);
		k <- ncol(X);
		if (k >= n)
			stop("Too many predictors!! The number of predictors 
			     is expected to be much less than the sample size!")
		qr.X <- qr(X);
		beta <- qr.coef(qr.X, Y);
		if (center) {
			mu <- beta[1, ];
			if (k == 1) {
				beta1 <- NULL;
			} else
				beta1 <- beta[-1, ];
		} else {
			mu <- NULL;
			beta1 <- beta;
		}
		Y <- qr.qty(qr.X, Y)[-(1:k), ];
		n <- n-k;
    } else {
        beta1 <- NULL;
        mu <- NULL;
    }
    
    # initializing Sigma
    Sigma <- apply(Y, 2, var); 

    if (r == 0)
        return(list(estSigma = Sigma, estU = NULL, estD = NULL, 
                    estV = NULL, estS = NULL, beta = beta1, mu = mu));
    if (r >= min(p, n)) {
        svd.Y <- ComputeSVD(Y, svd.method = svd.method);
        if (!is.null(X)) {
        	estU <- qr.qy(qr.X, rbind(matrix(rep(0, k * r), nrow = k), 
        						      svd.Y$u))
        } else
        	estU <- svd.Y$u
        return(list(estSigma = rep(0,p), estU = estU, 
                    estD = svd.Y$d / sqrt(n),
                    estV = svd.Y$v, estS = Y.ori, beta = beta1, mu = mu));   
    } 

    iter <- 0;
    while (1) {
        iter = iter + 1;
        if( iter > niter ){ 
            break;
        }
        Y.tilde <- sapply(1:p, function(i)
                            Y[, i] / sqrt(Sigma[i]));
        svd.Ytilde <- ComputeSVD(Y.tilde, rank = r);
        U.tilde <- svd.Ytilde$u;
        V.tilde <- svd.Ytilde$v;
        estU <- U.tilde;
        estD <- diag(svd.Ytilde$d[1:r], r, r)

        estV <- sqrt(Sigma) * V.tilde;     
        res <- Y - estU %*% estD %*% t(estV);
        Sigma <- apply(res, 2, function(x) sum(x^2)/length(x));
    }
    if (!is.null(X)) 
		estU <- rbind(matrix(rep(0, k * r), nrow = k), estU);
    estS <- estU %*% estD %*% t(estV);
    svd.Y <- ComputeSVD(estS, svd.method, r);
	if (!is.null(X)) {
		estU <- qr.qy(qr.X, svd.Y$u);
		estS <- estS + X %*% beta;
	} else
		estU <- svd.Y$u;
    estD <- svd.Y$d[1:r]/sqrt(n);
    return(list(estSigma = Sigma, estU = estU, estD = estD, 
                estV = svd.Y$v, estS = estS, beta = beta1, mu = mu));    
}


# A wrapper for computing SVD
ComputeSVD <- function(Y, svd.method = "fast", rank = NULL, kmax.r = 10) {
    if (is.null(rank)) {
        rank <- min(dim(Y));
    }
	if(svd.method == "propack") {
		svd.Y <- try(suppressWarnings(propack.svd(Y, neig = rank)));
		if (class(svd.Y) == "try-error" & (rank < min(dim(Y)))) {
			svd.Y <- suppressWarnings(propack.svd(Y, neig = rank + 1));
			if (!length(svd.Y$d) == (rank + 1)) {
				svd.Y <- propack.svd(Y, neig = rank + 1, 
									 opts = list(kmax = kmax.r * (rank + 1)));
			}
			svd.Y$u <- svd.Y$u[, 1:rank];
			svd.Y$v <- svd.Y$v[, 1:rank];
			svd.Y$d <- svd.Y$d[1:rank];
		} else if ((class(svd.Y) == "try-error") | (length(svd.Y$d) != rank)) {
			warning("PROPACK fails, used fast.svd to compute SVD instead!!");
			svd.Y <- ComputeSVD(Y, "fast", rank);
		}
	} else if (svd.method == "fast") {
        tol <- max(dim(Y)) * .Machine$double.eps;

        svd.Y <- fast.svd(Y, tol);

        if (rank < min(dim(Y))) {
            svd.Y$u <- matrix(svd.Y$u[, 1:rank], ncol = rank);
            svd.Y$v <- matrix(svd.Y$v[, 1:rank], ncol = rank);
            svd.Y$d <- svd.Y$d[1:rank];
        }
    } else    {
        svd.Y <- svd(Y, nu = rank, nv = rank);
    }
    return(svd.Y)
}
