# PF, 2008-09-18
# Computes principal factor analysis for compositional data
# Uniquenesses are nor longer of diagonal form



#' Factor analysis for compositional data
#' 
#' Computes the principal factor analysis of the input data which are
#' transformed and centered first.
#' 
#' The main difference to usual implementations is that uniquenesses are nor
#' longer of diagonal form. This kind of factor analysis is designed for
#' centered log-ratio transformed compositional data. However, if the
#' covariance is not specified, the covariance is estimated from isometric
#' log-ratio transformed data internally, but the data used for factor analysis
#' are backtransformed to the clr space (see Filzmoser et al., 2009).
#' 
#' @param x (robustly) scaled input data
#' @param factors number of factors
#' @param data default value is NULL
#' @param covmat (robustly) computed covariance or correlation matrix
#' @param n.obs number of observations
#' @param subset if a subset is used
#' @param na.action what to do with NA values
#' @param start starting values
#' @param scores which method should be used to calculate the scores
#' @param rotation if a rotation should be made
#' @param maxiter maximum number of iterations
#' @param control default value is NULL
#' @param \dots arguments for creating a list
#' @return \item{loadings }{A matrix of loadings, one column for each factor.
#' The factors are ordered in decreasing order of sums of squares of loadings.}
#' \item{uniquness }{uniquness} \item{correlation }{correlation matrix}
#' \item{criteria}{The results of the optimization: the value of the negativ
#' log-likelihood and information of the iterations used.} \item{factors }{the
#' factors } \item{dof }{degrees of freedom} \item{method }{\dQuote{principal}}
#' \item{n.obs }{number of observations if available, or NA} \item{call }{The
#' matched call.} \item{STATISTIC, PVAL }{The significance-test statistic and
#' p-value, if they can be computed}
#' @author Peter Filzmoser, Karel Hron, Matthias Templ
#' @references C. Reimann, P. Filzmoser, R.G. Garrett, and R. Dutter (2008):
#' Statistical Data Analysis Explained.  \emph{Applied Environmental Statistics
#' with R}.  John Wiley and Sons, Chichester, 2008.
#' 
#' P. Filzmoser, K. Hron, C. Reimann, R. Garrett (2009): Robust Factor Analysis
#' for Compositional Data.  \emph{Computers and Geosciences}, \bold{35} (9),
#' 1854--1861.
#' @keywords multivariate
#' @export
#' @importFrom MASS ginv
#' @examples
#' 
#' data(expenditures)
#' x <- expenditures
#' res0 <- pfa(x, factors=1)
#' 
#' ## the following produce always the same result:
#' res1 <- pfa(x, factors=1, covmat="covMcd")
#' res2 <- pfa(x, factors=1, covmat=covMcd(isomLR(x))$cov)
#' res3 <- pfa(x, factors=1, covmat=covMcd(isomLR(x)))
#' 
pfa <-
function (x, factors, data = NULL, covmat = NULL, n.obs = NA, 
    subset, na.action, start = NULL, scores = c("none", "regression", 
        "Bartlett"), rotation = "varimax", maxiter = 5, control = NULL, 
    ...) 
{
	z <- -isomLR(x)    #ilr transformed data
	## orthonormal basis:
	V <- matrix(0,nrow=ncol(x),ncol=ncol(x)-1)
	for (i in 1:ncol(V)){
		V[1:i,i] <- 1/i
		V[i+1,i] <- (-1)
		V[,i] <- V[,i]*sqrt(i/(i+1))
	}
	y <- z%*%t(V)  #clr transformed data
#	if(transformation=="ilr") x <- ilr(x) else if(transformation=="clr") x <- clr(x)$x.clr else stop("This transformation is not supported by pfa().")
	x <- scale(x, scale=FALSE)  ## TODO: check this line if needed! (not in Peters code)
	sortLoadings <- function(Lambda) {
	  cn <- colnames(Lambda)
	  Phi <- attr(Lambda, "covariance")
	  ssq <- apply(Lambda, 2, function(x) -sum(x^2))
	  Lambda <- Lambda[, order(ssq), drop = FALSE]
	  colnames(Lambda) <- cn
	  neg <- colSums(Lambda) < 0
	  Lambda[, neg] <- -Lambda[, neg]
	  if (!is.null(Phi)) {
	    unit <- ifelse(neg, -1, 1)
	    attr(Lambda, "covariance") <- unit %*% Phi[order(ssq), 
	                                               order(ssq)] %*% unit
	  }
	  Lambda
	}
	cl <- match.call()
	na.act <- NULL
	if (is.character(covmat) && length(covmat)==1){
	  if(covmat == "cov"){
	    cv <- cov(z)
	    n.obs <- nrow(x)
	  } else{
	    cv <- get(covmat)(z)$cov
	    n.obs <- nrow(x)
	  }
	} else if (is.list(covmat)) {
	  if (any(is.na(match(c("cov", "n.obs"), names(covmat))))) 
	    stop("covmat is not a valid covariance list")
	  cv <- covmat$cov
	  n.obs <- covmat$n.obs
	  have.x <- FALSE
	}
	else if (is.matrix(covmat)) {
	  cv <- covmat
	  have.x <- FALSE
	}
	else if (is.null(covmat)) {
	  if (missing(x)) 
	    stop("neither x nor covmat supplied")
	  have.x <- TRUE
	  if (inherits(x, "formula")) {
	    mt <- terms(x, data = data)
	    if (attr(mt, "response") > 0) 
	      stop("response not allowed in formula")
	    attr(mt, "intercept") <- 0
	    mf <- match.call(expand.dots = FALSE)
	    names(mf)[names(mf) == "x"] <- "formula"
	    mf$factors <- mf$covmat <- mf$scores <- mf$start <- mf$rotation <- mf$control <- mf$... <- NULL
	    mf[[1]] <- as.name("model.frame")
	    mf <- eval(mf, parent.frame())
	    na.act <- attr(mf, "na.action")
	    z <- model.matrix(mt, mf)
	  }
	  else {
	    z <- as.matrix(x)
	    if (!missing(subset)) 
	      z <- z[subset, , drop = FALSE]
	  }
	  covmat <- cov.wt(z)
	  cv <- covmat$cov
	  n.obs <- covmat$n.obs
	}
	else stop("covmat is of unknown type")
	scores <- match.arg(scores)
	if (scores != "none" && !have.x) 
	  z <- x
	sds <- sqrt(diag(cv))
	cv <- cv/(sds %o% sds)
	p <- ncol(cv)
	dof <- 0.5 * ((p - factors)^2 - p - factors)
	cn <- list(nstart = 1, trace = FALSE, lower = 0.005)
	cn[names(control)] <- control
	more <- list(...)[c("nstart", "trace", "lower", "opt", "rotate")]
	if (length(more)) 
	  cn[names(more)] <- more
	if (is.null(start)) {
	  start <- (1 - 0.5 * factors/p)/diag(solve(cv))
	}
	start <- as.matrix(start)
	if (nrow(start) != p) 
	  stop(paste("start must have", p, "rows"))
	nc <- ncol(start)
	if (nc < 1) 
	  stop("no starting values supplied")
	fit <- factanal.fit.principal1(cv, factors, p = p, start = start[, 
	                                                                 1], iter.max = maxiter)
	load <- fit$loadings
	if (rotation != "none") {
	  rot <- do.call(rotation, c(list(load), cn$rotate))
	  load <- if (is.list(rot)) 
	    rot$loadings
	  else rot
	}
	fit$loadings <- sortLoadings(load)
	class(fit$loadings) <- "loadings"
	fit$na.action <- na.act
	if (scores != "none") {
	  Lambda <- fit$loadings
	  zz <- z
	  switch(scores, regression = {
	    sc <- as.matrix(zz) %*% solve(cv, Lambda)
	    if (!is.null(Phi <- attr(Lambda, "covariance"))) 
	      sc <- sc %*% Phi
	  }, Bartlett = {
	    psiinv <- ginv(fit$psi)
	    sc <- t(ginv(t(Lambda)%*%psiinv%*%Lambda)%*%t(Lambda)%*%psiinv%*%t(zz))
	    ###            d <- 1/fit$uniquenesses
	    ###            tmp <- t(Lambda * d)
	    ###            sc <- t(solve(tmp %*% Lambda, tmp %*% t(zz)))
	  })
	  rownames(sc) <- rownames(z)
	  colnames(sc) <- colnames(Lambda)
	  if (!is.null(na.act)) 
	    sc <- napredict(na.act, sc)
	  fit$scores <- sc
	}
	if (!is.na(n.obs) && dof > 0) {
	  fit$STATISTIC <- (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3) * 
	    fit$criteria["objective"]
	  fit$PVAL <- pchisq(fit$STATISTIC, dof, lower.tail = FALSE)
	}
	fit$n.obs <- n.obs
	fit$call <- cl
	fit
}
