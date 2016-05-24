#  MIXEDMODEL.R

randomizedBlock <- mixedModel2 <- function(formula, random, weights=NULL, only.varcomp=FALSE, data=list(), subset=NULL, contrasts=NULL, tol=1e-6, maxit=50, trace=FALSE)
#	REML for mixed linear models with 2 variance components
#	Gordon Smyth, Walter and Eliza Hall Institute
#	28 Jan 2003.  Last revised 20 October 2005.
{
#	Extract model from formula
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	mf$only.varcomp <- mf$tol <- mf$tol <- mf$maxit <- NULL
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	xvars <- as.character(attr(mt, "variables"))[-1]
	if((yvar <- attr(mt,"response")) > 0) xvars <- xvars[-yvar]
	xlev <- if(length(xvars) > 0) {
		xlev <- lapply(mf[xvars], levels)
		xlev[!sapply(xlev, is.null)]
	}
	y <- model.response(mf, "numeric")
	w <- model.weights(mf)
	x <- model.matrix(mt, mf, contrasts)
	random <- mf[["(random)"]]

#	Missing values not allowed
	if(any(is.na(y)) || any(is.na(x)) || any(is.na(random))) stop("Missing values not allowed")
	if(!is.null(weights)) if(any(is.na(weights))) stop("Missing values not allowed")

#	Design matrix for random effects
	lev <- unique.default(random)
	z <- 0 + (matrix(random,length(random),length(lev)) == t(matrix(lev,length(lev),length(random))))

	mixedModel2Fit(y,x,z,w=w,only.varcomp=only.varcomp,tol=tol,maxit=maxit,trace=trace)
}

randomizedBlockFit <- mixedModel2Fit <- function(y,X,Z,w=NULL,only.varcomp=FALSE,tol=1e-6,maxit=50,trace=FALSE)
#	REML for mixed linear models with 2 variance components
#	Fits the model  Y = X*BETA + Z*U + E  where BETA is fixed
#	and U is random.
#
#	GAMMA holds the variance components.  The errors E and
#	random effects U are assumed to have covariance matrices
#	EYE*GAMMA(1) and EYE*GAMMA(2) respectively.

#	Gordon Smyth, Walter and Eliza Hall Institute
#	Matlab version 19 Feb 94.  Converted to R, 28 Jan 2003.
#	Last revised 20 Oct 2005
{
#  Prior weights
if(!is.null(w)) {
	sw <- sqrt(w)
	y <- sw * y
	X <- sw * X
}

#  Find null space Q of X
X <- as.matrix(X)
Z <- as.matrix(Z)
mx <- nrow(X)
nx <- ncol(X)
nz <- ncol(Z)
fit <- lm.fit(X,cbind(Z,y))
r <- fit$rank
QtZ <- fit$effects[(r+1):mx,1:nz]

#  Apply Q to Z and transform to independent observations
mq <- mx-r
if(mq == 0) return(list(varcomp=c(NA,NA)))
s <- La.svd(QtZ,nu=mq,nv=0)
uqy <- crossprod(s$u,fit$effects[(r+1):mx,nz+1])
d <- rep(0,mq)
d[1:length(s$d)] <- s$d^2
dx <- cbind(Residual=1,Block=d)
dy <- uqy^2

#  Try unweighted starting values
dfit <- lm.fit(dx,dy)
varcomp <- dfit$coefficients
dfitted.values <- dfit$fitted.values

#  Main fit
if(mq > 2 && sum(abs(d)>1e-15)>1 && var(d)>1e-15) {
	if(all(dfitted.values >= 0))
		start <- dfit$coefficients
	else
		start <- c(Residual=mean(dy),Block=0)
#	fit gamma glm identity link to dy with dx as covariates
	dfit <- glmgam.fit(dx,dy,coef.start=start,tol=tol,maxit=maxit,trace=trace)
	varcomp <- dfit$coefficients
	dfitted.values <- dfit$fitted.values
}
out <- list(varcomp=dfit$coef)
out$reml.residuals <- uqy/sqrt(dfitted.values)
if(only.varcomp) return(out)

#  Standard errors for variance components
dinfo <- crossprod(dx,vecmat(1/dfitted.values^2,dx))
out$se.varcomp=sqrt(2*diag(chol2inv(chol(dinfo))))

#  fixed effect estimates
s <- La.svd(Z,nu=mx,nv=0)
d <- rep(0,mx)
d[1:length(s$d)] <- s$d^2
v <- drop( cbind(Residual=1,Block=d) %*% varcomp )
mfit <- lm.wfit(x=crossprod(s$u,X),y=crossprod(s$u,y),w=1/v)
out <- c(out,mfit)
out$se.coefficients <- sqrt(diag(chol2inv(mfit$qr$qr)))

out
}

