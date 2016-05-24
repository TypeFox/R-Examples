# added GCV calculation 12-14-2011

## TODO: add formula interface

ridge <- function(y, ...) {
	UseMethod("ridge")
}

ridge.formula <-
		function(formula, data, lambda=0, df, svd=TRUE, ...){
	
	#code from MASS:::lm.ridge
	m <- match.call(expand.dots = FALSE)
	m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <-m$df <- NULL
	m[[1L]] <- as.name("model.frame")
	m <- eval.parent(m)
	Terms <- attr(m, "terms")
	Y <- model.response(m)
	X <- model.matrix(Terms, m, contrasts)
	n <- nrow(X)
	p <- ncol(X)
	offset <- model.offset(m)
	if (!is.null(offset)) 
		Y <- Y - offset
	if (Inter <- attr(Terms, "intercept")) {
		Xm <- colMeans(X[, -Inter])
		Ym <- mean(Y)
		p <- p - 1
		X <- X[, -Inter] - rep(Xm, rep(n, p))
		Y <- Y - Ym
	}
	ridge.default(Y, X, lambda=lambda, df=df, svd=svd)
}


ridge.default <-
		function(y, X, lambda=0, df, svd=TRUE, ...){
	#dimensions	
	n <- nrow(X)
	p <- ncol(X)
	#center X and y
	Xm <- colMeans(X)
	ym <- mean(y)
	X <- X - rep(Xm, rep(n, p))
	y <- y - ym
	#scale X, as in MASS::lm.ridge 
	Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
	X <- as.matrix(X/rep(Xscale, rep(n, p)))
	
	XPX <- crossprod(X)
	XPy <- crossprod(X,y)
	I <- diag(p)
	lmfit <- lm.fit(X, y)
	MSE <- sum(lmfit$residuals^2) / (n-p-1)
	HKB <- (p - 2) * MSE/sum(lmfit$coefficients^2)
	LW <- (p - 2) * MSE * n/sum(lmfit$fitted.values^2)
	
	# from ElemStatLearn:::simple.ridge
	svd.x <- svd(X, nu = p, nv = p)
	dd <- svd.x$d
	u <- svd.x$u
	v <- svd.x$v
	if (missing(df)) {
		df <- sapply(lambda, function(x) sum(dd^2/(dd^2 + x)))
	}
	else {
		fun <- function(df, lambda) df - sum(dd^2/(dd^2 + lambda))
		lambda <- sapply(df, FUN = function(df) uniroot(f = function(lambda) fun(df, 
										lambda), lower = 0, upper = 1000, maxiter = 10000)$root)
	}
	
	# prepare output    
	coef <- matrix(0, length(lambda), p)
	cov <- as.list(rep(0, length(lambda)))
	mse <- rep(0, length(lambda))
	
	# loop over lambdas
	for(i in seq(length(lambda))) {
		lam <- lambda[i]
		XPXr <- XPX + lam * I
		XPXI <- solve(XPXr)
		coef[i,] <- XPXI %*% XPy
		cov[[i]] <- MSE * XPXI %*% XPX %*% XPXI
		res <- y - X %*% coef[i,]
		mse[i] <- sum(res^2) / (n-p) 
		dimnames(cov[[i]]) <- list(colnames(X), colnames(X))
	}
	dimnames(coef) <- list(format(lambda), colnames(X))
	
	# calculate GCV, from MASS::lm.ridge
	dn <- length(dd)
	nl <- length(lambda)
	div <- dd^2 + rep(lambda, rep(dn, nl))
	GCV <- colSums((y - X %*% t(coef))^2)/(n - colSums(matrix(dd^2/div, dn)))^2
	k <- seq_along(GCV)[GCV == min(GCV)]
	kGCV <- lambda[k]
	
	result <- list(lambda=lambda, df=df, coef=coef, cov=cov, mse=mse, scales=Xscale, kHKB=HKB, kLW=LW,
			GCV=GCV, kGCV=kGCV)
	if (svd) {
		rownames(u) <- rownames(X)
		colnames(u) <- colnames(v) <- paste("dim", 1:p, sep="")
		rownames(v) <- colnames(X)
		result <- c(result, list(svd.D=dd, svd.U=u, svd.V=v))
	}
	class(result) <- "ridge"
	result
}



coef.ridge <-
function(object, ...) {
	object$coef
}

print.ridge <-
function(x, digits = max(5, getOption("digits") - 5),...) {
  if (length(coef(x))) {
      cat("Ridge Coefficients:\n")
      print.default(format(coef(x), digits = digits), print.gap = 2, 
          quote = FALSE)
  }
  invisible(x)
}

vcov.ridge <- function(object,  ...) {
	object$cov
}
