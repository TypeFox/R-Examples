###            Fit a linear quantile model (continuous and count reponses)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


##########################################################################################
# lqm functions (independent data)
# (negative) Laplace log-likelihood
"loglik.al" <- function(theta, x, y, tau) {

res <- (y - x %*% matrix(theta))
ind <- ifelse(res < 0, 1, 0)
val <- sum(res * (tau - ind))

return(val)

}


"switch_check" <- function(l0, l1, tol_ll, t0, t1, tol_theta, rule = "1"){

deltal <- abs(l1/l0 - 1)
deltat <- max(abs(t1/t0 - 1))

switch(rule,
	"1" = deltal < tol_ll,
	"2" = deltal < tol_ll && deltat < tol_theta
)

}


"lqmgsR" <- function(theta, x, y, weights, tau, control){

n <- sum(weights)
wx <- x*weights
wy <- y*weights

step <- control$loop_step
maxiter <- control$loop_max_iter

theta_0 <- theta
ll_0 <- loglik.al(theta_0, wx, wy, tau)
eps <- .Machine$double.eps^(1/2)

for(i in 1:maxiter) {
	if(control$verbose) cat(paste0("  (", i, ") logLik = ", round(ll_0,12), "\n"))
	# line search
	grad <- t(wx) %*% (tau - as.vector(wy < wx %*% theta_0))/n
	theta_1 <- theta_0 + grad*step
	ll_1 <- loglik.al(theta_1, wx, wy, tau)
	if(ll_1 < ll_0 + eps) {
		rule <- if(control$check_theta) "2" else "1"
		check <- switch_check(ll_0, ll_1, control$loop_tol_ll, theta_0, theta_1, control$loop_tol_theta, rule = rule)
		if(check) break
		theta_0 <- theta_1
			ll_0 <- ll_1
			step <- step*control$gamma
		} else {
			if(control$verbose) cat("  Decreasing step...\n")
			step <- step*control$beta
	}
}

list(theta = as.numeric(theta_1), grad = grad, optimum = ll_1, CONVERGE = if(i==maxiter) -1 else i)

}

"lqm.fit.gs" <- function(theta, x, y, weights, tau, control) {

n <- length(y)
p <- ncol(x)
if(is.null(p)) stop("x must be a matrix")
if(missing(theta)) theta <- lm.fit(as.matrix(x), y)$coefficients
if(is.null(control$loop_step)) control$loop_step <- sd(as.numeric(y))
if(length(tau) > 1) {tau <- tau[1]; warning("Length of tau is greater than 1. Only first value taken")}

if(control$method == "gs1"){
	wx <- x*weights
	wy <- y*weights
	fit <- .C("gradientSd_s", theta = as.double(theta), as.double(wx), as.double(wy), as.single(tau), as.integer(n), as.integer(p),
			  as.double(control$loop_step), as.double(control$beta), as.double(control$gamma), as.integer(control$reset_step), as.double(control$loop_tol_ll), as.double(control$loop_tol_theta), as.integer(control$check_theta), as.integer(control$loop_max_iter), as.integer(control$verbose), CONVERGE = integer(1), grad = double(p), optimum = double(1), PACKAGE = "lqmm")
} else if(control$method == "gs2") {
	fit <- lqmgsR(theta, x, y, weights, tau, control)
}

fit$residuals <- y - x%*%matrix(fit$theta)
fit$scale <- weighted.mean(fit$residuals * (tau - (fit$residuals < 0)), weights)
fit$logLik <- n*log(tau*(1-tau)/fit$scale) - 1/fit$scale * fit$optimum
OPTIMIZATION <- list(loop = fit$CONVERGE)

errorHandling(OPTIMIZATION$loop, "low", control$loop_max_iter, control$loop_tol_ll, "lqm")

list(theta = fit$theta, scale = fit$scale, gradient = fit$grad, logLik = fit$logLik, opt = OPTIMIZATION)

}

"lqmControl" <- function(method = "gs1", loop_tol_ll = 1e-5, loop_tol_theta = 1e-3, check_theta = FALSE, loop_step = NULL, beta = 0.5, gamma = 1.25, reset_step = FALSE, loop_max_iter = 1000, verbose = FALSE)
{
if(beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
if(gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
if(loop_max_iter < 0) stop("Number of iterations cannot be negative")

list(method = method, loop_tol_ll = loop_tol_ll, loop_tol_theta = loop_tol_theta, check_theta = check_theta, loop_step = loop_step, beta = beta, gamma = gamma, reset_step = reset_step, loop_max_iter = as.integer(loop_max_iter), verbose = verbose)

}

"lqm" <- function(formula, data, subset, na.action, weights = NULL, tau = 0.5, contrasts = NULL, control = list(), fit = TRUE){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
nq <- length(tau)

Call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")

y <-  model.response(mf, "numeric")
w <- as.vector(model.weights(mf))
if (!is.null(w) && !is.numeric(w)) 
  stop("'weights' must be a numeric vector")
if(is.null(w)) w <- rep(1, length(y))
x <- model.matrix(mt, mf, contrasts)
dim_theta <- ncol(x)

# Control

if(is.null(names(control))) control <- lqmControl()
    else {
    control_default <- lqmControl();
    control_names <- intersect(names(control), names(control_default));
    control_default[control_names] <- control[control_names];
    control <- control_default
    }
if(is.null(control$loop_step)) control$loop_step <- sd(as.numeric(y))
if(control$beta > 1 || control$beta < 0) stop("Beta must be a decreasing factor in (0,1)")
if(control$gamma < 1) stop("Beta must be a nondecreasing factor >= 1")


# Starting values

theta_0 <- lm.wfit(x = as.matrix(x), y = y, w = w)$coefficients

if(!fit) return(list(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau, control = control))

	   if(nq == 1){
		  fit <- lqm.fit.gs(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau, control = control)}
     else {
      fit <- vector("list", nq);
      names(fit) <- format(tau, digits = 4);
      for (i in 1:nq) fit[[i]] <- lqm.fit.gs(theta = theta_0, x = as.matrix(x), y = y, weights = w, tau = tau[i], control = control)
     }

term.labels <- colnames(x)

if(nq > 1) {
  fit$theta <- matrix(NA, dim_theta, nq);
  fit$scale <- rep(NA, nq);
 
 for(i in 1:nq){
    fit$theta[,i] <- fit[[i]]$theta;
    fit$scale[i] <- fit[[i]]$scale
    }
  rownames(fit$theta) <- term.labels;
  colnames(fit$theta) <- format(tau, digits = 4);
  }

class(fit) <- "lqm"

fit$call <- Call
fit$na.action <- attr(mf, "na.action")
fit$contrasts <- attr(x, "contrasts")
fit$term.labels <- term.labels
fit$terms <- mt
fit$nobs <- length(y)
fit$dim_theta <- fit$edf <- dim_theta
fit$rdf <- fit$nobs - fit$edf
fit$tau <- tau
fit$x <- as.matrix(x)
fit$y <- y
fit$weights <- w
fit$levels <- .getXlevels(mt, mf)
fit$InitialPar <- list(theta = theta_0)
fit$control <- control

fit
}

print.lqm <- function(x, digits = max(6, getOption("digits")), ...){

tau <- x$tau
nq <- length(tau)

if(nq == 1){
  theta <- x$theta
  names(theta) <- x$term.labels
  psi <- varAL(x$scale, tau)

  cat("Call: ")
  dput(x$call)
  cat("\n")
  cat(paste("Quantile", tau, "\n"))
  cat("Fixed effects:\n")
  print.default(format(theta, digits = digits), print.gap = 2, quote = FALSE)

  cat("\nDegrees of freedom:", x$nobs, "total;", x$rdf, "residual\n")
  cat(paste("Residual scale parameter: ",format(x$scale, digits = digits),
          " (standard deviation ",format(sqrt(psi), digits = digits),")","\n",sep=""))
  cat(paste("Log-likelihood (Laplace):", format(x$logLik, digits = digits),"\n"))
  } else {
  theta <- x$theta;
  rownames(theta) <- x$term.labels;
  colnames(theta) <- paste("tau = ", format(tau, digits = digits), sep ="")

  cat("Call: ")
  dput(x$call)
  cat("\n")
  cat("Fixed effects:\n");
  print.default(format(theta, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nDegrees of freedom:", x$nobs, "total;", x$rdf, "residual\n")
  }

invisible(x)

}

coef.lqm <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$theta

if(nq == 1){
  names(ans) <- object$term.labels
}

return(ans)

}

boot.lqm <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE){

set.seed(seed)
tau <- object$tau
nq <- length(tau)
obsS <- replicate(R, sample(1:object$nobs, replace = TRUE))
npars <- object$dim_theta

control <- object$control
control$verbose <- FALSE

if(nq == 1){
  bootmat <- matrix(NA, R, npars);
  colnames(bootmat) <- object$term.labels
  FIT_ARGS <- list(theta = object$theta, tau = tau, control = control)
  for(i in 1:R){
    a <- table(obsS[,i]);
    s <- as.numeric(names(a));
    test <- "try-error";
    while(test=="try-error"){
      FIT_ARGS$x <- as.matrix(object$x[s,])
      FIT_ARGS$y <- object$y[s]
      FIT_ARGS$weights <- as.numeric(a)
      FIT_ARGS$theta <- if(!startQR) lm.wfit(x = FIT_ARGS$x, y = FIT_ARGS$y, w = FIT_ARGS$weights)$coefficients
      fit <- try(do.call(lqm.fit.gs, FIT_ARGS));
      test <- class("fit")
      }
    bootmat[i,] <- fit$theta  
  }
} else {
  bootmat <- array(NA, dim = c(R, npars, nq), dimnames = list(NULL, object$term.labels, paste("tau = ", format(tau, digits = 4), sep ="")));
  FIT_ARGS <- list(control = control)
  for(i in 1:R){
    a <- table(obsS[,i]);
    s <- as.numeric(names(a));
    for (j in 1:nq){
      FIT_ARGS$x <- as.matrix(object$x[s,])
      FIT_ARGS$y <- object$y[s]
      FIT_ARGS$weights <- as.numeric(a)
      FIT_ARGS$theta <- if(startQR) object[[j]]$theta else lm.wfit(x = FIT_ARGS$x, y = FIT_ARGS$y, w = FIT_ARGS$weights)$coefficients
      FIT_ARGS$tau <- tau[j]
      fit <- try(do.call(lqm.fit.gs, FIT_ARGS))
      if(class(fit)!="try-error") bootmat[i,,j] <- fit$theta
    }
  }
}

class(bootmat) <- "boot.lqm"
attr(bootmat, "tau") <- tau
attr(bootmat, "estimated") <- object$theta
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "indices") <- obsS
attr(bootmat, "rdf") <- object$rdf

return(bootmat)

}

summary.boot.lqm <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), ...){

tau <- attr(object, "tau")
nq <- length(tau)

est <- attr(object, "estimated")
npars <- attr(object, "npars")
rdf <- attr(object, "rdf")
R <- attr(object, "R")


nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound", "Pr(>|t|)")
		
if(nq == 1){
  bias <- est - apply(as.matrix(object), 2, mean)
  Cov <- cov(as.matrix(object))
  stds <- sqrt(diag(Cov))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2 * pt(-abs(est/stds), R - 1)
  ans <- cbind(est, bias, stds, lower, upper, tP)
  colnames(ans) <- nn
  printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
}
else {
  bias <- est - apply(object, 3, colMeans)
  Cov <- apply(object, 3, function(x) cov(as.matrix(x)))
  if(npars == 1) Cov <- matrix(Cov, nrow = 1)
  stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2*pt(-abs(est/stds), R - 1)
  for(i in 1:nq){
    if(npars == 1){
    ans <- c(est[i], bias[i], stds[i], lower[i], upper[i], tP[i]);
    ans <- matrix(ans, nrow = 1)
    } else {ans <- cbind(est[,i], bias[,i], stds[,i], lower[,i], upper[,i], tP[,i])}
    rownames(ans) <- rownames(est)
    colnames(ans) <- nn;
    cat(paste("tau = ", tau[i], "\n", sep =""))
    printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
    cat("\n")
  }

}

}

summary.lqm <- function(object, method = "boot", alpha = 0.05, covariance = FALSE, ...){

tau <- object$tau
nq <- length(tau)
theta <- object$theta
x <- object$x
y <- object$y
n <- nrow(x)
npars <- object$dim_theta
rdf <- object$rdf

nn <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")

if(method == "boot"){
B <- boot.lqm(object, ...)
R <- attr(B, "R")
	if(nq == 1){
		Cov <- cov(as.matrix(B))
		stds <- sqrt(diag(Cov))
		tP <- 2*pt(-abs(theta/stds), R - 1)
		lower <- theta + qt(alpha/2, R - 1)*stds
		upper <- theta + qt(1 - alpha/2, R - 1)*stds
		ans <- cbind(theta, stds, lower, upper, tP)
		colnames(ans) <- nn
		} else {
		Cov <- apply(B, 3, function(x) cov(as.matrix(x)))
		if(npars == 1) Cov <- matrix(Cov, nrow = 1)
		stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
		tP <- 2*pt(-abs(theta/stds), R - 1)
		lower <- theta + qt(alpha/2, R - 1)*stds
		upper <- theta + qt(1 - alpha/2, R - 1)*stds
		ans <- vector("list", nq)
		Cov.array <- array(NA, dim = c(npars,npars,nq))
			for(i in 1:nq){
				if(npars == 1){
				ans[[i]] <- matrix(c(theta[i], stds[i], lower[i], upper[i], tP[i]), nrow = 1);
				rownames(ans[[i]]) <- rownames(theta)
				colnames(ans[[i]]) <- nn;
				} else {
				ans[[i]] <- cbind(theta[,i], stds[,i], lower[,i], upper[,i], tP[,i])
				rownames(ans[[i]]) <- rownames(theta)
				colnames(ans[[i]]) <- nn;
				}
				Cov.array[,,i] <- matrix(Cov[,i], npars, npars)
			}
		Cov <- Cov.array
		dimnames(Cov) <- list(rownames(theta), rownames(theta), format(tau, digits = 4))
		}
}
if(method == "nid"){
		eps <- .Machine$double.eps^(2/3)
		h <- bandwidth.rq(tau, n, hs = TRUE)
	if(nq == 1){
		if(tau + h > 1 || tau - h < 0) stop("bandwith is too large or tau is too close to 0 or 1")
		theta.up <- update(object, tau = tau + h, evaluate = TRUE)$theta
		#theta.up <- lqm.fit.gs(theta, x, y, object$weights, tau + h, object$control)$theta
		theta.low <- update(object, tau = tau - h, evaluate = TRUE)$theta
		#theta.low <- lqm.fit.gs(theta, x, y, object$weights, tau - h, object$control)$theta
		dyhat <- x %*% matrix(theta.up - theta.low, ncol = 1)
		dens <- pmax(0, (2 * h)/(dyhat - eps))
		fxxinv <- diag(npars)
        fxxinv <- backsolve(qr(sqrt(dens) * x)$qr[1:npars, 1:npars, drop = FALSE], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        Cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
		dimnames(Cov) <- list(colnames(x), colnames(x))
		stds <- sqrt(diag(Cov))
		tP <- 2*pt(-abs(theta/stds), object$rdf)
		lower <- theta + qt(alpha/2, object$rdf)*stds
		upper <- theta + qt(1 - alpha/2, object$rdf)*stds
		ans <- cbind(theta, stds, lower, upper, tP)
		colnames(ans) <- nn
	} else {
		if(any(tau + h > 1) || any(tau - h < 0)) stop("bandwith is too large or tau is too close to 0 or 1")
		theta.up <- update(object, tau = tau + h, evaluate = TRUE)$theta
		theta.low <- update(object, tau = tau - h, evaluate = TRUE)$theta
		dyhat <- x %*% (theta.up - theta.low)
		Cov.array <- array(NA, dim = c(npars,npars,nq))
		for(i in 1:nq){
			dens <- pmax(0, (2 * h[i])/(dyhat[,i] - eps))
			fxxinv <- diag(npars)
			fxxinv <- backsolve(qr(sqrt(dens) * x)$qr[1:npars, 1:npars, drop = FALSE], fxxinv)
			fxxinv <- fxxinv %*% t(fxxinv)
			Cov <- tau[i] * (1 - tau[i]) * fxxinv %*% crossprod(x) %*% fxxinv
			Cov.array[,,i] <- Cov
		}
		stds <- apply(Cov.array, 3, function(x) sqrt(diag(x)))
		tP <- 2*pt(-abs(theta/stds), object$rdf)
		lower <- theta + qt(alpha/2, object$rdf)*stds
		upper <- theta + qt(1 - alpha/2, object$rdf)*stds
		ans <- vector("list", nq)
		for(i in 1:nq){
			if(npars == 1){
			ans[[i]] <- matrix(c(theta[i], stds[i], lower[i], upper[i], tP[i]), nrow = 1);
			rownames(ans[[i]]) <- rownames(theta)
			colnames(ans[[i]]) <- nn;
			} else {
			ans[[i]] <- cbind(theta[,i], stds[,i], lower[,i], upper[,i], tP[,i])
			rownames(ans[[i]]) <- rownames(theta)
			colnames(ans[[i]]) <- nn;
			}
		}
		Cov <- Cov.array
		dimnames(Cov) <- list(rownames(theta), rownames(theta), format(tau, digits = 4))
	}


}

if(covariance) object$Cov <- Cov
object$tTable <- ans

class(object) <- "summary.lqm"
return(object)

}

print.summary.lqm <- function(x, ...){

tau <- x$tau
nq <- length(tau)

cat("Call: ")
dput(x$call)

  if(nq == 1){
    cat(paste("Quantile", tau, "\n"))
    printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
  } else {
    for(i in 1:nq){
      cat(paste("tau = ", tau[i], "\n", sep =""))
      printCoefmat(x$tTable[[i]], signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }  
  }

invisible(x)

}

predict.lqm <- function(object, newdata, interval = FALSE,
	level = 0.95, na.action = na.pass, ...) 
{

tau <- object$tau
nq <- length(tau)
lp <- (1 - level)/2
up <- 1 - lp

qrow <- function(x,p) apply(x, 1, function(y) quantile(y, p))

	if(missing(newdata)){
		X <- object$x 
		#yhat <- X%*%as.matrix(object$theta)
		}
	else {
		objt <- terms(object)
		Terms <- delete.response(objt)
		m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
		if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
		X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
		#yhat <- X%*%as.matrix(object$theta)
		}

	if(nq == 1){
		yhat <- drop(X %*% object$theta)
		if (interval) {
			lin_pred <- X %*% t(boot.lqm(object, ...))
			low <- qrow(lin_pred, lp)
			upp <- qrow(lin_pred, up)
			yhat <- cbind(yhat, low, upp)
			colnames(yhat) <- c("fitted", "lower", "higher")
		}
	} else {
		yhat <- matrix(X%*%as.matrix(object$theta), ncol = nq)
		if (interval) {
			B <- boot.lqm(object, ...)
			Rboot <- attr(B, "R")
			lin_pred <- matrix(apply(B, 3, function(b,x) x%*%as.matrix(t(b)), x = X), ncol = nq)
			low <- matrix(apply(lin_pred, 2, function(x, R, p) qrow(matrix(x, ncol = R),p), R = Rboot, p = lp), ncol = nq)
			upp <- matrix(apply(lin_pred, 2, function(x, R, p) qrow(matrix(x, ncol = R),p), R = Rboot, p = up), ncol = nq)
			res <- array(NA, dim = c(nrow(X),3,nq), dimnames = list(NULL, c("fitted", "lower", "higher"), format(tau, 4)))
			for(i in 1:nq) res[,,i] <- cbind(yhat[,i],low[,i],upp[,i])
			yhat <- res
		}
	}

return(yhat)
}

residuals.lqm <- function(object, ...){

ans <- as.numeric(object$y) - predict(object)
return(ans)

}

logLik.lqm <- function(object, ...){

tdf <- object$edf + 1
tau <- object$tau
nq <- length(tau)

  if(nq == 1){
  ans <- object$logLik
  } else {
    ans <- NULL
    for(i in 1:nq) ans <- c(ans, object[[i]]$logLik);
    names(ans) <- as.character(format(tau, digits = 4))
  }

attr(ans, "nobs") <- object$nobs
attr(ans, "df") <- tdf
attr(ans, "class") <- "logLik"

return(ans)
}

##########################################################################################
# QR for counts (Machado, Santos Silva)

addnoise <- function(x, centered = TRUE, B = 0.999) 
{

	n <- length(x)
    if (centered) 
        z <- x + runif(n, -B/2, B/2)
    else z <- x + runif(n, 0, B)
	
    return(z)
}

F.lqm <- function(x, cn){

xf <- floor(x)
df <- x - xf
if(df < cn & x >= 1){
	val <- xf - 0.5 + df/(2*cn)
}
if(any(cn <= df & df < (1 - cn), x < 1)){
	val <- xf
}

if(df >= (1 - cn)){
	val <- xf + 0.5 + (df - 1)/(2*cn)
}

return(val)
}

lqm.counts <- function (formula, data, weights = NULL, offset = NULL, contrasts = NULL, tau = 0.5, M = 50, zeta = 1e-5, B = 0.999, cn = NULL, alpha = 0.05, control = list()) 
{
    nq <- length(tau)
    if (nq > 1) 
        stop("One quantile at a time")
    
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")

    y <-  model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
      stop("'weights' must be a numeric vector")
    if(is.null(w))
      w <- rep(1, length(y))
    x <- model.matrix(mt, mf, contrasts)
    p <- ncol(x)
    n <- nrow(x)
	term.labels <- colnames(x)
	
    if (is.null(offset)) 
        offset <- rep(0, n)
    if (is.null(names(control))) 
        control <- lqmControl()
    else {
        control_default <- lqmControl()
        control_names <- intersect(names(control), names(control_default))
        control_default[control_names] <- control[control_names]
        control <- control_default
    }
    if (is.null(control$loop_step)) 
        control$loop_step <- sd(as.numeric(y))
    if (control$beta > 1 || control$beta < 0) 
        stop("Beta must be a decreasing factor in (0,1)")
    if (control$gamma < 1) 
        stop("Beta must be a nondecreasing factor >= 1")
	if (p == 1) 
        control$loop_tol_ll <- 0.005
    theta_0 <- glm.fit(x = as.matrix(x), y = y, weights = w, offset = offset, family = poisson())$coefficients
    
	# Add noise
	Z <- replicate(M, addnoise(y, centered = FALSE, B = B))
	# Transform Z
    TZ <- apply(Z, 2, function(x, off, tau, zeta) log(ifelse((x - 
        tau) > zeta, x - tau, zeta)) - off, off = offset, tau = tau, zeta = zeta)
	# Fit linear QR on TZ
    fit <- apply(TZ, 2, function(y, x, weights, tau, control, 
        theta) lqm.fit.gs(theta = theta, x = x, y = y, weights = weights, tau = tau, control = control), x = x, weights = w, tau = tau, control = control, theta = theta_0)
	# Trasform back
    yhat <- sapply(fit, function(obj, x) x %*% obj$theta, x = x)
	yhat <- as.matrix(yhat)
    eta <- sweep(yhat, 1, offset, "+")
    zhat <- tau + exp(eta)
	#
    Fvec <- Vectorize(F.lqm)
    if(is.null(cn)) cn <- 0.5 * log(log(n))/sqrt(n)
    F <- apply(zhat, 2, Fvec, cn = cn)
    Fp <- apply(zhat + 1, 2, Fvec, cn = cn)
    
    multiplier <- (tau - (TZ <= yhat))^2
    a <- array(NA, dim = c(p, p, M))
    for (i in 1:M) a[, , i] <- t(x * multiplier[, i]) %*% x/n
    
    multiplier <- tau^2 + (1 - 2 * tau) * (y <= (zhat - 1)) + 
        ((zhat - y) * (zhat - 1 < y & y <= zhat)) * (zhat - y - 
            2 * tau)
    b <- array(NA, dim = c(p, p, M))
    for (i in 1:M) b[, , i] <- t(x * multiplier[, i]) %*% x/n
    
    multiplier <- exp(eta) * (F <= Z & Z < Fp)
    d <- array(NA, dim = c(p, p, M))
    sel <- rep(TRUE, M)
    for (i in 1:M) {
        tmpInv <- try(solve(t(x * multiplier[, i]) %*% x/n), 
            silent = TRUE)
        if (class(tmpInv) != "try-error") 
            {d[, , i] <- tmpInv}
        else {sel[i] <- FALSE}
    }
    
    dad <- 0
    dbd <- 0
    for (i in (1:M)[sel]) {
        dad <- dad + d[, , i] %*% a[, , i] %*% d[, , i]
        dbd <- dbd + d[, , i] %*% b[, , i] %*% d[, , i]
    }
    
	m.n <- sum(sel)
    if (m.n != 0) {
		V <- dad/(m.n^2) + (1 - 1/m.n) * dbd * 1/m.n
		V <- V/n
		stds <- sqrt(diag(V))
		} else {
		stds <- NA
        warning("Standard error not available")
		}

    est <- sapply(fit, function(x) x$theta)
    est <- if (p == 1) mean(est) else rowMeans(est)

    qfit <- if (p == 1) {
        tau + exp(mean(eta[1, ]))
    } else {
        tau + exp(rowMeans(eta))
    }

	lower <- est + qt(alpha/2, n - p) * stds
	upper <- est + qt(1 - alpha/2, n - p) * stds
	tP <- 2 * pt(-abs(est/stds), n - p)

	ans <- cbind(est, stds, lower, upper, tP)
	colnames(ans) <- c("Value", "Std. Error", "lower bound", "upper bound", 
        "Pr(>|t|)")
	rownames(ans) <- names(est) <- term.labels
	
	fit <- list()
	fit$call <- call
	fit$na.action <- attr(mf, "na.action")
	fit$contrasts <- attr(x, "contrasts")
	fit$term.labels <- term.labels
	fit$terms <- mt

	fit$theta <- est
	fit$tau <- tau
	fit$nobs <- n
	fit$M <- M
	fit$Mn <- m.n
	fit$rdf <- n - p
	fit$x <- x
	fit$y <- y
	fit$fitted <- qfit
	fit$offset <- offset
	fit$Cov <- V
	fit$tTable <- ans
	fit$levels <- .getXlevels(mt, mf)
	fit$InitialPar <- list(theta = theta_0)
	fit$control <- control

	class(fit) <- "lqm.counts"
	
    return(fit)
}

coef.lqm.counts <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$theta

if(nq == 1){
  names(ans) <- object$term.labels
}

return(ans)

}

predict.lqm.counts <- function(object, newdata, na.action = na.pass, ...) 
{

tau <- object$tau

	if(missing(newdata)){
		yhat <- drop(object$x %*% object$theta)
	}
	else {
		objt <- terms(object)
		Terms <- delete.response(objt)
		m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
		if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
		x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
		yhat <- drop(x %*% object$theta)
		
	}

return(yhat)
}

residuals.lqm.counts <- function(object, ...){

ans <- as.numeric(object$y) - predict(object)
return(ans)

}

print.lqm.counts <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    tau <- x$tau
    nq <- length(tau)
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (nq == 1) {
        cat(paste("Quantile", tau, "\n"))
        cat("\n")
        cat("Fixed effects:\n")
        printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
    }
    else {
    NULL
	}
}




