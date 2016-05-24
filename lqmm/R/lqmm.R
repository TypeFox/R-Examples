###            Fit a linear quantile mixed model
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
# lqmm functions (hierarchical data)


"Tfun" <- function(n, type = "pdSymm") {

val <- 0;

	if(type == "pdIdent"){
		val <- matrix(diag(n), ncol = 1)
	}
	
	if(type == "pdDiag"){
		val <- sapply(1:n, function(x, n) {z <- matrix(0,n,n); z[x,x] <- 1; z}, n = n)
	}
	
	if(type == "pdSymm"){
		val <- apply(diag(n*(n+1)/2), 2, invTfun, n = n, type = type)
	}
	if(type == "pdCompSymm"){
		A <- matrix(1, n, n);
		diag(A) <- rep(0, n);
		val <- if(n > 1) cbind(as.vector(diag(n)), as.vector(A)) else 1
	}

return(val)
}

"invTfun" <- function(x, n, type = "pdSymm") 
{

val <- NULL

	if(type == "pdCompSymm"){
		val <- matrix(x[2], n, n);
		diag(val) <- rep(x[1], n)
	}

    if(type == "pdSymm"){
		dim(x) <- NULL
		val <- matrix(0, n, n)
		val[lower.tri(val, diag = TRUE)] <- x
		hold <- val
		hold[upper.tri(hold, diag = TRUE)] <- 0
		val <- val + t(hold)
	}

return(val)

}

## Use `SparseGrid' version 0.8.1 (Jelmer Ypma)

"createLaguerre" <- function(rule, q, k)
{

odd <- (k %% 2) > 0
QUAD <- NULL

if(rule == "product"){
	if(odd)
		{QUAD <- gauss.quad(n = (k - 1)/2, kind = "laguerre");
		QUAD$nodes <- c(-rev(QUAD$nodes),0,QUAD$nodes)}
	else
		{QUAD <- gauss.quad(n = k/2, kind = "laguerre");
		QUAD$nodes <- c(-rev(QUAD$nodes),QUAD$nodes)}
	QUAD$weights <- dal(QUAD$nodes, mu = 0, sigma = 1, tau = 0.5, log = FALSE)
	QUAD$nodes <- permutations(n = k, r = q, v = QUAD$nodes, set = FALSE, repeats.allowed = TRUE);
	QUAD$weights <- apply(permutations(n = k, r = q, v = QUAD$weights, set = FALSE, repeats.allowed = TRUE), 1, prod)
}	

if(rule == "sparse"){
	if(odd)
		{QUAD <- createSparseGrid(function(k) gauss.quad(n = k, kind = "laguerre"), dimension = q, k = (k - 1)/2, sym = FALSE);
		QUAD$nodes <- rbind(-QUAD$nodes,rep(0,q),QUAD$nodes);
		QUAD$weights <- c(QUAD$weights,dal(rep(0,q), mu = 0, sigma = 1, tau = 0.5, log = FALSE)^q,QUAD$weights)}
	else
		{QUAD <- createSparseGrid(function(k) gauss.quad(n = k, kind = "laguerre"), dimension = q, k = k/2, sym = FALSE);
		QUAD$nodes <- rbind(-QUAD$nodes,QUAD$nodes);
		QUAD$weights <- c(QUAD$weights,QUAD$weights)}
}

return(QUAD)

}


"quad" <- function(q, k, type = c("normal","robust"), rule = 1){

if(!(rule %in% 1:4)) {warning(paste("Rule ", rule, " not recognised. Rule 1 used (see details in '?lqmm'", sep = "")); rule <- 1}

if(rule == 1){

	odd <- (k %% 2) > 0
	if(type == "normal")
		{QUAD <- gauss.quad.prob(k, type);
		if(odd) QUAD$nodes[floor(k/2)+1] = 0; # ensure middle value is 0
		QUAD$nodes <- permutations(n = k, r = q, v = QUAD$nodes, set = FALSE, repeats.allowed = TRUE);
		QUAD$weights <- apply(permutations(n = k, r = q, v = QUAD$weights, set = FALSE, repeats.allowed = TRUE), 1, prod)
	}
	   
	if(type == "robust")
		QUAD <- createLaguerre(rule = "product", q = q, k = k)
	
} else if(rule == 2){

	if(k > 25) stop("This rule is for k < 25 only")

	if(type == "normal")
		QUAD <- createSparseGrid(type = "GQN", dimension = q, k = k)
	
	if(type == "robust")
		QUAD <- createLaguerre(rule = "sparse", q = q, k = k);
	
} else if(rule == 3){

	if(k > 25) stop("This rule is for k < 25 only")
	
	QUAD <- createSparseGrid(type = "KPN", dimension = q, k = k)
	
	if(type == "robust")
		warning("Nested rule for integral with Laguerre weights not implemented. Gaussian weights used")

} else {

	if(k > 25) stop("This rule is for k < 25 only")
	rule <- 2
	
	if(type == "normal"){
		QUAD <- createSparseGrid(type = "GQN", dimension = q, k = k);
		if (length(QUAD$weights) > k^q) {
			QUAD <- createProductRuleGrid(type = "GQN", dimension = q, k = k);
			rule <- 1;
		}
	}

	if(type == "robust"){
		QUAD <- createLaguerre(rule = "sparse", q = q, k = k)
		if (length(QUAD$weights) > k^q) {
			QUAD <- createLaguerre(rule = "product", q = q, k = k);
			rule <- 1;
		}
	}
}

attr(QUAD, "rule") <- rule
attr(QUAD, "dimension") <- q
attr(QUAD, "k") <- k

return(QUAD)

}

##

"loglik.t" <- function(theta, sigma, x, y, z, weights, Tq, V, W, tau, p, q, m, M, N, Kq, minn, maxn){

if(length(theta) != (p + m)) stop("Check length theta")

ans = .C("ll_h_R", theta = as.double(theta), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W),	as.double(sigma), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M),
		as.integer(N), as.integer(Kq), as.integer(minn-1), as.integer(maxn), loglik = double(1), PACKAGE = "lqmm")

ans$loglik
}


"loglik.s" <- function(sigma, theta, x, y, z, weights, Tq, V, W, tau, p, q, m, M, N, Kq, minn, maxn){

if(length(theta) != (p + m)) stop("Check length theta")

ans = .C("ll_h_R", theta = as.double(theta), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W),	as.double(sigma), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M),
		as.integer(N), as.integer(Kq), as.integer(minn-1), as.integer(maxn), loglik = double(1), PACKAGE = "lqmm")

ans$loglik
}

##

"lqmm.fit.df" <- function(theta_0, x, y, z, weights, cov_name, V, W, sigma_0, tau, group, control){

if(length(tau) > 1) {tau <- tau[1]; warning("Length of tau is greater than 1. Only first value taken")}

x <- as.matrix(x)
z <- as.matrix(z)
V <- as.matrix(V)
W <- as.matrix(W)

p <- ncol(x)
q <- ncol(z)
m <- theta.z.dim(cov_name, q)
#cov_type <- cov.sel(cov_name)
Tq <- Tfun(n = q, type = cov_name)

N <- nrow(x)
if(length(y) != N) stop("Check dim of x and y")
Kq <- nrow(V)

ns <- as.integer(table(group))
M <- length(ns)
minn <- c(1,cumsum(ns[-M])+1)
maxn <- cumsum(ns)

if(length(weights) != M) stop("Length of \"weights\" does not match number of groups")

UP_max_iter <- control$UP_max_iter
if(UP_max_iter == 0) break("Increase number of maximum iterations", " (", UP_max_iter,")", sep = "")
r <- 0

while(r < UP_max_iter){

if(control$verbose) cat(paste("Upper loop = ", r + 1, "\n",sep=""))

ans <- optim(par = theta_0, fn = loglik.t, sigma = sigma_0, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W,	tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn, control = list(maxit = control$LP_max_iter, abstol = control$LP_tol_ll, trace = control$verbose))

if(control$verbose) cat(paste("(", r+1, ")", " logLik = ", round(-ans$value,3), "\n",sep =""))

theta_1 <- ans$par

opt_s <- optimize(f = loglik.s, interval = c(.Machine$double.eps, 1e3*sigma_0), theta = theta_1, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn)

sigma_1 <- opt_s$minimum

delta <- abs(sigma_1 - sigma_0)

if (delta < control$UP_tol) {
	break
	} else {
	 r <- r + 1;
	 theta_0 <- theta_1;
	 sigma_0 <- sigma_1}

}

low_loop <- ans$convergence

if (r  < UP_max_iter) upp_loop <- r + 1
if (r  == UP_max_iter & UP_max_iter > 0) upp_loop <- -1
if (r  == UP_max_iter & UP_max_iter == 0) upp_loop <- -2
low_loop <- if (low_loop == 1)  -1 else as.numeric(ans$counts[1])
OPTIMIZATION <- list(low_loop = low_loop, upp_loop = upp_loop)

errorHandling(OPTIMIZATION$low_loop, "low", control$LP_max_iter, control$LP_tol_ll, "lqmm")
errorHandling(OPTIMIZATION$upp_loop, "upp", control$UP_max_iter, control$UP_tol, "lqmm")

list(theta = theta_1, scale = sigma_1, logLik = -ans$value, opt = OPTIMIZATION)

}

##

"lqmm.fit.gs" <- function(theta_0, x, y, z, weights, cov_name, V, W, sigma_0, tau, group, control){

if(length(tau) > 1) {tau <- tau[1]; warning("Length of tau is greater than 1. Only first value taken")}

x <- as.matrix(x)
z <- as.matrix(z)
V <- as.matrix(V)
W <- as.matrix(W)

p <- ncol(x)
q <- ncol(z)
m <- theta.z.dim(cov_name, q)
#cov_type <- cov.sel(cov_name)
Tq <- Tfun(n = q, type = cov_name)

N <- nrow(x)
if(length(y) != N) stop("Check dim of x and y")
Kq <- nrow(V)

ns <- as.integer(table(group))
M <- length(ns)
minn <- c(1,cumsum(ns[-M])+1)
maxn <- cumsum(ns)

if(length(weights) != M) stop("Length of \"weights\" does not match number of groups")

if(is.null(control$LP_step)) control$LP_step <- sd(as.numeric(y))
UP_max_iter <- control$UP_max_iter
if(UP_max_iter == 0) break("Increase number of maximum iterations", " (", UP_max_iter,")", sep = "")
r <- 0

while(r < UP_max_iter){

if(control$verbose) cat(paste("Upper loop = ", r + 1, "\n",sep=""))

ans <- .C("gradientSd_h", theta = as.double(theta_0), as.double(x), as.double(y), as.double(z), as.double(weights), as.double(Tq), as.double(V), as.double(W),
		as.double(sigma_0), as.single(tau), as.integer(p), as.integer(q), as.integer(m), as.integer(M), as.integer(N),	as.integer(Kq), as.integer(minn-1),
		as.integer(maxn), as.double(control$LP_step), as.double(control$beta), as.double(control$gamma), as.integer(control$reset_step),
		as.double(control$LP_tol_ll), as.double(control$LP_tol_theta), as.integer(control$check_theta), as.integer(control$LP_max_iter),
		as.integer(control$verbose), low_loop = integer(1), double(1), grad = double(p + m), opt_val = double(1), PACKAGE = "lqmm")

theta_1 <- ans$theta
grad <- ans$grad

opt_s <- optimize(f = loglik.s, interval = c(.Machine$double.eps, 10*sigma_0), theta = theta_1, x = x, y = y, z = z, weights = weights, Tq = Tq, V = V, W = W, tau = tau, p = p, q = q, m = m, M = M, N = N, Kq = Kq, minn = minn, maxn = maxn)

sigma_1 <- opt_s$minimum

delta <- abs(sigma_1 - sigma_0)

if (delta < control$UP_tol) {
	break
	} else {
	 r <- r + 1;
	 theta_0 <- theta_1;
	 sigma_0 <- sigma_1}

}

if (r  < UP_max_iter) upp_loop <- r + 1
if (r  == UP_max_iter & UP_max_iter > 0) upp_loop <- -1
if (r  == UP_max_iter & UP_max_iter == 0) upp_loop <- -2
OPTIMIZATION <- list(low_loop = ans$low_loop, upp_loop = upp_loop)

errorHandling(OPTIMIZATION$low_loop, "low", control$LP_max_iter, control$LP_tol_ll, "lqmm")
errorHandling(OPTIMIZATION$upp_loop, "upp", control$UP_max_iter, control$UP_tol, "lqmm")

list(theta = theta_1, scale = sigma_1, gradient = grad, logLik = -ans$opt_val, opt = OPTIMIZATION)

}

##

lqmmControl <- function(method = "gs", LP_tol_ll = 1e-5, LP_tol_theta = 1e-5, check_theta = FALSE, LP_step = NULL, beta = 0.5, gamma = 1, reset_step = FALSE, LP_max_iter = 500, UP_tol = 1e-4, UP_max_iter = 20, startQR = FALSE, verbose = FALSE){

if(beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
if(gamma < 1) stop("Beta must be a nondecreasing factor >= 1")
if(LP_max_iter < 0 || UP_max_iter < 0) stop("Number of iterations cannot be negative")

list(method = method, LP_tol_ll = LP_tol_ll, LP_tol_theta = LP_tol_theta, check_theta = check_theta, LP_step = LP_step, beta = beta, gamma = gamma, reset_step = reset_step, LP_max_iter = as.integer(LP_max_iter), UP_tol = UP_tol, UP_max_iter = as.integer(UP_max_iter), startQR = startQR, verbose = verbose)

}


errorHandling <- function(code, type, maxit, tol, fn){

txt <- switch(type, low = "Lower loop", upp = "Upper loop")

if(code == -1) warning(paste(txt, " did not converge in: ", fn, ". Try increasing max number of iterations ", "(", maxit, ") or tolerance (", tol,
")\n", sep =""))
if(code == -2) warning(paste(txt, " did not start in: ", fn, ". Check max number of iterations ", "(", maxit, ")\n", sep =""))

}

"lqmm" <- function(fixed, random, group, covariance = "pdDiag", tau = 0.5, nK = 7, type = "normal", rule = 1, data = sys.frame(sys.parent()), subset, weights, na.action = na.fail, control = list(), contrasts = NULL, fit = TRUE)
{
 
Call <- match.call()
  
if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
nq <- length(tau)

if(!is.data.frame(data)) stop("`data' must be a data frame")
if(!(type %in% c("normal","robust"))) stop("type must be either `normal' or `robust'")

  
# check arguments
if(!inherits(fixed, "formula") || length(fixed) != 3) {
	stop("\nFixed-effects model must be a formula of the form \"y ~ x\"")
}
if(!inherits(random, "formula") || length(random) != 2) {
	stop("\nRandom-effects model must be a formula of the form \" ~ x\"")
}

groupFormula <- asOneSidedFormula(Call[["group"]])
group <- groupFormula[[2]]
  
## extract data frame with all necessary information

mfArgs <- list(formula = asOneFormula(random, fixed, group), data = data, na.action = na.action)
if(!missing(subset)) {
	mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
}
if(!missing(weights)) {
	mfArgs[["weights"]] <- weights
}
mfArgs$drop.unused.levels <- TRUE
dataMix <- do.call("model.frame", mfArgs)
origOrder <- row.names(dataMix)	# preserve the original order
for(i in names(contrasts)) contrasts(dataMix[[i]]) = contrasts[[i]]

## sort the model.frame by groups
grp <- model.frame(groupFormula, dataMix)

## ordering data by groups
ord <- order(unlist(grp, use.names = FALSE))
grp <- grp[ord,,drop = TRUE]
dataMix <- dataMix[ord, ,drop = FALSE]
revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
ngroups <- length(unique(grp))

## obtaining model matrices and vectors
y <- eval(fixed[[2]], dataMix)
mmr <- model.frame(random, dataMix)
mmr <- model.matrix(random, data = mmr)
# Likelihood weights
if(!missing(weights)) weights <- model.weights(dataMix)[!duplicated(grp)]
if(!missing(weights) && is.null(weights)) weights <- rep(1,ngroups)
if(missing(weights))  weights <- rep(1,ngroups)
# keeping the contrasts for use in predict
contr <- attr(mmr, "contr")
mmf <- model.frame(fixed, dataMix)
Terms <- attr(mmf, "terms")
auxContr <- lapply(mmf, function(el)
	if (inherits(el, "factor") && length(levels(el)) > 1) contrasts(el))
	
contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
contr <- contr[!unlist(lapply(contr, is.null))]
mmf <- model.matrix(fixed, data = mmf)

## define dimensions
cov_name <- covariance
if(type == "robust" & !(cov_name %in% c("pdIdent","pdDiag"))) stop("Gauss-Laguerre quadrature available only for uncorrelated random effects.")
  
dim_theta <- integer(2)
dim_theta[1] <- ncol(mmf)
dim_theta[2] <- ncol(mmr)
dim_theta_z <- theta.z.dim(type = cov_name, n = dim_theta[2])

## Check if product rule quadrature is computationally heavy

if(rule == 1){
	if(dim_theta[2] > 4 && nK > 11) {
	warning(paste("For current value of \"nK\" the total number of quadrature knots is ", nK^dim_theta[2], sep = ""))
	}
}

## Quandrature nodes and weights
QUAD <- quad(q = dim_theta[2], k = nK, type = type, rule = rule)

## Control
if(is.null(names(control))) control <- lqmmControl()
    else {
    control_default <- lqmmControl();
    control_names <- intersect(names(control), names(control_default));
    control_default[control_names] <- control[control_names];
    control <- control_default
    }
if(is.null(control$LP_step)) control$LP_step <- sd(as.numeric(y))
method <- control$method
if(method == "gs" & cov_name == "pdCompSymm") {method <- "df"; cat("Switching to Nelder-Mead optimization \n")}

# Starting values
theta_z <- if(type == "normal") rep(1, dim_theta_z) else rep(invvarAL(1, 0.5), dim_theta_z)
lmfit <- lm.wfit(x = mmf, y = y, w = rep(weights, table(grp)))
theta_x <- lmfit$coefficients

	if(control$startQR){
		q_close <- if(nq == 1) tau else 0.5
		fit_rq <- lqm.fit.gs(theta = theta_x, x = as.matrix(mmf), y = y, weights = rep(weights, table(grp)), tau = q_close,
                control = lqmControl(loop_step = sd(as.numeric(y))));
		theta_x <- fit_rq$theta;
		sigma_0 <- fit_rq$scale
	} else {
		sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
	}
theta_0 <- c(theta_x, theta_z)

## Create list with all necessary arguments
FIT_ARGS <- list(theta_0 = theta_0, x = as.matrix(mmf), y = y, z = as.matrix(mmr), weights = weights, V = QUAD$nodes, W = QUAD$weights, sigma_0 = sigma_0, tau = tau, group = grp, cov_name = cov_name, control = control)

if(!fit) return(FIT_ARGS)

## Estimation
	if(method == "gs"){
	   if(nq == 1){
		  fit <- do.call(lqmm.fit.gs, FIT_ARGS)}
     else {
      fit <- vector("list", nq);
      names(fit) <- format(tau, digits = 4);
      for (i in 1:nq){
		    FIT_ARGS$tau <- tau[i];
			fit[[i]] <- do.call(lqmm.fit.gs, FIT_ARGS)
      }
     }
	}
	if(method == "df"){
		if(nq == 1){
		  fit <- do.call(lqmm.fit.df, FIT_ARGS)}
		else {
			fit <- vector("list", nq);
			names(fit) <- format(tau, digits = 4);
			for (i in 1:nq){
				FIT_ARGS$tau <- tau[i];
				fit[[i]] <- do.call(lqmm.fit.df, FIT_ARGS)
			}
		 }
	}

nn <- colnames(mmf)
mm <- colnames(mmr)

if(nq > 1) {
  fit$theta_x <- matrix(NA, dim_theta[1], nq);
  fit$theta_z <- matrix(NA, dim_theta_z, nq);
 for(i in 1:nq){
    fit$theta_x[,i] <- fit[[i]]$theta_x <- fit[[i]]$theta[1:dim_theta[1]];
    fit$theta_z[,i] <- fit[[i]]$theta_z <- fit[[i]]$theta[-(1:dim_theta[1])]
    }
  rownames(fit$theta_x) <- nn;
  colnames(fit$theta_x) <- colnames(fit$theta_z) <- format(tau, digits = 4);
  }
else {
  fit$theta_x <- fit$theta[1:dim_theta[1]];
  fit$theta_z <- fit$theta[-(1:dim_theta[1])]
  }


fit$call <- Call
fit$nn <- nn
fit$mm <- mm
fit$nobs <- length(y)
fit$dim_theta <- dim_theta
fit$dim_theta_z <- dim_theta_z
fit$edf <- fit$dim_theta[1] + fit$dim_theta_z
fit$rdf <- fit$nobs - fit$edf
fit$df <- dim_theta[1] +  dim_theta_z + 1
fit$tau <- tau
fit$mmf <- as.matrix(mmf)
fit$mmr <- as.matrix(mmr)
fit$y <- y
fit$revOrder <- revOrder
fit$weights <- weights
fit$contrasts <- contr
fit$group <- grp
fit$ngroups <- ngroups
fit$QUAD <- QUAD
fit$type <- type
fit$rule <- rule
fit$InitialPar <- list(theta = theta_0, sigma = sigma_0)
fit$control <- control
fit$cov_name <- cov_name
#attr(fit$cov_name, "cov_type") <- cov.sel(cov_name)
fit$mfArgs <- mfArgs

class(fit) <- "lqmm"
fit
}

theta.z.dim <- function(type, n){
	switch(type,
		"pdIdent" = 1,
		"pdDiag" = n,
		"pdCompSymm" = if(n == 1) 1 else 2,
		"pdSymm" = n*(n+1)/2)
}


covHandling <- function(theta, n, cov_name, quad_type){

	if(cov_name %in% c("pdIdent","pdDiag")){
		if(quad_type == "robust"){
			sigma <- theta;
			if(any(sigma < 0)){
				warning("Not positive-definite variance-covariance of random effects.");
				sigma[sigma < 0] <- .Machine$double.eps
			}
			sigma <- varAL(sigma, 0.5);
		} else {
			sigma <- theta;
			if(any(sigma < 0)){
				warning("Not positive-definite variance-covariance of random effects.");
				sigma[sigma < 0] <- .Machine$double.eps
			}
			sigma <- sigma^2;
		}
	}


	if(cov_name == "pdCompSymm"){
		if(quad_type == "robust"){
			stop("Not implemented yet: Gauss-Laguerre quadrature requires uncorrelated random effects.")
		} else {
			sigma <- as.matrix(invTfun(x = theta, n = n, type = cov_name));
			sigma <- sigma%*%sigma;
			if(!is.positive.definite(sigma)){
				warning("Not positive-definite variance-covariance of random effects.");
				sigma <- make.positive.definite(sigma)
			}
		}
	}

	if(cov_name == "pdSymm"){
		if(quad_type == "robust"){
			stop("Not implemented yet: Gauss-Laguerre quadrature requires uncorrelated random effects.")
		} else {
			sigma <- as.matrix(invTfun(x = theta, n = n, type = cov_name));
			sigma <- sigma%*%sigma;
			if(!is.positive.definite(sigma)){
				warning("Not positive-definite variance-covariance of random effects.");
				sigma <- make.positive.definite(sigma)
			}
		}
	}

return(sigma)	
}

VarCorr.lqmm <- function(x, sigma = NULL, ...){

tau <- x$tau
nq <- length(tau)

theta_z <- x$theta_z
dim_theta <- x$dim_theta
q <- dim_theta[2]
cov_name <- x$cov_name
type <- x$type
mm <- x$mm

	if(nq == 1){
		sigma <- covHandling(theta = theta_z, n = q, cov_name = cov_name, quad_type = type);
		if(cov_name == "pdIdent") {sigma <- rep(sigma, q); names(sigma) <- mm}
		if(cov_name == "pdDiag") {names(sigma) <- mm}
		if(cov_name == "pdCompSymm") {rownames(sigma) <- colnames(sigma) <- mm}
		if(cov_name == "pdSymm") {rownames(sigma) <- colnames(sigma) <- mm}
	} else {
		sigma <- vector("list", nq);
		names(sigma) <- format(tau, digits = 4);
		for(i in 1:nq){
			sigma[[i]] <- covHandling(theta = theta_z[,i], n = q, cov_name = cov_name, quad_type = type)
			if(cov_name == "pdIdent") {sigma[[i]] <- rep(sigma[[i]], q); names(sigma[[i]]) <- mm}
			if(cov_name == "pdDiag") {names(sigma[[i]]) <- mm}
			if(cov_name == "pdCompSymm") {rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm}
			if(cov_name == "pdSymm") {rownames(sigma[[i]]) <- colnames(sigma[[i]]) <- mm}
		}
	}

return(sigma)

}


print.lqmm <- function(x, digits = max(3, getOption("digits") - 3), ...){

tau <- x$tau
nq <- length(tau)

if(nq == 1){
  theta_x <- x$theta_x
  names(theta_x) <- x$nn
  sigma <- VarCorr(x)
  psi <- varAL(x$scale, tau)

  cat("Call: ")
  dput(x$call)
  cat("\n")
  cat(paste("Quantile", tau, "\n"))
  cat("\n")

  cat("Fixed effects:\n")
  print.default(format(theta_x, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Covariance matrix of the random effects:\n")
  print.default(format(sigma, digits = digits), quote = FALSE)

  cat("\n")
  cat(paste("Residual scale parameter: ",format(x$scale, digits = digits),
          " (standard deviation ",format(sqrt(psi), digits = digits),")","\n",sep=""))
  cat(paste("Log-likelihood:", format(x$logLik, digits = digits),"\n"))
  cat(paste("\nNumber of observations:", length(x$y), "\n"))
  cat(paste("Number of groups:", x$ngroups, "\n"))
  
} else {
  theta_x <- x$theta_x
  colnames(theta_x) <- paste("tau = ", format(tau, digits = digits), sep ="")
  rownames(theta_x) <- x$nn
  Scale <- sapply(x[1:nq], function(z) z$scale)
  psi <- varAL(sigma = Scale, tau = tau)
  sigma <- VarCorr(x)
  ll <- sapply(x[1:nq], function(z) z$logLik)

  cat("Call: ")
  dput(x$call)
  cat("\n")
  cat("Fixed effects:\n")
  print.default(format(theta_x, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Covariance matrix of the random effects:\n")
  for(i in 1:nq){
  cat(paste("tau = "), tau[i], "\n", sep = "")
  print.default(format(sigma[[i]], digits = digits), quote = FALSE)
  }
  
  cat("\n")
  cat("Residual scale parameter: ")
  cat(paste(format(Scale, digits = digits), " (tau = ", tau, ") ", sep = ""))
  cat("\n")
  cat("Log-likelihood: ")
  cat(paste(format(ll, digits = digits), " (tau = ", tau, ") ", sep = ""))
  cat("\n")
  cat(paste("\nNumber of observations:", length(x$y), "\n"))
  cat(paste("Number of groups:", x$ngroups, "\n"))

}

invisible(x)

}

coef.lqmm <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$theta_x

if(nq == 1){
names(ans) <- object$nn
}

return(ans)

}


ranef.lqmm <- function(object, ...){

tau <- object$tau
nq <- length(tau)
group <- object$group
M <- object$ngroups
BLPu <- vector("list", M)
q <- object$dim_theta[2]
mmr.l <- split(object$mmr, group)
sigma <- VarCorr(object)
cov_name <- object$cov_name

if(cov_name %in% c("pdIdent","pdDiag")){
  if(nq ==1) sigma <- diag(x = sigma, nrow = length(sigma))
    else for(j in 1:nq) sigma[[j]] <- diag(x = sigma[[j]], nrow = length(sigma[[j]]))
}

if(nq == 1){
  psi <- varAL(object$scale, tau)
  INV <- lapply(mmr.l, function(x, a, b, q) {x <- matrix(x, ncol = q); n <- nrow(x); y <- x%*%a%*%t(x) + diag(b, n); solve(y)},
                                      a = sigma, b = psi, q = q)
  GZ <- lapply(mmr.l, function(x, a, q) {x <- matrix(x, ncol = q); a%*%t(x)}, a = sigma, q = q)
  RES <- split(object$y - object$mmf%*%matrix(object$theta_x) - meanAL(0, object$scale, tau), group)
  for(i in 1:M){
    BLPu[[i]] <- GZ[[i]]%*%INV[[i]]%*%matrix(RES[[i]])
  }
ans <- data.frame(matrix(unlist(BLPu), ncol = q, byrow = TRUE));
rownames(ans) <- unique(group);
colnames(ans) <- object$mm;
} else {
  ans <- vector("list", nq)
  for(j in 1:nq){  
    tmp <- object[[j]];
    psi <- varAL(tmp$scale, tau[j])
    INV <- lapply(mmr.l, function(x, a, b, q) {x <- matrix(x, ncol = q); n <- nrow(x); y <- x%*%a%*%t(x) + diag(b, n); solve(y)},
                                      a = sigma[[j]], b = psi, q = q)
    GZ <- lapply(mmr.l, function(x, a, q) {x <- matrix(x, ncol = q); a%*%t(x)}, a = sigma[[j]], q = q)
    RES <- split(object$y - object$mmf%*%matrix(tmp$theta_x) - meanAL(0, tmp$scale, tau[j]), group)
      for(i in 1:M){
        BLPu[[i]] <- GZ[[i]]%*%INV[[i]]%*%matrix(RES[[i]])
      }
   ans[[j]] <- data.frame(matrix(unlist(BLPu), ncol = q, byrow = TRUE));
   rownames(ans[[j]]) <- unique(group);
   colnames(ans[[j]]) <- object$mm;
  }
 names(ans) <- format(tau, digits = 4)
 }
return(ans)
}


predict.lqmm <- function(object, level = 0, ...){

tau <- object$tau
nq <- length(tau)
group <- object$group
M <- object$ngroups
q <- object$dim_theta[2]

if(nq == 1){
  FXD <- object$mmf%*%matrix(object$theta_x)
} else {
  FXD <- object$mmf%*%object$theta_x
}

if(level == 1){
  RE <- ranef(object);
  mmr.l <- split(object$mmr, group)
  if(nq == 1){
    RE.l <- split(RE, unique(group));
    RND <- NULL;
    for(i in 1:M) {RND <- rbind(RND, matrix(as.numeric(mmr.l[[i]]), ncol = q)%*%matrix(as.numeric(RE.l[[i]]), nrow = q))}
    } else {
      RND <- matrix(NA, length(object$y), nq)
      for(j in 1:nq){
      RE.l <- split(RE[[j]], unique(group));
      tmp <- NULL;
      for(i in 1:M) {tmp <- rbind(tmp, matrix(as.numeric(mmr.l[[i]]), ncol = q)%*%matrix(as.numeric(RE.l[[i]]), nrow = q))}
      RND[,j] <- tmp
      }
    }
}

if(level == 0) {
	colnames(FXD) <- format(tau, digits = 4)
	ans <- FXD[object$revOrder,]
}
if(level == 1) {
	ans <- FXD + RND
	colnames(ans) <- format(tau, digits = 4)
	ans <- ans[object$revOrder,]
	}

return(ans)
}


predint.lqmm <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000))){

tau <- object$tau
nq <- length(object$tau)
p <- object$dim_theta[1]
m <- object$dim_theta_z

B <- boot(object, R = R, seed = seed)
tmp <- object

if(nq == 1){
	yhat <- matrix(NA, object$nobs, R)
	for(i in 1:R){
		tmp$theta <- B[i,1:(p+m)]
		tmp$theta_x <- B[i,1:p]
		tmp$theta_z <- B[i,(p+1):(p+m)]
		tmp$scale <- B[i,(p+m+1)]
		yhat[,i] <- predict(tmp, level = level)
	}
	LB <- apply(yhat, 1, quantile, probs = alpha/2)
	UB <- apply(yhat, 1, quantile, probs = 1 - alpha/2)
	ans <- data.frame(yhat = predict(object, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
} else {
	ans <- list()
	for(j in 1:nq){
		tmp$tau <- tau[j]
		yhat <- matrix(NA, object$nobs, R)
		for(i in 1:R){
			tmp$theta <- B[i,1:(p+m),j]
			tmp$theta_x <- B[i,1:p,j]
			tmp$theta_z <- B[i,(p+1):(p+m),j]
			tmp$scale <- B[i,(p+m+1),j]
			yhat[,i] <- predict(tmp, level = level)
		}
		LB <- apply(yhat, 1, quantile, probs = alpha/2)
		UB <- apply(yhat, 1, quantile, probs = 1 - alpha/2)
		ans[[j]] <- data.frame(yhat = predict(tmp, level = level), lower = LB, upper = UB, SE = apply(yhat, 1, sd))
	}
	names(ans) <- format(tau, digits = 4)
}

return(ans)

}


residuals.lqmm <- function(object, level = 0, ...){

object$y[object$revOrder] - predict(object, level = level)

}


logLik.lqmm <- function(object, ...){

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


summary.lqmm <- function(object, method = "boot", alpha = 0.05, covariance = FALSE, ...){

tau <- object$tau
nq <- length(tau)
object$logLik <- logLik(object)
object$aic <- AIC(object)
est <- extractAll(object)
nn <- object$nn
rdf <- object$rdf
ddf <- object$dim_theta[1] - 1
npars <- object$dim_theta[1] + object$dim_theta_z + 1

coln <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")

  if(ddf > 0){
	newF <- update.formula(as.formula(object$call[['fixed']]), as.formula(. ~ 1))
	FITNULL <- update(object, fixed = as.formula(newF), evaluate = TRUE)
    LR <- -2*(logLik(FITNULL) - object$logLik);
    if(any(LR < 0)) {LR[LR < 0] <- 0; warning("Negative LR test (set to zero).")};
    LRp <- 1 - pchisq(LR, df = ddf);
    names(LR) <- tau;
    attr(LR, "pvalue") <- LRp;
    attr(LR, "df") <- ddf;
    object$LRtest <- LR;
    object$null_model <- FITNULL
  }

  if(method == "boot"){
    B <- boot.lqmm(object, ...)
	R <- attr(B, "R")
    if(nq == 1){
		Cov <- cov(as.matrix(B))
		stds <- sqrt(diag(Cov))
		tP <- 2*pt(-abs(est/stds), R - 1)
		lower <- est + qt(alpha/2, R - 1)*stds
		upper <- est + qt(1 - alpha/2, R - 1)*stds
		ans <- cbind(est, stds, lower, upper, tP)
		colnames(ans) <- coln
		ans <- ans[c(nn),,drop=FALSE]
    } else {
		Cov <- apply(B, 3, function(x) cov(as.matrix(x)))
		stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
		tP <- 2*pt(-abs(est/stds), R - 1)
		lower <- est + qt(alpha/2, R - 1)*stds
		upper <- est + qt(1 - alpha/2, R - 1)*stds
		ans <- vector("list", nq)
		names(ans) <- tau
		Cov.array <- array(NA, dim = c(object$df,object$df,nq))
		  for(i in 1:nq){
			ans[[i]] <- cbind(est[,i], stds[,i], lower[,i], upper[,i], tP[,i]);
			rownames(ans[[i]]) <- rownames(est)
			colnames(ans[[i]]) <- coln;
			ans[[i]] <- ans[[i]][c(nn),,drop=FALSE]
			Cov.array[,,i] <- matrix(Cov[,i], object$df, object$df)
			}
		Cov <- Cov.array
		dimnames(Cov) <- list(rownames(attr(B, "estimated")), rownames(attr(B, "estimated")), format(tau, digits = 4))
		}
	}
	if(method == "nid"){
	
	
	
	}
if(covariance) object$Cov <- Cov
object$tTable <- ans
object$B <- B
class(object) <- "summary.lqmm"
return(object)

}


print.summary.lqmm <- function(x, digits = max(3, getOption("digits") - 3), ...){

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
  } else {
    for (i in 1:nq) {
      cat(paste("tau = ", tau[i], "\n", sep = ""))
	  cat("\n")
      cat("Fixed effects:\n")

      printCoefmat(x$tTable[[i]], signif.stars = TRUE, P.values = TRUE)
      cat("\n")
    }
  }



LRp <- attr(x$LRtest, "pvalue")
cat("Null model (likelihood ratio):\n")
print.default(paste(format(x$LRtest, digits = digits), " (p = ", format(LRp, digits = digits), ")", sep =""), quote = FALSE)
cat("AIC:\n")
  print.default(
    paste(format(x$aic, digits = digits), " (df = ", attr(x$logLik, "df"), ")", sep =""),
    quote = FALSE
  )
}

boot.lqmm <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE){

if(startQR) warning("Standard errors may be underestimated when 'startQR = TRUE'")

set.seed(seed)
tau <- object$tau
nq <- length(tau)

ngroups <- object$ngroups
group_all <- object$group
group_unique <- unique(group_all)
obsS <- replicate(R, sample(group_unique, replace = TRUE))
dim_theta_z <- object$dim_theta_z

npars <- object$dim_theta[1] + dim_theta_z + 1
dimn <- c(object$nn, paste("reStruct", 1:dim_theta_z, sep=""), "scale")

control <- object$control
control$verbose <- FALSE

FIT_ARGS <- list(x = as.matrix(object$mmf), y = object$y, z = as.matrix(object$mmr), cov_name = object$cov_name,
            V = object$QUAD$nodes, W = object$QUAD$weights, group = group_all, control = control)

if(nq == 1){

  FIT_ARGS$theta_0 <- object$theta;
  FIT_ARGS$sigma_0 <- object$scale;
  FIT_ARGS$tau <- object$tau;

  bootmat <- matrix(NA, R, npars);
  colnames(bootmat) <- dimn

  for(i in 1:R){
    group_freq <- table(obsS[,i]);
    sel_unique <- group_unique%in%names(group_freq);
    w <- rep(0, ngroups); w[sel_unique] <- group_freq;
    FIT_ARGS$weights <- w;
	
	if(!startQR){
		lmfit <- lm.wfit(x = as.matrix(object$mmf), y = object$y, w = rep(w, table(group_all)))
		theta_x <- lmfit$coefficients
		theta_z <- if (object$type == "normal")
			rep(1, dim_theta_z)
				else rep(invvarAL(1, 0.5), dim_theta_z)
		FIT_ARGS$theta_0 <- c(theta_x, theta_z)
		FIT_ARGS$sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
	}

    test <- "try-error";

    while(test=="try-error"){
		  if(control$method == "gs") fit <- try(do.call(lqmm.fit.gs, FIT_ARGS));
		  if(control$method == "df") fit <- try(do.call(lqmm.fit.df, FIT_ARGS));
      test <- class("fit")
      }
    bootmat[i,] <- c(fit$theta , fit$scale)
  }
} else {

  bootmat <- array(NA, dim = c(R, npars, nq), dimnames = list(NULL, dimn, paste("tau = ", format(tau, digits = 4), sep ="")));

  for(i in 1:R){
    group_freq <- table(obsS[,i]);
    sel_unique <- group_unique%in%names(group_freq);
    w <- rep(0, ngroups); w[sel_unique] <- group_freq;
    FIT_ARGS$weights <- w;
    for (j in 1:nq){
		
		if(startQR){
			FIT_ARGS$theta_0 <- object[[j]]$theta;  
			FIT_ARGS$sigma_0 <- object[[j]]$scale
		} else {
			lmfit <- lm.wfit(x = as.matrix(object$mmf), y = object$y, w = rep(w, table(group_all)))
			theta_x <- lmfit$coefficients
			theta_z <- if(object$type == "normal")
				rep(1, dim_theta_z) else rep(invvarAL(1, 0.5), dim_theta_z)
			FIT_ARGS$theta_0 <- c(theta_x, theta_z)
			FIT_ARGS$sigma_0 <- invvarAL(mean(lmfit$residuals^2), 0.5)
		}

		FIT_ARGS$tau <- object$tau[j];

		if(control$method == "gs") fit <- try(do.call(lqmm.fit.gs, FIT_ARGS));
		if(control$method == "df") fit <- try(do.call(lqmm.fit.df, FIT_ARGS));
		if(class(fit)!="try-error") bootmat[i,,j] <- c(fit$theta , fit$scale)
    }
  }
}

class(bootmat) <- "boot.lqmm"
attr(bootmat, "tau") <- tau
attr(bootmat, "estimated") <- extractAll(object)
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "nn") <- object$nn
attr(bootmat, "npars") <- npars
attr(bootmat, "indices") <- obsS
attr(bootmat, "dim_theta") <- object$dim_theta
attr(bootmat, "dim_theta_z") <- object$dim_theta_z

return(bootmat)

}


extractBoot.boot.lqmm <- function(object, which = "fixed"){

tau <- attr(object, "tau")
nq <- length(tau)
nn <- attr(object, "nn")
dim_theta <- attr(object, "dim_theta")
dim_theta_z <- attr(object, "dim_theta_z")

ind.f <- 1:dim_theta[1]
ind.r <- (dim_theta[1] + 1):(dim_theta[1] + dim_theta_z)
ind.s <- dim_theta[1] + dim_theta_z + 1

  if(which == "fixed"){
    if(nq == 1){
      ans <- as.matrix(object[,c(ind.f,ind.s)])
    } else {
      ans <- object[,c(ind.f,ind.s),]
    }
  }
  
  if(which == "random"){
    if(nq == 1){
      ans <- as.matrix(object[,ind.r])
    } else {
      ans <- object[,ind.r,]
    }
  }
  
return(ans)


}


summary.boot.lqmm <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), ...){

est <- attr(object, "estimated")
R <- attr(object, "R")
tau <- attr(object, "tau")
nq <- length(tau)
nn <- attr(object, "nn")
npars <- attr(object, "npars")
coln <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")
		
if(nq == 1){
  Cov <- cov(as.matrix(object))
  stds <- sqrt(diag(Cov))
  tP <- 2*pt(-abs(est/stds), R - 1)
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  ans <- cbind(est, stds, lower, upper, tP)
  colnames(ans) <- coln
  ans <- ans[c(nn,"scale"),]
  cat(paste("Quantile", tau, "\n"))
  printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
}
else {
  Cov <- apply(object, 3, function(x) cov(as.matrix(x)))
  stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = npars))
  tP <- 2*pt(-abs(est/stds), R - 1)
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  for(i in 1:nq){
    ans <- cbind(est[,i], stds[,i], lower[,i], upper[,i], tP[,i]);
    rownames(ans) <- rownames(est)
    colnames(ans) <- coln;
    ans <- ans[c(nn,"scale"),]
    cat(paste("tau = ", tau[i], "\n", sep =""))
    printCoefmat(ans, signif.stars = TRUE, P.values = TRUE)
    cat("\n")
  }

}

}


extractAll <- function(object){

nq <- length(object$tau)

dim_theta_z <- object$dim_theta_z
dimn <- c(object$nn, paste("reStruct", 1:dim_theta_z, sep=""), "scale")
  
  if(nq == 1){
  ans <- c(object$theta, object$scale);
  names(ans) <- dimn
  } else {
  val <- NULL
  for(i in 1:nq) val <- c(val, object[[i]]$scale)
  ans <- rbind(object$theta_x, object$theta_z, val);
  rownames(ans) <- dimn
  }
return(ans)

}

