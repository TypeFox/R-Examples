# Copyright (C) 2014 Stephan Dlugosz, Limin Peng and Ruosha Li
# distributed under the terms of the GNU public license
crrQR <-
# function for quantile regression modeling of subdistribution functions
# arguments:
#  ftime = vector of failure/censoring times (log-transformed)
#  fstatus = vector with a unique code for each failure type and a
#     separate code for censored observations
#  X = (nobs x ncovs) matrix of fixed covariates
#  failcode = code of fstatus that denotes the failure type of interest
#  cencode = code of fstatus that denotes censored observations (default: 0)
#  tau.range = subset of (0,1) denoting the range of qwuantiles to be computed
#  tau.step = stepsize to produce a grid for covering the tau.range
#  subset = logical vector length(ftime) indicating which cases to include
#  na.action = function defining action to take for cases that have NA for
#     any of ftime, fstatus, X, or subset.
#  variance = should coefficient estimator variances be computed?
#

function(ftime, fstatus, X,
	failcode = 1, cencode = 0, tau.range = c(0.01,0.99), tau.step = 0.01,
	subset, na.action = na.omit, rq.method = "br", variance = TRUE) {

	tau.range <- sort(tau.range)
	if (length(tau.range) != 2 || tau.range[1] < 0 || tau.range[2] > 1) stop('invalid range for tau')
	if (!missing(subset) && length(subset) != length(ftime)) stop('dimensions of subset and ftime do not match')
	if (nrow(X) != length(ftime) | length(ftime) != length(fstatus)) stop('dimensions of ftime, fstatus and X do not match')
	n <- length(ftime)

	check <- data.frame(ftime = ftime, fstatus = fstatus)
	if (!missing(subset)) check <- cbind(check, subset)
	check <- cbind(check, X)

	check <- na.action(check)

	ftime <- check$ftime
	fstatus <- check$fstatus
	if (!missing(subset)) {
		subset <- check$subset
		X <- check[,-(1:3), drop = FALSE]
	} else {
		X <- check[,-(1:2), drop = FALSE]
	}

	rm(check)

	if (!missing(subset)) {
		ftime <- ftime[subset]
		fstatus <- fstatus[subset]
		X <- X[subset,, drop = FALSE]
	}
	X <- as.matrix(X)

	ret <- crrQR.int(ftime, fstatus, X, tau.L = tau.range[1], tau.U = tau.range[2], tau.step = tau.step,
	                 outcome = failcode, cencode = cencode, rq.method = rq.method, orig.num = n, variance = variance)

	return(ret)
}


crrQR.int <- function(ftime, fstatus, X, tau.L, tau.U, tau.step,
                      outcome = 1, cencode = 0, noconst = FALSE, variance = TRUE, offset=0, max.lmt = 10^8,
                      rq.method = "br", orig.num = nrow(X), retry = TRUE, ...) {

	eddcmp <- function(M) {
		ev <- eigen(M)
		return(list(evectors = ev$vectors, evalues = ev$values, im.evalues = NULL))
	}

	if (!noconst) {
		X <- cbind(rep.int(1, dim(X)[1]), X)
		dimnames(X)[[2]][1] <- "const"
	}

	o <- order(ftime, decreasing = TRUE)

	ftime <- ftime[o]
	fstatus <- fstatus[o]
	X <- X[o,, drop = FALSE]

	cvt.length <- ncol(X)  ### number of covariates

	tau.seq <- seq(tau.L, tau.U, tau.step)
	L.tau.seq <- length(tau.seq)

	num <- nrow(X)

  FT0 <- as.numeric(fstatus == cencode)
	cens <- (sum(FT0) > 0)
  FT1 <- as.numeric(fstatus == outcome)

  smallest <- .Machine$double.eps ^ 0.5

  csurv.x <- 1.0
  if (cens) {
		#cat(' with censoring\n')
		csurv.fit <- survfit(Surv(ftime, FT0) ~ 1)
		#L.tpt <- length(csurv.fit$time)

		csurv.indx <- unlist(lapply(ftime,
		                            function(et) length(which(et > csurv.fit$time | (abs(et - csurv.fit$time) < smallest)))))
		csurv.x <- csurv.fit$surv[csurv.indx] * FT1 + (1 - FT1)

		rm(csurv.indx)
		gc()
  }

	FT1.csurv.x <- FT1 / csurv.x

  pseudo.resp <- c((ftime - offset) * FT1.csurv.x, max.lmt, max.lmt)
  n.cvt.1 <- crossprod(-X, FT1.csurv.x)

  est.beta.seq <- NULL
	if (variance) {
  	est.var.seq <- NULL
		inf.est.func <- list()
		ind.conv.var <- rep(1, cvt.length)
		est.var <- rep(0, cvt.length)
		curr.inf <- array(rep(0, num*(cvt.length)), c(num, cvt.length))
		compare.vec <- rank(-ftime, ties.method = "max")
	}
  L.bsq <- 1
  ind.conv <- 1

	X.FT1.csurv.x <- X * FT1.csurv.x

  #cat('computing')
	goon <- TRUE
  while (goon) {
  	est.beta <- rep(0, cvt.length)
  	curr.tau <- tau.seq[L.bsq]
  	n.cvt.2 <- apply(2 * X * curr.tau, 2, sum)
		pseudo.cvt <- rbind(X.FT1.csurv.x, as.vector(n.cvt.1), as.vector(n.cvt.2))
  	est.beta.obj <- try(switch(rq.method, fn = rq.fit.fnb(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				fnb = rq.fit.fnb(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				fnc = rq.fit.fnc(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				pfn = rq.fit.pfn(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				br = rq.fit.br(pseudo.cvt, pseudo.resp, tau = curr.tau, ci = FALSE, ...),
				lasso = rq.fit.lasso(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				scad = rq.fit.scad(pseudo.cvt, pseudo.resp, tau = curr.tau, ...),
				{
					what <- paste("rq.fit.", rq.method, sep = "")
					if (exists(what, mode = "function")) (get(what, mode = "function"))(pseudo.cvt,
								pseudo.resp, ...) else stop(paste("unimplemented method:",
										rq.method))
				}))
		if (inherits(est.beta.obj, "try-error")) {
			if (retry) {
				cat('retry with method br ')
				est.beta.obj <- rq.fit.br(pseudo.cvt, pseudo.resp, tau = curr.tau, ci = FALSE, ...)
			} else {
				stop('estimation not successful, try different method')
			}
		}

		est.beta <- est.beta.obj$coef
	  ind.conv <- as.numeric(abs(est.beta.obj$residuals[num + 1L]) > 1) * as.numeric(abs(est.beta.obj$residuals[num + 2L]) > 1)

    if (ind.conv == 1) {
    	est.beta.seq <- rbind(est.beta.seq, est.beta)
			if (variance) {
				helper <- as.numeric((ftime) < X %*% as.vector(est.beta)) * FT1.csurv.x

      	tmp.1 <- ((helper - curr.tau) * X)
      	var.matrix.cmp.1 <- crossprod(tmp.1, tmp.1) / num

      	ZI.mat <- X * helper

	  		tmp.2.1 <- .Call(comptmp21, num, ZI.mat, compare.vec)

      	tmp.2 <- FT0*(tmp.2.1/compare.vec)

      	var.matrix.cmp.2 <- crossprod(tmp.2, tmp.2) / num
      	var.matrix <- var.matrix.cmp.1 - var.matrix.cmp.2
	  		var.dcmp <- eddcmp(var.matrix)
	  		sigma.sqrt <- (var.dcmp$evector %*% diag(sqrt(var.dcmp$evalues),
	  			                                       length(var.dcmp$evalues),
	  			                                       length(var.dcmp$evalues)) %*% solve(var.dcmp$evector))  #eigen decomposition

				est.beta.var <- NULL
    		ind.conv.var <- NULL
    		for (k in 1:cvt.length) {
    			pseudo.cvt.var <- rbind(X.FT1.csurv.x, as.vector(n.cvt.1), as.vector(n.cvt.2) + 2 * sigma.sqrt[,k] * sqrt(num))
    			est.beta.var.obj <- try(switch(rq.method, fn = rq.fit.fnb(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
		  		fnb = rq.fit.fnb(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
					fnc = rq.fit.fnc(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
					pfn = rq.fit.pfn(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
					br = rq.fit.br(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ci = FALSE, ...),
					lasso = rq.fit.lasso(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
					scad = rq.fit.scad(pseudo.cvt.var, pseudo.resp, tau = curr.tau, ...),
					  {
							what <- paste("rq.fit.", rq.method, sep = "")
							if (exists(what, mode = "function")) (get(what, mode = "function"))(pseudo.cvt,
										pseudo.resp, ...) else stop(paste("unimplemented method:",
												rq.method))
						}))
					if (inherits(est.beta.var.obj, "try-error")) {
						if (retry) {
							cat('retry with method br ')
							est.beta.var.obj <- rq.fit.br(pseudo.cvt, pseudo.resp, tau = curr.tau, ci = FALSE, ...)
						} else {
							stop('estimation not successful, try different method')
						}
					}

        	ind.conv.var <- c(ind.conv.var, as.numeric(abs(est.beta.var.obj$residuals[num + 1]) > 1) *
        	                    as.numeric(abs(est.beta.var.obj$residuals[num + 2]) > 1))
        	est.beta.var <- cbind(est.beta.var, est.beta.var.obj$coef - est.beta)
      	}

      	if (sum(ind.conv.var) == cvt.length) {
        	est.var <- diag(tcrossprod(est.beta.var, est.beta.var))
        	est.inv.deriv.matrx <- est.beta.var %*% solve(sigma.sqrt)*sqrt(num)
        	est.inv.deriv.matrx <- (est.inv.deriv.matrx + t(est.inv.deriv.matrx))/2
        	curr.inf <- tcrossprod(tmp.1 - tmp.2, est.inv.deriv.matrx)
      	} else {
				  goon <- FALSE
				}
				inf.est.func[[L.bsq]] <- curr.inf
      	est.var.seq <- rbind(est.var.seq, est.var)
			}
      L.bsq <- L.bsq + 1
		} else {
			goon <- FALSE
		}

		goon <- (L.bsq <= L.tau.seq & goon)

		cat('.')
  }
  cat('done.\n')

	if (L.bsq < L.tau.seq) {
		cat('\n')
		cat('stopping at tau=')
		cat(tau.seq[L.bsq])
		cat('\n')
		L.tau.seq <- L.bsq - 1
		if (L.bsq <= 1) return(list(tau.seq = NULL, L.bsq = 0))
		tau.seq <- tau.seq[1:L.tau.seq]
		tau.U <- tau.seq[L.tau.seq]
	}

	dimnames(est.beta.seq)[[1]] <- tau.seq
	if (variance) dimnames(est.var.seq)[[1]] <- tau.seq

  # return results
	ret <- list(call = match.call(),
				beta.seq = est.beta.seq,
				tau.seq = tau.seq,
				cvt.length = cvt.length,
				n = num,
				n.missing = orig.num - num)
	if (variance) {
		ret$var.seq = est.var.seq
		ret$inf.func = inf.est.func
	}
	class(ret) <- c("crrQR")

  return(ret)
}

summary.crrQR <-
# compute confidence intervals and tests
# as well as some summary statistics
# for a crrQR object
function(object, conf.int = 0.95, digits = max(options()$digits - 5, 2), ...) {
	if (is.null(object$var.seq)) stop('please re-run crrQR with variance=TRUE')
	L.tau.seq <- length(object$tau.seq)
	tau.L <- object$tau.seq[1]
	tau.U <- object$tau.seq[L.tau.seq]
	tau.step <- {tau.U - tau.L}/L.tau.seq

	tau.w1 <- (as.numeric(object$tau.seq < (tau.L + tau.U) / 2))[1:(L.tau.seq - 1)] ### weight for constancy test

	inf.ave <- NULL
	est.ave <- NULL
	inf.cnst.test.w1 <- NULL
	cnst.test.w1 <- NULL
	for (k in 1:(object$cvt.length)) {
		infs <- object$inf.func
		infs[[length(infs)]] <- NULL
		infs <- do.call("cbind", lapply(infs, function(x) x[,k,drop = FALSE]))
		inf.ave <- cbind(inf.ave, apply(infs*tau.step, 1, sum) / (tau.U - tau.L))
		est.ave <- c(est.ave, sum(object$beta.seq[1:(L.tau.seq - 1), k, drop = FALSE]*tau.step) / (tau.U - tau.L))
		diff <- infs - matrix(rep.int(inf.ave[,k],ncol(infs)), ncol = (ncol(infs)), byrow = FALSE)
		inf.cnst.test.w1 <- cbind(inf.cnst.test.w1, apply(t(t(diff)*tau.w1)*tau.step, 1, sum) / sum(tau.w1 * tau.step))
		cnst.test.w1 <- c(cnst.test.w1,
		                  sum((object$beta.seq[1:(L.tau.seq - 1), k, drop = FALSE] - est.ave[k]) * tau.w1*tau.step) /
		                    sum(tau.w1*tau.step))
	}
	names(est.ave) <- dimnames(object$beta.seq)[[2]]
	var.ave <- apply(inf.ave^2 / object$n, 2, mean)
	names(var.ave) <- names(est.ave)
	var.cnst.test.w1 <- apply(inf.cnst.test.w1^2 / object$n, 2, mean)
	names(cnst.test.w1) <- names(est.ave)
	names(var.cnst.test.w1) <- names(est.ave)

	ret <- list(call = object$call,
			        ave.eff = est.ave,
			        var.ave.eff = var.ave,
				      p.signf.test = signif(1 - pchisq(est.ave^2/var.ave, 1), digits),
				      cnst.test = cnst.test.w1,
				      var.cnst.test = var.cnst.test.w1,
				      p.cnst.test = signif(1 - pchisq(cnst.test.w1^2/var.cnst.test.w1, 1), digits),
				      n = object$n,
				      n.missing = object$n.missing)

	class(ret) <- "summary.crrQR"

	return(ret)
}

print.summary.crrQR <-
# prints the summary object to a crrQR model
function(x, digits=max(options()$digits - 4, 3), ...) {
	cat("Competing Risks Quantile Regression\n\n")
	if (!is.null(x$call))
	{ cat("Call:\n")
		dput(x$call)
		cat("\n")
	}
	savedig <- options(digits = digits)
	on.exit(options(savedig))

	out.bcr <- data.frame(x$ave.eff, sqrt(x$var.ave.eff), x$p.signf.test, x$cnst.test, x$var.cnst.test, x$p.cnst.test)
	names(out.bcr) <- c("ave.eff", "se.ave.eff", "p.signf.test", "cnst.test", "var.cnst.test", "p.cnst.test")

	print(out.bcr)
	cat("\n")
	cat("Num. cases =", x$n)
	if (x$n.missing > 0) cat(" (", x$n.missing, " cases omitted due to missing values)", sep = "")
	cat("\n")

	invisible()
}

predict.crrQR <-
# for a crrQR object x, estimates subdistributions at covariate
# combinations given by rows of X. The terms in X
# must correspond exactly to the corresponding call to crrQR.
# add rearrangement=TRUE for using the rearrangement method to ensure increasing CiC predictions
# (cf. Chernozhukov and Fernandez-Val and Galichon, Quantile and probability curves without crossing, Econometrica 78, 1093--1125)
function(object, x, ...) {
	rearrangement <- FALSE
	if (!missing(rearrangement)) rearrangement <- rearrangement
	betas <- object$beta.seq
	if (inherits(x, "data.frame")) x <- as.matrix(x)
	if (inherits(x, "matrix")) {
		if (ncol(x) != ncol(betas) - 1) stop('dimensions of X and estimated betas do not match')
		x <- x[,dimnames(betas)[[2]][-1], drop = FALSE]
		x <- cbind(1, x)
	} else {
		if (length(x) != ncol(betas) - 1) stop('dimensions of X and estimated betas do not match')
		x <- x[dimnames(betas)[[2]][-1]]
		x <- matrix(c(1, x), nrow = 1)
	}

	xb <- tcrossprod(betas, x)

	exb <- exp(xb)
	if (rearrangement) exb <- sort(exb)

	lhat <- cbind(exb, object$tau.seq)
	class(lhat) <- 'predict.crr'
	return(lhat)
}

plot.predict.crrQR <-
# plots estimated subdistributions from predict.crrQR
function(x, lty=1:(ncol(x) - 1), color = 1, ylim = c(0,max(x[,-1])), xmin = 0, xmax = max(x[,1]), ...) {
	if (length(lty) < ncol(x) - 1) lty <- rep(lty[1], ncol(x) - 1)
	if (length(color) < ncol(x) - 1) color <- rep(color[1], ncol(x) - 1)
	if (xmax < max(x[,1])) x <- x[x[,1] < xmax,]
	times <- c(xmin,rep(x[,1], rep(2,nrow(x))),xmax)
	plot(c(xmin,xmax), ylim, type = "n", ...)
	for (j in 2:ncol(x)) lines(times, c(0, 0, rep(x[,j], rep(2, nrow(x)))), lty = lty[j - 1], col = color[j - 1])
}

print.crrQR <-
# prints the essential parts of a crrQR Object
function(x, ...) {
	cat('coefficients:\n')
	print(x$beta.seq)
	v <- sqrt(x$var.seq)
	cat('standard errors:\n')
	print(signif(v,4),...)
	v <- 2 * (1 - pnorm(abs(x$beta.seq) / v))
	cat('two-sided p-values:\n')
	print(signif(v, 2),...)
	invisible()
}

plot.crrQR <-
# plots the profile graph of the coefficients
# subset defines the subset of betas to plot
function(x, subset=NULL, main=NULL, ...) {
	if (is.null(subset)) subset <- 1:(x$cvt.length)
	if (is.null(main)) main <- ""

	CI.lb.seq <- x$beta.seq - 1.96 * sqrt(x$var.seq)
	CI.ub.seq <- x$beta.seq + 1.96 * sqrt(x$var.seq)

	if (length(subset) > 1) op <- par(mfrow = c(ceiling((x$cvt.length + 1) / 2), 2))

	#grp.frm.1<-apply(results[[risk]]$CI.lb.seq, 2, min)*0.99
	#grp.frm.2<-apply(results[[risk]]$CI.ub.seq, 2, max)*1.01

	grp.frm.1 <- apply(CI.lb.seq, 2, function(x) { quantile(x, probs = c(0.1)) } )*0.99
	grp.frm.2 <- apply(CI.ub.seq, 2, function(x) { quantile(x, probs = c(0.9)) } )*1.01

	char.cvt <- dimnames(x$beta.seq)[[2]]

	for (k in subset) {
		plot(x$tau.seq, x$beta.seq[, k], ylim = c(grp.frm.1[k],
			grp.frm.2[k]), type = "s", xlab = "tau", ylab = "Coefficient",
	        main = paste(main, char.cvt[k], sep = ":"), ...)
		lines(x$tau.seq, CI.lb.seq[,k], type = "s", lty = 3)
		lines(x$tau.seq, CI.ub.seq[,k], type = "s", lty = 3)
		lines(x$tau.seq, rep(0, length(x$tau.seq)), lty = 2)
		lines(x$tau.seq, rep(x$est.ave[k], length(x$tau.seq)), col = 2)
	}
	if (length(subset) > 1) par(op)
}
