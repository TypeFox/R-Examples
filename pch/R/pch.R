#' @importFrom stats model.matrix model.response delete.response model.frame terms prcomp
#' @importFrom stats model.weights logLik coef nobs vcov .getXlevels quantile runif approxfun sd
#' @importFrom survival Surv is.Surv
#' @importFrom splines ns bs


#' @export
pchreg <- function(formula, breaks, data, weights, splinex = NULL){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	mf$formula <- formula
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	if(!survival::is.Surv(zyd <- model.response(mf)))
		{stop("the model response must be created with Surv()")}
	if((n <- nrow(zyd)) == 0){stop("zero non-NA cases")}
	type <- attributes(zyd)$type
	if(type == "right"){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]}
	else if(type == "counting"){z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]}
	else{stop("only type = 'right' and type = 'counting' are supported")}
	if(!(any(d == 1))){stop("all observation are censored")}
	
	x <- model.matrix(mt,mf); ax <- NULL
	if(!is.null(splinex)){
	  x <- build.splinex(x, splinex$method, splinex$df, splinex$degree, splinex$v)
	  ax <- attributes(x)
	 }
	w <- model.weights(mf)
	if(is.null(w)){w <- rep.int(1,n)}

	###

	fit <- pch.fit(z,y,d,x,w,breaks)

	###
  
	attr(mf, "rangex") <- fit$rangex
	attr(mf, "splinex") <- ax
	attr(mf, "u") <- fit$u
	attr(mf, "y") <- fit$y
	attr(mf, "type") <- type
	attr(mf, "n") <- length(y)
	attr(mf, "n.events") <- sum(d)

	beta <- fit$beta
	r <- dim(beta)
	colnames(beta) <- 1:r[2]
	vn <- paste(rep.int(rownames(beta), r[2]), rep(1:r[2], each = r[1]), sep = "_")
	dimnames(fit$covar) <- list(vn, vn)

	Hy <- predF.pch(fit, x, y)
	Hz <- (if(type == "counting") predF.pch(fit,x,z)[,"Haz"] else 0)
	logLik <- sum(d*log(Hy[,"haz"]), na.rm = TRUE) - sum(Hy[,"Haz"] - Hz)
	# note: 'haz' can be 0, I set 0*log(0) = 0.
	attr(logLik, "df") <- sum(fit$beta != 0)
	
	fit <- list(call = cl, beta = beta, breaks = fit$breaks, 
		covar = fit$covar, logLik = logLik, 
		lambda = fit$lambda, Lambda = fit$Lambda, mf = mf, x = x, conv.status = fit$conv.status)
	class(fit) <- "pch"
	fit
}


#' @export
predict.pch <- function(object, type = c("distr", "quantile", "sim"), 
	newdata, p, sim.method = c("quantile", "sample"), ...){
  
	if(is.na(match(type <- type[1], c("distr", "d", "quantile", "q", "sim", "s"))))
		{stop("invalid 'type'")}
	type <- strsplit(type, split = "")[[1]][1]
	if(type == "s"){
		if(is.na(match(method <- sim.method[1], c("quantile", "q", "sample", "s"))))
		{stop("invalid 'method'")}
		method <- strsplit(method, split = "")[[1]][1]
	}
	if(type == "q"){
	  if(missing(p)){stop("please indicate 'p'")}
	  if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}	
	}
	
	mf <- object$mf
	mt <- terms(mf)
	xlev <- .getXlevels(mt, mf)
	fittype <- attr(mf, "type")
	splinex <- attr(mf, "splinex")
	
	obj <- list(beta = object$beta, breaks = object$breaks, y = attr(mf, "y"), 
		lambda = object$lambda, Lambda = object$Lambda, u = attr(mf, "u"), 
		rangex = attr(mf, "rangex"))

	if((nodata <- missing(newdata))){
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	}
	else{
		if(type == "d"){
			yn <- as.character(if(fittype == "counting") mt[[2]][[3]] else mt[[2]][[2]])
			if(is.na(ind <- match(yn, colnames(newdata))))
			{stop("for 'type = distr', 'newdata' must contain the y-variable")}
			if(fittype == "right" && length(mt[[2]]) == 3)
			  {newdata[,as.character(mt[[2]][[3]])] <- 1}
			if(fittype == "counting"){
				newdata[,as.character(mt[[2]][[4]])] <- 1
				newdata[,as.character(mt[[2]][[2]])] <- -Inf
				newdata[,yn] <- pmax(newdata[,yn], obj$breaks[1] - 1)
			}
		}
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
		{stop("'newdata' must contain all x-variables")}
		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){stop("zero non-missing values in the supplied newdata")}

		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
    
		x <- model.matrix(mt, mf)
		if(!is.null(splinex)){
		  x <- predict.splinex(x,splinex)
		}
	}

	n <- length(miss) + length(nomiss)
	if(type == "d"){
		if(!nodata){y <- cbind(model.response(mf))[,1 + (fittype == "counting")]}
		out <- matrix(, n, 4)
		pred <- (if(nodata) predF.pch(obj) else predF.pch(obj,x,y))
		out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
  	else if(type == "q"){
		out <- matrix(, n, length(p))
		pred <- (if(nodata) predQ.pch(obj, p = p) else predQ.pch(obj,x,p))
        	out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
	else{
		tsim <- rep.int(NA,n)
		tsim[nomiss] <- (if(nodata) sim.pch(obj, method = method) else sim.pch(obj,x,method))
		return(tsim)
	}
}

# print and summary method
#' @export
print.pch <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
summary.pch <- function(object, ...){
  out <- list(call = object$call, 
    n = attr(object$mf, "n"), n.events = attr(object$mf, "n.events"), 
		n.free.par = attr(object$logLik, "df"), logLik = object$logLik, conv.status = object$conv.status)
  class(out) <- "summary.pch"
  out
}

#' @export
print.summary.pch <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of obs: ", paste(deparse(round(x$n)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of events: ", paste(deparse(round(x$n.events)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of free parameters: ", paste(deparse(round(x$n.free.par)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", paste(deparse(round(as.numeric(x$logLik),1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("convergence status: ", paste(deparse(round(x$conv.status)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
logLik.pch <- function(object, ...){object$logLik}
#' @export
nobs.pch <- function(object, ...){nrow(object$mf)}
#' @export
coef.pch <- function(object, ...){object$beta}
#' @export
vcov.pch <- function(object, ...){object$covar}

CumSum <- function(x){
	if((q <- ncol(x)) > 1){
		for(j in 2:q){x[,j] <- x[,j] + x[,j-1]}
	}
	x
}


makebreaks <- function(y,d,x, breaks){
  
	n <- length(y)
	n1 <- sum(d)
	q <- ncol(x)
	r <- range(y)
	y <- y[d == 1]
	
	if(missing(breaks)){breaks <- max(5, min(50, ceiling(n1/q/5)))}

	if(length(breaks) > 1){
		breaks <- sort(unique(breaks))
		k <- length(breaks) - 1
		if(r[1] < breaks[1] | r[2] > breaks[k + 1])
		  {stop("all y values must be within the breaks")}
		breaks <- breaks[1:which(breaks >= r[2])[1]]
		k <- length(breaks) - 1
	}
	else{
		k <- breaks
		breaks <- quantile(y, (0:k)/k)
		breaks[1] <- r[1]; breaks[k + 1] <- r[2]
		a <- duplicated(breaks)
		if(any(a)){
			for(j in which(a)){
				if(a[j - 1]){breaks[j] <- NA}
				else{
					h <- abs(y - breaks[j])
					h <- min(h[h > 0])
					h <- min(h/2, (r[2] - r[1])/n/100)
					breaks[j-1] <- breaks[j-1] - h
					breaks[j] <- breaks[j] + h
				}
			}
			breaks <- breaks[!is.na(breaks)]
			k <- length(breaks)
			h <- c(Inf, breaks[2:k] - breaks[1:(k - 1)])
			breaks[h <= 0] <- NA
			breaks <- breaks[!is.na(breaks)]
			k <- length(breaks) - 1
		}
	}
	breaks[1] <- breaks[1] - (breaks[2] - breaks[1])/n1
	names(breaks) <- NULL
	list(breaks = breaks, k = k)	
}




pois.newton <- function(start, f, tol = 1e-5, maxit = 200, safeit = 0, ...){

	f0 <- f(start, ..., deriv = 2)
	g <- attr(f0, "gradient")
	h <- attr(f0, "hessian")
	conv <- FALSE
	eps <- 1
	alg <- "gs"

	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####
		
	  	if(i > safeit){
		  	H1 <- try(chol(h), silent = TRUE)
		  	if(class(H1) != "try-error"){
		  		if(alg == "gs"){alg <- "nr"; eps <- 1}
		  		delta <- chol2inv(H1)%*%g
		  	}
		  	else{
		  		if(alg == "nr"){alg <- "gs"; eps <- 1}
		  		delta <- g
		  	}
	  	}
		else{delta <- g}
	  
		####

		f1 <- Inf
		while(f1 > f0){
			new.start <- start - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			f1 <- try(f(new.start, ..., deriv = 0), silent = TRUE)
			eps <- eps*0.5
			if(class(f1) == "try-error" || is.na(f1)){f1 <- Inf}
		}

		if(conv | f0 - f1 < tol){break}
		f1 <- f(new.start, ..., deriv = 2)
		g1 <- attr(f1, "gradient")
		h1 <- attr(f1, "hessian")

		start <- new.start; f0 <- f1; g <- g1; h <- h1
		eps <- min(eps*10,1)
	}

	list(estimate = start, n.it = i, minimum = as.numeric(f0), gradient = g, hessian = h)
}




pois.loglik <- function(beta,d,x,w,off, zeror, deriv = 0){

	log.lambda <- x%*%cbind(beta) - zeror*1e+10
	lambda <- exp(log.lambda)
	a <- lambda*off
	l <- sum(w*(d*log.lambda - a))
	if(deriv == 0){return(-l)}

	s <- x*c(w*(d - a))
	if(deriv == 1){return(s)}

	s <- colSums(s)
	h <- -t(x*(c(w*a)))%*%x
	out <- -l
	attr(out, "gradient") <- -s
	attr(out, "hessian") <- -h 
	out
}



myapply <- function(x,fun){
  out <- NULL
  for(j in 1:ncol(x)){out <- rbind(out, fun(x[,j]))}
  if(ncol(out) == 1){out <- c(out)}
  out
}

poisfit <- function(d,x,w,off){

	cn <- colnames(x)
	if(!any(d != 0)){
		sx <- myapply(x,sd)
		const <- (if(nrow(x) > 1){sx == 0} else {colnames(x) == "(Intercept)"})
		if(int <- any(const)){
			beta <- rep(NA, ncol(x))
			beta[const] <- -Inf
			names(beta) <- cn
			return(list(beta = beta, vcov = 0))
		}
		else{stop("zero-risk can not be fitted unless an intercept is included")}
	}

	# Handling zero-risk regions

	r <- myapply(x[d == 1,, drop = FALSE], range)
	delta <- r[,2] - r[,1]
	r[,1] <- r[,1] - 0.1*delta
	r[,2] <- r[,2] + 0.1*delta
	zeror <- rep.int(FALSE, length(d))
	for(j in 1:ncol(x)){
		out.l <- (x[,j] < r[j,1])
		out.r <- (x[,j] > r[j,2])
		if(mean(out.l) > mean(out.r)){outx <- out.l; r[j,2] <- Inf}
		else{outx <- out.r; r[j,1] <- -Inf}
		zeror <- (zeror | outx)
	}
	x <- x*(!zeror)

  
	# Scaling
	
	xx <- qr(x)
	sel <- xx$pivot[1:xx$rank]
	x <- x[,sel, drop = FALSE]
	q <- ncol(x)
	X <- list(x = x, off = off)

	mx <- colMeans(x)
	sx <- myapply(x,sd)
	M <- 10/max(off)

	const <- (if(nrow(x) > 1){sx == 0} else {colnames(x) == "(Intercept)"})
	if(int <- any(const)){
		const <- which(const)
		mx[const] <- 0
		sx[const] <- x[1,const]
		off <- off*M
	}
	else{
		const <- integer(0)
		mx <- rep.int(0, q)
		M <- 1
	}
	vars <- which(sx > 0)
	x <- scale(x, center = mx, scale = sx)

	# Fit
	
	conv <- TRUE
	beta0 <- rep.int(0,q)
	if(int){beta0[const] <- max(-10, log(mean(d[!zeror])))}
	safeit <- 0; fit.ok <- FALSE; count <- 0
	while(!fit.ok){
	  fit <- pois.newton(beta0, pois.loglik, tol = 1e-5, maxit = 10*(1 + q), safeit = safeit,
	                     d = d, x = x, w = w, off = off, zeror = zeror)

	  fit.ok <- all(abs(fit$gradient) < 1)
	  count <- count + 1; safeit <- safeit + 2
	  if(count == 20){conv <- FALSE; break}
	}

	beta <- fit$estimate
	beta[vars] <- beta[vars]/sx[vars]
	beta[const] <- beta[const] - sum(beta[vars]*mx[vars]) + log(M)
	loglik <- pois.loglik(beta,d,X$x,w,X$off, zeror, deriv = 2)

	Beta <- rep(NA, length(xx$pivot))
	Beta[sel] <- beta
	names(Beta) <- cn

	# covariance matrix

	h <- attr(loglik, "hessian")
	V <- 0*h
	sel <- 1:ncol(h)
	v <- try(chol2inv(chol(h)), silent = TRUE)
	while(class(v) == "try-error"){
	  qrh <- qr(h[sel,sel,drop = FALSE])
	  sel <- sel[qrh$pivot[1:qrh$rank]]
	  v <- try(chol2inv(chol(h[sel,sel, drop = FALSE])), silent = TRUE)
	}
	V[sel,sel] <- v
	list(beta = Beta, vcov = V, r = r, converged = conv)
}






pch.fit <- function(z,y,d,x,w,breaks){ 

	n <- length(y)
	q <- ncol(x)
	if(missing(z)){z <- rep.int(-Inf,n)}
	if(missing(d)){d <- rep.int(1,n)}
	if(missing(w)){w <- rep.int(1,n)}
	Breaks <- makebreaks(y,d,x,breaks)

	k <- Breaks$k
	Breaks <- Breaks$breaks
	z[z < min(y)] <- Breaks[1]

	zcut <- cut(z,Breaks, include.lowest = TRUE)
	ycut <- cut(y,Breaks, include.lowest = TRUE)
	tab.ycut <- table(ycut)
	end.z <- match(zcut, levels(zcut))
	end.y <- match(ycut, levels(ycut))
	h <- Breaks[2:(k + 1)] - Breaks[1:k]

	conv <- TRUE
	A <- cbind(z,end.z,y,end.y,d,w,x)
	beta <- NULL; rangex <- list()
	V <- matrix(0,k*q, k*q)
	for(j in 1:k){

		uy <- (A[,"end.y"] == j)
		uz <- (A[,"end.z"] == j)
		vz <- (A[,"end.z"] < j)

		zj <- yj <- rep.int(h[j],n); dj <- rep.int(0,n)
		zj[uz] <- A[uz, "z"] - Breaks[j]; zj[vz]<- 0
		yj[uy] <- A[uy, "y"] - Breaks[j]
		yz_zj <- yj - zj
		dj[uy] <- A[uy, "d"]; wj <- A[,"w"]
		xj <- A[,7:ncol(A),drop = FALSE]

		r <- (yz_zj != 0)
		dj <- dj[r]; yj_zj <- yz_zj[r]
		xj <- xj[r,, drop = FALSE]; wj <- wj[r]
		modj <- poisfit(dj, xj, wj, yj_zj)

		beta <- cbind(beta, modj$beta)
		rangex[[j]] <- modj$r
		v <- matrix(0,q,q)
		sel <- which(!is.na(modj$beta))
		ind <- ((j - 1)*q + 1):(j*q)
		v[sel,sel] <- modj$vcov
		V[ind,ind] <- v

		n <- nrow(A <- A[!uy,, drop = FALSE])
		conv <- (conv & modj$converged)
	}

	beta[is.na(beta)] <- 0	
	lambda <- cleanlambda(exp(x%*%beta),x,rangex)
	Lambda <- CumSum(t(t(lambda)*h))
	colnames(lambda) <- colnames(Lambda) <- 1:k

	# y and corresponding interval
	
	y <- cbind(y, end.y + 1, sort(y))
	colnames(y) <- c("y","interval","sort.y")
	attr(y, "tab.ycut") <- tab.ycut
	attr(Breaks, "h") <- h; attr(Breaks, "k") <- k
	
	# approxfun for quick prediction of the interval
	
	br <- c(Breaks[1] - 1, Breaks)
	u <- approxfun(br, c(1:k, k + 1, k + 1), rule = 2, method = "constant")  
	
	# check1: convergence
	# check2: there should not be many survival = 1
	
	conv.status <- 0
	if(!conv){conv.status <- 1; warning("the algorithm did not converge", call. = FALSE)}
	else{
	  surv <- exp(-cbind(0,Lambda)[cbind(1:nrow(Lambda), end.y)]) # approx survival
	  if(mean(surv <= 1e-3) > 0.01){
	    conv.status <- 2
	    warning(
	       "numerically zero survivals fitted: the solution may not be well-behaved
	       (or the model is misspecified)", 
	       call. = FALSE)
	  }
	}

	list(beta = beta, lambda = lambda, Lambda = Lambda,
		covar = V, breaks = Breaks, y = y, u = u, rangex = rangex, conv.status = conv.status)
}

cleanlambda <- function(lambda,x,rangex){
  for(h in 1:ncol(lambda)){
    r <- rangex[[h]]
    for(j in 1:ncol(x)){
      e <- which(x[,j] < r[j,1] | x[,j] > r[j,2])
      lambda[e,h] <- 0
    }
  }
  pmin(lambda, 1e+6)
}


predF.pch <- function(obj,x,y){

  Breaks <- obj$breaks
  k <- attr(Breaks, "k")
  h <- attr(Breaks, "h")
  beta <- obj$beta
  u <- obj$u
  rangex <- obj$rangex

  if(missing(x)){
    lambda <- obj$lambda
    Lambda <- obj$Lambda
  }
  else{
    lambda <- cleanlambda(exp(x%*%beta), x, rangex)
    Lambda <- CumSum(t(t(lambda)*h))
  }
  
  if(missing(y)){y <- obj$y[,1]; end.y <- obj$y[,2]}
  else{end.y <- u(y); y <- pmax(y, Breaks[1] - 1)}
  
  n <- length(y)
  t <- y - c(0,Breaks)[end.y]
  ind <- cbind(1:n,end.y)
  lambda <- cbind(0,lambda)[ind]
  Lambda <- cbind(0,0,Lambda)[ind] + lambda*t
  SF <- exp(-Lambda)
  PDF <- lambda*SF
  out <- cbind(haz = lambda, Haz = Lambda, Surv = SF, f = PDF)
  rownames(out) <- NULL
  out
}

predQ.pch <- function(obj,x,p){ # p can be a vector (to simulate)

	Breaks <- obj$breaks
	h <- attr(Breaks, "h")
	k <- attr(Breaks, "k")
	beta <- obj$beta
	rangex <- obj$rangex

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- cleanlambda(exp(x%*%beta), x, rangex)
		Lambda <- CumSum(t(t(lambda)*h))
	}
	
	n <- nrow(lambda)
	if(missing(p)){p <- runif(n); sim <- TRUE; r <- 1}
	else{sim <- FALSE; r <- length(p)}
	pp <- -log(1 - p)

	out <- NULL
	for(j in 1:r){
		tau <- (if(sim) pp else pp[j])
		ind <- .rowSums(Lambda < tau, n, ncol(Lambda)) + 1
		Lambdap <- cbind(0,Lambda)[cbind(1:n, ind)]
		c <- Breaks[ind]
		t <- (tau - Lambdap)/lambda[cbind(1:n, pmin(ind,k))]
		# I use the k-th hazard after the last break
		out <- cbind(out, c + t)
	}
	if(sim){return(c(out))}
	colnames(out) <- paste("p",p, sep = "")
	rownames(out) <- NULL
	out
}

# Note: method = "sample" will only work with "many" breaks.
# On the other hand, it will be more robust to model assumptions,
# and it will correcly reproduce a mass. Moreover, it will work well
# if the data are not compact, while "empty spaces" in the support
# will be filled with method = "quantile".

sim.pch <- function(obj,x,method = c("q","s")){

	if((method <- method[1]) == "q")
		{return(predQ.pch(obj,x))}
	
	Breaks <- obj$breaks
	h <- attr(Breaks, "h")
	k <- attr(Breaks, "k")
	beta <- obj$beta
	rangex <- obj$rangex
	y <- obj$y
	tab.ycut <- attr(y, "tab.ycut")
	y <- y[,"sort.y"]	

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- cleanlambda(exp(x%*%beta), x, rangex)
		Lambda <- CumSum(t(t(lambda)*h))
	}
	
	n <- nrow(lambda)
	p <- -log(1 - runif(n))
	ind <- .rowSums(Lambda < p, n, ncol(Lambda)) + 1 # can be k + 1
	i1 <- 0
	t <- rep.int(NA,n)
	
	for(j in 1:k){
		uj <- which(ind == j)
		i2 <- i1 + tab.ycut[j]
		t[uj] <- sample(y[(i1 + 1):i2], size = length(uj), replace = TRUE)
		i1 <- i2
	}
	if(any(u.out <- (is.na(t)))){ # where ind = k + 1
		u.out <- which(u.out)
		p.out <- p[u.out]
		lambda.out <- lambda[u.out,k]
		Lambda.out <- Lambda[u.out,k]
		t[u.out] <- Breaks[k + 1] + (p.out - Lambda.out)/lambda.out
	}
	t
}


# splinex functions

#' @export
splinex <- function(method = c("ns", "bs"), df = 2, degree = 2, v = 0.98, ...){
  if(is.na(match((method <- method[1]), c("ns", "bs")))){stop("invalid 'method'")}
  if((df <- round(df)) < 1){stop("invalid 'df'")}
  if((degree <- round(degree)) < 1){stop("invalid 'degree'")}
  if(method == "bs" && df < degree){df <- degree; warning(paste("'df' was too small; have used", df))}
  if(v <= 0 | v > 1){stop("0 < v <= 1 is required")}
  list(method = method, df = df, degree = degree, v = v)
}

build.splinex <- function(x, method, df, degree, v){

  m <- (if(method == "ns") function(x, df, deg) splines::ns(x, df = df) 
        else function(x, df, deg) splines::bs(x, df = df, degree = deg))
  x <- cbind("(Intercept)" = 1, as.matrix(x))
  qx <- qr(x)
  sel1 <- qx$pivot[1:qx$rank]
  x <- x[,sel1, drop = FALSE]
  cnx <- colnames(x)
  
  if(ncol(x) == 1){
    attr(x, "onlyint") <- TRUE
    return(x)
  }
  if(df == 1){
    X <- x
    df <- rep.int(1, ncol(X))
    sel2 <- 1:ncol(X)
    knots <- bknots <- NULL
  }
  else{
    u <- myapply(x, function(a){length(unique(a))})
    X <- cbind("(Intercept)" = rep.int(1, nrow(x)))
    dfX <- 1; cnX <- "(Intercept)"
    knots <- bknots <- list()
    for(j in 2:ncol(x)){
      dfj <- dfX[j] <- min(df, u[j] - 1)
      xx <- (if(dfj > 1) m(x[,j], df = dfj, deg = degree) else x[,j])
      X <- cbind(X, xx)
      cnX <- c(cnX, paste(cnx[j],1:dfj, sep = "."))
      knots[[j]] <- attr(xx, "knots")
      bknots[[j]] <- attr(xx, "Boundary.knots")
    }
    colnames(X) <- cnX
    df <- dfX
    qX <- qr(X)
    sel2 <- qX$pivot[1:qX$rank]
    X <- X[,sel2, drop = FALSE]
  }
  if(v == 1){rot <- diag(ncol(X) - 1); center <- scale <- FALSE; sel3 <- 1:(ncol(X) - 1)}
  else{
    XX <- prcomp(X[,-1, drop = FALSE], center = TRUE, scale. = TRUE)
    rot <- XX$rotation; center <- XX$center; scale <- XX$scale
    vx <- cumsum(XX$sdev^2); vx <- vx/vx[length(vx)]
    sel3 <- 1:(which(vx >= v)[1])
    X <- cbind(X[,1,drop = FALSE],XX$x[,sel3])
    colnames(X) <- c("(Intercept)", paste("pc", 1:(ncol(X) - 1), sep = ""))
  }
  
  attr(X, "onlyint") <- FALSE
  attr(X, "method") <- method[1]
  attr(X, "sel1") <- sel1
  attr(X, "df") <- df
  attr(X, "degree") <- degree
  attr(X, "knots") <- knots
  attr(X, "bknots") <- bknots
  attr(X, "sel2") <- sel2
  attr(X, "center") <- center
  attr(X, "scale") <- scale
  attr(X, "rot") <- rot
  attr(X, "sel3") <- sel3
  X
}

predict.splinex <- function(x, a){
  
  if(a$onlyint){return(x)}

  m <- (if(a$method == "ns") function(x, df, deg, kn, bkn) splines::ns(x, df = df, knots = kn, Boundary.knots = bkn) 
        else function(x, df, deg, kn, bkn) splines::bs(x, df = df, degree = deg, knots = kn, Boundary.knots = bkn))
  
  x <- cbind(1,x)
  x <- x[,a$sel1, drop = FALSE]
  X <- 1
  for(j in 2:ncol(x)){
    xx <- (if(a$df[j] == 1) x[,j] else m(x[,j], df = a$df[j], deg = a$degree, 
      kn = a$knots[[j]], bkn = a$bknots[[j]]))
    X <- cbind(X, xx)
  }
  X <- X[,a$sel2, drop = FALSE]
  X <- X[,-1,drop = FALSE]
  X <- scale(X, center = a$center, scale = a$scale)
  X <- X%*%a$rot
  X <- X[,a$sel3, drop = FALSE]
  cbind(1,X)
}


