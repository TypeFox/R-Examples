#' @importFrom stats integrate splinefun model.response model.weights model.matrix terms model.frame delete.response coef pnorm qexp
#' @importFrom stats sd prcomp lm.wfit pchisq weighted.mean printCoefmat .getXlevels pchisq runif ks.test vcov nobs
#' @importFrom survival Surv is.Surv survfit coxph
#' @importFrom graphics plot points abline
#' @importFrom utils menu
#' @import pch

#' @export
iqr <- function(formula, formula.p = ~ slp(p,3), weights, data, s, tol = 1e-5, maxit){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	ctiqr.internal(mf = mf,cl = cl, formula.p = formula.p, tol = tol, maxit = maxit, type = "iqr", s = s)
}


ctiqr <- function(formula, formula.p = ~ slp(p,3), weights, data, s, tol = 1e-5, maxit){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	ctiqr.internal(mf = mf,cl = cl, formula.p = formula.p, tol = tol, maxit = maxit, type = "ctiqr", s = s)
}


check.in <- function(mf, formula.p, type, s){

	if(!missing(s) && all(s == 0)){stop("'s' cannot be all zero")}
	explore.s <- function(s, dim){
		if(dim == 2){s <- t(s)}
		out <- 1
		if((r <- nrow(s)) > 1){
			for(j in 2:r){
				done <- FALSE; rj <- s[j,]
				for(h in 1:(j - 1)){
					if(all(rj == s[h,])){out[j] <- out[h]; done <- TRUE}
				}
				if(!done){out[j] <- max(out) + 1}
			}
		}
		out
	}

	# y,z,d, weights	

	zyd <- model.response(mf)
	if(type == "ctiqr" && !is.Surv(zyd))
	  {stop("the model response must be a Surv() object")}
	type <- attributes(zyd)$type
	if((n <- nrow(zyd <- cbind(zyd))) == 0){stop("zero non-NA cases")}
	if(is.null(type)){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- rep.int(1,n); type <- "iqr"}
	else if(type == "right"){
	  y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
	  type <- (if(any(d == 0)) "ciqr" else "iqr")
	}
	else if(type == "counting"){
	  z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- "ctiqr"
	  if(all(z < min(y))){type <- (if(any(d == 0)) "ciqr" else "iqr")}
	}
	if(!(any(d == 1))){stop("all data are censored")}
	if(any((weights <- model.weights(mf)) < 0)){stop("negative 'weights'")}
	weights <- (if(is.null(weights)) rep.int(1,n) else c(weights)/sum(weights)*n)


	# x and b(p)

	X <- model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
	termlabelsX <- attr(attr(mf, "terms"), "term.labels")
	assignX <- attr(X, "assign")
	coefnamesX <- colnames(X)
	p10 <- seq.int(1/1024, 1023/1024, length.out = 1023)
	p <- c(0.1^(6:4), p10, 1 - 0.1^(4:6))

	if((use.slp <- is.slp(formula.p))){
		k <- attr(use.slp, "k")
		intercept <- attr(use.slp, "intercept") 	# slp(0) = 0?
		intB <- attr(use.slp, "intB")			# b(p) includes 1?
		assignB <- (1 - intB):k
		termlabelsB <- paste("slp", 1:k, sep = "")
		coefnamesB <- (if(intB) c("(Intercept)", termlabelsB) else termlabelsB)
		k <- k + intB
	}
	else{
		k <- ncol(B <- model.matrix(formula.p, data = data.frame(p = p)))
		assignB <- attr(B, "assign")
		termlabelsB <- attr(terms(formula.p), "term.labels")
		coefnamesB <- colnames(B)
	}
	if(missing(s)){s <- matrix(1,q,k)}
	else{
		if(any(dim(s) != c(q,k))){stop("wrong size of 's'")}
		if(any(s != 0 & s != 1)){stop("'s' can only contain 0 and 1")}
	}
	
	# x singularities (set s = 0 where singularities occur)
	# x is dropped as in a linear model, irrespective of s.

	vx <- qr(X); selx <- vx$pivot[1:vx$rank]
	if(vx$rank < q){s[-selx,] <- 0}

	# b(p) singularities. Dropped row by row, based on s

	if(!use.slp && qr(B)$rank < k){
		u <- explore.s(s,1)
		for(j in unique(u)){
			sel <- which(s[which(u == j)[1],] == 1)
			if(length(sel) > 1){
				vbj <- qr(B[,sel, drop = FALSE])
				if((rj <- vbj$rank) < length(sel)){
					s[u == j, sel[-vbj$pivot[1:rj]]] <- 0
				}
			}
		}
	}

	# location-scale statistics for x, b(p), and y

	ry <- range(y); my <- ry[1]; My <- ry[2]

	sX <- apply(X,2,sd); mX <- colMeans(X)
	intX <- (length((constX <- which(sX == 0 & mX != 0))) > 0)
	varsX <- which(sX > 0); zeroX <- which(sX == 0 & mX == 0)
	sX[constX] <- X[1,constX]; mX[constX] <- 0; sX[zeroX] <- 1
	if(length(constX) > 1){zeroX <- c(zeroX, constX[-1]); constX <- constX[1]}

	if(!use.slp){
		sB <- apply(B,2,sd); mB <- colMeans(B)
		intB <- (length((constB <- which(sB == 0 & mB != 0))) > 0); varsB <- which(sB > 0)
		if(length(varsB) == 0){stop("the quantile function must depend on p")}
		if(length(constB) > 1){stop("remove multiple constant functions from 'formula.p'")}
		if(any(sB == 0 & mB == 0)){stop("remove zero functions from 'formula.p'")}
		sB[constB] <- B[1,constB]; mB[constB] <- 0
	}
	else{
		sB <- rep.int(1, k); mB <- rep.int(0, k)
		if(intB){constB <- 1; varsB <- 2:k}
		else{constB <- integer(0); varsB <- 1:k}
	}

	if(all(s[, varsB] == 0)){stop("the quantile function must depend on p (wrong specification of 's')")}
	if(!(theta00 <- ((intX & intB) && s[constX, constB] == 1)))
		{my <- 0; My <- sd(y)*5; mX <- rep.int(0,q)}
	else{for(j in varsX){if(any(s[j,] > s[constX,])){mX[j] <- 0}}}
	if(!intB | (intB && any(s[,constB] == 0))){mB <- rep.int(0,k)}


	# first scaling of x, b(p), y

	U <- list(X = X, y = y, z = z)
	X <- scale(X, center = mX, scale = sX)
	y <- (y - my)/(My - my)*10
	z <- z0 <- (z - my)/(My - my)*10
	z[z < min(y)] <- -Inf
	if(!use.slp){B <- scale(B, center = mB, scale = sB)}

	# principal component rotations that I can apply to x and b(p); second scaling

	rotX <- (diag(1,q))
	MX <- rep.int(0,q); SX <- rep.int(1,q)
	if(length(varsX) > 0){
		uX <- explore.s(s,1)
		X_in <- rowSums(s)	
		for(j in unique(uX)){
			sel <- which(uX == j)
			if(intX){sel <- sel[sel != constX]}
			if(length(sel) > 1 && X_in[sel[1]] != 0){
				PC <- prcomp(X[,sel], center = FALSE, scale. = FALSE)
				X[,sel] <- PC$x
				rotX[sel,sel] <- PC$rotation
			}
		}
		MX <- colMeans(X); MX[mX == 0] <- 0
		SX <- apply(X,2,sd); SX[constX] <- 1; SX[zeroX] <- 1
		X <- scale(X, center = MX, scale = SX)
	}

	rotB <- (diag(1,k))
	MB <- rep.int(0,k); SB <- rep.int(1,k)
	if(!use.slp){
		uB <- explore.s(s,2)
		B_in <- colSums(s)	
		for(j in unique(uB)){
			sel <- which(uB == j)
			if(intB){sel <- sel[sel != constB]}
			if(length(sel) > 1 && B_in[sel[1]] != 0){
				PC <- prcomp(B[,sel], center = FALSE, scale. = FALSE)
				B[,sel] <- PC$x
				rotB[sel,sel] <- PC$rotation
			}
		}
		MB <- colMeans(B); MB[mB == 0] <- 0
		SB <- apply(B,2,sd); SB[constB] <- 1
		B <- scale(B, center = MB, scale = SB)
	}

	# Create bfun

	if(!use.slp){
		bfun <- list()
		if(intB){
			bfun[[constB]] <- list(
				bfun = function(p, deriv = 0){rep.int(1 - deriv, length(p))},
				Bfun = function(p){p}, BBfun = function(p){0.5*p^2}, BB1 = 0.5
			)
		}
		for(j in varsB){bfun[[j]] <- make.bfun(p,B[,j], type)}
		names(bfun) <- coefnamesB
		bfun$BB1 <- sapply(bfun, function(x){x$BB1})
		bfun$BB1 <- matrix(rep(bfun$BB1, each = n), n)
		attr(bfun, "k") <- k
		bfun10 <- B[4:1026,, drop = FALSE]
	}
	else{
		k <- k - intB
		bfun <- slp.basis(k, intercept)
		bfun$BB1 <- colSums(bfun$AA)
		if(!intB){bfun$a[1,1] <- bfun$A[1,1] <- bfun$AA[1,1] <- 0; bfun$BB1 <- bfun$BB1[-1]}
		bfun$BB1 <- matrix(rep(bfun$BB1, each = n), n)
		attr(bfun, "intB") <- intB

		pp <- cbind(1,p10)
		if(k > 1){for(j in 2:k){pp <- cbind(pp, pp[,j]*p10)}}
		bfun10 <- pp%*%bfun$a
		if(!intB){bfun10 <- bfun10[,-1, drop = FALSE]}
	}
	attr(bfun, "p10") <- p10
	attr(bfun, "bfun10") <- bfun10

	# output. U = the original variables. V = the scaled/rotated variables.
	# stats.B, stats.X, stats.y = lists with the values use to scale/rotate

	stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
		intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
	stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
		intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
	stats.y <- list(m = my, M = My)

	V <- list(X = X, y = y, z = z, d = d, weights = weights)
	if(type == "ctiqr"){V$z0 <- z0}
	list(U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y, bfun = bfun, s = s, type = type)
}


check.out <- function(theta, S, covar){

	blockdiag <- function(A, d, type = 1){
		h <- nrow(A); g <- d/h
		if(type == 1){
			out <- diag(1,d)
			for(j in 1:g){ind <- (j*h - h  + 1):(j*h); out[ind,ind] <- A}
		}
		else{
			out <- matrix(0,d,d)
			for(i1 in 1:h){
				for(i2 in 1:h){
					ind1 <- (i1*g - g  + 1):(i1*g)
					ind2 <- (i2*g - g  + 1):(i2*g)
					out[ind1, ind2] <- diag(A[i1,i2],g)
				}
			}
			out <- t(out)
		}
		out
	}

	mydiag <- function(x){
		if(length(x) > 1){return(diag(x))}
		else{matrix(x,1,1)}
	}

	th <- cbind(c(theta))
	q <- nrow(theta)
	k <- ncol(theta)
	g <- q*k
	aX <- S$X; ay <- S$y; aB <- S$B
	cX <- aX$const; cB <- aB$const

	##########################

	A <- blockdiag(mydiag(1/aX$S), g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$M; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$S),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$M; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(aX$rot,g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	A <- blockdiag(t(aB$rot),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	##########################

	A <- blockdiag(mydiag(1/aX$s),g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$m/aX$s[cX]; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$s),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$m/aB$s[cB]; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	v <- (ay$M - ay$m)/10
	th <- th*v
	covar <- covar*(v^2)
	theta <- matrix(th,q,k)
	theta[cX,cB] <- theta[cX,cB] + ay$m/aB$s[cB]/aX$s[cX]
	
	list(theta = theta, covar = covar)
}


# integrated quantile regression

ctiqr.internal <- function(mf,cl, formula.p, tol = 1e-6, maxit, type, s){

	n <- nrow(mf)
	A <- check.in(mf, formula.p, type, s); V <- A$V; U <- A$U; s <- A$s; type <- A$type
	S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
	bfun <- A$bfun; attributes(bfun) <- c(attributes(bfun), S$B)

	if(missing(maxit)){maxit <- 10 + 10*attr(bfun, "k")*ncol(V$X)}
	else{maxit <- pmax(10, maxit)}

	if((q <- length(S$X$vars)) > 0){
		xx <- list(V$X, cbind(V$X, V$X^2))
		DF <- c(10,10,20,20,10,10,20,20)
		wX <- c(1,1,1,1,2,2,2,2)
		wy <- c(1,2,1,2,1,2,1,2)
	}
	else{xx <- list(V$X); wX <- c(1,1,1,1,1,1); DF <- c(5,10,20,5,10,20); wy <- c(1,1,1,2,2,2)}
	if(type != "iqr" | q > 0){
		yy <- list(V$y, qexp((rank(V$y) - 0.5)/n))
		if(type == "ctiqr"){zz <- list(V$z, splinefun(yy[[1]],yy[[2]], method = "monoH.FC")(V$z))}
		else{zz <- list(V$z, V$z)}
	}
	else{yy <- zz <- NULL; DF <- wX <- wy <- 1}

	Fit <- NULL
	for(i in 1:length(DF)){

		theta0 <- start.theta(V$y, V$z, V$d, V$X, V$weights, bfun, df = DF[i], 
			xx = xx[[wX[i]]], yy = yy[[wy[i]]], zz = zz[[wy[i]]], s = s)

		fit <- iqr.newton(theta0, V$y, V$z, V$d, V$X, V$weights, 
			bfun, type = type, tol = tol, maxit = maxit, s = s)

		if(fit.ok <- (fit$rank == ncol(fit$jacobian) & mean(abs(fit$ee)) < 1)){
			covar <- try(cov.theta(fit$coefficients, V$y, V$z, V$d, V$X,
				V$weights, bfun, fit$CDF.y, fit$CDF.z, type, s = s), silent = TRUE)
			fit.ok <- (class(covar) != "try-error")
		}

		covar.ok <- (if(fit.ok){(qr(covar$Q)$rank == fit$rank)} else FALSE)
		if(fit.ok & covar.ok){break}
		else if(fit.ok){Fit <- fit; Cov <- covar}
	}

	if(!fit.ok && is.null(Fit)){stop("unable to fit the model: this can be due to severe misspecification")}
	if(!covar.ok){warning("the estimated covariance matrix is deemed to be singular")}
	if(!fit.ok){fit <- Fit; covar <- Cov}
	if(!fit$converged){warning("the algorithm did not converge")}

	# minimized loss function

	if(type == "iqr"){
		v <- (S$y$M - S$y$m)/10
		fit$obj.function <- iobjfun(fit$coef, V$y,V$X,V$weights, bfun, fit$CDF.y)*v
	}

	# fitted CDFs

	CDFs <- data.frame(CDF.y = fit$CDF.y, CDF.t = fit$CDF.y,
                     CDF.z = (if(type == "ctiqr") fit$CDF.z else NA))
	if(type != "iqr"){
		kmfit <- suppressWarnings(km(fit$CDF.z, fit$CDF.y, V$d, V$weights, type))
		CDFs$CDF.t <- kmfit$Q((1:n)/(n + 1))
		attr(CDFs, "km") <- kmfit
	}

	# output
	
	bfun0 <- bfun # save it for test.fit
	bfun10 <- attr(bfun, "bfun10")
	for(j in 1:ncol(bfun10)){bfun10[,j] <- S$B$M[j] + bfun10[,j]*S$B$S[j]}
	bfun10 <- bfun10%*%qr.solve(qr(S$B$rot))
	for(j in 1:ncol(bfun10)){bfun10[,j] <- S$B$m[j] + bfun10[,j]*S$B$s[j]}
	attr(bfun, "bfun10") <- bfun10
  
	attr(mf, "assign") <- S$X$assign
	attr(mf, "stats") <- S
	attr(mf, "all.vars") <- V
	attr(mf, "CDFs") <- CDFs
	attr(mf, "bfun0") <- bfun0
	attr(mf, "bfun") <- bfun
	attr(mf, "theta") <- fit$coefficients
	attr(mf, "type") <- type

	out <- check.out(fit$coefficients, S, covar = covar$Q)
	fit <- list(coefficients = out$theta, call = cl, 
		converged = fit$converged, n.it = fit$n.it,
		obj.function = fit$obj.function,  
		covar = out$covar, mf = mf, s = s, CDF = fit$CDF.y)
	jnames <- c(sapply(attr(bfun, "coef.names"), 
		function(x,y){paste(x,y, sep = ":")}, y = S$X$coef.names))
	dimnames(fit$covar) <- list(jnames, jnames)
	dimnames(fit$coefficients) <- dimnames(fit$s) <- list(S$X$coef.names, S$B$coef.names)


	# PDF

	b1 <- apply_bfun(bfun, fit$CDF, "b1fun", convert = TRUE)
	fit$PDF <- 1/c(rowSums((U$X%*%fit$coef)*b1))
	fit$PDF[attr(fit$CDF, "out")] <- 0
	attributes(fit$CDF) <- attributes(fit$PDF) <- list(names = rownames(mf))
	if(any(fit$PDF < 0)){warning("quantile crossing detected (PDF < 0 at some y)")}

	# finish

	class(fit) <- "iqr"
	fit
}

#' @export
print.iqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	invisible(x)
}


apply_bfun <- function(bfun, p, fun, convert = FALSE){

	k <- attr(bfun, "k")
	if(class(bfun) == "slp.basis"){

		# for "bfun" (only used in p.bisec), the powers of p must be computed
		# otherwise, p is already a matrix (1,p,p^2,p^3,...)

		intB <- attr(bfun, "intB")
		if(fun == "bfun"){
			pp <- cbind(1,p)
			if(k > 1){for(j in 2:k){pp <- cbind(pp, pp[,j]*p)}}
			out <- pp%*%bfun$a
			if(!intB){out <- out[,-1, drop = FALSE]}
			attr(out, "p^") <- pp
			return(out)
		}
		p <- attr(p, "p^")
		if(fun == "b1fun"){out <- cbind(0,p[,1:k, drop = FALSE]%*%bfun$a1[-1,-1, drop = FALSE])}
		else if(fun == "Bfun"){out <- p[,2:(k + 2)]%*%bfun$A}
		else{out <- p[,3:(k + 3)]%*%bfun$AA}
		if(!intB){out <- out[,-1, drop = FALSE]}
		return(out)
	}
	else{
		out <- matrix(NA_real_)
		n <- length(p)
		length(out) <- n*k
		dim(out) <- c(n,k)
		if(fun == "b1fun"){for(j in 1:k){out[,j] <- bfun[[j]]$bfun(p, deriv = 1)}}
		else{for(j in 1:k){out[,j] <- bfun[[j]][[fun]](p)}}
		if(convert){
			ap <- attributes(bfun)
			if(fun == "bfun"){phi <- 1}
			else if(fun == "b1fun"){phi <- 0}
			else if(fun == "Bfun"){phi <- p}
			else{phi <- 0.5*p^2}

			for(j in 1:k){out[,j] <- ap$M[j]*phi + out[,j]*ap$S[j]}
			out <- out%*%qr.solve(qr(ap$rot))
			for(j in 1:k){out[,j] <- ap$m[j]*phi + out[,j]*ap$s[j]}
		}
		return(out)
	}
}


p.bisec <- function(theta, y,X, bfun, n.it = 17, convert = FALSE){

	p10 <- attr(bfun, "p10")
	bfun10 <- attr(bfun, "bfun10")%*%t(theta)
	n <- length(y); q <- ncol(X)
	m <- rep.int(512,n)
	for(i in 2:10){
		delta.m <- y - .rowSums(X*bfun10[m,, drop = FALSE], n, q)
		m <- m + sign(delta.m)*(2^(10 - i))
	}
	m <- p10[m]
	for(i in 11:n.it){
		bp <- apply_bfun(bfun, m, "bfun", convert = convert)
		delta.m <- y - .rowSums(X*(bp%*%t(theta)), n,q)
		m <- m + sign(delta.m)/2^i
	}
	m <- c(m)
	out.l <- which(m == 1/2^n.it)
	out.r <- which(m == 1 - 1/2^n.it)
	m[out.l] <- 0; m[out.r] <- 1
	attr(m, "out") <- c(out.l, out.r)
	attr(m, "out.r") <- out.r

	bp <- apply_bfun(bfun, m, "bfun", convert = convert)
	attr(m, "bfun") <- bp[,, drop = FALSE]

	if(class(bfun) == "slp.basis"){
		k <- attr(bfun, "k")
		pp <- attr(bp, "p^")
		pp <- cbind(pp, pp[,k + 1]*m, pp[,k + 1]*m^2)
		attr(m, "p^") <- pp
	}
	m
}


make.bfun <- function(p,x, type){
	bfun <- splinefun(p,x, method = "fmm")
	Bfun <- numint(bfun)
	if(type != "iqr"){
		BBfun <- numint(Bfun)
		BB1 <- BBfun(1)
	}
	else{
		BBfun <- NULL
		BB1 <- integrate(Bfun, lower = 0, upper = 1, subdivisions = 10000)$value
	}
	list(bfun = bfun, Bfun = Bfun, BBfun = BBfun, BB1 = BB1)
}

numint <- function(f){

	n <- 1000
	delta <- 1/n
	ul <- c(0,(0.1)^(6:4), 1/n)
	ur <- sort(1 - ul)
	u <- seq.int(delta, 1 - delta, length.out = n - 1)

	ff <- f(u)
	FF <- (ff[1:(n - 2)] + ff[2:(n - 1)])/2*delta

	# left tail

	FF_l <- 0
	for(j in 1:(length(ul) - 1)){
		FF_l[j + 1] <- integrate(f, lower = ul[j], upper = ul[j + 1], 
		subdivisions = 1000)$value
	}

	# right tail

	FF_r <- NULL
	for(j in 1:(length(ur) - 1)){
		FF_r[j] <- integrate(f, lower = ur[j], upper = ur[j + 1], 
		subdivisions = 1000)$value
	}

	## end

	u <- c(ul,u[-1],ur[-1])
	FF <- cumsum(c(FF_l, FF, FF_r))
	splinefun(u, FF, method = "fmm")
}

#' @export
plf <- function(p, knots){ # basis of piecewise linear function
	if(is.null(knots)){return(cbind(b1 = p))}
	k <- length(knots <- sort(knots))
	ind1 <- 1
	ind2 <- NULL
	for(j in knots){ind1 <- cbind(ind1, (p > j))}
	for(j in k:1){ind2 <- cbind(ind2, ind1[,j] - ind1[,j+1])}
	ind2 <- cbind(ind1[,k+1], ind2)[,(k + 1):1]
	ind1 <- cbind(ind1,0)
	knots <- c(0,knots,1)
	a <- NULL
	for(j in 1:(k + 1)){
		a <- cbind(a, (p - knots[j])*ind2[,j] + (knots[j + 1] - knots[j])*ind1[,j + 1])
	}
	colnames(a) <- paste("b", 1:(k + 1), sep = "")
	attr(a, "knots") <- knots[2:(k+1)]
	a
}

slp.basis <- function(k, intercept){ # shifted Legendre polynomials basis

	K <- k + 1

	# matrix a such that P%*%a is an orthogonal polynomial, P = (1, p, p^2, p^3, ...)

	a <- matrix(0, K,K)
	for(i1 in 0:k){
		for(i2 in 0:i1){
			a[i2 + 1, i1 + 1] <- choose(i1,i2)*choose(i1 + i2, i2)
		}
	}
	a[,seq.int(2, K, 2)] <- -a[,seq.int(2, K, 2)]
	a[seq.int(2, K, 2),] <- -a[seq.int(2, K, 2),]

	# a1 = first derivatives to be applied to P' = (0,1, p, p^2, ...)
	# A = first integral to be applied to PP = (p, p^2, p^3, p^4, ...)
	# AA = second integral to be applied to PPP = (p^2, p^3, p^4, p^5, ...)


	a1 <- A <- AA <- matrix(,K,K)
	for(j in 0:k){
		a1[j + 1,] <- a[j + 1,]*j
		A[j + 1,] <- a[j + 1,]/(j + 1)
		AA[j + 1,] <- A[j + 1,]/(j + 2)
	}

	if(!intercept){a[1,-1] <- A[1,-1] <- AA[1, -1] <- 0}
	out <- list(a = a, a1 = a1, A = A, AA = AA)
	attr(out, "k") <- k
	class(out) <- "slp.basis"
	out
}


#' @export
slp <- function(p, k = 3, intercept = FALSE){
	if((k <- round(k)) < 1){stop("k >= 1 is required")}
	P <- cbind(1, outer(p, seq_len(k), "^"))
	B <- P%*%slp.basis(k, intercept)$a
	colnames(B) <- paste("slp", 0:k, sep = "")
	B <- B[,-1, drop = FALSE]
	attr(B, "k") <- k
	class(B) <- "slp"
	B
}


is.slp <- function(f){
	test.p <- seq(0,1,0.1)
	B <- model.matrix(f, data = data.frame(p = test.p))
	if(nrow(B) == 0){return(FALSE)}
	a <- attr(B, "assign")
	if(any(a > 1)){return(FALSE)}
	B <- B[,a == 1, drop = FALSE]
	k <- ncol(B)
	intercept <- FALSE
	if(any(B != slp(test.p, k = k, intercept))){
		intercept <- TRUE
		if(any(B != slp(test.p, k = k, intercept))){
			return(FALSE)
		}
	}
	out <- TRUE
	attr(out, "k") <- k
	attr(out, "intercept") <- intercept
	attr(out, "intB") <- any(a == 0)
	out
}


iobjfun <- function(theta, y,X,weights, bfun, p.star){
	s <- NULL
	B <- apply_bfun(bfun, p.star, "Bfun")
	for(j in 1:ncol(B)){s <- cbind(s, bfun$BB1[1,j] - B[,j])}
	sum(y*weights*(p.star - 0.5)) + sum(((X*weights)%*%theta)*s)
}

iqr.ee <- function(theta, y,z,d,X,weights, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE){

	k <- ncol(theta)
	BB1 <- bfun$BB1
	if(missing(G)){
		B <- apply_bfun(bfun, p.star.y, "Bfun")
		S1 <- BB1 - B
		if(!i){g <- c(t(X*weights)%*%S1)}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B <- G$B; g <- G$g}
 
	if(J){
		b <- attr(p.star.y, "bfun")
		b1 <- apply_bfun(bfun, p.star.y, "b1fun")
		A1 <- 1/c(rowSums((X%*%theta)*b1))
		A1 <- pmax(0,A1)
		A1[attr(p.star.y, "out")] <- 0
		Xw <- X*weights*A1

		J <- NULL
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){h.temp <- cbind(h.temp, t(Xw)%*%(X*(b[,i2]*b[,i1])))}
			J <- rbind(J, h.temp)
		}
	}
	
	list(g = g, J = J, B = B)
}

ciqr.ee <- function(theta, y,z,d,X,weights, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE){

	k <- ncol(theta)
	BB1 <- bfun$BB1

	if(missing(G)){
		B <- apply_bfun(bfun, p.star.y, "Bfun")
		BB <- apply_bfun(bfun, p.star.y, "BBfun")
		a <- (1 - d)/(1 - p.star.y); a[attr(p.star.y, "out.r")] <- 0
		S1 <- a*(BB - BB1)
		a <- d; a[attr(p.star.y, "out.r")] <- 1
		S1 <- BB1 - a*B + S1
		if(!i){g <- c(t(X*weights)%*%S1)}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B <- G$B; BB <- G$BB; g <- G$g}

	if(J){
		b1 <- apply_bfun(bfun, p.star.y, "b1fun")
		b <- attr(p.star.y, "bfun")
		A1 <- 1/c(rowSums((X%*%theta)*b1))
		A1 <- pmax(0,A1)
		A1[attr(p.star.y, "out")] <- 0
		a <- 1 - p.star.y
		a[attr(p.star.y, "out.r")] <- 1 
		# the vale 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1)
		A2 <- d*b + (1 - d)/(a^2)*(BB1 - B*a - BB)
		Xw <- X*weights*A1

		J <- NULL
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){h.temp <- cbind(h.temp, t(Xw)%*%(X*(b[,i2]*A2[,i1])))}
			J <- rbind(J, h.temp)
		}
	}
	list(g = g, J = J, B = B, BB = BB)
}




ctiqr.ee <- function(theta, y,z,d,X,weights, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE){

	k <- ncol(theta)
	BB1 <- bfun$BB1
	if(missing(G)){
		B.y <- apply_bfun(bfun, p.star.y, "Bfun")	
		BB.y <- apply_bfun(bfun, p.star.y, "BBfun")
		B.z <- apply_bfun(bfun, p.star.z, "Bfun")
		BB.z <- apply_bfun(bfun, p.star.z, "BBfun")

		out.y <- attr(p.star.y, "out.r")
		out.z <- attr(p.star.z, "out.r")
		a <- d

		S1.y <- (1 - d)/(1 - p.star.y)*(BB.y - BB1)
		S1.z <- 1/(1 - p.star.z)*(BB.z - BB1)
		S1.y[out.y,] <- 0; a[out.y] <- 1
		S1.y[out.z,] <- S1.z[out.z,] <- a[out.z] <- 0

		S1 <- -a*B.y + S1.y - S1.z
		if(!i){g <- c(t(X*weights)%*%S1)}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B.y <- G$B.y; BB.y <- G$BB.y; B.z <- G$B.z; BB.z <- G$BB.z; g <- G$g}

	if(J){
		b1.y <- apply_bfun(bfun, p.star.y, "b1fun")
		b.y <- attr(p.star.y, "bfun")
		b1.z <- apply_bfun(bfun, p.star.z, "b1fun")
		b.z <- attr(p.star.z, "bfun")

		A1.y <- 1/c(rowSums((X%*%theta)*b1.y))
		A1.z <- 1/c(rowSums((X%*%theta)*b1.z))
		A1.y <- pmax(A1.y,0)
		A1.z <- pmax(A1.z,0)
		A1.y[attr(p.star.y, "out")] <- 0
		A1.z[attr(p.star.z, "out")] <- 0

		ay <- 1 - p.star.y; az <- 1 - p.star.z
		ay[attr(p.star.y, "out.r")] <- az[attr(p.star.z, "out.r")] <- 1 
		# the value 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1.y and A1.z)

		A2.y <- d*b.y - (1 - d)/(ay^2)*(B.y*ay + BB.y - BB1)
		A2.z <- 1/(az^2)*(B.z*az + BB.z - BB1)
		H.y <- A2.y*A1.y 
		H.z <- A2.z*A1.z
		Xw <- X*weights

		J <- NULL
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){h.temp <- cbind(h.temp, t(Xw)%*%(X*(b.y[,i2]*H.y[,i1] + b.z[,i2]*H.z[,i1])))}
			J <- rbind(J, h.temp)
		}
	}
	list(g = g, J = J, B.y = B.y, BB.y = BB.y, B.z = B.z, BB.z = BB.z)
}




start.theta <- function(y,z,d, x, weights, bfun, df, xx, yy, zz, s){

	if(is.null(yy)){p.star <- (rank(y) - 0.5)/length(y)}
	else{
	  m0 <- suppressWarnings(pch:::pch.fit(z = zz, y = yy, d = d, 
		x = cbind(1,xx), w = weights, breaks = df))
	  p.star <- 1 - pch:::predF.pch(m0)[,3]
	}

	b.star <- apply_bfun(bfun, p.star, "bfun")
	X <- model.matrix(~ -1 + x:b.star)
	X <- X[, c(s) == 1, drop = FALSE]
	start.ok <- FALSE
	while(!start.ok){
		m <- lm.wfit(X, y, weights)
		res <- m$residuals
		start.ok <- all(w <- (abs(res)/sd(res) < 4))
		if(start.ok | sum(w) < 0.5*length(y)){break}
		X <- X[w,, drop = FALSE]
		y <- y[w]
		weights <- weights[w]
	}
	out <- rep.int(0, length(s))
	out[c(s) == 1] <- m$coef
	out <- matrix(out, ncol(x))
	out[is.na(out)] <- 0
	out
}


# Note: if s has some zeroes, the covariance matrix Q will contain some zero-columns and rows,
# while the gradient and jacobian will just omit the parameters that are not estimated

cov.theta <- function(theta, y,z,d,X,weights, bfun, 
		p.star.y, p.star.z, type, s){

	if(type == "iqr"){ee <- iqr.ee}
	else if(type == "ciqr"){ee <- ciqr.ee}
	else{ee <- ctiqr.ee}
	s <- c(s == 1)
	
	G.i <- ee(theta, y,z,d,X,weights, bfun, p.star.y, p.star.z, J = TRUE, i = TRUE)
	s.i <- G.i$g[,s, drop = FALSE]
	G <- G.i$J[s,s, drop = FALSE]

	Omega <- chol2inv(chol(t(s.i*weights)%*%s.i))
	Q <- t(G)%*%Omega%*%G
	Q <- chol2inv(chol(Q))
	U <- matrix(0, length(s), length(s))
	U[s,s] <- Q

	list(Q = U, jacobian = G, ee = colSums(s.i*weights), Omega = Omega, s.i = s.i)
}



# Observed = n. of CDF.y in each interval, no matter if type = "iqr" or not.
# If type = "iqr", expected = n/k.
# Otherwise, "expected" = n. of cases in which phi_t = km(CDF.z, CDF.y, d)$F(CDF.y) 
# falls in the same interval. If phi_t = U(0,1), and a value of CDF.y
# is in a given interval, then its phi_t should also be in the same interval.
# I quote "expected" because both observed and expected are estimated quantities in this case.

chitest <- function(y, w, k = 10, type, Fy){
  
	if(k < 3){stop("k >= 3 is required")}
	n <- length(y)
	u <- seq(0,1,length = k + 1)
	u[1] <- -Inf # to include the zeroes
  
	# expected
  
	if(type == "iqr"){p <- 1/k; e <- n/k}
	else{
		p <- NULL
		for(j in 1:k){
			if(is.null(w)){p[j] <- mean(Fy > u[j] & Fy <= u[j + 1])}
			else{p[j] <- weighted.mean(Fy > u[j] & Fy <= u[j + 1], w)}
		}
		e <- p*n
	}
	v <- n*p*(1 - p)

	# observed
  
	o <- NULL
	for(j in 1:k){
		if(is.null(w)){o[j] <- mean(y > u[j] & y <= u[j + 1])}
		else{o[j] <- weighted.mean(y > u[j] & y <= u[j + 1], w)}
	}
	if(type != "iqr"){v <- v + n*o*(1 - o)}
	o <- o*n
  
	# finish

	chi2 <- sum(((o - e)^2)/v)
	rbind(c('chi^2' = chi2, df = k - 2, 
	'p-value' = pchisq(chi2, df = k - 2, lower.tail = FALSE)))
}

#' @export
summary.iqr <- function(object, p, cov = FALSE, k = 10, ...){

	if(missing(p)){
		mf <- object$mf
		theta <- object$coefficients
		CDF.y <- attr(mf, "CDFs")$CDF.y
		w <- attr(mf, "all.vars")$weights
		type <- attr(mf, "type")
		if(type != "iqr"){
			km <- attr(attr(mf, "CDFs"), "km")
			Fy <- km$F(CDF.y)
		}
  
		u <- sqrt(diag(object$covar))
		u <- matrix(u, q <- nrow(theta), r <- ncol(theta))
		dimnames(u) <- dimnames(theta)
		test <- (if(q*r == 1) NULL else iqr.waldtest(object))
		out <- list(converged = object$converged, n.it = object$n.it,
			coefficients = theta, se = u, 
			test.x = test$test.x, test.p = test$test.p)

		out$obj.function <- object$obj.function
		testfit <- NULL; k <- unique(round(k))

		for(j in k){testfit <- rbind(testfit, chitest(CDF.y, w, j, type, Fy))}
		out$test.fit <- testfit
		rownames(out$test.fit) <- paste("k =", k)
		
		out$n <- nrow(object$mf)
		out$free.par <- sum(theta != 0)
	}
	else{
		out <- list()
		for(i in 1:length(p)){out[[i]] <- extract.p(object, p[i], cov)}
		names(out) <- paste("p =", p)
		attr(out, "nacoef") <- which(apply(object$coefficients,1, function(v){all(v == 0)}))
	}
	out$call <- object$call
	class(out) <- "summary.iqr"
	out	
}

#' @export
print.summary.iqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if(!is.null(x$coef)){

		nacoef <- which(x$coef == 0)
		x$coef[nacoef] <- x$se[nacoef] <- NA

		cat("converged:", x$converged, "\n")
		cat("n. of iterations:", x$n.it, "\n")
		cat("n. of observations:", x$n, "\n")
		cat("n. of free parameters:", x$free.par, "\n\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Coefficients:\n")
		print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("Standard errors:\n")
		print.default(format(x$se, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		if(!is.null(x$test.x)){
			cat("Wald test for x:\n")

			printCoefmat(x$test.x, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n\n")
		}

		if(!is.null(x$test.p)){
			cat("Wald test for b(p):\n")
			printCoefmat(x$test.p, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}

		if(!is.null(x$obj.function)){
			cat("\n")
			cat("Minimized loss function:", x$obj.function)
		}
		cat("\n\n")
		cat("Chi-squared test for goodness-of-fit:\n")
		printCoefmat(x$test.fit, digits = digits, signif.stars = TRUE, tst.ind = 1, 
			signif.legend = FALSE, P.values = TRUE, has.Pvalue = TRUE)
	}

	else{
		nacoef <- attr(x, "nacoef")
		for(j in 1:(length(x) - 1)){
			cat(paste(names(x)[j], "\n"))
			cat("\n")
			cat("Coefficients:\n")
			coe <- x[[j]]$coef; coe[nacoef,] <- NA
			printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n")

			if(!is.null(x[[j]]$cov)){
				cat("Covar:\n")
				print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
			}
			cat("\n\n")
		}
	}

	invisible(x)
}



extract.p <- function(model,p, cov = FALSE){

	theta <- model$coefficients
	v <- model$covar
	q <- nrow(theta)
	k <- ncol(theta)

	bfun <- attr(model$mf, "bfun")
	pred.p <- apply_bfun(bfun, p, "bfun", convert = TRUE)
	beta <- c(pred.p%*%t(theta))

	cov.beta <- matrix(NA,q,q)
	for(j1 in 1:q){
		w1 <- seq.int(j1,q*k,q)
		for(j2 in j1:q){
			w2 <- seq.int(j2,q*k,q)
			cc <- v[w1,w2, drop = FALSE]
			cov.beta[j1,j2] <- cov.beta[j2,j1] <- pred.p%*%cc%*%t(pred.p)
		}
	}
	se <- sqrt(diag(cov.beta))
	z <- beta/se
	out <- cbind(beta, se, z, 2*pnorm(-abs(z)))
	colnames(out) <- c("Estimate", "std.err", "z value", "p(>|z|))")
	rownames(out) <- colnames(cov.beta) <- rownames(cov.beta) <- rownames(theta)
	if(cov){list(coef = out, cov = cov.beta)}
	else{list(coef = out)}
}

#' @export
plot.iqr <- function(x, conf.int = TRUE, which = NULL, ask = TRUE, ...){

	plot.iqr.int <- function(p,u,j,conf.int,L){
		beta <- u[[j]]$beta
		if(is.null(L$ylim)){
			if(conf.int){y1 <- min(u[[j]]$low); y2 <- max(u[[j]]$up)}
			else{y1 <- min(beta); y2 <- max(beta)}
			L$ylim <- c(y1,y2)
		}
		plot(p, u[[j]]$beta, xlab = L$xlab, ylab = L$ylab, main = L$labels[j], 
		type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col)
		if(conf.int){
		  points(p, u[[j]]$low, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		  points(p, u[[j]]$up, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		}
	}

	L <- list(...)
	if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
	if(is.null(L$lwd)){L$lwd <- 2}
	if(is.null(L$col)){L$col <- "black"}
	if(is.null(L$xlab)){L$xlab <- "p"}
	if(is.null(L$ylab)){L$ylab <- "beta(p)"}
	L$labels <- rownames(x$coefficients)
	q <- length(L$labels)
	L$labels <- c(L$labels, "qqplot")

	p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
	u <- predict.iqr(x, p = p, type = "beta", se = conf.int)

	if(!is.null(which) | !ask){
		if(is.null(which)){which <- 1:q}
		for(j in which){plot.iqr.int(p,u,j,conf.int,L)}
	}
	else{
		pick <- 1
		while(pick > 0 && pick <= q + 1){
			pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
			if(pick > 0 && pick <= q){plot.iqr.int(p,u,pick,conf.int,L)}
			else if(pick == q + 1){
				n <- length(FF <- sort(attr(x$mf, "CDFs")$CDF.t))
				plot((1:n)/(n + 1), FF, type = "l", lwd = 2, xlab = "U(0,1) quantiles", ylab = "fitted CDF quantiles")
				abline(0,1)
			}
		}
	}
}


# predict function.
# p: default to percentiles for type = "beta". No default for "fitted". Ignored for "CDF".
# se: ignored for type = "CDF"
# x: only for type = "CDF" or type = "fitted"
# y: only for type = "CDF"


#' @export
predict.iqr <- function(object, type = c("beta", "CDF", "QF", "sim"), newdata, p, se = TRUE, ...){

	if(is.na(match(type <- type[1], c("beta", "CDF", "QF", "sim")))){stop("invalid 'type'")}
	if(type == "beta"){
		if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
		return(pred.beta(object, p, se))
	}

	mf <- object$mf
	mt <- terms(mf)
	miss <- attr(mf, "na.action")
	fittype <- attr(mf, "type")
	nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	xlev <- .getXlevels(mt, mf)

	if(!missing(newdata)){

	  if(type == "CDF"){
	    yn <- as.character(if(fittype == "ctiqr") mt[[2]][[3]] 
              else if(fittype == "ciqr") mt[[2]][[2]] else mt[[2]])
	    if(is.na(ind <- match(yn, colnames(newdata))))
	    {stop("for 'type = CDF', 'newdata' must contain the y-variable")}
	    if(fittype == "ciqr"){newdata[,as.character(mt[[2]][[3]])] <- 1}
	    if(fittype == "ctiqr"){newdata[,as.character(mt[[2]][[4]])] <- 1
              newdata[,as.character(mt[[2]][[2]])] <- -Inf}
	  }
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}

		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){
			nr <- nrow(newdata)
			if(type == "CDF"){
				out <- data.frame(matrix(NA,nr,3))
				colnames(out) <- c("log.f", "log.F", "log.S")
				rownames(out) <- rownames(newdata)
			}
			else if(type == "QF"){
				out <- data.frame(matrix(NA,nr,length(p)))
				colnames(out) <- paste("p",p, sep = "")
				rownames(out) <- rownames(newdata)
			}
			else{out <- rep.int(NA, nr)}
			return(out)
		}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
	}

	x <- model.matrix(mt, mf)

	if(type == "CDF"){
		bfun <- attr(object$mf, "bfun")
		y <- cbind(model.response(mf))[,1 + (fittype == "ctiqr")]
		Fy <- p.bisec(object$coefficients, y,x, bfun, convert = TRUE)
		b1 <- apply_bfun(bfun, Fy, "b1fun", convert = TRUE)
		fy <- 1/c(rowSums((x%*%object$coefficients)*b1))
		fy[attr(Fy, "out")] <- 0
		if(any(fy < 0)){warning("some PDF values are negative (quantile crossing)")}
		CDF <- PDF <- NULL
		CDF[nomiss] <- Fy
		PDF[nomiss] <- fy
		CDF[miss] <- PDF[miss] <- NA
		out <- data.frame(CDF = CDF, PDF = PDF)
		rownames(out)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
		return(out)
	}

	else if(type == "QF"){
		if(missing(p)){stop("please indicate the value(s) of 'p' to compute x*beta(p)")}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}

		fit <- se.fit <- matrix(, length(c(miss,nomiss)), length(p))
		for(j in 1:length(p)){
			fit.beta <- extract.p(object,p[j], cov = se)
			fit[nomiss,j] <- x%*%cbind(fit.beta$coef[,1])
			if(se){se.fit[nomiss,j] <- sqrt(diag(x%*%fit.beta$cov%*%t(x)))}
		}
		fit <- data.frame(fit)
		colnames(fit) <- paste("p",p, sep = "")
		rownames(fit)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(fit)[miss] <- names(miss)}
		if(se){
			se.fit <- data.frame(se.fit)
			colnames(se.fit) <- paste("p",p, sep = "")
			rownames(se.fit)[nomiss] <- rownames(mf)
			if(!is.null(miss)){rownames(se.fit)[miss] <- names(miss)}
			return(list(fit = fit, se.fit = se.fit))
		}
		else{return(fit)}
	}	
	else{
		p <- runif(nrow(x))
		beta <- apply_bfun(attr(object$mf, "bfun"), p, "bfun", convert = TRUE)%*%t(object$coefficients)
		y <- NULL; y[nomiss] <- rowSums(beta*x); y[miss] <- NA
		return(y)
	}
}




pred.beta <- function(model, p, se = FALSE){

	if(se){
		Beta <- NULL
		SE <- NULL
		for(j in p){
			b <- extract.p(model,j)$coef
			Beta <- rbind(Beta, b[,1])
			SE <- rbind(SE, b[,2])
		}
		out <- list()
		for(j in 1:ncol(Beta)){
			low <- Beta[,j] - 1.96*SE[,j]
			up <- Beta[,j] + 1.96*SE[,j]
			out[[j]] <- data.frame(p = p, beta = Beta[,j], se = SE[,j], low = low, up = up)
		}
		names(out) <- rownames(model$coefficients)
		return(out)
	}
	else{
		theta <- model$coefficients
		beta <- apply_bfun(attr(model$mf, "bfun"), p, "bfun", convert = TRUE)%*%t(theta)
		out <- list()
		for(j in 1:nrow(theta)){out[[j]] <- data.frame(p = p, beta = beta[,j])}
		names(out) <- rownames(theta)
		return(out)
	}
}



iqr.newton <- function(theta, y,z,d,X,weights, bfun, type, tol = 1e-5, maxit = 200, s){

	if(type == "iqr"){ee <- iqr.ee}
	else if(type == "ciqr"){ee <- ciqr.ee}
	else{ee <- ctiqr.ee}

	q <- nrow(theta)
	k <- ncol(theta)
	s <- c(s == 1)

	p.star.y <- p.bisec(theta, y,X, bfun)
	p.star.z <- (if(type == "ctiqr") pmin(p.bisec(theta, z,X, bfun), p.star.y) else NULL)
	G <- ee(theta, y,z,d,X,weights, bfun, p.star.y, p.star.z)
	g <- G$g[s]
	h <- G$J[s,s, drop = FALSE]
	conv <- FALSE
	eps <- .1
	alg <- "nr"

	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####

		if(type == "iqr"){
			H1 <- try(chol(h), silent = TRUE)
			err <- (class(H1) == "try-error")
		}
		else{
			H1 <- qr(h)
			r <- H1$rank
			err <- (r != ncol(h))
		}
		if(!err){
			if(alg == "gs"){alg <- "nr"; eps <- 1}
			delta <- (if(type == "iqr") chol2inv(H1)%*%g else qr.solve(H1)%*%g)
		}
		else{
			if(alg == "nr"){alg <- "gs"; eps <- 1}
			delta <- g
		}

		u <- rep.int(0, q*k)
		u[s] <- delta
		delta <- matrix(u, q,k)	
		delta[is.na(delta)] <- 0
		cond <- FALSE
		while(!cond){
			new.theta <- theta - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			p.star.y <- p.bisec(new.theta, y,X, bfun)
			if(type == "ctiqr"){p.star.z <- pmin(p.bisec(new.theta, z,X, bfun), p.star.y)}
			G1 <- ee(new.theta, y,z,d,X,weights, bfun, p.star.y, p.star.z, J = FALSE)
			g1 <- G1$g[s]
			cond <- (sum(g1^2) < sum(g^2))
			eps <- eps*0.5
		}

		if(conv){break}
		g <- g1
		G <- G1
		theta <- new.theta
		h <- ee(theta, y,z,d,X,weights, bfun, p.star.y, p.star.z, J = TRUE, G = G)$J[s,s, drop = FALSE]
		if(i > 1){eps <- min(eps*10,1)}
		else{eps <- min(eps*10,0.1)}
	}

	p.star.y <- p.bisec(theta, y,X, bfun)
	p.star.z <- (if(type == "ctiqr") pmin(p.bisec(theta, z,X, bfun), p.star.y - 1e-6) else NULL)
	list(coefficients = matrix(theta, q, k),
		converged = (i < maxit), n.it = i, CDF.y = p.star.y, CDF.z = p.star.z,
		ee = g, jacobian = h, rank = (alg == "nr")*sum(s))
}


iqr.waldtest <- function(obj){
	bfun <- attr(obj$mf, "bfun")
	ax <- attr(obj$mf, "assign")
	ap <- attr(bfun, "assign")
	theta <- obj$coef
	q <- nrow(theta)
	k <- ncol(theta)
	s <- obj$s
	cc <- obj$covar
	ind.x <- rep(ax,k)
	ind.p <- sort(rep.int(ap,q))
	K <- tapply(rowSums(s), ax, sum)
	Q <- tapply(colSums(s), ap, sum)

	testx <- testp <- NULL

	if(q > 1){
		for(i in unique(ax)){
			theta0 <- c(theta[which(ax == i),])
			w <- which(ind.x == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]
			
			if(length(theta0) == 0){tx <- NA}
			else{tx <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testx <- c(testx, tx)
		}
		testx <- cbind(testx, df = K, pchisq(testx, df = K, lower.tail = FALSE))
		colnames(testx) <- c("chi-square", "df", "P(> chi)")

		nx <- attr(attr(obj$mf, "terms"), "term.labels")
		if(attr(attr(obj$mf, "terms"), "intercept") == 1){nx <- c("(Intercept)", nx)}
		rownames(testx) <- nx
	}
	if(k > 1){
		for(i in unique(ap)){
			theta0 <- c(theta[,which(ap == i)])
			w <- which(ind.p == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]

			if(length(theta0) == 0){tp <- NA}
			else{tp <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testp <- c(testp, tp)
		}
		testp <- cbind(testp, df = Q, pchisq(testp, df = Q, lower.tail = FALSE))
		colnames(testp) <- c("chi-square", "df", "P(> chi)")
		np <- attr(bfun, "term.labels")
		if(any(ap == 0)){np <- c("(Intercept)", np)}
		rownames(testp) <- np
	}

	list(test.x = testx, test.p = testp)
}


# kaplan-meier estimator for ct data. Assumes 0 <= y <= 1
# and returns an interpolating spline that computes the quantile function.

km <- function(z,y,d,w, type = c("ciqr", "ctiqr")){

	if(type == "ciqr"){m <- survfit(Surv(y,d) ~ 1, weights = w)}
	else{m <- survfit(coxph(Surv(z,y,d) ~ 1), type = "kaplan-meier", weights = w)}
	surv <- c(1,m$surv,0)
	time <- c(0,m$time,1)
	list(
    		Q = splinefun(1 - surv, time, method = "hyman"),
		F = splinefun(time, 1 - surv, method = "hyman")
	)
}




#' @export
test.fit <- function(object, R = 100){


	# ks and cvm statistic for uniformity with (possibly) censored and truncated data
	test.unif.ct <- function(z,y,d,w, type = c("iqr", "ciqr", "ctiqr")){
		n <- length(y)
		o <- order(y)
		y <- y[o]
		if(type != "iqr"){d <- d[o]}
		if(type == "ctiqr"){z <- z[o]}
		if(missing(w) || is.null(w)){w <- rep.int(1,n)}
		else{w <- w[o]}
		if(type == "iqr"){
			W <- cumsum(w)
			hat.Fy <- W/W[n]
			Fy <- y # = punif(y)
		}
		else{
			fit <- (if(type == "ctiqr") 
				survfit(coxph(Surv(z, y, d) ~ 1), type = "kaplan-meier", weights = w)
				else
				survfit(Surv(y, d) ~ 1, weights = w)
			)
			hat.Fy <- 1 - fit$surv
			Fy <- fit$time # = punif(fit$time)
			n <- length(Fy)
		}
	
		# kolmogorov - smirnov
	
		DD <- Fy - hat.Fy
		ks <- max(abs(DD))
	
		# cramer - von mises

		Fy <- c(0, Fy, 1)
		hat.Fy <- c(0, hat.Fy, 1)
		y <- c(0,y,1)
		U <- (hat.Fy - Fy)^2
		n <- n + 2
		h <- y[2:n] - y[1:(n-1)]
		b1 <- U[1:(n-1)]
		b2 <- U[2:n]
		A <- (b1 + b2)*h/2
		cvm <- sum(A)

		###
	
		c(ks = ks, cvm = cvm)	
	}


	# main function

	s <- object$s
	n <- N <- nrow(mf <- object$mf)
	rho <- rep.int(1,n)
	type <- attr(mf, "type")
	bfun <- bfun.b <- attr(mf, "bfun0")
	theta <- attr(mf, "theta")
	maxit <- 10 + 2*sum(!is.na(theta))
	V <- attr(mf, "all.vars")  
	x <- V$X; y <- V$y; z <- V$z0; d <- V$d; w <- V$weights
	CDFs <- attr(mf, "CDFs")
	br <- max(15, min(25, ceiling(n/50)))


	if(type == "ciqr"){mc <- suppressWarnings(pch:::pch.fit(y = y, 
        d = 1 - d, x = cbind(1,x), w = w, breaks = br))}
	else if(type == "ctiqr"){
		z <- pmax(z, min(z[z != - Inf]))
		mz <- suppressWarnings(pch:::pch.fit(z = -y, y = -z, 
        d = rep.int(1,n), x = cbind(1,x), w = w, breaks = br))
		mc <- suppressWarnings(pch:::pch.fit(z = z, y = y, 
        d = 1 - d, x = cbind(1,x,z), w = w, breaks = br))
		alphax <- 0
		B <- 250
		for(i in 1:B){
			beta <- apply_bfun(bfun, runif(n), "bfun", convert = FALSE)%*%t(theta)
			t.b <- rowSums(beta*x)
			z.b <- -pch:::sim.pch(mz, method = "s")
			c.b <- pch:::sim.pch(mc, x = cbind(1,x,z.b), method = "s")
			y.b <- pmin(t.b, c.b)
			alphax <- alphax + (y.b > z.b)
		}

		alphax <- pmax(alphax/B, 0.025) # pr(y > z | x), bounded to avoid too big N
		N <- round(n/mean(alphax)) # n. of observations to sample in order to obtain n final obs.
		rho <- 1/alphax # inclusion probability of each obs.
	}

	test0 <- test.unif.ct(z = CDFs$CDF.z, y = CDFs$CDF.y, d = d, w = w, type = type)
	test <- matrix(,R,2)

	for(b in 1:R){

		# x

		id <- sample.int(n, size = N, replace = TRUE, prob = rho)
		xb <- x[id,, drop = FALSE]
		wb <- w[id]

		# t
		beta <- apply_bfun(bfun, runif(N), "bfun", convert = FALSE)%*%t(theta)
		tb <- yb <- rowSums(beta*xb)

		# c,y,d
		if(type == "ciqr"){
			cb <- pch:::sim.pch(mc, x = cbind(1,xb), method = "s")
			yb <- pmin(cb, tb)
			db <- (tb <= cb)
		}
		if(type == "ctiqr"){
			zb <- -pch:::sim.pch(mz, x = cbind(1,xb), method = "s")
			cb <- pch:::sim.pch(mc, x = cbind(1,xb,zb), method = "s")
			yb <- pmin(cb, tb)
			db <- (tb <= cb)
			nb <- length(obs <- which(yb > zb))
			yb <- yb[obs]; zb <- zb[obs]; wb <- wb[obs]
			db <- db[obs]; xb <- xb[obs,, drop = FALSE]
			if(class(bfun) == "slp.basis"){
				bfun.b$BB1 <- matrix(rep(bfun$BB1[1,], each = nb), nb)
			}
			else{
			  bfun.b <- bfun[-length(bfun)]
			  bfun.b$BB1 <- sapply(bfun.b, function(x){x$BB1})
			  bfun.b$BB1 <- matrix(rep(bfun.b$BB1, each = nb), nb)
			  attributes(bfun.b) <- attributes(bfun)
			}
		}

		# fit the model

		fit.b <- iqr.newton(theta = theta, y = yb, z = zb, d = db, X = xb, 
			weights = wb, bfun = bfun.b, type = type, s = s, maxit = maxit, tol = 1e-5)
		test[b,] <- test.unif.ct(z = fit.b$CDF.z, y = fit.b$CDF.y, d = db, w = wb, type = type)

	}
  
	out <- cbind(test0*c(1,n), c(mean(test[,1] >= test0[1]), mean(test[,2] >= test0[2])))
	rownames(out) <- c("Kolmogorov-Smirnov", "Cramer-Von Mises")
	colnames(out) <- c("statistic", "p-value")
	out
}




#' @export
terms.iqr <- function(x, ...){attr(x$mf, "terms")}
#' @export
model.matrix.iqr <- function(object, ...){
  mf <- object$mf
  mt <- terms(mf)
  model.matrix(mt, mf)
}
#' @export
vcov.iqr <- function(object, ...){object$covar}
#' @export
nobs.iqr <- function(object, ...){nrow(object$mf)}

