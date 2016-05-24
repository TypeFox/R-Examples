nmk <-
function(par, fn, control=list(), ...) {
    ctrl <- list(tol=1.e-06, maxfeval = min(5000, max(1500, 20*length(par)^2)), regsimp=TRUE, maximize=FALSE, restarts.max=3, trace=FALSE)
	namc <- match.arg(names(control), choices=names(ctrl), several.ok=TRUE)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if (!is.null(names(control))) ctrl[namc] <- control
    ftol <- ctrl$tol
    maxfeval <- ctrl$maxfeval
    regsimp <- ctrl$regsimp
    restarts.max <- ctrl$restarts.max
    maximize <- ctrl$maximize
    trace <- ctrl$trace

	if (maximize) fnm <- function(par, ...) -fn(par, ...) else fnm <- function(par, ...) fn(par, ...) 

	x0 <- par
	n <- length(par)
	if (n == 1) stop(call. = FALSE, "Use `optimize' for univariate optimization")
 	if (n > 30) warning("Nelder-Mead should not be used for high-dimensional optimization")

	V <- cbind(rep(0, n), diag(n))
	f <- rep(0, n+1)
	f[1] <- fnm(x0, ...)
	V[, 1] <- x0
	scale <- max(1, sqrt(sum(x0^2)))

	if (regsimp) {
		alpha <- scale / (n * sqrt(2)) * c(sqrt(n+1) + n - 1, sqrt(n+1) -1)
		V[, -1] <- (x0 + alpha[2])
		diag(V[, -1]) <- x0[1:n] + alpha[1]
		for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...) 
	} else {
		V[, -1] <- x0 + scale * V[, -1] 
		for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...) 
	}

	f[is.nan(f)] <- Inf

	nf <- n + 1
	ord <- order(f)
	f <- f[ord]
	V <- V[, ord]
	
	rho <- 1
	gamma <- 0.5
	chi <- 2
	sigma <- 0.5
	conv <- 1
	oshrink <- 0
	restarts <- 0
	orth <- 0
	dist <- f[n+1] - f[1]
	
	v <- V[, -1] - V[, 1]
	delf <- f[-1] - f[1]
	diam <- sqrt(colSums(v^2))
	sgrad <- c(crossprod(t(v), delf))
	alpha <- 1.e-04 * max(diam) / sqrt(sum(sgrad^2))
	simplex.size <- sum(abs(V[, -1] - V[, 1])) / max(1, sum(abs(V[, 1])))

	itc <- 0
	conv <- 0
	message <- "Succesful convergence"

	while (nf < maxfeval & restarts < restarts.max & dist > ftol & simplex.size > 1.e-06) {

		fbc <- mean(f)
		happy <- 0
		itc <- itc + 1
		xbar <- rowMeans(V[, 1:n])
		xr <- (1 + rho) * xbar - rho * V[, n+1]
		fr <- fnm(xr, ...)
		nf <- nf + 1
		if(is.nan(fr)) fr <- Inf

		if (fr >= f[1] & fr < f[n]) {
			happy <- 1
			xnew <- xr
			fnew <- fr
		} else if (fr < f[1]) {
			xe <- (1 + rho * chi) * xbar - rho * chi * V[, n+1]
			fe <- fnm(xe, ...)
		if(is.nan(fe)) fe <- Inf
			nf <- nf + 1
			if (fe < fr) {
				xnew <- xe
				fnew <- fe
				happy <- 1
			} else {
				xnew <- xr
				fnew <- fr
				happy <- 1
			}
		} else if (fr >= f[n] & fr < f[n+1]) {
			xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, n+1]
			fc <- fnm(xc, ...)
		if(is.nan(fc)) fc <- Inf
			nf <- nf + 1
			if (fc <= fr) {
				xnew <- xc
				fnew <- fc
				happy <- 1
			}
		} else if (fr >= f[n+1]) {
			xc <- (1 - gamma) * xbar + gamma * V[, n+1]
			fc <- fnm(xc, ...)
		if(is.nan(fc)) fc <- Inf
			nf <- nf + 1
			if (fc < f[n+1]) {
				xnew <- xc
				fnew <- fc
				happy <- 1
			}
		}

		if (happy == 1 & oshrink == 1) {
			fbt <- mean(c(f[1:n], fnew))
			delfb <- fbt - fbc
			armtst <- alpha * sum(sgrad^2)
			if (delfb > - armtst/n) {
			if (trace) cat("Trouble - restarting: \n")
				restarts <- restarts + 1
				orth <- 1
				diams <- min(diam)
				sx <- sign(0.5 * sign(sgrad))
				happy <- 0
				V[, -1] <- V[, 1]
				diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
			}
		}

		if (happy == 1) {
			V[, n+1] <- xnew
			f[n+1] <- fnew
			ord <- order(f)
			V <- V[, ord]
			f <- f[ord]
		} else if (happy == 0 & restarts < restarts.max) {
			if (orth == 0) orth <- 1
			V[, -1] <- V[, 1] - sigma * (V [, -1] - V[, 1])
			for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...)  ## kmm change
			nf <- nf + n
			ord <- order(f)
			V <- V[, ord]
			f <- f[ord]
		}

		v <- V[, -1] - V[, 1]
		delf <- f[-1] - f[1]
		diam <- sqrt(colSums(v^2))
		simplex.size <- sum(abs(v)) / max(1, sum(abs(V[, 1])))

		f[is.nan(f)] <- Inf

		dist <- f[n+1] - f[1]
		sgrad <- c(crossprod(t(v), delf))
		if (trace & !(itc %% 2)) cat("iter: ", itc, "\n", "value: ", f[1], "\n")
	}

	if (dist <= ftol | simplex.size <= 1.e-06) {
		conv <- 0
		message <- "Successful convergence"
		} else if (nf >= maxfeval) {
		conv <- 1
		message <- "Maximum number of fevals exceeded"
		} else if (restarts >= restarts.max) {
		conv <- 2
		message <- "Stagnation in Nelder-Mead"
		}	

	return(list(par = V[, 1], value=f[1]*(-1)^maximize, feval=nf, restarts=restarts, convergence=conv, message=message))
}

