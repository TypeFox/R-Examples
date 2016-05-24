rlm <-
	function(formula, weights, acc = 1e-2)
{
	options("na.action"=stats::na.pass)

	mf <- match.call(expand.dots = FALSE)
	mf$acc <- NULL
	mf[[1L]] <- quote(stats::model.frame)
	mf <- eval.parent(mf)
	mt <- attr(mf, "terms")

	y <- stats::model.response(mf)
	y <- stats::na.pass(y)

	if(is.vector(y)) {
		N <- length(y)
		M <- 1
	} else if (is.matrix(y)){	
		nrow_y <- dim(y)[1]
		ncol_y <- dim(y)[2]
		N <- nrow_y
		M <- ncol_y
	} else {
		stop("y is neither vector nor matrix")	
	}


	if(N > 80) {
		stop("number of observations is too big")		
	}


	x <- stats::model.matrix(mt, mf, stats::contrasts)
	nrow_x <- dim(x)[1]
	ncol_x <- dim(x)[2]
	K <- ncol_x
	if(K > 4) {
		stop("number of regressors is too big")		
	}

	w <- stats::model.weights(mf)
	nrow_w <- length(w)


	est <- matrix(data = 0, nrow = M, ncol = K)

	out <-.C("rlm_cpu", NAOK = TRUE, as.double(y), as.double(x), as.double(w), est = as.double(est), 
		as.integer(N), as.integer(K), as.integer(M), as.double(acc))
	est = matrix (out$est, nrow = M, ncol = K, byrow = TRUE)

	est
}