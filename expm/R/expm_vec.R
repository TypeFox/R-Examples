#### Originally by  Roger B. Sidje (rbs@maths.uq.edu.au)
####	EXPOKIT: Software Package for Computing Matrix Exponentials.
####   ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

##' Performs exp(A t) %*% v  directly  w/o  evaluating exp(A)
##' Originally by  Roger B. Sidje
##'    EXPOKIT: Software Package for Computing Matrix Exponentials.
##'    ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
##' @title Compute   exp(A t) %*% v   directly
##' @param A n x n matrix
##' @param v n - vector
##' @param t number (scalar) ___ FIXME allow vector ? ___
##' @param tol
##' @param btol
##' @param m.max integer constants you should only change if you know what you're doing
##' @param mxrej
##' @param verbose flag indicating if the algorithm should be verbose..
##' @return a list with components
##' @author Ravi Varadhan, Johns Hopkins University; Martin Maechler (cosmetic)
expAtv <- function(A, v, t=1,
                   method = "Sidje98",
                   ## currently only one method, with these arguments:
                   ## FIXME argMeth=list( ... )
                   rescaleBelow = 1e-6,
                   tol=1e-7, btol = 1e-7, m.max = 30, mxrej = 10,
		   verbose = getOption("verbose"))
{
    ## R translation:  Ravi Varadhan, Johns Hopkins University
    ##		   "cosmetic", apply to sparse A: Martin Maechler, ETH Zurich
    if(length(d <- dim(A)) != 2) stop("'A' is not a matrix") # <- also for sparseMatrix
    stopifnot(length(v) == (n <- d[1]), m.max >= 2)
    if(n <= 1) {
	if(n == 1) return(list(eAtv = exp(A*t)*v, error = 0, nstep = 0L, n.reject = 0L))
	stop("nrow(A) must be >= 1")
    }
    method <- match.arg(method)
    m <- min(n, m.max)# >= 2
    ##-<FIXME> these are function arguments as well :
    gamma <- 0.9
    delta <- 1.2
    ##-</FIXME>
    nA <- norm(A, "I")
    if(nA < rescaleBelow) { ## rescaling, by MMaechler, needed for small norms
	A <- A/nA
	t <- t*nA
	nA <- 1
    }
    rndoff <- nA * .Machine$double.eps

    t_1 <- abs(t)
    sgn <- sign(t)
    t_now <- 0
    s_error <- 0
    k1 <- 2
    mb <- m
    xm <- 1/m
    beta <- sqrt(sum(v*v))# = norm(v) = |\ v ||
    if(beta == 0) ## border case: v is all 0, and the result is too
	return(list(eAtv = v, error = 0L, nstep = 0L, n.reject = 0L))
    fact <- (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1))
    myRound <- function(tt) {
	s <- 10^(floor(log10(tt)) - 1)
	ceiling(tt/s)*s
    }
    t_new <- myRound( (1/nA)*(fact*tol/(4*beta*nA))^xm )

    V <- matrix(0, n, m+1)
    H <- matrix(0, m+2, m+2)
    nstep <- n.rej <- 0L
    w <- v
    while (t_now < t_1) {
	nstep <- nstep + 1L
	t_step <- min(t_1 - t_now, t_new)
	if(verbose) cat(sprintf("while(t_now = %g < ..): nstep=%d, t_step=%g\n",
				t_now, nstep, t_step))
	V[,1] <- (1/beta)*w
	for (j in 1:m) {
	    p <- as.vector(A %*% V[,j])
	    for (i in 1:j) {
		H[i,j] <- s <- sum(V[,i] *  p)
		p <- p - s * V[,i]
	    }
	    s <- sqrt(sum(p*p))
	    if (s < btol) {
		k1 <- 0
		mb <- j
		t_step <- t_1 - t_now
		break
	    }
	    H[j+1, j] <- s
	    V[, j+1] <- p / s
	} ## j-loop complete
	if (k1 != 0) {
	    H[m+2, m+1] <- 1
	    av <- A %*% V[, m+1]
	    avnorm <- sqrt(sum(av * av))
	}
	i.rej <- 0L
	while (i.rej <= mxrej) {
	    mx <- mb + k1; imx <- seq_len(mx) # = 1:mx
	    if(verbose) cat(sprintf("	inner while: k1=%d -> mx=%d\n",
				    k1, mx))
	    F <- expm(sgn * t_step * H[imx,imx, drop=FALSE])
	    if (k1 == 0) {
		err_loc <- btol
		break
	    } else {
		phi1 <- abs(beta * F[m+1,1])
		phi2 <- abs(beta * F[m+2,1] * avnorm)
		if(is.nan(phi1) || is.nan(phi2))
		    stop("NaN phi values; probably overflow in expm()")
		if (phi1 > 10*phi2) {
		    err_loc <- phi2
		    xm <- 1/m
		} else if (phi1 > phi2) {
		    err_loc <- (phi1 * phi2)/(phi1 - phi2)
		    xm <- 1/m
		} else {
		    err_loc <- phi1
		    xm <- 1/(m-1)
		}
	    }
	    if (err_loc <= delta * t_step*tol) break
	    else {
		if (i.rej == mxrej)
		    stop(gettextf('The requested tolerance (tol=%g) is too small for mxrej=%d.',
				  tol, mxrej))
		t_step <- gamma * t_step * (t_step * tol / err_loc)^xm
		s <- 10^(floor(log10(t_step))-1)
		t_step <- s * ceiling(t_step / s)
		i.rej <- i.rej + 1L
	    }
	}## end{ while (i.rej < mx..) }
	n.rej <- n.rej + i.rej
	mx <- mb + max(0, k1-1); imx <- seq_len(mx) # = 1:mx
	w <- as.vector(V[, imx] %*% (beta*F[imx,1, drop=FALSE]))
	beta <- sqrt(sum(w*w))
	t_now <- t_now + t_step
	t_new <- myRound(gamma * t_step * (t_step*tol/err_loc)^xm)
	err_loc <- max(err_loc, rndoff)
	s_error <- s_error + err_loc
    }# end{ while }
    list(eAtv = w, error = s_error, nstep = nstep, n.reject = n.rej)
}
