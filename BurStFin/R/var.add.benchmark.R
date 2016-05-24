"var.add.benchmark" <-
function (variance, benchmark.weights, name="benchmark", sum.to.one=TRUE) 
{
        fun.copyright <- "Placed in the public domain 2009-2014 by  Burns Statistics Ltd."
	fun.version <- "var.add.benchmark 005"

	subfun.varadd <- function(varmat, lwt, p, pseq)
	{
		bcov <- drop(varmat %*% lwt)
		sf.ans <- array(NA, dim(varmat) + 1)
		sf.ans[pseq, pseq] <- varmat
		sf.ans[pseq, p+1] <- bcov
		sf.ans[p+1, pseq] <- bcov
		sf.ans[p+1, p+1] <- sum(lwt * bcov)
		sf.ans
	}

	#
	# start of main function
	#

	vnam <- dimnames(variance)[[1]]
	bnam <- names(benchmark.weights)
	if(!length(vnam)) {
		if(!length(bnam)) {
			stop(paste("need asset names for both",
				"'variance' and 'benchmark.weights'"))
		}
		stop("no asset names for 'variance'")
	} else if(!length(bnam)) {
		stop("no asset names for 'benchmark.weights'")
	}
	if(any(nchar(c(vnam, bnam)) == 0)) 
		stop("no asset name may be missing")
	if(any(duplicated(bnam))) {
		stop("duplicate names in 'benchmark.weights'")
	}
	if(length(unique(intersect(bnam, vnam))) < length(bnam)) {
		nmiss <- length(bnam) - length(unique(intersect(bnam, vnam)))
		stop(paste(nmiss, "asset(s) in 'benchmark.weights' are",
			"not in 'variance'"))
	}
	if(any(is.na(benchmark.weights))) {
		stop(paste(sum(is.na(benchmark.weights)),
			"missing value(s) in 'benchmark.weights'"))
	}
	if(sum.to.one && abs(sum(abs(benchmark.weights)) - 1) > 1e-10) {
		wsum <- sum(abs(benchmark.weights))
		warning(paste("absolute of 'benchmark.weights' sums to", wsum,
			" adjusting so it sums to 1",
			"-- use 'sum.to.one=FALSE' to avoid adjustment"))
		benchmark.weights <- benchmark.weights / wsum
	}

	dv <- dim(variance)
	ldv <- length(dv)
	if(ldv != 2 && ldv != 3) {
		stop("'variance' must be a matrix or 3D array")
	}
	p <- dv[1]
	if(dv[2] != p || any(dimnames(variance)[[2]] != vnam)) {
		stop(paste("second dimension of 'variance' does not match",
			"the first"))
	}
	if(any(is.na(variance))) {
		stop(paste(sum(is.na(variance)),
			"missing value(s) in variance"))
	}

	pseq <- 1:p
	lwt <- rep(0, p)
	names(lwt) <- vnam
	lwt[bnam] <- benchmark.weights
	vnamp <- c(vnam, name)

	if(ldv == 2) {
		ans <- subfun.varadd(variance, lwt, p, pseq)
		dimnames(ans) <- list(vnamp, vnamp)
	} else {
		ans <- array(NA, dv + c(1,1,0))
		for(i in 1:dv[3]) {
			ans[, , i] <- subfun.varadd(variance[,, i], lwt, p,
				pseq)
		}
		dimnames(ans) <- list(vnamp, vnamp, dimnames(variance)[[3]])
	}

	ans
}

