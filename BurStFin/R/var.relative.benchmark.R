"var.relative.benchmark" <-
function (variance, benchmark)
{
        fun.copyright <- "Placed in the public domain 2009-2012  Burns Statistics Ltd."
	fun.version <- "var.relative.benchmark 004"

	subfun.varrel <- function(varmat, b.ind) {
		sf.ans <- varmat[-b.ind, -b.ind, drop=FALSE] + 
			varmat[b.ind, b.ind]
		bencor <- varmat[-b.ind, b.ind]
		sf.ans <- t(sf.ans - bencor) - bencor
		sf.ans
	}

	#
	# start of main function
	#

	dv <- dim(variance)
	ldv <- length(dv)

	if(ldv != 2 && ldv != 3) {
		stop("'variance' must be a matrix or 3D array")
	}
	if(dv[1] != dv[2]) stop("bad dimensions for 'variance'")

	vnam <- dimnames(variance)[[1]]
	if(!length(vnam)) {
		stop("no asset names for 'variance'")
	}
	if(any(dimnames(variance)[[2]] != vnam)) {
		stop("mismatch of 'variance' names in 1st and 2nd dimensions")
	}
	if(length(benchmark) != 1 || !is.character(benchmark) ||
			nchar(benchmark) == 0) {
		stop(paste("'benchmark' must be a single non-empty",
			"character string -- given has mode",
			mode(benchmark), "and length", length(benchmark)))
	}

	b.ind <- match(benchmark, vnam, nomatch=NA)
	if(is.na(b.ind)) {
		stop(paste("benchmark (", benchmark, 
			") not an asset in 'variance'", sep=""))
	}

	if(ldv == 2) {
		ans <- subfun.varrel(variance, b.ind)
	} else {
		ans <- array(NA, dv - c(1, 1, 0), list(vnam[-b.ind],
			vnam[-b.ind], dimnames(variance)[[3]]))
		for(i in 1:dv[3]) {
			ans[, , i] <- subfun.varrel(variance[, , i, drop=TRUE],
				b.ind)
		}
	}

	attr(ans, "call") <- match.call()
	ans
}

