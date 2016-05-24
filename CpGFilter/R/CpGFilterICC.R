
CpGFilterICC <-
function(dat, rep.design, logit.transform=TRUE, verbose=TRUE) {
	# Computes the ICC for each probe based on fast LMM
	#
	# Args:
	#   dat: a matrix of CpG beta-values, row - CpG, column - sample
	#   rep.design: a vector indicating the replicate desgin, it could be factor, character or numeric vectors
	#      Example - c(1, 2, 3, 4, 4, 4, 5, 5) OR c('S1', 'S2', 'S2', 'S2', 'S1')
	#   logit.transform: If TRUE, beta-value will be converted into M-value;  Default is TRUE.
	#   verbose: If TRUE, print run information
	#
	# Returns:
	#   The ICCs for all probes

	ptm <- proc.time()

	if (verbose == TRUE) 	cat("Fast LMM for calculating ICC based on replicates\n")
    if (nrow(dat) < ncol(dat)) {
		warning("Less rows than columns. Rows should be CpGs! \n")
	}
	if (ncol(dat) != length(rep.design)) {
		stop("Number of columns in 'dat' is not equal to the length of 'rep.design'!\n")
	}
	
	if (sum(is.na(dat)) != 0) {
		warning("NAs detected in your data! NAs will be replaced by row means!\n")
		dat[is.na(dat)] <- rowMeans(dat, na.rm=T)[which(is.na(dat), arr.ind=T)[, 1]]
	}
	if (logit.transform == TRUE) {
		dat.range <- range(dat)
		if (dat.range[1] < 0 | dat.range[2] > 1) {
			warning("Data value range is not between 0 and 1! Logit transform will not be performed!\n")
		} else {
			dat <- B2M(dat)
		}
    }
	rep.design <- as.numeric(factor(rep.design))
	m <- length(unique(rep.design))
	n <- length(rep.design)
	
	if (verbose == TRUE) 	cat("The data has ", nrow(dat), " CpGs, ",  m, " unique samples and ", n-m, "replicates.\n")
	if (m < 40 | (n-m) < 10) {
		warning("Recommended sample size and replicate number are 40 and 8 at least!\n")
	}
	
	# order the samples
	dat <- dat[, order(rep.design), drop=F]
	rep.design <- sort(rep.design)
	
	# replicate number
	rep.no <- table(rep.design)
	rep.no <- rep.no[paste(1:m)]
	
	# replicates list
	rep.id.list <- lapply(1:m, function(i) which(rep.design == i))
		
	# Estimate sd and mean using all biological replicates
	indep.id <- sapply(rep.id.list, function(x) x[1])
	dat.temp <- t(dat[, indep.id, drop=F])
	dat.temp <- scale(dat.temp)
#	cpg.sd <- rowSds(dat[, indep.id, drop=F])
#	cpg.m <- rowMeans(dat[, indep.id, drop=F])
	
    cpg.sd <- attr(dat.temp, 'scaled:scale')
    cpg.m <- attr(dat.temp, 'scaled:center')
	
	gc(dat.temp)
	
	# retain samples with replicates
	dat <- dat[, unlist(rep.id.list[rep.no != 1]), drop=F]
	rep.design <- rep.design[unlist(rep.id.list[rep.no != 1])]
	rep.no <- rep.no[rep.no != 1]
	rep.id.list <- lapply(unique(rep.design), function(i) which(rep.design == i))
	m2 <- length(rep.no)
	
	# center the data
	dat <- (dat - cpg.m) / cpg.sd
	
	# create vectors for computation
	rep.no.cum <- c(0, cumsum(rep.no))
	exp.row <- unlist(lapply(1:m2, function(i) 
						expand.grid(1:rep.no[i], 1:rep.no[i])[, 1] + rep.no.cum[i]))
	exp.col <- unlist(lapply(1:m2, function(i) 
						expand.grid(1:rep.no[i], 1:rep.no[i])[, 2] + rep.no.cum[i]))
	exp.ind <- rep(1:m2, rep.no^2)
	
	temp <- c(0, cumsum(rep.no^2))
	ind <- unlist(lapply(1:m2, function(i) (1:rep.no[i]-1)*rep.no[i] + 1:rep.no[i] + temp[i]))
	mask.c0 <- rep(FALSE, length(exp.ind))
	mask.c0[ind] <- TRUE
	mask.c1 <- !(mask.c0)
	
	# Score function - vectorization to speed up computation
	f1 <- function(rho, y) {
		b0 <- (1 + (rep.no-2)*rho) / (1 - rho) / (1 + (rep.no-1)*rho)
		b1 <- -rho / (1 - rho) / (1 + (rep.no-1)*rho)
		c0 <- (rep.no - 1) * (rep.no - 2) * b1^2 + 2 * (rep.no - 1) * b0 * b1
		c1 <- ((rep.no - 1) * (rep.no - 2) + 1) * b1^2 + 2 * (rep.no - 2) * b0 * b1 + b0^2
		sum(rep.no * (rep.no - 1) * b1) / 2 -
				sum(y[exp.row] * (c0[exp.ind]*mask.c0 + c1[exp.ind]*mask.c1) * y[exp.col]) / 2

	}
	
	f2 <- function(y) {
		i <<- i + 1
		if (i %% 5000 == 0 & verbose == TRUE) 	cat(i, "\n")
	    a1 <- f1(0, y)
		if (a1 >= 0) {
			rho <- 0
		} else {
			# Root finding
			rho <- uniroot(f1, c(0, 0.9999), y, tol=1e-4)$root
		}
		rho
	}
	
	if (verbose == TRUE) 	cat("Start to fit fast LMM ...\n")
	i <- 0
	res <- apply(dat, 1, f2)
	
	if (verbose == TRUE) 	cat("Complete. ", (proc.time() - ptm)[3], "seconds used.\n")
	
	return(res)
}
