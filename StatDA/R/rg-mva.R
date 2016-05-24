"rg.mva" <- 
function(x, main = deparse(substitute(x)))
{
	# Procedure to undertake non-robust multivariate data analyses; the object
	# generated is identical to that of rg.robmva so that the saved list may be
	# passed to other rotation and display functions.  Thus weights are set to
	# 1, and other variables are set to appropriate defaults.  The estimation
	# of Mahalanobis distances is only undertaken if x is non-singular, i.e., 
	# the lowest eigenvalue is > 10e-4.
	#
	# Note this procedure uses svd() rather than the classic solve().
	#
	# PCA output may be plotted with rg.rqpca.plot() and rg.rqpca.screeplot(), 
	# and Mahalanobis distances may be plotted with rg.md.plot().
	#
	# Determine the length of the vectors
	if(!is.matrix(x)) stop("Not a Matrix")
	n <- length(x[, 1])
	p <- length(x[1,  ])
	matnames <- dimnames(x)
	wts <- numeric(n)
	wts[1:n] <- 1
	nc <- n
	cat("  n =", n, "\tnc =", n, "\tp =", p, "\t\tnc/p =", round(nc/p, 2), "\n")
	if(nc <= 5 * p)
		cat("  *** Proceed with Care, n is < 5p ***\n")
	if(nc <= 3 * p)
		cat("  *** Proceed With Great Care, n = ", n, ", which is < 3p ***\n")
	# Compute means & SDs, and standardize the data set.  Note cov.wt() is used in order
	# to have cov and r entries for the saved object
	save <- cov.wt(x, wt = wts, cor = TRUE)
	xmean <- save$center
	xsd <- sqrt(diag(save$cov))
	# Compute SNDs
	temp <- sweep(x, 2, xmean, "-")
	snd <- sweep(temp, 2, xsd, "/")
	# The following procedure duplicates that in rg.rqpca() for iwght = "r"
	# Standardize centred x for R-mode PCA and compute
	xsd2 <- sqrt(n) * xsd
	w <- sweep(temp, 2, xsd2, "/")
	wt <- t(as.matrix(w))
	a <- wt %*% as.matrix(w)
	b <- svd(a)
	cat("  Eigenvalues:", signif(b$d, 4), "\n")
	sumc <- sum(b$d)
	econtrib <- 100 * (b$d/sumc)
	rqscore <- w %*% b$v
###	vcontrib <- colVars(rqscore)
        vcontrib <- apply(rqscore,2,var)
	sumv <- sum(vcontrib)
	pvcontrib <- (100 * vcontrib)/sumv
	cpvcontrib <- cumsum(pvcontrib)
	b1 <- b$v * 0
	diag(b1) <- sqrt(b$d)
	rload <- b$v %*% b1
	rcr <- rload[,  ] * 0
	rcr1 <- apply(rload^2, 1, sum)
	rcr <- 100 * sweep(rload^2, 1, rcr1, "/")
	# Test for non-singularity and compute Mahalanobis distances
	if(b$d[p] > 0.001) {
		md <- mahalanobis(x, save$center, save$cov)
		temp <- (nc - p)/(p * (nc + 1))
		ppm <- 1 - pf(temp * md, p, nc - p)
		epm <- 1 - pchisq(md, p)
	}
	else {
		cat("  Lowest eigenvalue < 10^-4, Mahalanobis distances not computed\n")
		md <- NULL
		ppm <- NULL
		epm <- NULL
	}
	invisible(list(main = main, input = deparse(substitute(x)), proc = "cov", n = n, nc = nc,
		p = p, matnames = matnames, wts = wts, mean = xmean, cov = save$cov, sd = xsd,
		snd = snd, r = save$cor, eigenvalues = b$d, econtrib = econtrib, eigenvectors = 
		b$v, rload = rload, rcr = rcr, rqscore = rqscore, vcontrib = vcontrib, pvcontrib
		 = pvcontrib, cpvcontrib = cpvcontrib, md = md, ppm = ppm, epm = epm, nr = NULL)
		)
}

