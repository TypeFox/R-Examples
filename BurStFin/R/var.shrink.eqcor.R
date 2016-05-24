"var.shrink.eqcor" <-
function (x, weights=seq(0.5, 1.5, length=nt), shrink=NULL, center=TRUE, 
	vol.shrink=0, sd.min=20, quan.sd=0.9, tol=1e-4, compatible=FALSE,
	verbose=2) 
{
	fun.copyright <- "Placed in the public domain 2009-2014 by Burns Statistics Ltd."
	fun.version <- "var.shrink.eqcor 004"

	x <- as.matrix(x)
	nassets <- ncol(x)
	if(nassets < 2) stop("'x' must have at least 2 columns")
	nt <- nrow(x)
	if(nt < 2) stop("'x' must have at least 2 rows")

	# for use in finance, try to check if prices rather than returns
	if(verbose >= 1 && min(x, na.rm=TRUE) >= 0) {
		warning(paste("minimum value in 'x' is",
			min(x, na.rm=TRUE), "are you giving a price",
			"matrix rather than a return matrix?",
			"(warning suppressed with verbose < 1)"))
	}

	if(is.null(weights)) weights <- rep(1, nt)
	if(any(weights < 0)) {
		stop(paste(sum(weights < 0), "negative value(s) in 'weights'"))
	}
	if(length(weights) != nt) {
		if(length(weights) == 1 && weights > 0) {
			weights <- rep(1, nt)
		} else {
			stop(paste("bad length for 'weights' argument",
				"-- length must be the number of rows of 'x'",
				"or a single positive number (meaning",
				"equal weighting)"))
		}
	}
	if(any(weights == 0)) {
		x <- x[weights > 0, , drop=FALSE]
		nt <- nrow(x)
		weights <- weights[weights > 0]
	}

	if(!is.numeric(sd.min) || length(sd.min) != 1) {
		stop(paste("'sd.min' should be a single number -- given",
			"has mode", mode(sd.min), "and length", 
			length(sd.min)))
	} else if(sd.min < 2) {
		stop("'sd.min' must be at least 2")
	}
	if(!is.numeric(quan.sd) || length(quan.sd) != 1) {
		stop(paste("'quan.sd' should be a single number -- given",
			"has mode", mode(quan.sd), "and length", 
			length(quan.sd)))
	}
	nobs <- nt - colSums(is.na(x))
	good <- nobs >= min(sd.min, nt)
	if(sum(good) < 2) {
		stop(paste("not enough columns with enough data",
			"-- fewer than 2 columns with at least", sd.min,
			"non-missing observations"))
	}
	weights <- weights / mean(weights)

	if(!is.numeric(tol) || length(tol) != 1) {
		stop(paste("bad value for 'tol' which should be a",
			"single numeric value -- given has mode",
			mode(tol), "and length", length(tol)))
	}

	if(is.logical(center)) {
		if(center) {
			center <- colMeans(x * weights, na.rm=TRUE)
		} else {
			center <- rep(0, nassets)
		}
	} else if(length(center) != nassets || !is.numeric(center)) {
		stop(paste("bad value for 'center' which should be",
			"either a single logical value or a numeric vector",
			"of length", nassets, "-- given has mode",
			mode(center), "and length", length(center)))
	}
	x <- sweep(x, 2, center, "-")
	svar <- var(x * sqrt(weights), use='pairwise')
	if(compatible) {
		# to match Ledoit-Wolf code, do rescaling as in next line
		svar <- svar * (nt - 1) / nt
		if(diff(range(weights)) > 1e-7) {
			warning("in compatible mode but time weighting")
		}
	}
	svar.sup <- svar[good, good, drop=FALSE]
	sdiag <- sqrt(diag(svar))
	sdiag[!good] <- quantile(sdiag[good], quan.sd)
	if(vol.shrink > 0) {
		if(vol.shrink > 1) vol.shrink <- 1
		meanvol <- mean(sdiag)
		sdiag <- vol.shrink * meanvol + (1 - vol.shrink) * sdiag
	}
	meancor <- mean(cov2cor(svar.sup)[which(lower.tri(svar.sup, 
		diag=FALSE))])
	prior <- svar
	prior[] <- meancor
	diag(prior) <- 1
	prior <- sdiag * prior * rep(sdiag, each=nassets)
	svar.sup <- which(is.na(svar))
	svar[svar.sup] <- prior[svar.sup]
	if(length(shrink) != 1) {
		if(length(shrink)) {
			if(!is.numeric(shrink)) {
				stop(paste("'shrink' should be either NULL",
					"or a single numeric value -- given",
					"has mode", mode(shrink), 
					"and length", length(shrink)))
			}
			warning(paste("'shrink' ignored, it has length", 
				length(shrink), "-- should be either NULL",
				"or a single numeric value"))
		}
		gamma <- sum((svar - prior)^2)
		x[which(is.na(x))] <- 0
		pi.mat <- theta1.mat <- theta2.mat <- array(0, 
			c(nassets, nassets))
		for(i in 1:nt) {
			this.wt <- weights[i]
			this.cross <- crossprod(x[i, , drop=FALSE]) - svar
			pi.mat <- pi.mat + this.wt * this.cross^2
			theta1.mat <- theta1.mat + this.wt * 
				diag(this.cross) * this.cross
			theta2.mat <- theta2.mat + rep(diag(this.cross), 
				each=nassets) * this.cross * this.wt
		}
		pi.hat <- sum(pi.mat)/nt
		theta1.mat <- theta1.mat / nt
		theta2.mat <- theta2.mat / nt
		diag(theta1.mat) <- 0
		diag(theta2.mat) <- 0
		theta1.mat <- rep(sdiag, each=nassets) * theta1.mat / sdiag
		theta2.mat <- rep(1/sdiag, each=nassets) * theta2.mat * sdiag
		rho <- sum(diag(pi.mat)) / nt + meancor * 0.5 * 
			(sum(theta1.mat) + sum(theta2.mat))
		shrink <- (pi.hat - rho) / gamma / nt
		# allow memory clean up
		pi.mat <- theta1.mat <- theta2.mat <- this.cross <- NULL
	} else {
		if(!is.numeric(shrink)) {
			stop(paste("'shrink' should be either NULL or",
				"a single numeric value -- given has mode",
				mode(shrink), "and length", length(shrink)))
		}
		if(is.na(shrink)) {
			stop("missing value for 'shrink' argument")
		}
	}
	shrink <- min(1, max(0, shrink))
	
	# overwrite svar to save space
	svar <- shrink * prior + (1 - shrink) * svar

	if(tol > 0 || any(is.na(x))) {
		svar.sup <- eigen(svar, symmetric=TRUE)
		tol <- tol * max(svar.sup$values)
		if(min(svar.sup$values) < tol) {
			vals <- svar.sup$values
			vals[which(vals < tol)] <- tol
			svar[] <- svar.sup$vectors %*% (vals * t(svar.sup$vectors))
		}
	}
	attr(svar, "shrink") <- shrink
	attr(svar, "timestamp") <- date()
	svar
}

