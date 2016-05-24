"factor.model.stat" <-
function (x, weights=seq(0.5, 1.5, length.out=nobs), output="full", center=TRUE, 
	frac.var=.5, iter.max=1, nfac.miss=1, full.min=20, reg.min=40, 
	sd.min=20, quan.sd=.90, tol=1e-3, zero.load=FALSE, 
	range.factors=c(0, Inf), constant.returns.okay=FALSE,
	specific.floor=0.1, floor.type="quantile", verbose=2)
{
        fun.copyright <- "Placed in the public domain 2006-2014 by Burns Statistics Ltd."
	fun.version <- "factor.model.stat 014"

	subfun.ssd <- function(z, weights, sd.min) {
		nas <- is.na(z)
		if(any(nas)) {
			if(sum(!nas) < sd.min) return(NA)
			sum(weights[!nas] * z[!nas]^2) / sum(weights[!nas])
		} else {
			sum(weights * z^2)
		}
	}

	#
	# start of main function
	#

	x <- as.matrix(x)
	if(!is.numeric(x)) stop("'x' needs to be numeric")
	x[!is.finite(x)] <- NA

	# for use in finance, try to check it is returns and not prices
	if(verbose >= 1 && min(x, na.rm=TRUE) >= 0) {
		warning(paste("minimum of values in 'x' is",
			min(x, na.rm=TRUE), "are you giving a price",
			"matrix rather than a return matrix?",
			"(warning suppressed if verbose < 1)"))
	}

	xna <- is.na(x)
	allmis <- rowSums(xna) == ncol(x)
	if(any(allmis)) {
		x <- x[!allmis, , drop=FALSE]
		xna <- is.na(x)
	}
	num.mis <- colSums(xna)

	if(any(num.mis > 0)) {
		if(sum(num.mis == 0) < full.min) 
			stop("not enough columns without missing values")
		if(!length(dimnames(x)[[2]])) 
			stop("'x' needs column names when missing values exist")

		max.miss <- max(num.mis)
		lnfm <- length(nfac.miss)
		if(lnfm == 0) stop("'nfac.miss' must have positive length")
		nfac.miss <- round(nfac.miss)
		if(any(nfac.miss < 0)) 
			stop("negative values in 'nfac.miss'")
		if(lnfm < max.miss) {
			nfac.miss <- c(nfac.miss, rep(nfac.miss[lnfm],
				max.miss - lnfm))
		}
	}

	if(!is.character(output) || length(output) != 1) {
                stop(paste("'output' should be a single character string",
                        "-- given has mode", mode(output), "and length",
                        length(output)))
        }
        output.menu <- c("full", "factor", "systematic", "specific")
        output.num <- pmatch(output, output.menu, nomatch=0)
        if(output.num == 0) {
                stop(paste("unknown or ambiguous input for 'output'",
                        "-- the allowed choices are:",
                        paste(output.menu, collapse=", ")))
        }
        output <- output.menu[output.num]

	nassets <- ncol(x)
	nobs <- nrow(x)
	if(is.null(weights)) {
		weights <- rep(1, nobs)
	} else if(!is.numeric(weights)) {
		stop(paste("'weights' must be numeric -- given has mode",
			mode(weights), "and length", length(weights)))
	}
	if(length(weights) != nobs) {
		if(length(weights) == nobs + sum(allmis)) {
			weights <- weights[!allmis]
		} else if(length(weights) == 1 && weights > 0) {
			weights <- rep(1, nobs)
		} else {
			stop(paste("bad value for 'weights'",
				"-- must be a single positive number",
				"(meaning equal weighting) or have length",
				"equal to the number of observations"))
		}
	}
	if(any(weights < 0)) {
		stop(paste(sum(weights < 0), "negative value(s) in 'weights'"))
	}
	weights <- weights / sum(weights)

	if(is.logical(center)) {
		if(center) {
			center <- colSums(x * weights, na.rm=TRUE)
		} else {
			center <- rep(0, nassets)
		}
	} else if(length(center) != nassets) stop("wrong length for 'center'")
	x <- sweep(x, 2, center, "-")

	sdev <- sqrt(apply(x, 2, subfun.ssd, weights=weights, sd.min=sd.min))
	sdzero.names <- NULL
	sdzero <- FALSE
	if(any(sdev <= 0, na.rm=TRUE)) {
		sdzero <- !is.na(sdev) & sdev <= 0
		sdzero.names <- dimnames(x)[[2]][sdzero]
		if(constant.returns.okay) {
		    sdev[which(sdzero)] <- 1e-16
		    if(verbose >= 1) {
			warning(paste(sum(sdzero),
				"asset(s) with constant returns:",
				paste(sdzero.names, collapse=", "),
				"(warning suppressed with verbose < 1)"))
		    }
		} else {
			stop(paste(sum(sdzero),
				"asset(s) with constant returns:",
				paste(dimnames(x)[[2]][sdzero], collapse=", ")))
		}
	}
	if(any(is.na(sdev))) {
		sdev[is.na(sdev)] <- quantile(sdev, quan.sd, na.rm=TRUE)
	}
	x <- scale(x, scale=sdev, center=FALSE)

	x <- sqrt(weights) * x   # x is now weighted
	fullcolnames <- dimnames(x)[[2]][num.mis == 0]
	fullcols <- which(num.mis == 0)
	decomp <- try(svd(x[, fullcols, drop=FALSE], nu=0), silent=TRUE)
	if(inherits(decomp, "try-error")) {
		# presumably Lapack error, failing to converge
		rever <- nrow(x):1
		decomp <- svd(x[rever, fullcols, drop=FALSE])
		decomp$v <- decomp$v[rever, , drop=FALSE]
	}
	svdcheck <- colSums(decomp$v^2)
	if(any(abs(svdcheck - 1) > 1e-3)) {
		stop("bad result from 'svd'")
	}
	cumvar <- cumsum(decomp$d^2) / sum(decomp$d^2)
	nfac <- sum(cumvar < frac.var) + 1
	if(nfac > max(range.factors)) {
		nfac <- max(range.factors)
	} else if(nfac < min(range.factors)) {
		nfac <- min(range.factors)
	}
	if(nfac > length(cumvar)) nfac <- length(cumvar)
	fseq <- 1:nfac
	loadings <- scale(decomp$v[, fseq, drop=FALSE], 
		scale=1/decomp$d[fseq], center=FALSE)
	svd.d <- decomp$d

	if(iter.max > 0) {
		cmat <- crossprod(x[, fullcols, drop=FALSE]) 
                uniqueness <- 1 - rowSums(loadings^2)
                uniqueness[which(uniqueness < 0)] <- 0
                uniqueness[which(uniqueness > 1)] <- 1
		start <- uniqueness
		converged <- FALSE
                for(i in 1:iter.max) {
                        cor.red <- cmat
                        diag(cor.red) <- diag(cor.red) - uniqueness
                        decomp <- eigen(cor.red, symmetric=TRUE)
                        t.val <- decomp$value[fseq]
                        t.val[t.val < 0] <- 0
                        loadings <- scale(decomp$vector[, fseq, drop=FALSE], 
				center=FALSE, scale=1/sqrt(t.val))
                        uniqueness <- 1 - rowSums(loadings^2)
                        uniqueness[which(uniqueness < 0)] <- 0
                        uniqueness[which(uniqueness > 1)] <- 1
                        if(all(abs(uniqueness - start) < tol)) {
                                converged <- TRUE
                                break
                        }
                        start <- uniqueness
                }
	}
	dimnames(loadings) <- list(fullcolnames, NULL)

	if(any(num.mis > 0)) {
		# calculate loadings for columns with NAs
		floadings <- loadings
		if(zero.load) {
			loadings <- array(0, c(nassets, nfac))
		} else {
			meanload <- colMeans(floadings)
			loadings <- t(array(meanload, c(nfac, nassets)))
		}
		dimnames(loadings) <- list(dimnames(x)[[2]], NULL)
		loadings[dimnames(floadings)[[1]], ] <- floadings
		scores <- x[, fullcols, drop=FALSE] %*% floadings
		dsquare <- svd.d[1:nfac]^2
		nfac.miss[nfac.miss > nfac] <- nfac
		
		xna <- is.na(x)
		for(i in (1:nassets)[num.mis > 0 & nobs - num.mis > reg.min]) {
			t.nfac <- nfac.miss[ num.mis[i] ]
			if(t.nfac == 0) next
			t.okay <- which(!xna[, i])
			t.seq <- 1:t.nfac
			t.load <- lsfit(x[t.okay, i], scores[t.okay, t.seq], 
				intercept=FALSE)$coef / dsquare[t.seq]
			loadings[i, t.seq] <- t.load
			NULL
		}
	}

	comm <- rowSums(loadings^2)
	if(any(comm > 1)) {
		# adjust loadings where communalities too large
		toobig <- comm > 1
		if(verbose >= 2 && sum(reallytoobig <- comm > 1+1e-5)) {
			anam <- dimnames(loadings)[[1]]
			if(!length(anam)) {
				anam <- paste("V", 1:nrow(loadings), sep="")
			}
			warning(paste(sum(reallytoobig), 
				"asset(s) being adjusted",
				"from negative specific variance",
				"-- the assets are:", paste(anam[reallytoobig],
				collapse=", "), 
				"(warning suppressed with verbose < 2)"))
		}
		loadings[toobig,] <- loadings[toobig,] / sqrt(comm[toobig])
		comm[toobig] <- 1
	}

	uniqueness <- 1 - comm
	if(!is.numeric(specific.floor) || length(specific.floor) != 1) {
		stop(paste("'specific.floor' should be a single number",
			"-- given has mode", mode(specific.floor),
			"and length", length(specific.floor)))
	}
	if(specific.floor > 0) {
		if(!is.character(floor.type) || length(floor.type) != 1) {
			stop(paste("'floor.type' must be a single character",
				"string -- given has mode", mode(floor.type),
				"and length", length(floor.type)))
		}
		floor.menu <- c("quantile", "fraction")
		floor.num <- pmatch(floor.type, floor.menu, nomatch=0)
		if(floor.num == 0) {
			stop(paste("unknown or ambiguous input for",
				"'floor.type' -- valid choices are:",
				paste(floor.menu, collapse=", ")))
		}
		floor.type <- floor.menu[floor.num]
		switch(floor.type,
			"quantile" = {
				uf <- quantile(uniqueness, specific.floor)
				uniqueness[which(uniqueness < uf)] <- uf
			},
			"fraction" = {
				uniqueness[which(uniqueness < 
					specific.floor)] <- specific.floor
			}
		)
	}
	
	sdev[which(sdzero)] <- 0
	switch(output, 
		full= {
			cmat <- loadings %*% t(loadings)
			cmat <- t(sdev * cmat) * sdev
			diag(cmat) <- diag(cmat) + uniqueness * sdev^2
			attr(cmat, "number.of.factors") <- ncol(loadings)
			attr(cmat, "timestamp") <- date()
		},
		systematic=,
		specific=,
		factor={
			cmat <- list(loadings=loadings, 
				uniquenesses= uniqueness, sdev=sdev, 
				constant.names=sdzero.names,
				cumulative.variance.fraction=cumvar,
				timestamp=date(), call=match.call())
			class(cmat) <- "statfacmodBurSt"
		})
	if(output == "systematic" || output == "specific") {
		fitted(cmat, output=output)
	} else {
		cmat
	}
}

