"lo.wam" <-
function(x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-7, trace = FALSE,
	se = TRUE, ...)
{
	if(is.data.frame(smooth.frame)) {
		first <- TRUE
		# first call to wam; set up some things
		#on first entry, smooth.frame is a data frame with elements the terms to be
		#smoothed in which
		data <- smooth.frame[, names(which), drop = FALSE]
		smooth.frame <- gam.match(data)
		dx <- as.integer(dim(x))
		oldClass(data) <- NULL
		atts <- lapply(data, attributes)
		span <- sapply(atts, "[[", "span")
		degree <- sapply(atts, "[[", "degree")
		nvars <- sapply(atts, "[[", "ncols")
		ndim <- sapply(atts, "[[", "dim")[2.,  ]
		npetc <- as.integer(c(dx, length(which), 0., maxit, 0.))
		nef <- smooth.frame$nef
		nvmax <- 200. + 300. * (1. - 1./log(apply(cbind(nef - 200.,
			3.), 1., max)))
		nspar <- (nef * span + 1.)
		liv <- 50. + (2.^nvars + 4.) * nvmax + 2. * nef
		lv <- 50. + (3. * nvars + 3.) * nvmax + nef + (ifelse(degree ==
			2., ((nvars + 2.) * (nvars + 1.))/2., nvars + 1.) +
			2.) * nspar
		LL <- nspar * nvmax
		liv <- liv + LL
		lv <- lv + (nvars + 1.) * LL
		which <- sapply(which, "[", 1.)
		wddnfl <- cbind(unlist(which), nvars, ndim, degree, nef, liv,
			lv, nvmax)
		storage.mode(wddnfl) <- "integer"
		spatol <- as.double(c(span, tol))
		nwork <- 9. * dx[1.] + sum(nef * (nvars + ndim + 4.) + 5. +
			3. * ndim)
		liv <- sum(liv)
		lv <- sum(lv)
		smooth.frame <- c(smooth.frame, list(npetc = npetc, wddnfl = 
			wddnfl, spatol = spatol,niwork=2*sum(nvars), nwork = nwork, liv = liv,
			lv = lv))
	}
	else first <- FALSE
	storage.mode(y) <- "double"
	storage.mode(w) <- "double"
	n <- smooth.frame$npetc[1.]
	p <- smooth.frame$npetc[2.]
	q <- smooth.frame$npetc[3.]
	fit <- .Fortran("baklo",
		x,
		y = y,
		w = w,
		npetc = smooth.frame$npetc,
		smooth.frame$wddnfl,
		smooth.frame$spatol,
		smooth.frame$o,
		etal = double(n),
		s = s,
		eta = double(n),
		beta = double(p),
		var = s,
		df = double(q),
		qr = x,
		qraux = double(p),
		qpivot = as.integer(1.:p),
                effects=double(n),
		integer(smooth.frame$liv),
		double(smooth.frame$lv),
                integer(smooth.frame$niwork),
		double(smooth.frame$nwork),
                        PACKAGE="gam")

	nit <- fit$npetc[4.]
	qrank <- fit$npetc[6.]
	if((nit == maxit) & maxit > 1.)
		warning(paste("lo.wam convergence not obtained in ", maxit,
			" iterations"))
	names(fit$df) <- dimnames(s)[[2]]
	names(fit$beta) <- labels(x)[[2]]
                qrx <- structure(list(qr = fit$qr,qraux = fit$qraux,
                     rank = qrank, pivot = fit$qpivot,tol=1e-7),class="qr")
        effects<-fit$effects
        r1 <- seq(len = qrx$rank)
        dn <- colnames(x)
        if (is.null(dn)) 
          dn <- paste("x", 1:p, sep = "")
        names(effects) <- c(dn[qrx$pivot[r1]], rep.int("", n - qrx$rank))
	rl <- list(coefficients = fit$beta, residuals = fit$y - fit$eta, 
                   fitted.values = fit$eta,
                   effects=effects, weights=w, rank=qrank,
                   assign=attr(x,"assign"),
                   qr=qrx,
                   smooth = fit$s,
                   nl.df = fit$df
                   )
	rl$df.residual <- n - qrank - sum(rl$nl.df) - sum(fit$w == 0.)
	if(se)
		rl <- c(rl, list(var = fit$var))
	c(list(smooth.frame = smooth.frame), rl)
}
