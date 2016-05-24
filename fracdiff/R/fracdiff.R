
### Original file:
### copyright 1991 Department of Statistics, Univeristy of Washington

### Patched by Friedrich.Leisch, for use with R, 22.1.1997
### fixed & changed by Martin Maechler, since Dec 2003

if(getRversion() < "2.15")
paste0 <- function(...) paste(..., sep="")

.fdcov <- function(x, d, h, nar, nma, hess, fdf.work)
{
    npq <- as.integer(nar + nma)
    npq1 <- npq + 1L # integer, too
    stopifnot(length(di <- dim(hess)) == 2, di == c(npq1, npq1))
    fdc <- .C(fdcov, ## --> ../src/fdhess.c
	      x,
	      d,
	      h = as.double(if(missing(h)) -1 else h),
	      hd = double(npq1),
	      cov = hess, npq1,
	      cor = hess, npq1,
	      se = double(npq1),
	      fdf.work,
	      info = integer(1))[c("h","hd", "cov","cor", "se", "info")]

    f.msg <-
	if(fdc$info) {
	    msg <-
		switch(fdc$info,
		       "fdcov problem in gamma function",	# 1
		       "singular Hessian",			# 2
		       ## FIXME improve: different reasons for info = 3 :
		       "unable to compute correlation matrix; maybe change 'h'",
								# 3
		       stop("error in gamma function"))		# 4
	    warning(msg, call. = FALSE)
	    msg
	} else "ok"
    se.ok <- fdc$info %in% 0:2
    nam <- "d"
    if(nar) nam <- c(nam, paste0("ar", 1:nar))
    if(nma) nam <- c(nam, paste0("ma", 1:nma))

    dimnames(fdc$cov) <- dn <- list(nam, nam)
    if(se.ok) dimnames(fdc$cor) <- dn
    list(msg = f.msg, d = d, nam = nam,
	 h = fdc$h, hd = fdc$hd, se.ok = se.ok,
	 covariance.dpq = fdc$cov,
	 stderror.dpq	= if(se.ok) fdc$se, # else NULL
	 correlation.dpq= if(se.ok) fdc$cor)
}## end{.fdcov}

fracdiff <- function(x, nar = 0, nma = 0,
                     ar = rep(NA, max(nar, 1)), ma = rep(NA, max(nma, 1)),
                     dtol = NULL, drange = c(0, 0.5), h, M = 100, trace = 0)
{
    ## #########################################################################
    ##
    ##   x      - time series for the ARIMA model
    ##   nar    - number of autoregressive parameters
    ##   nma    - number of moving average parameters
    ##   ar     - initial autoregressive parameters
    ##   ma     - initial moving average parameters
    ##   dtol   - desired accurcay for d
    ##            by default (and if negative), (4th root of machine precision)
    ##		  is used.  dtol will be changed internally if necessary
    ##   drange - interval over which the likelihood function is to be maximized
    ##            as a function of d
    ##   h      - finite difference interval
    ##   M      - number of terms in the likelihood approximation
    ##
    ##           (see Haslett and Raftery 1989)
    ##
    ## ########################################################################

    cl <- match.call()
    if(any(is.na(x)))
        stop("missing values not allowed in time series")
    if(is.matrix(x) && ncol(x) > 2)
        stop("multivariate time series not allowed")
    n <- length(x)
    if(round(nar) != nar || nar < 0 || round(nma) != nma || nma < 0)
        stop("'nar' and 'nma' must be non-negative integer numbers")
    npq <- as.integer(nar + nma)
    npq1 <- npq + 1L # integer, too
    lenw <- max(npq + 2*(n + M),
                3*n + (n+6)*npq + npq %/% 2 + 1,
                31 * 12, ## << added because checked in ../src/fdcore.f
                (3 + 2*npq1) * npq1 + 1)## << this is *not* checked (there)
    lenw <- as.integer(lenw)
    ar[is.na(ar)] <- 0
    ma[is.na(ma)] <- 0
    if(is.null(dtol))
        dtol <- .Machine$double.eps^0.25 # ~ 1.22e-4
    ## if dtol < 0: the fortran code will choose defaults
    x <- as.double(x)

    ## this also initializes "common blocks" that are used in .C(.) calls :
    fdf <- .C(fracdf,
	      x,
	      n,
	      as.integer(M),
	      as.integer(nar),
	      as.integer(nma),
	      dtol = as.double(dtol),
	      drange = as.double(drange),
	      hood.etc = double(3),
	      d = double(1),
	      ar = as.double(ar),
	      ma = as.double(ma),
	      w = double(lenw),
	      lenw = lenw,
	      iw = integer(npq), ## <<< new int-work array
	      info = as.integer(trace > 0),## <- "verbose" [input]
	      .Machine$double.xmin,
	      .Machine$double.xmax,
	      .Machine$double.neg.eps,
	      .Machine$double.eps)[c("dtol","drange","hood.etc",
	      "d", "ar", "ma", "w", "lenw", "info")]

    fd.msg <-
        if(fdf$info) {
	    msg <-
                switch(fdf$info,
                       stop("insufficient workspace; need ", fdf$lenw,
                            " instead of just ", lenw),		# 1
                       stop("error in gamma function"),		# 2
                       stop("invalid MINPACK input"),		# 3
                       "warning in gamma function",		# 4
                       "C fracdf() optimization failure",	# 5
                       "C fracdf() optimization limit reached") # 6
            ## otherwise
            ## stop("unknown .C(fracdf, *) info -- should not happen")
	    warning(msg, call. = FALSE, immediate. = TRUE)
	    msg
	} else "ok"

    if(nar == 0) fdf$ar <- numeric(0)
    if(nma == 0) fdf$ma <- numeric(0)

    hess <- .C(fdhpq,
               hess = matrix(double(1), npq1, npq1),
               npq1,
               fdf$w)$hess

    ## NOTA BENE: The above  hess[.,.]  is further "transformed",
    ##            well, added to  and inverted  in fdcov :
    ## Cov == (-H)^{-1} == solve(-H)

    ## Note that the following can be "redone" using fracdiff.var() :

    fdc <- .fdcov(x, fdf$d, h,
                  nar=nar, nma=nma, hess=hess, fdf.work = fdf$w)

    dimnames(hess) <- dimnames(fdc$covariance.dpq)
    hess[1, ] <- fdc$hd
    hess[row(hess) > col(hess)] <- hess[row(hess) < col(hess)]

    hstat <- fdf[["hood.etc"]]
    var.WN <- hstat[3]
    structure(list(log.likelihood = hstat[1],
                   n = n,
		   msg = c(fracdf = fd.msg, fdcov = fdc$msg),
		   d = fdf$d, ar = fdf$ar, ma = fdf$ma,
		   covariance.dpq = fdc$covariance.dpq,
		   fnormMin = hstat[2], sigma = sqrt(var.WN),
		   stderror.dpq	  = if(fdc$se.ok) fdc$stderror.dpq, # else NULL
		   correlation.dpq= if(fdc$se.ok) fdc$correlation.dpq,
		   h = fdc$h, d.tol = fdf$dtol, M = M, hessian.dpq = hess,
		   length.w = lenw, call = cl),
	      class = "fracdiff")
}

### FIXME [modularity]: a lot of this is "cut & paste" also in fracdiff() itself
### ----- NOTABLY, now use  .fdcov() !

fracdiff.var <- function(x, fracdiff.out, h)
{
    if(!is.numeric(h))
        stop("h must be numeric")
    if(!is.list(fracdiff.out) || !is.numeric(M <- fracdiff.out$M))
        stop("invalid ", sQuote("fracdiff.out"))
    p <- length(fracdiff.out$ar)
    q <- length(fracdiff.out$ma)
    n <- length(x)
    npq <- p + q
    npq1 <- npq + 1
    lwork <- max(npq + 2 * (n + M),
                 3 * n + (n + 6) * npq + npq %/% 2 + 1,
                 (3 + 2 * npq1) * npq1 + 1)
    ## Initialize
    .C(fdcom,
       n,
       as.integer(M),
       (p),
       (q),
       as.double(fracdiff.out$log.likelihood),
       .Machine$double.xmin,
       .Machine$double.xmax,
       .Machine$double.neg.eps,
       .Machine$double.eps)
    ## Re compute Covariance Matrix:
    fdc <- .C(fdcov,
              as.double(x),
              as.double(fracdiff.out$d),
              h = as.double(h),
              hd = double(npq1),
              cov = as.double(fracdiff.out$hessian.dpq),
              as.integer(npq1),
              cor = as.double(fracdiff.out$hessian.dpq),
              as.integer(npq1),
              se = double(npq1),
              as.double(c(fracdiff.out$ma,
                          fracdiff.out$ar,
                          rep(0, lwork))),
              info = integer(1))
## FIXME: should be *automatically* same messages as inside fracdiff() above!
    fracdiff.out$msg <-
        if(fdc$info) {
            msg <-
                switch(fdc$info,
                       "warning in gamma function",
                       "singular Hessian",
                       "unable to compute correlation matrix",
                       stop("error in gamma function"))
            warning(msg)
        } else "ok"
    se.ok <- fdc$info != 0 || fdc$info < 3 ## << FIXME -- illogical!!
    nam <- "d"
    if(p) nam <- c(nam, paste0("ar", 1:p))
    if(q) nam <- c(nam, paste0("ma", 1:q))
    fracdiff.out$h <- fdc$h
    fracdiff.out$covariance.dpq <- array(fdc$cov, c(npq1,npq1), list(nam,nam))
    fracdiff.out$stderror.dpq    <- if(se.ok) fdc$se # else NULL
    fracdiff.out$correlation.dpq <- if(se.ok) array(fdc$cor, c(npq1, npq1))
    fracdiff.out$hessian.dpq[1, ] <- fdc$hd
    fracdiff.out$hessian.dpq[, 1] <- fdc$hd
    fracdiff.out
}## end{ fracdiff.var() }

## MM:  Added things for more  arima.sim() compatibility.
##      really, 'mu' is nonsense since can be done separately (or via 'innov').
fracdiff.sim <- function(n, ar = NULL, ma = NULL, d, rand.gen = rnorm,
                         innov = rand.gen(n+q, ...), n.start = NA,
			 backComp = TRUE, allow.0.nstart = FALSE, # <- for back-compatibility
                         start.innov = rand.gen(n.start, ...), ..., mu = 0)
{
    p <- length(ar)
    q <- length(ma)
    if(p) {
        minroots <- min(Mod(polyroot(c(1, -ar))))
        if(minroots <= 1) {
            warning("'ar' part of fracdiff model is not stationary!!")
            minroots <- 1.01 # -> n.start= 603 by default
        }
    }
    if(is.na(n.start))
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)), 0)
    if(n.start < p + q && !allow.0.nstart)
        stop("burn-in 'n.start' must be as long as 'ar + ma'")
    if(missing(start.innov)) {
        if(!backComp) force(start.innov)
    } else if(length(start.innov) < n.start)
        stop(gettextf("'start.innov' is too short: need %d points",
                      n.start), domain = NA)
    if(length(innov) < n+q) stop("'innov' must have length >= n + q")
    y <- c(start.innov[seq_len(n.start)], innov[1:(n+q)])
    stopifnot(is.double(y), length(y) == n + q + n.start)
    if(d < -1/2 || d > 1/2)
	stop("'d' must be in [-1/2, 1/2].  Consider using cumsum(.) or diff(.)
 for additional integration or differentiation")
    ii <- n.start - (if(backComp) 0L else q) + 1:n
    y <- .C(fdsim,
            as.integer(n + n.start),
            (p),
            (q),
            as.double(ar),
            as.double(ma),
            as.double(d),
            as.double(mu),
            y = y,
            s = double(length(y)),
            .Machine$double.xmin,
            .Machine$double.xmax,
            .Machine$double.neg.eps,
            .Machine$double.eps)[["s"]][ii]
    list(series = y, ar = ar, ma = ma, d = d, mu = mu, n.start = n.start)
}
