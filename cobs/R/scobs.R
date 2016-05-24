#### $Id: scobs.R,v 1.46 2011/04/29 15:16:14 maechler Exp $

.onAttach <- function(lib, pkg) {
    ## now have NAMESPACE library.dynam("cobs", pkg, lib)
    if(interactive() || getOption("verbose")) # not in test scripts
	packageStartupMessage(sprintf(
		"Package %s (%s) attached.  To cite, see citation(\"%s\")\n",
		pkg, utils::packageDescription(pkg)$Version, pkg))
}

## S+ does not allow "cut(*, labels = FALSE)" -- use cut00() for compatibility:
if(is.R()) {
    cut00 <- function(x, breaks)
	cut.default(x, breaks, labels = FALSE, include.lowest = TRUE)
} else { ## S-plus  (tested only with S+ 6.0):
    cut00 <- function(x, breaks)
	as.integer(cut.default(x, breaks, include.lowest = TRUE))
}

mk.pt.constr <- function(pointwise)
{
    ## Purpose: produce the 'ptConstr' list from the 'pointwise' 3-col. matrix
    ## Author: Martin Maechler, 29 Jul 2006

    if(is.null(pointwise)) {
	list(n.equal = 0, n.greater = 0, n.smaller = 0, n.gradient = 0)
    }
    else { ## 'pointwise'
	if(!is.matrix(pointwise)|| dim(pointwise)[2] != 3)
	    stop("	Argument 'pointwise' has to be a three-column matrix.")
	kind <- pointwise[,1]		# .Alias
	equal	<- pointwise[kind ==  0, , drop = FALSE]
	greater <- pointwise[kind ==  1, , drop = FALSE]
	smaller <- pointwise[kind == -1, , drop = FALSE]
	gradient <-pointwise[kind ==  2, , drop = FALSE]
	list(equal   = equal,
	     greater = greater,
	     smaller = smaller,
	     gradient= gradient,
	     n.equal   = nrow(equal),# can be 0 ..
	     n.greater = nrow(greater),
	     n.smaller = nrow(smaller),
	     n.gradient= nrow(gradient))
    }
}

cobs <-
function(x, y,
	 constraint = c("none", "increase", "decrease",
			"convex", "concave", "periodic"),
	 w = rep(1,n),
	 knots, nknots = if(lambda == 0) 6 else 20,
         method = "quantile", degree = 2, tau = 0.5, lambda = 0,
         ic = c("AIC", "SIC", "BIC", "aic", "sic", "bic"),
	 knots.add = FALSE, repeat.delete.add = FALSE, pointwise = NULL,
         keep.data = TRUE, keep.x.ps = TRUE,
	 print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
         lambdaSet = exp(seq(log(lambda.lo), log(lambda.hi), length= lambda.length)),
	 lambda.lo = f.lambda*1e-4, lambda.hi = f.lambda*1e3, lambda.length = 25,
	 maxiter = 100, rq.tol = 1e-8, toler.kn = 1e-6, tol.0res = 1e-6, nk.start = 2,
### old "back-compatibility-only" arguments:
	 eps, n.sub, coef, lstart, factor)

{
    ## preamble
    ##
    cl <- match.call()
    constraint <- match.arg(constraint, several.ok = TRUE)
    if(length(constraint) == 0 || any(constraint == "none"))
	constraint <- "none"

    if(any(oldN <- names(cl) %in% c("eps", "n.sub", "coef", "lstart", "factor"))) {
	n <- sum(oldN)
	warning(paste("The use of", ngettext(n,"argument", "arguments"),
		      paste(sQuote(names(cl)[oldN]), collapse=", "),
		      "is deprecated.\n ", ngettext(n,"It is","They are"),
		      "now unused and will be removed in the future."))
    }
    ## FIXME: add "proper" na.action (as lm(), ..)
    na.idx <- is.na(x) | is.na(y)
    x <- x[!na.idx]
    y <- y[!na.idx]
    minx <- min(x)
    maxx <- max(x)
    n <- nrq <- length(x)
    ox <- order(x)
    xo <- x[ox]

    select.lambda <- (lambda < 0)
    if(select.lambda) {
        ## To make things scale-equivariant, the default lambdaSet
        ## must scale according scales of x and y :
        ## oops: does NOT scale with 'y'!   f.lambda <- sd(y) * sd(x)^degree
        f.lambda <- sd(x)^degree
	if(lambda.lo >= lambda.hi)
	    stop("The argument 'lambda.hi' has to be greater than 'lambda.lo'.")
	if(lambda.lo <= 0) stop("The argument 'lambda.lo' has to be positive.")
    }
    if(length(unique(y)) == 2 && print.warn)
	## warn(7)
	cat("\n It looks like you are fitting a binary choice model.",
	    "We recommend pre-smoothing your data using smooth.spline(),",
	    "loess() or ksmooth() before calling COBS\n", sep="\n ")

    ##
    ## generate default knots sequence
    ##
    lux <- length(ux <- unique(xo)) # needed in both cases
    if(missing(knots)) {
	select.knots <- TRUE
	nk.max <- if(degree == 1) lux else lux - 1
	if(nknots > nk.max) {
	    if(!missing(nknots)) ## 'nknots' specified explicitly
		warning("You chose nknots = ", nknots,
			". It has to be no greater than ", nk.max,
			" for degree = ", degree, ".\n",
			" cobs() has automatically reduced it to ", nk.max, ".\n")
	    nknots <- nk.max
	}

	knots <-
	    if(method == "quantile") {
		if(degree == 1 && lux == nknots)
		    ux
		else ux[seq(1, lux, len = nknots)] ## ``rounded'' quantiles
	    }
	    else ## "equidistant"
		seq(xo[1], xo[n], len = nknots)
    }
    else {
	names(knots) <- NULL
	lKn <- length(knots <- sort(knots))
	select.knots <-
	    if(!missing(nknots) && nknots != lKn) {
		warning("you specified 'nknots' and 'knots' with  nknots != length(knots).\n",
			" Hence will disregard 'nknots' and select from 'knots'")
		TRUE
	    } else missing(nknots)
	## => select.knots is FALSE  iff  nknots == length(knots), both specified
	nknots <- lKn
	if(degree == 2 && nknots == lux) {
	    ## ensure that  max #{knots} is n-1 for quadratic spline
            ## If we don't do this we get "singularity problem .." warnings
            ## in some cases
	    warning("The number of knots can't be equal to the number of unique x for degree = 2.\n",
		    "'cobs' has automatically deleted the middle knot.")
	    nknots <- length(knots <- knots[-round(nknots/2)])
	}
	if(knots[1] > minx || knots[nknots] < maxx)
	    stop("  The range of knots should cover the range of x.")
    }
    if(nknots < 2) stop("  A minimum of two knots is needed.")
    if(nknots == 2) {
	if(lambda == 0 && select.knots)
	    stop("  Can't perform automatic knot selection with nknots == 2.",
		 if(degree == 2) " Try 'degree=1'.")
	else if(degree == 1)
	    stop("  You need at least 3 knots when lambda!=0 and degree==1")
    }
    ## make sure that there is at least one observation between any pair of
    ## adjacent knots
    if(length(unique(cut00(x, knots))) != nknots - 1)
	stop("There is at least one pair of adjacent knots that contains no observation.")

    ptConstr <- mk.pt.constr(pointwise)

    ##
    ## set up proper dimension for the pseudo design matrix
    ##
    dim.o <- getdim2(degree,nknots,constraint)
    ks <- dim.o$ks
    neqc <- dim.o$n.eqc + ptConstr$n.equal + ptConstr$n.gradient
    niqc <- dim.o$n.iqc + ptConstr$n.greater + ptConstr$n.smaller
    if(lambda == 0) { ## quantile *regression* splines (no penalty)
	##
	## compute B-spline coefficients for quantile *regression* B-spline with
	## stepwise knots selection or with fixed knots
	##
        ic <- toupper(match.arg(ic))
	rr <- qbsks2(x, y, w, pw = 0, knots, nknots, degree, Tlambda = lambda,
		     constraint, ptConstr, maxiter, trace,
		     nrq, nl1 = 0, neqc, tau, select.lambda = select.lambda,
		     ks, do.select = select.knots, knots.add, repeat.delete.add, ic,
		     print.mesg = print.mesg,
                     give.pseudo.x = keep.x.ps,
		     rq.tol = rq.tol, tol.kn = toler.kn, tol.0res = tol.0res,
		     print.warn = print.warn, nk.start = nk.start)
	knots <- rr$knots
	nknots <- rr$nknots
    }
    else { ## lambda !=0 : quantile *smoothing* B-Splines with penalty

	if(degree == 1) {
	    nl1 <-  nknots - 2
	    nvar <- dim.o$nvar
	}
	else { ## degree == 2
	    nl1 <- 1
	    nvar <- dim.o$nvar + 1 # one more parameter for sigma
	    niqc <- niqc + 2 * (nknots - 1)
	}
        pw <- rep(1, nknots + degree-3)

        ## => lambda is chosen by ic (SIC / AIC )
        ## and lambdaSet := grid of lambdas in log scale [ == sfsmisc::lseq() ]

	##
	## compute B-spline coefficients for quantile smoothing B-spline
	##
	if(select.lambda && print.mesg)
	    cat("\n Searching for optimal lambda. This may take a while.\n",
		"  While you are waiting, here is something you can consider\n",
		"  to speed up the process:\n",
		"      (a) Use a smaller number of knots;\n",
		"      (b) Set lambda==0 to exclude the penalty term;\n",
		"      (c) Use a coarser grid by reducing the argument\n",
		"	   'lambda.length' from the default value of 25.\n") # 3

        ## shift the first and last knot a tiny bit "outside" {same as in qbsks2()}:
        rk <- diff(range(knots))
        knots[1] <- knots[1] - toler.kn*rk
        knots[nknots] <- knots[nknots] + toler.kn*rk

	rr <- drqssbc2(x, y, w, pw = pw, knots = knots, degree = degree,
		       Tlambda = if(select.lambda) lambdaSet else lambda,
		       constraint = constraint, ptConstr = ptConstr,
		       maxiter = maxiter, trace = trace-1,
		       nrq, nl1, neqc, niqc, nvar,
		       tau = tau, select.lambda = select.lambda,
		       give.pseudo.x = keep.x.ps,
		       rq.tol = rq.tol, tol.0res = tol.0res,
		       print.warn = print.warn)
    }
    if(any(rr$icyc >= maxiter))
	warning("The algorithm has not converged after ", maxiter, " iterations",
		if(select.lambda) " for at least one lambda", if(!is.null(pointwise)|
		constraint!="none") "\nCheck to make sure that your constraint is feasible",
		immediate. = TRUE)
    if(!all(rr$ifl == 1)) ## had problem
	warning("Check 'ifl'")

    nvar <- rr$nvar
    if(length(rr$coef) != nvar) message("length(rr$coef) != nvar -- and MM thought it should")
    Tcoef <- rr$coef[1:nvar]
    ##
    ## compute the residual
    ##
    y.hat <- .splValue(degree, knots, Tcoef, xo)
    y.hat <- y.hat[order(ox)]		# original (unsorted) ordering

    r <- list(call = cl,
	      tau = tau, degree = degree, constraint = constraint,
	      ic = if(lambda == 0) ic, pointwise = pointwise,
	      select.knots = select.knots, select.lambda = select.lambda,
	      x = if(keep.data) x,
	      y = if(keep.data) y,
	      resid = y - y.hat, fitted = y.hat,
	      coef = Tcoef, knots = knots,
	      k0 = rr$k,
	      k	 = if(select.lambda) rr$kk else rr$k, ##was: min(rr$k, nknots-2+ks),
	      x.ps = rr$pseudo.x,
	      SSy = sum((y - mean(y))^2),
	      lambda = rr$lambda,
	      icyc   = rr$icyc,
	      ifl    = rr$ifl,
	      pp.lambda = if(select.lambda) rr$pp.lambda,
	      pp.sic	= if(select.lambda) rr$sic,
	      i.mask	= if(select.lambda) rr$i.mask)

    class(r) <- "cobs"
    r
} ## cobs()

knots.cobs <- function(Fn, ...) Fn$knots

coef.cobs  <- function(object, ...) object$coef

.print.part <- function(x, nMin, namX = deparse(substitute(x)),
                        digits = getOption("digits"))
{
    ## Purpose: print (only part) of (numeric) vector 'x'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 7 Aug 2006
    force(namX)
    stopifnot(length(nMin) == 1, nMin == round(nMin))
    nx <- length(x)
    cat(namX, "[1:", nx,"]: ", sep = "")
    chk <- format(x[if(nx <= nMin) 1:nx else c(1:(nMin-1), nx)], digits=digits)
    if(nx > nMin) chk[nMin-1] <- "... "
    cat(chk, sep = ", "); cat("\n")
}

print.cobs <- function(x, digits = getOption("digits"), ...) {
    if(!is.numeric(lam <- x$lambda))
	stop("'x' is not a valid \"cobs\" object")
    cat("COBS ", if(lam == 0) "regression" else "smoothing",
	" spline (degree = ", x$degree, ") from call:\n	 ", deparse(x$call),
	if(any(x$ifl != 1)) {
	    if(all(x$ifl != 1))
		sprintf("\n\n **** ERROR in algorithm: ifl = %d\n\n", x$ifl[1])
	    else "\n * Warning in algorithm: some ifl != 1\n"
	},
	"\n{tau=",format(x$tau,digits),"}-quantile",
	";  dimensionality of fit: ",x$k," from {",
	paste(unique(x$k0),collapse=","), "}\n", sep="")

    .print.part(x$knots, nMin = 5, digits = digits)
    if(lam != 0) {
	cat("lambda =", format(lam, digits = digits))
	if((nlam <- length(x$pp.lambda))) {
	    cat(", selected via ", "SIC",# x$ic,
		", out of ", nlam, " ones.", sep='')
	}
	cat("\n")
    }
    invisible(x)
} # print

summary.cobs <- function(object, digits = getOption("digits"), ...) {
    if(!is.numeric(object$lambda))
	stop("'object' is not a valid \"cobs\" object")
    print(object, digits = digits, ...)# includes knots
    if(!is.null(pw <- object$pointwise)) {
	cat("with",nrow(pw),"pointwise constraints\n")
    }
    .print.part(object$coef, nMin = 7, namX = "coef", digits = digits)
    tau <- object$tau
    if(abs(tau - 0.50) < 1e-6)
	cat("R^2 = ", round(100 * (1 - sum(object$resid^2) / object$SSy), 2),
	    "% ;  ", sep="")
    k <- sum((r <- resid(object)) <= 0)
    n <- length(r)
    cat("empirical tau (over all): ",k,"/",n," = ", format(k/n,digits),
	" (target tau= ",tau,")\n", sep="")

    ## add more -- maybe finally *return* an object and define
    ## print.summary.cobs <- function(x, ...)

} # summary

residuals.cobs <- function (object, ...) object$resid
fitted.cobs <- function (object, ...) object$fitted

predict.cobs <-
    function(object, z, minz = knots[1], maxz = knots[nknots], nz = 100,
	     interval = c("none", "confidence", "simultaneous", "both"),
	     level = 0.95, ...)
{
    if(is.null(knots <- object$knots)	||
       is.null(coef  <- object$coef)	||
       is.null(degree<- object$degree)	||
       is.null(tau   <- object$tau)) stop("not a valid 'cobs' object")

    interval <- match.arg(interval)

    big	       <- if(is.R()) 3.4028234663852886e+38 else .Machine$single.xmax
    ##IN
    single.eps <- if(is.R())1.1920928955078125e-07 else .Machine$single.eps

    nknots <- length(knots)
    ord <- as.integer(degree + 1)
    nvar   <- length(coef)

    ##DBG cat("pr..cobs(): (ord, nknots, nvar) =",ord, nknots, nvar, "\n")

    backOrder <- function(m) m
    ##
    ## compute fitted value at z
    ##
    ## MM: why should z be *inside* (even strictly) the knots interval? _FIXME_
    if(missing(z)) {
	if(minz >= maxz) stop("minz >= maxz")
	##NOT YET (for "R CMD check" compatibility):
	## z <- seq(minz, maxz, len = nz)
	z <- seq(max(minz,knots[1]     + single.eps),
                 min(maxz,knots[nknots]- single.eps), len = nz)
    }
    else {
        ## Careful:  predict(*, x) should predict at 'x', not sort(x) !!
        ## Note this is needed, since .splValue() requires *increasing* z
        notOrd <- is.unsorted(z)
        if (notOrd) {
            iz.ord <- order(z)
            z <- z[iz.ord]
            i.rev <- order(iz.ord)
            backOrder <- function(m) m[ i.rev , , drop = FALSE]
        }
	##IN z <- z[z > knots[1] & z < knots[nknots]]
	nz <- length(z)
    }

    fit <- .splValue(degree, knots, coef, z)

    if(interval != "none") {
	##
	## compute confidence bands
	## both (pointwise and simultaneous : as cheap as only one !
	##
	z3 <- .splBasis(ord = ord, knots, ncoef = nknots + degree - 1, xo = z)
	idx <- cbind(rep(1:nz, rep(ord, nz)),
		     c(outer(1:ord, z3$offsets, "+")))
	X <- matrix(0, nz, nvar)
	X[idx] <- z3$design
        ## This is a residuum from cobs99 to make the old Bartel and Conns'
        ## algorithm more stable (and will probably never be used now):
	if(any(ibig <- abs(X) > big)) {
	    X[ibig] <- X[ibig] / big^0.5
	    warning("re-scaling ", sum(ibig), "values of spline basis 'X'")
	}

        if(is.null(object$x.ps))
	    stop("no 'x.ps' pseudo.x component; must use cobs(*, keep.x.ps = TRUE)")
        ##MM: We should have crossprod() and qr() working for sparse matrices,
        ## 	at least '%*%' and t() work; [hmm, but there is  slm() in 'SparseM' !]
	## Tqr <- qr(crossprod(object$x.ps))
        Tqr <- qr(as.matrix(t(object$x.ps) %*% object$x.ps))
	if(Tqr$rank != dim(object$x.ps)[2])
	    stop("The pseudo design matrix is singular; this can most likely be solved by using a smaller lambda")# when obtaining qsbs.out
	## Improved way of computing xQx = diag(X %*% solve(Tqr) %*% t(X)) :
	tX <- t(X)
	xQx <- colSums(qr.coef(Tqr, tX) * tX)

	res <- object$resid
	n <- length(res)

	s <- shat(res, tau, 1 - level, hs = TRUE)
	cn <- sqrt(xQx * tau * (1 - tau))
	sde <- cn * s
	an <- sqrt(qchisq(level, object$k))   * sde
	bn <- qt((1 + level)/2, n - object$k) * sde

	backOrder(cbind(z = z, fit = fit,
                        cb.lo = fit - an, cb.up = fit + an,
                        ci.lo = fit - bn, ci.up = fit + bn))
    }# interval
    else
	backOrder(cbind(z = z, fit = fit))
} # predict



##   "2" : 2nd ("scobs") cobs version
getdim2 <- function(degree, nknots,
		   constraint = c("none", "increase", "decrease",
		   "convex", "concave", "periodic"))
{
    ##=########################################################################
    ##
    ## Compute the appropriate dimension for the pseudo design matrix
    ##
    ##=########################################################################
    ## Note: "periodic" is here treated differently than in "cobs99" cobs()
    ##	   because the new FN algorithm doesn't handle separate equality
    ##     constraints. Hence, all equality constraints are processed through
    ##	   pairs of inequality contraints.
    ## Note 2: now also works for  *multiple* constraints
    ##	      n.iqc will be a vector of the same length as 'constraint'
    if(degree == 1)	ks <- 2
    else if(degree == 2)ks <- 3
    else stop("degree has to be either 1 or 2")
    deg1 <- as.integer(degree == 1)
    nvar <- nknots - 2 + ks
    n.eqc <- # the number of EQuality Constraints
        0 ## if(constraint == "periodic") 2 else 0
    n.iqc <- # the number of IneQuality Constraints,  will be summed up :
	sapply(match.arg(constraint, several.ok = TRUE),
	       function(constr) {
		   if(constr == "increase" || constr == "decrease")
		       nknots - deg1
		   else if(constr == "concave" || constr == "convex")
		       nknots - 1 - deg1
		   else if(constr == "periodic" )
		       4
		   else if(constr == "none")
		       0
	       })
    list(n.iqc = n.iqc, n.eqc = n.eqc, ks = ks, nvar = nvar)
}

### These are (only) used for confidence intervals :

shat <- function(residual, tau, alpha, hs)
{
    ##=########################################################################
    ##
    ## sparsity estimate from empirical quantile function using residuals
    ##
    ##=########################################################################
    residual <- sort(residual)
    n <- length(residual)
    residual <- c(residual[1], residual, residual[n])
    grid <- c(0, seq(0.5/n, 1 - 0.5/n, 1/n), 1)
    hn <- dn(tau, n, hs = hs, alpha)
    ## for small n,  tau +/- hn might be outside [0,1]
    bound <- pmax(0, pmin(1, c(tau - hn, tau + hn)))
    idx <- cut00(bound, grid)
    lambda <- bound * n - (idx - 1) + 0.5
    return(diff(lambda * residual[idx + 1] +
		(1 - lambda) * residual[idx])/(2 * hn))
}

dn <- function(p, n, hs = FALSE, alpha)
{
    ##=########################################################################
    ##
    ## compute window width for sparsity estimator
    ## at quantile p, level = 1-alpha,	n observations
    ## according to
    ##	 Hall and Sheather (1988),  <<-	 hs=TRUE,   or
    ##	 Bofinger (1975),	    <<-	 hs=FALSE
    ##=########################################################################
    x0 <- qnorm(p)
    f0 <- dnorm(x0)
    stopifnot(is.logical(hs))
    if(hs)
	n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
	    ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3)
    else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
}

plot.cobs <-
    function(x, which = if(x$select.lambda) 1:2 else 2, show.k = TRUE,
	     col = par("col"), l.col = c("red","pink"), k.col = gray(c(0.6, 0.8)),
             lwd = 2, cex = 0.4, ylim = NULL,
	     xlab = deparse(x$call[[2]]),
	     ylab = deparse(x$call[[3]]),
	     main = paste(deparse(x$call, width.cutoff = 100), collapse="\n"),
	     subtit= c("choosing lambda", "data & spline curve") , ...)
{
    stopifnot(all((which <- sort(unique(which))) %in% 1:2))
    both.plots <- all(which == 1:2)

    if(any(which == 1)) { ## --- SIC ~ lambda --------------
	stopifnot(x$select.lambda) # have no $pp.lambda . . .

	if(both.plots) {
	    op <- par(mfcol = 1:2, mgp = c(1.6, 0.5, 0),
		      mar = c(3.1, 3.1, 4.1, 0.6))
	}
	i.good <- x$ifl == 1
	i.bad  <- !i.good
        i.mask <- x$i.mask
        i.show <- !i.mask

	plot(x$pp.lambda, x$pp.sic, type = "n", log = "x",
	     xlab = expression(lambda),
	     ylab = "SIC",# paste("ic =", x$ic))
	     main = if(!both.plots) main) # no, col.axis="red")
	mtext(subtit[1], side=3, font = par("font.main"))
	lines (x$pp.lambda[i.good], x$pp.sic[i.good], type = 'o',
	       col = l.col[2], pch = 16)
	lines (x$pp.lambda[i.good & i.show], x$pp.sic[i.good & i.show], type = 'o',
	       col = l.col[1], pch = 16)
	lines(x$lambda, x$pp.sic[x$pp.lambda == x$lambda],# <- the chosen lambda
	      type = "h", lty = 3, col = l.col[1])
        mtext(substitute(hat(lambda) == LAM, list(LAM = formatC(x$lambda, digits=3))),
              side = 3, line = -1, col = l.col[1])
	if(any(i.bad) || any(i.mask)) {
	    all.are.11 <- all(x$ifl[i.bad] == 11)
	    all.are.18 <- all(x$ifl[i.bad] == 18) # = 1+{17 : the new code}
	    ch <- if(all.are.11) "ifl = 10+1"
	    else  if(all.are.18) "ifl = 17+1"
	    else paste("ifl in {",
		       paste(unique(x$ifl[i.bad]), collapse= ","), "}", sep='')
	    ## TODO(?) could use even more legend categories ..
	    points(x$pp.lambda [i.bad], x$pp.sic[i.bad], col = l.col[2], pch = 1)
	    points(x$pp.lambda [i.bad & i.show], x$pp.sic[i.bad & i.show], col = l.col[1], pch = 1)
	    legend("right", ## MM had "topleft"
		   c(if(any(i.bad)) c("ifl = 1, good fits",
                                      paste(ch,", bad fits", sep='')),
                     if(any(i.mask)) c("SIC when k <= sqrt(n)",
                                       "SIC when k > sqrt(n)")),
		   pch = c(if(any(i.bad)) c(16, 1), if(any(i.mask)) c(16,16)),
                   col = c(if(any(i.bad)) rep(l.col[1],2), if(any(i.mask)) l.col))

	}
	if(show.k) {
	    par(new = TRUE)
            ## plot(x$pp.lambda[i.good], x$k0[i.good], axes = FALSE, ann = FALSE,
            ##     type = 'l', log = "x", col = col.k, lty = 2)
	    plot(x$pp.lambda, x$k0, log = "x", type = "n",
                 axes = FALSE, ann = FALSE)
            lines (x$pp.lambda[i.good], x$k0[i.good], type = 'o',
                   col = k.col[2], pch = 16)
            lines (x$pp.lambda[i.good & i.show], x$k0[i.good & i.show], type = 'o',
                   col = k.col[1], pch = 16)
	    axis(4, at = unique(c(0, x$k0[i.good])),
		 las = 2, col = k.col[1], col.axis = k.col[1])
	    mtext("k", side = 4, line = .5, at = par("usr")[4],
                  las = 2, col = k.col[1])
            if(any(i.bad)) {
                points(x$pp.lambda[i.bad], x$k0[i.bad], col = gray(.9), pch = 1)
                points(x$pp.lambda[i.bad & i.show], x$k0[i.bad & i.show],
                       col = gray(.8), pch = 1)
            }
	}
    }

    ## warning("case of lambda = 0  not fully implemented")

    ## FIXME -- do something (?) for 'lambda = 0', i.e. select.knots

    if(any(which == 2)) { ## ---------- (x,y) + smooth cobs curve ----
	if(is.null(x$x)) {
	    ## MM: "works"	only if original (x,y) objects are still "there":
	    x$x <- eval.parent(x$call[[2]])
	    x$y <- eval.parent(x$call[[3]])
	}
        px <- predict(x)
	if(is.null(ylim)) # make sure ylim fits (y & predict(.)):
            ylim <- range(x$y, px[,"fit"], finite=TRUE)
        plot(x$x, x$y, xlab = xlab, ylab = ylab, ylim= ylim,
             main = if(!both.plots) main, col = col, cex = cex, ...)
	lines(px, col = l.col[1], lwd = lwd, ...)
	if(both.plots) {
	    mtext(subtit[2], side = 3, font = par("font.main"))
	    par(op)		# <- back to 1x1 figure (we assume ..)
	    mtext(main, side=3, line=1,
		  cex = par("cex.main"), font = par("font.main"))
	}
    }
}
