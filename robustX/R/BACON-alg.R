#### From  Atanassios KONDYLIS  @ uni ne . ch
####  15 Dec 2005     --> see ./BACON-alg_email.txt
####

#################### MAIN PROGRAM BACON

BACON <- function(x, y = NULL, intercept = TRUE,
                  m = min(collect * p, n * 0.5),
                  init.sel = c("Mahalanobis", "dUniMedian", "random", "manual"),
                  man.sel, init.fraction = 0, collect = 4, alpha = 0.95,
		  maxsteps = 100, verbose = TRUE)
{
    ## This S-Plus function performs an outlier identification
    ## algorithm to the data in the x array [n x p] and y vector [n]
    ## followin the lines described by Hadi et al. for their
    ## BACON outlier procedure.
    ##
    ## Ref:  Billor, N., Hadi, A. S., and Velleman , P. F. (2000),
    ##	     "BACON: Blocked Adaptive Computationally-Efficient
    ##	     Outlier Nominators,"
    ##	     Computational Statistics & Data Analysis, (in press)
    ##
    ## x:	 a multivariate matrix of dimension [n x p]
    ##		 considered as containing no missing values
    ## y:	 the response array [n] for the x data
    ## intercept: a logical indicating if an intercept has to be used
    ##		 for the regression
    ## init.sel: the initial selection mode; implemented modes are:
    ##	     "Mah" -> based on Mahalanobis distance (default)
    ##	     "dis" -> based on the distances from the medians
    ##	     "ran" -> based on a random selection
    ##	     "man" -> based on manual selection
    ##		 in this case the vector 'man.sel' which contains
    ##		 the indices of the selected observations
    ##		 must be given.
    ##	     "Mah" and "dis" are proposed by Hadi while "ran" and "man"
    ##	     were implemented in order to study the behaviour of BACON.
    ## init.fraction: if this parameter is > 0 then the tedious steps
    ##		 of selecting the initial subset are skipped and
    ##		 an initial subset of size n * init.fraction is
    ##		 chosen (with smallest dis)
    ## collect:	 a factor chosen by the user to define the size of
    ##		 the initial subset (p * collect)
    ## alpha:	 significance level
    ## maxsteps: the maximal number of iteration steps
    ##		 (to prevent infinite loops)
    ## verbose: prints messages on the standard output in order
    ##		 to inform the user about the progress
    ##
    ## Written 25.05.01 by U. Oetliker for S-Plus 5.1 (Unix).
    ## Modified 28.05.01 / 31.05.01 / 01.06.01 / 09.06.01
    ##	 10.06.01 / 17.06.01

    ## 'U. Oetliker' = Ueli Oetliker --- SFSO (Swiss Federal Statistical Office)
    ##
########## BEGIN BACON

    if(is.vector(x))
	x <- as.matrix(x)

    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n > 0, p > 0, n > p, n >= m, m > 0,
              0 < alpha, alpha < 1)
    if(length(y) == 0)
	return(mvBACON(x, collect=collect, m=m,
                       alpha=alpha, init.sel=init.sel, man.sel=man.sel,
                       maxsteps=maxsteps, verbose=verbose))
    else {
	mvB <- mvBACON(x, collect=collect, m=m,
                       alpha=alpha, init.sel=init.sel, man.sel=man.sel,
                       maxsteps=maxsteps, verbose=verbose)
	regB <- .lmBACON (x, y, intercept,
				 mvB$dis, init.fraction, collect, alpha,
				 maxsteps, verbose)
	return(list(subset=  regB$subset,
		    tis=     regB$tis,
		    mv.subset= mvB$subset))
    }
}

########## END BACON

#################### FUNCTION mvBACON

mvBACON <-
    function(x, collect = 4, m = min(collect * p, n * 0.5), alpha = 0.95,
             init.sel = c("Mahalanobis", "dUniMedian", "random", "manual"),
             man.sel, maxsteps = 100, allowSingular = FALSE, verbose = TRUE)
{
    ## docu: --> ../man/mvBACON.Rd
    ##                  ~~~~~~~~~~  help file

    ## Written 25.05.01 by U. Oetliker for S-Plus 5.1 (Unix).
    ## Modified 28.05.01 / 31.05.01 / 01.06.01 / 09.06.01 / 10.06.01

    ## more Modifications in porting to R : Martin Maechler, 2009-05 ++

    trace1 <- function(i, r, n)
    {
	cat("MV-BACON (subset no. ", i, "): ", r,
	    " of ", n, " (", round(r/n*100, digits = 2), " %)",
	    sep = "", fill = TRUE)
    }

    chiCrit <- function(n, p, r, alpha)
    {
	h <- (n + p + 1) / 2
	chr <- max(0, (h - r) / (h + r))
	cnp <- 1 + (p + 1) / (n - p) + 2 / (n - 1 - 3*p)
	cnpr <- cnp + chr
	chisq <- qchisq(alpha / n, p, lower.tail = FALSE)
	return(sqrt(chisq) * cnpr)
    }

    init.sel <- match.arg(init.sel)

    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n > 0, p > 0, n > p, n >= m, m > 0,
              0 < alpha, alpha < 1)

    ordered.indices <-
        switch(init.sel,
               "Mahalanobis" = {
                   center <- colMeans(x)
                   cov <- var(x)
                   order(mahalanobis(x, center, cov))
               },
               "random" = sample(n),
               "manual" = {
                   m <- length(man.sel)
                   stopifnot(is.numeric(man.sel), m > 0,
                             1 <= man.sel, man.sel <= n,
                             man.sel == as.integer(man.sel))
                   c(man.sel, c(1:n)[-man.sel])
               },
               "dUniMedian" = {
                   x.centr <- sweep(x, 2, apply(x, 2, median))
                   order(mahalanobis(x.centr, 0, cov(x.centr)))
               },
               ## otherwise:
               stop("invalid 'init.sel' -- should not happen; please report!")
               )

    ordered.x <- x[ordered.indices, , drop = FALSE]
    ## FIXME{MM}: qr(.) is used to determine rank -- *not* state of the art
    ## see Matrix::rankMatrix
    while (m < n && p > (rnk <- qr(var(ordered.x[1:m, , drop = FALSE]))$rank))
	m <- m + 1

    if(rnk < p && !allowSingular)
        stop("matrix-rank ( x[1:m,] ) < p  for all m <= n")

    subset <- 1:n %in% ordered.indices[1:m]

    presubset <- rep(FALSE, n)
    converged <- FALSE
    steps <- 1
    repeat {
        r <- sum(subset)
        if(verbose) trace1(steps, r, n)
	center <- colMeans(x[subset, , drop = FALSE])
	cov <- var(x[subset, , drop = FALSE])
	dis <- sqrt(mahalanobis(x, center, cov))

        converged <- !any(xor(presubset, subset))
        if(converged)
            break
	##print(dis)
	presubset <- subset
	limit <- chiCrit(n, p, r, alpha)
	subset <- dis < limit
	steps <- steps + 1
        if (steps > maxsteps)
            break
    }
    if(steps > maxsteps)
        warning("basic subset has not converged in ", maxsteps, " steps")

    list(dis = dis, subset = subset, center = center, cov = cov,
         steps = steps, limit = limit, converged = converged)
}

########## END mvBACON




#################### FUNCTION .lmBACON

.lmBACON <- function(x, y, intercept = TRUE, init.dis,
                    init.fraction = 0, collect = 4, alpha = 0.95,
                    maxsteps = 100, verbose = TRUE)
{
    ## This function performs the regression part of the
    ## outlier identification algorithm described by Hadi et al.
    ## for their BACON outlier procedure.
    ##
    ## Ref:  Billor, N., Hadi, A. S., and Velleman , P. F. (2000),
    ##	     "BACON: Blocked Adaptive Computationally-Efficient
    ##	     Outlier Nominators,"
    ##	     Computational Statistics & Data Analysis, (in press)
    ##
    ## x:	 a multivariate matrix of dimension [n x p]
    ##		 considered as containing no missing values
    ## y:	 the response array [n] for the x data
    ## intercept: a logical indicating if an intercept has to be used
    ##		 for the regression
    ## init.dis: the distances of the x matrix used for the initial
    ##		 subset determined by the mvBACON
    ## init.fraction: if this parameter is > 0 then the tedious steps
    ##		 of selecting the initial subset are skipped and
    ##		 an initial subset of size n * init.fraction is
    ##		 chosen (with smallest dis)
    ## collect:	 a factor chosen by the user to define the size of
    ##		 the initial subset (p * collect)
    ## alpha:	 significance level
    ## maxsteps: the maximal number of iteration steps
    ##		 (to prevent infinite loops)
    ## verbose: prints messages on the standard output in order
    ##		 to inform the user about the progress
    ##
    ## Written 09.06.01 by U. Oetliker for S-Plus 5.1 (Unix).
    ## Modified 10.06.01
    ##
########## INTERNAL FUNCTIONS OF .lmBACON

    trace1 <- function(i, r, n,
                       init.steps = FALSE, skip.init = FALSE) {
        cat("Reg-BACON (",if(init.steps)"init ", "subset no. ", i,
            if(skip.init) " after skipping init", "): ",
            r, " of ", n, " (", round(r/n*100, digits = 2), " %)",
            sep = "", fill = TRUE)
    }

    GiveTis <- function(x, y, subset, ordered.x, ordered.y, m = 0, check.rank = FALSE)
    {
	## This function calculates the t(i)s
	## for each observation of the data set y ~ x.
	## tis = resid / [sigmahatm * sqrt(1 - PP)]  for tis \in the subset
	## tis = resid / [sigmahatm * sqrt(1 + PP)]  for tis \in the rest
	## with PP = xiT(xmTxm)-1xiT
	## resid = yi - xiTbetahatm
	##
	## betahatm = (xTx)-1xTy
	##
	## sigmahatm = sqrt(s2)
	##   s2 = SSE/(n-p)
	##   SSE = eTe
	##   e = (In - PP)y = y - PPy = residuals
	##   PP = x(xTx)-1xT
	##
	## subset is a logical vector of length n
	## with T indicating belonging to the subset
	##
	## Written 09.06.01 by U. Oetliker for S-Plus 5.1 (Unix).
	##

	n <- nrow(x)
	p <- ncol(x)

	if(check.rank)
	{
	    while (
		   m < n
		   &&
		   p > qr(crossprod(ordered.x[1:m,,drop = FALSE]))$rank
		   )
		m <- m + 1

	    xm <- ordered.x[1:m,, drop = FALSE]
	    ym <- ordered.y[1:m]
	}
	else
	{
	    xm <- ordered.x
	    ym <- ordered.y
	}
	nm <- nrow(xm)

	xm  <- as.matrix(xm)
	fit.m <- lm(ym ~ xm[,1:ncol(xm)] -1)

	betahatm <- fit.m$coef
	x  <- as.matrix(x)
	resid <- y - as.vector(x %*% betahatm)

	em <- fit.m$residuals
	sigmahatm <- sqrt(sum(em * em) / (nm - p))

	hiim  <- hatvalues(fit.m)

	sqroot <- 1 + ifelse(subset, -1, 1) * hiim

	tis <- resid / (rep(sigmahatm, n) * sqrt(sqroot))
	list(m = m, tis = tis)
    }


########## BEGIN .lmBACON

    ordered.indices <- order(init.dis)
    if(intercept)
	x <- cbind(1, x)
    n <- nrow(x)
    p <- ncol(x)
    steps <- 0

    ordered.x <- x[ordered.indices, , drop = FALSE]
    ordered.y <- y[ordered.indices]
    if(init.fraction > sqrt(.Machine$double.eps)) {
	m <- round(init.fraction * n)

	tis <- GiveTis(x, y, subset, ordered.x, ordered.y, m, TRUE)
	subset <- is.element(1:n, ordered.indices[1:m])
	tis <- GiveTis(x, y, subset, ordered.x, ordered.y, m, TRUE)
	m <- tis$m
	tis <- abs(tis$tis)
	if(verbose) trace1(steps, m, n, skip.init = TRUE)

	ordered.indices <- order(tis)
	r <- m
    }
    else {
	m <- collect * p

	subset <- is.element(1:n, ordered.indices[1:m])
	tis <- GiveTis(x, y, subset, ordered.x, ordered.y, m, TRUE)
	m <- tis$m
	if(verbose) trace1(steps, m, n, init.steps = TRUE)
	tis <- abs(tis$tis)

	r <- p + 1
	if(verbose) trace1(steps, r, n, init.steps = TRUE)
	while (r < n && r < m) {
	    ordered.indices <- order(tis)
	    ordered.x <- x[ordered.indices, , drop = FALSE]
	    ordered.y <- y[ordered.indices]

	    subset <- is.element(1:n, ordered.indices[1:r])
	    tis <- GiveTis(x, y, subset, ordered.x, ordered.y, r, TRUE)
	    r <- tis$m
	    tis <- abs(tis$tis)

	    r <- r + 1
	    steps <- steps + 1
	    if(verbose) trace1(steps, r, n, init.steps = TRUE)
	}
	ordered.indices <- order(tis)
    }

    subset <- is.element(1:n, ordered.indices[1:r])
    presubset <- rep(FALSE, n)
    pre.r <- sum(presubset)
    prepre.r <- pre.r
    steps <- 0
    while(steps <= maxsteps &&
          !(pre.r == r && (!any(xor(presubset, subset)) || prepre.r == r))) {
	presubset <- subset
	prepre.r <- pre.r
	pre.r <- r
	tis <- abs(GiveTis(x, y, subset,
			   x[subset, , drop = FALSE],
			   y[subset])$tis)
	limit <- qt(alpha / (2*(r + 1)), r - p, lower.tail=FALSE)
	subset <- tis < limit
	r <- sum(subset)
	steps <- steps + 1
	if(verbose) trace1(steps, r, n)
    }
    list(subset = subset, tis = tis)
}

########## END .lmBACON
