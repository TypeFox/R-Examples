# $Id: mvt.R 327 2016-01-29 21:10:37Z bbnkmp $

##' Do we have a correlation matrix?
##' @param x typically a matrix
chkcorr <- function(x) {
    if (!is.matrix(x) || (d <- dim(x))[1] != d[2])
        return(FALSE)
    rownames(x) <- colnames(x) <- NULL
    storage.mode(x) <- "numeric"
    ONE <- 1 + sqrt(.Machine$double.eps)

    ## return
    -ONE <= min(x) && max(x) <= ONE && isTRUE(all.equal(diag(x), rep(1, d[1])))
}

checkmvArgs <- function(lower, upper, mean, corr, sigma)
{
    if (!is.numeric(lower) || !is.vector(lower))
        stop(sQuote("lower"), " is not a numeric vector")
    if (!is.numeric(upper) || !is.vector(upper))
        stop(sQuote("upper"), " is not a numeric vector")
    if (!is.numeric(mean) || !is.vector(mean))
        stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(lower) || any(is.na(lower)))
        stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper)))
        stop(sQuote("upper"), " not specified or contains NA")
    rec <- cbind(lower, upper, mean)# <--> recycling to same length
    lower <- rec[,"lower"]
    upper <- rec[,"upper"]
    if (!all(lower <= upper))
        stop("at least one element of ", sQuote("lower"), " is larger than ",
             sQuote("upper"))
    mean <- rec[,"mean"]
    if (any(is.na(mean)))
        stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
        corr <- diag(length(lower))
        # warning("both ", sQuote("corr"), " and ", sQuote("sigma"),
        # " not specified: using sigma=diag(length(lower))")
    }
    if (!is.null(corr) && !is.null(sigma)) {
        sigma <- NULL
        warning("both ", sQuote("corr"), " and ", sQuote("sigma"),
                " specified: ignoring ", sQuote("sigma"))
    }
    UNI <- FALSE
    if (!is.null(corr)) {
         if (!is.numeric(corr))
             stop(sQuote("corr"), " is not numeric")
         if (!is.matrix(corr)) {
             if (length(corr) == 1)
                UNI <- TRUE
             if (length(corr) != length(lower))
               stop(sQuote("diag(corr)"), " and ", sQuote("lower"),
                    " are of different length")
         } else {
             if (length(corr) == 1) {
                 UNI <- TRUE
                 corr <- corr[1,1]
                 if (length(lower) != 1)
                   stop(sQuote("corr"), " and ", sQuote("lower"),
                        " are of different length")
             } else {
                 if (length(diag(corr)) != length(lower))
                     stop(sQuote("diag(corr)"), " and ", sQuote("lower"),
                          " are of different length")
                 if (!chkcorr(corr))
                     stop(sQuote("corr"), " is not a correlation matrix")
             }
         }
    }
    if (!is.null(sigma)) {
         if (!is.numeric(sigma))
             stop(sQuote("sigma"), " is not numeric")
         if (!is.matrix(sigma)) {
            if (length(sigma) == 1)
                UNI <- TRUE
            if (length(sigma) != length(lower))
               stop(sQuote("diag(sigma)"), " and ", sQuote("lower"),
                    " are of different length")
         } else {
            if (length(sigma) == 1) {
                UNI <- TRUE
                sigma <- sigma[1,1]
                if (length(lower) != 1)
                  stop(sQuote("sigma"), " and ", sQuote("lower"),
                       " are of different length")
            } else {
              if (length(diag(sigma)) != length(lower))
                 stop(sQuote("diag(sigma)"), " and ", sQuote("lower"),
                      " are of different length")
              if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 0))
                 stop(sQuote("sigma"), " is not a covariance matrix")
            }
         }
    }
    list(lower=lower, upper=upper, mean=mean, corr=corr, sigma=sigma, uni=UNI)
}


pmvnorm <- function(lower=-Inf, upper=Inf, mean=rep(0, length(lower)), corr=NULL, sigma=NULL,
                    algorithm = GenzBretz(), ...)
{
    carg <- checkmvArgs(lower=lower, upper=upper, mean=mean, corr=corr,
                        sigma=sigma)
    if (!is.null(carg$corr)) {
      corr <- carg$corr
      if (carg$uni) {
          stop(sQuote("sigma"), " not specified: cannot compute pnorm")
      } else {
          lower <- carg$lower - carg$mean
          upper <- carg$upper - carg$mean
          mean <- rep(0, length(lower))
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     algorithm = algorithm, ...)
      }
    } else {
      if (carg$uni) {
        RET <- list(value = pnorm(carg$upper, mean=carg$mean, sd=sqrt(carg$sigma)) -
                            pnorm(carg$lower, mean=carg$mean, sd=sqrt(carg$sigma)),
                    error = 0, msg="univariate: using pnorm")
      } else {
          lower <- (carg$lower - carg$mean)/sqrt(diag(carg$sigma))
          upper <- (carg$upper - carg$mean)/sqrt(diag(carg$sigma))
          mean <- rep(0, length(lower))
          corr <- cov2cor(carg$sigma)
          RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                     algorithm = algorithm, ...)
      }
    }
    ## return
    structure(RET$value, "error" = RET$error, "msg" = RET$msg)
}

pmvt <- function(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
                 df=1, corr=NULL, sigma=NULL,
                 algorithm = GenzBretz(),
                 type = c("Kshirsagar", "shifted"), ...)
{
    type <- match.arg(type)
    carg <- checkmvArgs(lower=lower, upper=upper, mean=delta, corr=corr,
                        sigma=sigma)
    if (type == "shifted") { # can be handled by integrating over central t
      if(!is.null(carg$corr)){ # using transformed integration bounds
        d <- 1
      } else {
        if(!is.null(carg$sigma)){
          d <- sqrt(diag(carg$sigma))
          carg$corr <- cov2cor(carg$sigma)
        }
      }
      carg$lower <- (carg$lower - carg$mean)/d
      carg$upper <- (carg$upper - carg$mean)/d
      carg$mean <- rep(0, length(carg$mean))
    }

    if (is.null(df))
        stop(sQuote("df"), " not specified")
    if (df < 0) # MH: was any(..)
        stop("cannot compute multivariate t distribution with ",
             sQuote("df"), " < 0")
    if(is.finite(df) && (df != as.integer(df))) # MH: was !isTRUE(all.equal(as.integer(df), df))
        stop(sQuote("df"), " is not an integer")
    if (carg$uni) {
        if (df > 0) # df = Inf is taken care of by pt()
            RET <- list(value = pt(carg$upper, df=df, ncp=carg$mean) -
                                pt(carg$lower, df=df, ncp=carg$mean),
                       error = 0, msg="univariate: using pt")
        else
            RET <- list(value = pnorm(carg$upper, mean = carg$mean) -
                                pnorm(carg$lower, mean=carg$mean),
                       error = 0, msg="univariate: using pnorm")
    } else { # mvt() takes care of df = 0 || df = Inf
        if (!is.null(carg$corr)) {
            RET <- mvt(lower=carg$lower, upper=carg$upper, df=df, corr=carg$corr,
                       delta=carg$mean,  algorithm = algorithm, ...)
        } else { # need to transform integration bounds and delta
            d <- sqrt(diag(carg$sigma))
            lower <- carg$lower/d
            upper <- carg$upper/d
            corr <- cov2cor(carg$sigma)
            RET <- mvt(lower=lower, upper=upper, df=df, corr=corr,
                       delta=carg$mean/d, algorithm = algorithm, ...)
        }
    }
    attr(RET$value, "error") <- RET$error
    attr(RET$value, "msg") <- RET$msg
    return(RET$value)
}

## identical(., Inf) would be faster but not vectorized
isInf <- function(x) x > 0 & is.infinite(x) # check for  Inf
isNInf <- function(x) x < 0 & is.infinite(x) # check for -Inf

mvt <- function(lower, upper, df, corr, delta, algorithm = GenzBretz(), ...)
{

    ### only for compatibility with older versions
    addargs <- list(...)
    if (length(addargs) > 0)
        algorithm <- GenzBretz(...)
    else if (is.function(algorithm) || is.character(algorithm))
        algorithm <- do.call(algorithm, list())

    ### handle cases where the support is the empty set
    ##  Note: checkmvArgs() has been called ==> lower, upper are *not* NA
    if (any(abs(d <- lower - upper) < sqrt(.Machine$double.eps)*(abs(lower)+abs(upper)) |
            lower == upper)) ## e.g. Inf == Inf
	return(list(value = 0, error = 0, msg = "lower == upper"))

    n <- ncol(corr)
    if (is.null(n) || n < 2) stop("dimension less then n = 2")

    if (length(lower) != n) stop("wrong dimensions")
    if (length(upper) != n) stop("wrong dimensions")

    if (n > 1000) stop("only dimensions 1 <= n <= 1000 allowed")

    infin <- rep(2, n)
    infin[ isInf(upper)] <- 1
    infin[isNInf(lower)] <- 0
    infin[isNInf(lower) & isInf(upper)] <- -1

    ### fix for Miwa algo:
    ### pmvnorm(lower=c(-Inf, 0, 0), upper=c(0, Inf, Inf),
    ###         mean=c(0, 0, 0), sigma=S, algorithm = Miwa())
    ###         returned NA

    if (inherits(algorithm, "Miwa")) {
        if (n >= 3 && any(infin == -1)) {
            WhereBothInfIs <- which(infin == -1)
            n <- n - length(WhereBothInfIs)
            corr <- corr[-WhereBothInfIs, -WhereBothInfIs]
            upper <- upper[-WhereBothInfIs]
            lower <- lower[-WhereBothInfIs]
        }

        if (n >= 2 && any(infin == 0)) {
            WhereNegativInfIs <- which(infin==0)
            inversecorr <- rep(1, n)
            inversecorr[WhereNegativInfIs] <- -1
            corr <- diag(inversecorr) %*% corr %*% diag(inversecorr) ## MM_FIXME
            infin[WhereNegativInfIs] <- 1

            tempsaveupper <- upper[WhereNegativInfIs]
            upper[WhereNegativInfIs] <- -lower[WhereNegativInfIs]
            lower[WhereNegativInfIs] <- -tempsaveupper
        }
    }

    ### this is a bug in `mvtdst' not yet fixed
    if (all(infin < 0))
        return(list(value = 1, error = 0, msg = "Normal Completion"))

    if (n > 1) {
        corrF <- matrix(as.vector(corr), ncol=n, byrow=TRUE)
        corrF <- corrF[upper.tri(corrF)]
    } else corrF <- corr


    ret <- probval(algorithm, n, df, lower, upper, infin, corr, corrF, delta)
    inform <- ret$inform
    msg <-
        if (inform == 0) "Normal Completion"
    else if (inform == 1) "Completion with error > abseps"
    else if (inform == 2) "N greater 1000 or N < 1"
    else if (inform == 3) "Covariance matrix not positive semidefinite"
    else inform

    ## return including error est. and msg:
    list(value = ret$value, error = ret$error, msg = msg)
}

rmvt <- function(n, sigma = diag(2), df = 1,
                 delta = rep(0, nrow(sigma)),
                 type = c("shifted", "Kshirsagar"), ...)
{
    if (length(delta) != nrow(sigma))
        stop("delta and sigma have non-conforming size")
    if (hasArg(mean)) # MH: normal mean variance mixture != t distribution (!)
        stop("Providing 'mean' does *not* sample from a multivariate t distribution!")
    if (df == 0 || isInf(df)) # MH: now (also) properly allow df = Inf
        return(rmvnorm(n, mean = delta, sigma = sigma, ...))
    type <- match.arg(type)
    switch(type,
           "Kshirsagar" = {
               return(rmvnorm(n, mean = delta, sigma = sigma, ...)/
                      sqrt(rchisq(n, df)/df))
           },
           "shifted" = {
               sims <- rmvnorm(n, sigma = sigma, ...)/sqrt(rchisq(n, df)/df)
               return(sweep(sims, 2, delta, "+"))
           },
           stop("wrong 'type'"))
}

dmvt <- function(x, delta = rep(0, p), sigma = diag(p), df = 1,
                 log = TRUE, type = "shifted")
{
    if (is.vector(x))
        x <- matrix(x, ncol = length(x))
    p <- ncol(x)
    if (df == 0 || isInf(df)) # MH: now (also) properly allow df = Inf
        return(dmvnorm(x, mean = delta, sigma = sigma, log = log))
    if(!missing(delta)) {
	if(!is.null(dim(delta))) dim(delta) <- NULL
	if (length(delta) != p)
	    stop("delta and sigma have non-conforming size")
    }
    if(!missing(sigma)) {
	if (p != ncol(sigma))
	    stop("x and sigma have non-conforming size")
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
			 check.attributes = FALSE))
	    stop("sigma must be a symmetric matrix")
    }
    type <- match.arg(type)

    dec <- tryCatch(chol(sigma), error=function(e)e)
    if (inherits(dec, "error")) {
	x.is.d <- colSums(t(x) != delta) == 0
	logretval <- rep.int(-Inf, nrow(x))
	logretval[x.is.d] <- Inf # and all other f(.) == 0
    } else {
	R.x_m <- backsolve(dec, t(x) - delta, transpose = TRUE)
	rss <- colSums(R.x_m ^ 2)
	logretval <- lgamma((p + df)/2) -
	    (lgamma(df / 2) + sum(log(diag(dec))) + p/2 * log(pi * df)) -
		0.5 * (df + p) * log1p(rss / df)
    }
    names(logretval) <- rownames(x)
    if (log) logretval else exp(logretval)
}

## get start interval for root-finder used in qmvnorm and qmvt
getInt <- function(p, delta, sigma, tail,
                   type = c("Kshirsagar", "shifted"), df){
  type <- match.arg(type)
  sds <- sqrt(diag(sigma))
  if(df == 0 | df == Inf){
    df <- Inf
    cdf <- function(x, ...)
      pnorm(x, delta, sds, ...)
  } else {
    if(type == "shifted"){
      cdf <- function(x, ...)
        pt((x-delta)/sds, df=df, ...)
    }
    if(type == "Kshirsagar"){
      cdf <- function(x, ...)
        pt(x, ncp=delta/sds, df=df, ...)
    }
  }
  switch(tail, both.tails = {
    interval <- c(0,10)
    func <- function(x, delta, sds)
      prod(cdf(x)-cdf(-x))
    UB <- max(abs(delta+sds*qt(1-(1-p)/2, df=df)))
  }, upper.tail = {
    interval <- c(-10,10)
    func <- function(x, delta, sds)
      prod(cdf(x, lower.tail=FALSE))
    UB <- min(delta+sds*qt(1-p, df=df))
  }, lower.tail = {
    interval <- c(-10,10)
    func <- function(x, delta, sds)
      prod(cdf(x))
    UB <- max(delta+sds*qt(p, df=df))
  }, )
  LB <- uniroot(function(x)
                func(x, delta=delta, sds=sds)-p,
                interval, extendInt = "yes")$root
  sort(c(LB, UB))
}


qmvnorm <- function(p, interval = NULL,
                    tail = c("lower.tail", "upper.tail", "both.tails"),
                    mean = 0, corr = NULL, sigma = NULL, algorithm =
                    GenzBretz(),
                    ptol = 0.001, maxiter = 500, trace = FALSE, ...)
{
    if (length(p) != 1 || p < 0 || p > 1)
        stop(sQuote("p"), " is not a double between zero and one")

    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm) && !is.null(algorithm))
        algorithm <- dots$algorithm

    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5)
        stop("cannot compute two-sided quantile for p < 0.5")
    dim <- length(mean)
    if (is.matrix(corr)) dim <- nrow(corr)
    if (is.matrix(sigma)) dim <- nrow(sigma)
    lower <- upper <- rep.int(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    if (args$uni) {
        if (is.null(args$sigma))
          stop(sQuote("sigma"), " not specified: cannot compute qnorm")
        if (tail == "both.tails") p <- ifelse(p < 0.5, p / 2, 1 - (1 - p)/2)
        q <- qnorm(p, mean = args$mean, sd = args$sigma,
                   lower.tail = (tail != "upper.tail"))
        return( list(quantile = q, f.quantile = p) )
    }
    dim <- length(args$mean)

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)

    pfct <- function(q) {
        ### use the same seed for different values of q
        assign(".Random.seed", R.seed, envir = .GlobalEnv)
        switch(tail, "both.tails" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep( abs(q), dim)
           }, "upper.tail" = {
                  low <- rep(      q, dim)
                  upp <- rep(    Inf, dim)
           }, "lower.tail" = {
                  low <- rep(   -Inf, dim)
                  upp <- rep(      q, dim)
           },)
           ret <- pmvnorm(lower = low, upper = upp, mean = args$mean,
                          corr = args$corr, sigma = args$sigma,
                          algorithm = algorithm)
           if(tail == "upper.tail") ## get_quant_loclin assumes an increasing function
             ret <- 1-ret
           return(ret)
    }
    if(is.null(interval)){
      if(is.null(args$sigma)){
        sig <- args$corr
      } else {
        sig <- args$sigma
      }
      interval <- getInt(p=p, delta=args$mean, sigma=sig,
                         tail=tail, df=Inf)
      dif <- diff(interval)
      interval <- interval+c(-1,1)*0.2*max(dif,0.1) ## extend range slightly
    } 
    if(tail == "upper.tail") ## get_quant_loclin assumes an increasing function
      p <- 1-p
    qroot <- get_quant_loclin(pfct, p, interval=interval,
                              link="probit",
                              ptol=ptol, maxiter=maxiter, verbose=trace)
    qroot$f.quantile <- qroot$f.quantile - p
    qroot
}

qmvt <- function(p, interval = NULL,
                 tail = c("lower.tail", "upper.tail", "both.tails"),
                 df = 1, delta = 0, corr = NULL, sigma = NULL,
                 algorithm = GenzBretz(),
                 type = c("Kshirsagar", "shifted"),
                 ptol = 0.001, maxiter = 500, trace = FALSE, ...) {

    if (length(p) != 1 || (p <= 0 || p >= 1))
        stop(sQuote("p"), " is not a double between zero and one")

    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm)  && !is.null(algorithm))
        algorithm <- dots$algorithm
    type <- match.arg(type)

    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5)
        stop("cannot compute two-sided quantile for p < 0.5")
    dim <- 1
    if (!is.null(corr)) dim <- NROW(corr)
    if (!is.null(sigma)) dim <- NROW(sigma)
    lower <- upper <- rep.int(0, dim)
    args <- checkmvArgs(lower, upper, delta, corr, sigma)
    if (args$uni) {
        if (tail == "both.tails") p <- ifelse(p < 0.5, p / 2, 1 - (1 - p)/2)
        if (df == 0 || isInf(df)) { # MH: now (also) properly allow df = Inf
            q <- qnorm(p, mean = args$mean, lower.tail = (tail != "upper.tail"))
        } else {
            q <- qt(p, df = df, ncp = args$mean, lower.tail = (tail != "upper.tail"))
        }
        qroot <- list(quantile = q, f.quantile = p)
        return(qroot)
    }

    dim <- length(args$mean)

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)

    pfct <- function(q) {
        ### use the same seed for different values of q
        assign(".Random.seed", R.seed, envir = .GlobalEnv)
        switch(tail, "both.tails" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep( abs(q), dim)
           }, "upper.tail" = {
                  low <- rep(      q, dim)
                  upp <- rep(    Inf, dim)
           }, "lower.tail" = {
                  low <- rep(   -Inf, dim)
                  upp <- rep(      q, dim)
           },)
           ret <- pmvt(lower = low, upper = upp, df = df, delta = args$mean,
                       corr = args$corr, sigma = args$sigma,
                       algorithm = algorithm, type = type)
        if(tail == "upper.tail") ## get_quant_loclin assumes an increasing function
          ret <- 1-ret
        return(ret)
    }
    if(is.null(interval)){
      if(is.null(args$sigma)){
        sig <- args$corr
      } else {
        sig <- args$sigma
      }
      interval <- getInt(p=p, delta=args$mean, sigma=sig,
                         tail=tail, type=type, df=df)
      dif <- diff(interval)
      interval <- interval+c(-1,1)*0.2*max(dif,0.1) ## extend range slightly
    } 
    if(tail == "upper.tail") ## get_quant_loclin assumes an increasing function
      p <- 1-p
    link <- ifelse(df <= 7 & df > 0, "cauchit", "probit")
    qroot <- get_quant_loclin(pfct, p, interval=interval,
                              link=link,
                              ptol=ptol, maxiter=maxiter, verbose=trace)
    qroot$f.quantile <- qroot$f.quantile - p
    qroot
}


GenzBretz <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    structure(list(maxpts = maxpts, abseps = abseps, releps = releps),
              class = "GenzBretz")
}

Miwa <- function(steps = 128) {
    if (steps > 4098) stop("maximum number of steps is 4098")
    structure(list(steps = steps), class = "Miwa")
}

probval <- function(x, ...)
    UseMethod("probval")

probval.GenzBretz <- function(x, n, df, lower, upper, infin, corr, corrF, delta) {

    if(isInf(df)) df <- 0 # MH: deal with df=Inf (internally requires df=0!)

    lower[isNInf(lower)] <- 0
    upper[ isInf(upper)] <- 0

    error <- 0; value <- 0; inform <- 0
    .C("C_mvtdst",
       N = as.integer(n),
       NU = as.integer(df),
       LOWER = as.double(lower),
       UPPER = as.double(upper),
       INFIN = as.integer(infin),
       CORREL = as.double(corrF),
       DELTA = as.double(delta),
       MAXPTS = as.integer(x$maxpts),
       ABSEPS = as.double(x$abseps),
       RELEPS = as.double(x$releps),
       error = as.double(error),
       value = as.double(value),
       inform = as.integer(inform),
       RND = as.integer(1)) ### init RNG
}

probval.Miwa <- function(x, n, df, lower, upper, infin, corr, corrF, delta) {

    if (!( df==0 || isInf(df) ))
        stop("Miwa algorithm cannot compute t-probabilities")

    if (n > 20)
        stop("Miwa algorithm cannot compute probabilities for dimension n > 20")

    sc <- try(solve(corr))
    if (inherits(sc, "try-error"))
        stop("Miwa algorithm cannot compute probabilities for singular problems")

    p <- .Call("C_miwa", steps = as.integer(x$steps),
                         corr = as.double(corr),
                         upper = as.double(upper),
                         lower = as.double(lower),
                         infin = as.integer(infin))
    list(value = p, inform = 0, error = NA)
}

dots2GenzBretz <- function(...) {
    addargs <- list(...)
    fm1 <- sapply(names(addargs), function(x) length(grep(x, names(formals(GenzBretz)))) == 1)
    fm2 <- sapply(names(addargs), function(x) length(grep(x, names(formals(uniroot)))) == 1)
    algorithm <- NULL
    uniroot <- NULL
    if (any(fm1))
        algorithm <- do.call("GenzBretz", addargs[fm1])
    if (any(fm2))
        uniroot <- addargs[fm2]
    list(algorithm = algorithm, uniroot = uniroot)
}
