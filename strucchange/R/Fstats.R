Fstats <- function(formula, from = 0.15, to = NULL, data = list(), vcov. = NULL)
{
  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- length(y)
  e <- lm.fit(X,y)$residuals

  ## check if tsp are available and may be used
  ## (potentially not if NAs were removed)
  ytsp <- NULL
  orig.y <- NULL
  tsp_ok <- FALSE
  if(is.ts(data)){
      if(NROW(data) == n) {
          ytime <- time(data)
          ytsp <- tsp(data)
	  tsp_ok <- TRUE
      }
  } else {
      env <- environment(formula)
      if(missing(data)) data <- env
      orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
      if(is.ts(orig.y) & (NROW(orig.y) == n)){
	    ytime <- time(orig.y)
	    ytsp <- tsp(orig.y)
  	    tsp_ok <- TRUE
      }
  }

  ts.eps <- getOption("ts.eps")

  if(length(from) > 1) {
    if(!is.null(ytsp) && from[2] <= ytsp[3]) {
      from <- which(abs(ytime-(from[1]+(from[2]-1)/ytsp[3])) < ts.eps)
      if(!is.null(to)) to <- which(abs(ytime-(to[1]+(to[2]-1)/ytsp[3])) < ts.eps)
    } else {
      stop(paste(sQuote("from"), "does not specify a valid time point"))
    }
  }
  else if(from < 1)
  {
    from <- floor(from*n)
    if(!is.null(to)) to <- floor(to*n)
  }
  if(is.null(to)) to <- n - from
  if(from < (k+1))
  {
    from <- k+1
    warning("'from' changed (was too small)")
  }
  if(to > (n-k-1))
  {
    to <- n-k-1
    warning("'to' changed (was too large)")
  }
  if(from <= to)
    point <- (from:to)
  else
    stop("inadmissable change points: 'from' is larger than 'to'")

  sume2 <- sum(e^2)
  lambda <- ((n-from)*to)/(from*(n-to))
  np <- length(point)
  stats <- rep(0,np)
  for(i in 1:np)
  {
    X1 <- as.matrix(X[(1:point[i]),])
    X2 <- as.matrix(X[((point[i]+1):n),])

    if(is.null(vcov.)) {
      fm1 <- lm.fit(X1,y[1:point[i]])
      fm2 <- lm.fit(X2,y[((point[i]+1):n)])
      u <- c(fm1$residuals, fm2$residuals)
      sigma2 <- (sum(u^2))/(n-2*k)
      stats[i] <- (sume2-sum(u^2))/sigma2
    }
    else {
      allX <- cbind(X1, matrix(rep(0, point[i]*k), ncol=k))
      allX <- rbind(allX, cbind(X2, X2))
      fm2 <- lm(y ~ 0 + allX)
      beta2 <- coef(fm2)[-(1:k)]
      V <- vcov.(fm2)
      stats[i] <- as.vector(t(beta2) %*% chol2inv(chol(V[-(1:k),-(1:k)])) %*% beta2)
     }
  }

  sup.point <- which.max(stats) + from - 1
  if(is.null(vcov.))
    min.RSS <- sume2/(1 + max(stats)/(n - 2*k))
  else
    min.RSS <- NA

  if(is.ts(data) & tsp_ok){
      stats <- ts(stats, start = time(data)[from], frequency = frequency(data))
      datatsp <- tsp(data)
  }
  else if(!is.null(orig.y) & tsp_ok) {
      stats <- ts(stats, start = time(orig.y)[from], frequency = frequency(orig.y))
      datatsp <- tsp(orig.y)
  }
  else{
      stats <- ts(stats, start = from/n, frequency = n)
      datatsp <- c(0, 1, n)
  }

  retval <- list(Fstats = stats,
                 nreg = k,
                 nobs = n,
                 par = lambda,
                 call = match.call(),
                 formula = formula,
		 breakpoint = sup.point,
		 RSS = min.RSS,
                 datatsp = datatsp)

  class(retval) <- "Fstats"
  return(retval)
}

print.Fstats <- function(x, ...)
{
    cat("\nF statistics \n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
}

sctest.Fstats <- function(x, type = c("supF", "aveF", "expF"), asymptotic = FALSE, ...)
{
    dname <- paste(deparse(substitute(x)))
    type <- match.arg(type)
    switch(type,
           supF = {
               STATISTIC <- max(x$Fstats)
               names(STATISTIC) <- "sup.F"
               METHOD <- "supF test"
           },
           aveF = {
               STATISTIC <- mean(x$Fstats)
               names(STATISTIC) <- "ave.F"
               METHOD <- "aveF test"
           },
           expF = {
               STATISTIC <- log(mean(exp(0.5*x$Fstats)))
               names(STATISTIC) <- "exp.F"
               METHOD <- "expF test"
           })
    if((x$par == 1) & !(type == "expF") & !asymptotic)
    {
        METHOD <- "Chow test"
        PVAL <- 1 - pf(STATISTIC, x$nreg, (x$nobs-2*x$nreg))
    }
    else
        PVAL <- pvalue.Fstats(STATISTIC, type = type,
                              k=x$nreg, lambda=x$par)

    RVAL <- list(statistic = STATISTIC, p.value = PVAL,
                 method = METHOD, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

plot.Fstats <- function(x, pval = FALSE, asymptotic = FALSE,
                        alpha = 0.05, boundary = TRUE, aveF = FALSE,
                        xlab = "Time", ylab = NULL,
                        ylim = NULL, ...)
{
    k <- x$nreg
    n <- x$nobs
    bound <- boundary(x, alpha = alpha, pval = pval, aveF = aveF, asymptotic =
                      asymptotic)
    x <- x$Fstats

    if(pval)
    {
        if(asymptotic)
            x <- 1 - pchisq(x, k)
        else
            x <- 1 - pf(x, k, (n-2*k))
        if(is.null(ylab)) ylab <- "p values"
    }
    else
        if(is.null(ylab)) ylab <- "F statistics"
    if(is.null(ylim)) ylim <- c(0, max(c(x,bound)))
    plot(x, xlab = xlab, ylab = ylab, ylim = ylim, ...)
    abline(0,0)
    if(boundary)
    {
        lines(bound, col=2)
        if(aveF) lines(ts(rep(mean(x),length(x)),start=start(x),
                       frequency = frequency(x)),lty=2)
    }
}

boundary.Fstats <- function(x, alpha = 0.05, pval = FALSE, aveF =
                            FALSE, asymptotic = FALSE, ...)
{
    if(aveF)
    {
      myfun <-  function(y) {pvalue.Fstats(y, type="ave", x$nreg, x$par) - alpha}
      upper <- 40
    }
    else
    {
      myfun <-  function(y) {pvalue.Fstats(y, type="sup", x$nreg, x$par) - alpha}
      upper <- 80
    }
    bound <- uniroot(myfun, c(0,upper))$root
    if(pval)
    {
        if(asymptotic)
            bound <- 1 - pchisq(bound, x$nreg)
        else
            bound <- 1 - pf(bound, x$nreg, (x$nobs-2*x$nreg))
    }
    bound <- ts(bound,
                start = start(x$Fstats),
                end = end(x$Fstats),
                frequency = frequency(x$Fstats))
    return(bound)
}

lines.Fstats <- function(x, ...)
{
    lines(x$Fstats, ...)
}

