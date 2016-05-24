
### Based on polr() function in MASS (originally developed for categorical wage data in Social Capital Benchmark Survey)
intReg <- function(formula, start, boundaries,
                   ...,
                 contrasts = NULL, Hess = FALSE,
                 model = TRUE,
                 method = c("probit", "logistic", "cloglog", "cauchit", "model.frame"),
                     print.level=0,
                   data, subset, weights, na.action,
                   iterlim=100)
{
   ## beta     parameters for the x-s (linear model)
   ## nBeta    number of x-s (including constant)
    logit <- function(p) log(p/(1 - p))
    ## log likelihood
    loglik <- function(theta) {
       beta <- theta[iBeta]
        zeta <- theta[iBoundaries]
        sd <- theta[iStd]
        eta <- offset
        if (nBeta > 0)
            eta <- eta + drop(x %*% beta)
       Pr <- pfun((zeta[boundaryInterval + 1] - eta)/sd) - pfun((zeta[boundaryInterval] - eta)/sd)
       if(any(Pr <= 0))
           return(NA)
       names(Pr) <- NULL
                                        # individual probability to be in interval (zeta[y+1], zeta[y]])
       ll <- wt * log(Pr)
       (ll)
    }
    ## gradient
    gradlik <- function(theta)
    {
        jacobian <- function(theta) { ## dzeta by dtheta matrix
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0 , k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
         }
        beta <- theta[iBeta]
        zeta <- theta[iBoundaries]
        sigma <- theta[iStd]
        eta <- offset
        if(nBeta > 0)
            eta <- eta + drop(x %*% beta)
        normArg1 <- (zeta[boundaryInterval +1] - eta)/sigma
        normArg2 <- (zeta[boundaryInterval] - eta)/sigma
        Pr <- pfun(normArg1) - pfun(normArg2)
        p1 <- dfun(normArg1)
        p2 <- dfun(normArg2)
        g1 <-
                                        # d loglik / d beta
            if(nBeta > 0)
                x * (wt*(p2 - p1)/Pr/sigma)
            else
                matrix(0, 0, nBeta)
        xx <- .polrY1*p1 - .polrY2*p2
        dg.dzeta <- xx * (wt/Pr/sigma)
##         g2 <- - t(xx) %*% (wt/Pr/sigma)
##         g2 <- t(g2) %*% jacobian(zeta)
##         if(all(Pr > 0))
##             -c(g1, g2)
##         else
##             rep(NA, nBeta+nInterval)
        s1 <- p1*normArg1
        s1[is.infinite(normArg1)] <- 0
                                        # if a boundary is Inf, we get 0*inf type of NaN
        s2 <- p2*normArg2
        s2[is.infinite(normArg2)] <- 0
        dg.dsd <- -(s1 - s2)*(wt/Pr/sigma)
        grad <- cbind(g1, dg.dzeta, dg.dsd)
        grad
     }
    ## ---------- main function -------------
    cl <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    if(is.matrix(eval.parent(cl$data)))
        cl$data <- as.data.frame(data)
    cl$start <- cl$Hess <- cl$method <- cl$model <- cl$boundaries <- cl$... <- cl$print.level <- NULL
    cl[[1]] <- as.name("model.frame")
    mf <- eval.parent(cl)
    if (method == "model.frame")
        return(mf)
    mt <- attr(mf, "terms")
    ## Select the correct model
    pfun <- switch(method, logistic = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    x <- model.matrix(mt, mf, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    nObs <- nrow(x)
    nBeta <- ncol(x)
    cons <- attr(x, "contrasts") # will get dropped by subsetting
##     if(xint > 0) {
##        x <- x[, -xint, drop=FALSE]
##     }
##     else
##         warning("an intercept is needed and assumed")
    wt <- model.weights(mf)
    if(!length(wt)) wt <- rep(1, nObs)
    offset <- model.offset(mf)
    if(length(offset) <= 1) offset <- rep(0, nObs)
    y <- model.response(mf)
    if(is.matrix(y)) {
       ## Use the intervals, given for each observation.  We have interval regression and the interval boundaries are fixed.
       ## Save boundaries as a sequence of pairs of L,B boundaries
       ordered <- FALSE
       dimnames(y) <- NULL
       lowerBound <- y[,1]
       upperBound <- y[,2]
       ## in case of interval regression, we have to construct a set of intervals and pack them correctly to the
       ## parameter vector
       intervals <- sets::set()
       for(i in 1:length(lowerBound)) {
          intervals <- intervals | c(lowerBound[i], upperBound[i])
       }
       ## Now make the set to a list to have ordering
       intervals <- as.list(intervals)
       ## Now find which interval each observation falls into
       y <- boundaryInterval <- numeric(length(lowerBound))
                                        # which interval observation falls to 
                                        # y:  in terms of ordered intervals
                                        # boundaryInterval         in terms of ordered boundaries
       for(i in seq(along=intervals)) {
          j <- lowerBound == intervals[[i]][1] & upperBound == intervals[[i]][2]
          boundaryInterval[j] <- 1 + 2*(i - 1)
          y[j] <- i
       }
       boundaries <- unlist(intervals)
                                        # boundaries as a vector (note the joint boundaries are twice
       names(boundaries) <- paste(c("L", "U"), rep(seq(along=intervals), each=2))
       nInterval <- length(boundaries) - 1
    }
    else {
       ## response is given as an ordered factor, boundaries must be given separately.
       ## Save them as a vector of boundaries, all numbers (except the first, last) represent the upper boundary of
       ## the smaller interval and the lower boundary of the upper interval at the same time.
       lev <- levels(y)
       if(length(lev) <= 2)
           stop("response must have 3 or more levels")
       if(!is.factor(y))
           stop("response must be a factor or Nx2 matrix of boundaries")
       y <- unclass(y)
       nInterval <- length(lev)
       if(missing(boundaries))
           ordered <- TRUE
       else {
           ordered <- FALSE
           intervals <- vector("list", length(boundaries) - 1)
           for(i in seq(length=length(boundaries) - 1)) {
              intervals[[i]] <- c(boundaries[i], boundaries[i+1])
           }
           boundaryInterval <- y
                                        # y falls inbetween boundaries 'boundaryInterval' and 'boundaryInterval + 1'
        }
       if(is.null(names(boundaries)))
           names(boundaries) <- paste("Boundary", seq(along=boundaries))
    }
    Y <- matrix(0, nObs, nInterval + 1)
    .polrY1 <- col(Y) == boundaryInterval + 1
    .polrY2 <- col(Y) == boundaryInterval 
                                        # .polr are markers for which interval the boundaryInterval falls to
    ## starting values
    iBeta <- seq(length=ncol(x))
                           # coefficients
    iBoundaries <- nBeta + seq(along=boundaries)
                           # boundaries
    iStd <- max(iBoundaries) + 1
                           # standard deviation
    if(missing(start)) {
       start <- numeric(max(iBeta, iBoundaries, iStd))
                           # +1 for the error variance 'sigma'
       activePar <- logical(length(start))
       if(ordered) {
          ## try logistic/probit regression on 'middle' cut
          q1 <- nInterval %/% 2
          y1 <- (y > q1)
          fit <-
              switch(method,
                     "logistic"= glm.fit(x, y1, wt, family = binomial(), offset = offset),
                     "probit" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                     ## this is deliberate, a better starting point
                     "cloglog" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                     "cauchit" = glm.fit(x, y1, wt, family = binomial("cauchit"), offset = offset))
          if(!fit$converged)
              warning("attempt to find suitable starting values did not converge")
          coefs <- fit$coefficients
          if(any(is.na(coefs))) {
             warning("design appears to be rank-deficient, so dropping some coefs")
             keep <- names(coefs)[!is.na(coefs)]
             coefs <- coefs[keep]
#          x <- x[, keep[-1], drop = FALSE]
             ## note: we keep the intercept
             nBeta <- ncol(x)
          }
          spacing <- logit((1:nInterval)/(nInterval+1)) # just a guess
          if(method != "logit") spacing <- spacing/1.7
          zetas <- -coefs[2] + spacing - spacing[q1]
          coefs[1] <- 0
          activePar <- c(FALSE, rep(TRUE, length(coef) - 1 + length(zetas) + 1))
                                        # intercept is fixed to 0
          sigma <- 1
       }
       else {
          ## not ordered: estimate OLS on interval middle points
          means <- sapply(intervals, mean)
                                        # we have to put a reasonable value to infinite intervals.
                                        # Pick the average width of the interval and use it as the meanpoint
          widths <- sapply(intervals, function(x) x[2] - x[1])
          meanWidth <- mean(widths[!is.infinite(widths)])
          negInf <- is.infinite(means) & means < 0
          if(any(negInf)) {
                                        # if none is true, sapply returns 'list()' and transforms means to a list
             means[negInf] <- sapply(intervals[negInf], function(x) x[2] - meanWidth)
          }
          posInf <- is.infinite(means) & means > 0
          if(any(posInf)) {
             means[posInf] <- sapply(intervals[posInf], function(x) x[1] + meanWidth)
          }
          yMean <- means[y]
          fit <- lm(yMean ~ x - 1)
          xCoefs <- coef(fit)
          if(any(is.na(xCoefs))) {
             cat("Suggested initial values:\n")
             print(xCoefs)
             stop("NA in the initial values")
          }
          names(xCoefs) <- gsub("^x", "", names(xCoefs))
          sigma <- sqrt(var(fit$residuals))
       }
       start[iBeta] <- xCoefs
       names(start)[iBeta] <- names(xCoefs)
       start[iBoundaries] <- boundaries
       names(start)[iBoundaries] <- names(boundaries)
       start[iStd] <- sigma
       names(start)[iStd] <- "sigma"
    }
    else
        if(length(start) != iStd)
            stop("'start' is of wrong length:\n",
                 "The current model includes ", nBeta,
                 " explanatory variables plus\n",
                 length(iBeta), " interval boundaries ",
                 "plus 1 disturbance standard deviation\n",
                 "(", iStd, " in total).\n",
                 "However, 'start' is of length ",
                 length(start))
    if(print.level > 0) {
       cat("Initial values:\n")
       print(start)
    }
    if(!ordered) {
       ## Not ordered model: fix the fixed parameters
       activePar <- logical(length(start))
       activePar[iBeta] <- TRUE
       activePar[iBoundaries] <- FALSE
       activePar[iStd] <- TRUE
    }
##     compareDerivatives(loglik, gradlik, t0=start)
##     stop()
    estimation <- maxLik(loglik, gradlik, start=start,
                  method="BHHH", activePar=activePar, iterlim=iterlim, ...)
    res <- c(estimation,
             param=list(list(ordered=ordered,
                      boundaries=boundaries,
                      index=list(beta=iBeta, boundary=iBoundaries, std=iStd),
                      df=nObs - sum(activePar),
                      nObs=nObs
                      )),
             call=cl,
             terms=mt,
             method=method,
             na.action=list(attr(mf, "na.action"))
             )
    class(res) <- c("intReg", class(estimation))
    return(res)
}
