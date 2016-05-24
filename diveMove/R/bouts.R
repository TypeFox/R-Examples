
"logit" <- function(p) log(p / (1 - p))

"unLogit" <- function(logit) exp(logit) / (exp(logit) + 1)

"boutfreqs" <- function(x, bw, method=c("standard", "seq.diff"),
                        plot=TRUE, ...)
{
    ## Value: data frame with log frequencies and bin mid-points
    ## --------------------------------------------------------------------
    ## Arguments: x=numeric vector, bw=bin width for histogram,
    ## method=method used to construct the histogram, plot=logical whether
    ## to plot or not; ...=arguments passed to hist (must exclude 'breaks'
    ## and 'include.lowest')
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    method <- match.arg(method)
    switch(method,
           standard = {upper <- max(x, na.rm=TRUE)
                       brks <- seq(min(x, na.rm=TRUE), upper, bw)
                       if (brks[length(brks)] < upper) {
                           brks <- c(brks, brks[length(brks)] + bw)
                       }
                       h <- hist(x, breaks=brks, include.lowest=TRUE,
                                 plot=plot, ...)},
           seq.diff = {diff.x <- abs(diff(x))
                       upper <- max(diff.x, na.rm=TRUE)
                       brks <- seq(0, upper, bw)
                       if (brks[length(brks)] < upper) {
                           brks <- c(brks, brks[length(brks)] + bw)
                       }
                       h <- hist(diff.x, breaks=brks, include.lowest=TRUE,
                                 plot=plot, ...)})
    ok <- which(h$counts > 0)
    freq.adj <- h$counts[ok] / diff(c(0, ok))
    data.frame(lnfreq=log(freq.adj), x=h$mids[ok])
}

"boutinit" <- function(lnfreq, x.break, plot=TRUE, ...)
{
    ## Value: list with starting values for nls bout function
    ## --------------------------------------------------------------------
    ## Arguments: lnfreq=data frame with 'lnfreq' (log frequencies) and 'x'
    ## (midpoints); x.break=vector of length 1 or 2 with x value(s)
    ## defining the break(s) point(s) for broken stick model, such that x <
    ## x.break[1] is 1st process, and x >= x.break[1] & x < x.break[2] is
    ## 2nd one, and x >= x.break[2] is 3rd one; plot=logical whether to
    ## plot or not; ... arguments passed to lattice's xyplot()
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    nproc <- length(x.break)
    if (nproc > 2) stop ("x.break must be of length 1 or 2")
    procf <- cut(lnfreq$x, breaks=c(min(lnfreq$x), x.break, max(lnfreq$x)),
                 include.lowest=TRUE, right=TRUE)
    coefs <- by(lnfreq, procf, function(k) {coef(lm(lnfreq ~ x, k))})
    pars <- lapply(coefs, function(p) {
        lambda <- as.vector(-p[2])
        a <- as.vector(exp(p[1]) / lambda)
        c(a=a, lambda=lambda)
    })
    if (plot) {
        requireNamespace("lattice", quietly=TRUE) ||
            stop("lattice package is not available")
        pp <- lattice::xyplot(lnfreq ~ x, lnfreq, groups=procf,
                              pars=pars, panel=function(x, y, ...,
                                pars=pars, ab=coefs) {
                                  lattice::panel.xyplot(x, y, ...)
                                  a1 <- pars[[1]][1]
                                  lambda1 <- pars[[1]][2]
                                  a2 <- pars[[2]][1]
                                  lambda2 <- pars[[2]][2]
                                  if (length(pars) < 3) {
                                      "procFun2" <- function(x) {
                                          log(a1 * lambda1 * exp(-lambda1 * x) +
                                              a2 * lambda2 * exp(-lambda2 * x))
                                      }
                                      lattice::panel.curve(procFun2,
                                                           min(x), max(x),
                                                           add=TRUE)
                                      lattice::panel.abline(ab[[1]], lty=2)
                                      lattice::panel.abline(ab[[2]], lty=3)
                                  } else {
                                      a3 <- pars[[3]][1]
                                      lambda3 <- pars[[3]][2]
                                      "procFun3" <- function(x) {
                                          log(a1 * lambda1 * exp(-lambda1 * x) +
                                              a2 * lambda2 * exp(-lambda2 * x) +
                                              a3 * lambda3 * exp(-lambda3 * x))
                                      }
                                      lattice::panel.curve(procFun3,
                                                           min(x), max(x),
                                                           add=TRUE)
                                      lattice::panel.abline(ab[[1]], lty=2)
                                      lattice::panel.abline(ab[[2]], lty=3)
                                      lattice::panel.abline(ab[[3]], lty=4)
                                  }
                              }, ...)
        print(pp)
    }
    pars
}

"bouts2.nlsFUN" <- function(x, a1, lambda1, a2, lambda2) {
    log(a1 * lambda1 * exp(-lambda1 * x) + a2 * lambda2 * exp(-lambda2 * x))
}

"bouts2.nls" <- function(lnfreq, start, maxiter)
{
    ## Value: list with non linear fitted model and bout ending criterion
    ## --------------------------------------------------------------------
    ## Arguments: lnfreq=data frame with 'lnfreq' (log frequencies) and 'x'
    ## (midpoints), start, maxiter=arguments for nls.
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    fit.nls <- nls(lnfreq ~ bouts2.nlsFUN(x, a1, lambda1, a2, lambda2),
                   data=lnfreq, start=start,
                   control=nls.control(maxiter=maxiter))
    fit.nls
}

"bouts2.nlsBEC" <- function(fit)
{
    ## Value: Numeric with bout ending criterion
    ## --------------------------------------------------------------------
    ## Arguments: list with nls fit
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    if (length(coefs) != 4) {
        stop("fit must have 4 coefficients in a 2-process model")
    }
    a1_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    a2_hat <- as.vector(coefs[3])
    lambda2_hat <- as.vector(coefs[4])
    log((a1_hat * lambda1_hat) / (a2_hat * lambda2_hat)) /
        (lambda1_hat - lambda2_hat)
}

"plotBouts2.nls" <- function(fit, lnfreq, bec.lty=2, ...)
{
    ## Value: plot of fitted model of log frequencies on x, with bec line.
    ## --------------------------------------------------------------------
    ## Arguments: fit=nls list, lnfreq=data frame with named objects lnfreq
    ## and x, bec.lty=line type for arrow; ...=arguments passed to
    ## plot()
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    bec <- bouts2.nlsBEC(fit)
    a1_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    a2_hat <- as.vector(coefs[3])
    lambda2_hat <- as.vector(coefs[4])
    plot(lnfreq ~ x, lnfreq, type="n", ...)
    curve(log(a1_hat * lambda1_hat * exp(-lambda1_hat * x) +
              a2_hat * lambda2_hat * exp(-lambda2_hat * x)),
          min(lnfreq$x), max(lnfreq$x), add=TRUE)
    points(lnfreq ~ x, lnfreq, pch=21, bg="white")
    becy <- predict(fit, list(x=bec))
    usr <- par("usr")
    arrows(bec, becy, bec, usr[3], code=0, lty=bec.lty)
    legend(bec, usr[3] + ((usr[4] - usr[3]) * 0.08),
           paste("bec = ", round(bec, 2), sep=""), bty="n", cex=0.8)
    a1_hat <- round(a1_hat, 2)
    a2_hat <- round(a2_hat, 3)
    lambda1_hat <- round(lambda1_hat, 3)
    lambda2_hat <- round(lambda2_hat, 4)
    legend("topright",
           legend=bquote(y == log(.(a1_hat) %.% .(lambda1_hat) %.%
                             e^(- .(lambda1_hat) * x) +
                             .(a2_hat) %.% .(lambda2_hat) %.%
                             e^(- .(lambda2_hat) * x))),
           bty="n", cex=0.8, adj=c(0, 1))
}

"bouts3.nlsFUN" <- function(x, a1, lambda1, a2, lambda2, a3, lambda3) {
    log(a1 * lambda1 * exp(-lambda1 * x) +
        a2 * lambda2 * exp(-lambda2 * x) +
        a3 * lambda3 * exp(-lambda3 * x))
}

"bouts3.nls" <- function(lnfreq, start, maxiter)
{
    ## Value: list with non linear fitted model and bout ending criterion
    ## --------------------------------------------------------------------
    ## Arguments: lnfreq=data frame with 'lnfreq' (log frequencies) and 'x'
    ## (midpoints), start, maxiter=arguments for nls.
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    fit.nls <- nls(lnfreq ~ bouts3.nlsFUN(x, a1, lambda1, a2, lambda2,
                                          a3, lambda3), data=lnfreq,
                   start=start, control=nls.control(maxiter=maxiter))
    fit.nls
}

"bouts3.nlsBEC" <- function(fit)
{
    ## Value: Numeric with bout ending criterion
    ## --------------------------------------------------------------------
    ## Arguments: list with nls fit
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    if (length(coefs) != 6) {
        stop("fit must have 6 coefficients in a 3-process model")
    }
    a1_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    a2_hat <- as.vector(coefs[3])
    lambda2_hat <- as.vector(coefs[4])
    a3_hat <- as.vector(coefs[5])
    lambda3_hat <- as.vector(coefs[6])
    b1 <- log((a1_hat * lambda1_hat) / (a2_hat * lambda2_hat)) /
        (lambda1_hat - lambda2_hat)
    b2 <- log((a2_hat * lambda2_hat) / (a3_hat * lambda3_hat)) /
        (lambda2_hat - lambda3_hat)
    c(bec1=b1, bec2=b2)
}

"plotBouts3.nls" <- function(fit, lnfreq, bec.lty=2, ...)
{
    ## Value: plot of fitted model of log frequencies on x, with bec lines.
    ## --------------------------------------------------------------------
    ## Arguments: fit=nls list, lnfreq=data frame with named objects lnfreq
    ## and x, bec.lty=line type for arrow; ...=arguments passed to
    ## plot()
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    bec <- bouts3.nlsBEC(fit)
    a1_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    a2_hat <- as.vector(coefs[3])
    lambda2_hat <- as.vector(coefs[4])
    a3_hat <- as.vector(coefs[5])
    lambda3_hat <- as.vector(coefs[6])
    plot(lnfreq ~ x, lnfreq, type="n", ...)
    curve(log(a1_hat * lambda1_hat * exp(-lambda1_hat * x) +
              a2_hat * lambda2_hat * exp(-lambda2_hat * x) +
              a3_hat * lambda3_hat * exp(-lambda3_hat * x)),
          min(lnfreq$x), max(lnfreq$x), add=TRUE)
    points(lnfreq ~ x, lnfreq, pch=21, bg="white")
    becy <- predict(fit, list(x=bec))
    usr <- par("usr")
    arrows(bec, becy, bec, usr[3], code=0, lty=bec.lty)
    legend(bec[1], usr[3] + ((usr[4] - usr[3]) * 0.12),
           paste("bec1 = ", round(bec[1], 2), sep=""), bty="n", cex=0.8)
    legend(bec[2], usr[3] + ((usr[4] - usr[3]) * 0.08),
           paste("bec2 = ", round(bec[2], 2), sep=""), bty="n", cex=0.8)
    a1_hat <- round(a1_hat, 2)
    a2_hat <- round(a2_hat, 3)
    a3_hat <- round(a3_hat, 3)
    lambda1_hat <- round(lambda1_hat, 3)
    lambda2_hat <- round(lambda2_hat, 4)
    lambda3_hat <- round(lambda3_hat, 4)
    legend("topright",
           legend=bquote(y == log(.(a1_hat) %.% .(lambda1_hat) %.%
                           e^(- .(lambda1_hat) * x) +
                           .(a2_hat) %.% .(lambda2_hat) %.%
                           e^(- .(lambda2_hat) * x) +
                           .(a3_hat) %.% .(lambda3_hat) %.%
                           e^(- .(lambda3_hat) * x))),
           bty="n", cex=0.8, adj=c(0, 1))
}

"labelBouts" <- function(x, bec, bec.method=c("standard", "seq.diff"))
{
    ## Value: a numeric vector labelling each row in x with a unique,
    ## sequential bout number
    ## --------------------------------------------------------------------
    ## Arguments: x=numeric vector or matrix with variable or variables,
    ## respectively, to use for splitting bouts; bec=vector or matrix with
    ## corresponding bout ending criterion (i.e. each element/column of x
    ## is compared against the element in bec at the same index),
    ## bec.method=what method was used to identify bouts
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    if (!is(x, "matrix")) x <- as.matrix(x)
    if (!is(bec, "matrix")) bec <- as.matrix(bec)
    if (!identical(dim(x), dim(bec)))
        stop(paste("x and bec must have the same", "dimensions"))
    bec.method <- match.arg(bec.method)
    switch(bec.method,
           standard = {xx <- x[-1, ]},
           seq.diff = {xx <- apply(x, 2, function(k) abs(diff(k)))})
    testfun <- function(xi, beci) ifelse(xi > beci, 1, 2)
    bectest <- mapply(testfun, xx, bec[-1, ])
    dim(bectest) <- dim(xx)
    bectest.full <- rbind(1, bectest)
    bectest.any <- apply(bectest.full, 1, function(k) any(k < 2))
    chgbout <- which(bectest.any)
    boutno <- seq(along=chgbout)
    reps <- diff(c(chgbout, nrow(x) + 1))
    rep(boutno, reps)
}

"bouts2.mleFUN" <- function(x, p, lambda1, lambda2)
{
    log(p * lambda1 * exp(-lambda1 * x) +
        (1 - p) * lambda2 * exp(-lambda2 * x))
}

"bouts2.ll" <- function(x) # 2-process Poisson
{
    function(p, lambda1, lambda2) {
        -sum(diveMove::bouts2.mleFUN(x, p, lambda1, lambda2))
    }
}

"bouts2.LL" <- function(x) # 2-process Poisson; transformed model
{
    function(p, lambda1, lambda2) {
        p <- unLogit(p)
        lambda1 <- exp(lambda1)
        lambda2 <- exp(lambda2)
        -sum(diveMove::bouts2.mleFUN(x, p, lambda1, lambda2))
    }
}

## "bouts3.mleFUN" <- function(x, p1, lambda1, p2, lambda2, lambda3)
## {
##     log(p1 * lambda1 * exp(-lambda1 * x) +
##         p2 * lambda2 * exp(-lambda2 * x) +
##         (1 - (p1 + p2)) * exp(-lambda3 * x))
## }

## "bouts3.ll" <- function(x) # 3-process Poisson
## {
##     function(p1, lambda1, p2, lambda2, lambda3) {
##         -sum(diveMove::bouts3.mleFUN(x, p1, lambda1,
##                                      p2, lambda2, lambda3))
##     }
## }

## "bouts3.LL" <- function(x) # 3-process Poisson; transformed model
## {
##     function(p1, lambda1, p2, lambda2, lambda3) {
##         p1 <- unLogit(p1); p2 <- unLogit(p2)
##         lambda1 <- exp(lambda1)
##         lambda2 <- exp(lambda2)
##         lambda3 <- exp(lambda3)
##         -sum(diveMove::bouts3.mleFUN(x, p1, lambda1,
##                                      p2, lambda2, lambda3))
##     }
## }

"bouts.mle" <- function(ll.fun, start, x, ...)
{
    ## Value: An mle object with fitted parameters
    ## --------------------------------------------------------------------
    ## Arguments: loglik.fun=string naming the function to fit; start=named
    ## list with starting values (exactly as given in ll.fun2, i.e. the
    ## reparameterized versions); x=numeric vector with variable to model;
    ## ...=passed to mle
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    loglik.fun <- ll.fun(x)
    fit.mle <- stats4::mle(loglik.fun, start=start, ...)
    fit.mle
}

"bouts2.mleBEC" <- function(fit)
{
    ## Value: Numeric with bout ending criterion
    ## --------------------------------------------------------------------
    ## Arguments: fit=mle object with fitted 2-process model
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    if (length(coefs) != 3) {
        stop("fit must have 3 coefficients in a 2-process model")
    }
    p_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    lambda2_hat <- as.vector(coefs[3])
    log((p_hat * lambda1_hat) / ((1 - p_hat) * lambda2_hat)) /
        (lambda1_hat - lambda2_hat)
}

"plotBouts2.mle" <- function(fit, x, xlab="x", ylab="Log Frequency",
                             bec.lty=2, ...)
{
    ## Value: plot
    ## --------------------------------------------------------------------
    ## Arguments: fit=mle object with fitted 2-process model, x=numeric
    ## vector with observed data; xlab=ylab=strings for titles;
    ## bec.lty=line type for bec; ...=args to curve().
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    coefs <- coef(fit)
    p_hat <- as.vector(coefs[1])
    lambda1_hat <- as.vector(coefs[2])
    lambda2_hat <- as.vector(coefs[3])
    bec <- bouts2.mleBEC(fit)
    range.x <- range(x, na.rm=TRUE)
    curve(log(p_hat * lambda1_hat * exp(-lambda1_hat * x) +
              (1 - p_hat) * lambda2_hat * exp(-lambda2_hat * x)),
          from=range.x[1], to=range.x[2], xlab=xlab, ylab=ylab,
          xaxs="i", yaxs="i", ...)
    rug(jitter(x), side=3, ticksize=0.015, quiet=TRUE)
    becy <- bouts2.mleFUN(bec, p=p_hat, lambda1=lambda1_hat,
                          lambda2=lambda2_hat)
    usr <- par("usr")
    arrows(bec, becy, bec, usr[3], code=0, lty=bec.lty)
    legend(bec, usr[3] + ((usr[4] - usr[3]) * 0.08),
           paste("bec = ", round(bec, 2), sep=""), bty="n", cex=0.8)
    p_hat <- round(p_hat, 2)
    lambda1_hat <- round(lambda1_hat, 3)
    lambda2_hat <- round(lambda2_hat, 4)
    legend("topright",
           legend=bquote(y == log(.(p_hat) %.% .(lambda1_hat) %.%
                             e^(- .(lambda1_hat) * x) +
                             .(1 - p_hat) %.% .(lambda2_hat) %.%
                             e^(- .(lambda2_hat) * x))),
           bty="n", cex=0.8, adj=c(0, 1))
}

"plotBouts2.cdf" <- function(fit, x, draw.bec=FALSE, bec.lty=2, ...)
{
    ## Value: plot
    ## --------------------------------------------------------------------
    ## Arguments: fit=mle object with fitted 2-process model, x=numeric
    ## vector with observed data, draw.bec=logical; whether to draw the bec;
    ## bec.lty=line type for the bec reference line; ...=passed to plot()
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    cdf.fun <- function(x, p, lambda1, lambda2) {
        1 - p*exp(-lambda1*x) - (1 - p)*exp(-lambda2*x)
    }
    coefs <- coef(fit)
    if (!length(coefs) %in% c(3, 4))
        stop("Number of coefficients in 'fit' must be 3 or 4")
    if (length(coefs) == 4) {
        p <- coefs[1] / (coefs[1] + coefs[3])
        coefs <- c(p, coefs[2], coefs[4])
    }
    x <- log1p(x)
    x.ecdf <- ecdf(x)
    plot(x.ecdf, las=1, cex.p=0.5, pch=19, xaxt="n", ...)
    xorig.pretty <- axTicks(1, axp=c(range(exp(x)), 1), log=TRUE)
    xat <- c(0, log1p(xorig.pretty[-1]))
    axis(1, at=xat, labels=c(0, xorig.pretty[-1]))
    plot(function(x) {
        x <- expm1(x)
        cdf.fun(x, coefs[1], coefs[2], coefs[3])
    }, 0, max(x), add=TRUE)
    if (draw.bec) {
        bec <- bec2(fit)
        becy <- cdf.fun(bec, coefs[1], coefs[2], coefs[3])
        arrows(log1p(bec), becy, log1p(bec), 0, code=0, lty=bec.lty)
        legend(log1p(bec), 0.1, paste("bec = ", round(bec, 2), sep=""),
               bty="n", cex=0.8)
    }
}


## We set these here to avoid collating problems
setMethod("bec2", signature(fit="nls"), bouts2.nlsBEC)
setMethod("bec2", signature(fit="mle"), bouts2.mleBEC)
setMethod("bec3", signature(fit="nls"), bouts3.nlsBEC)

## Declare global variables, if needed
if (getRversion() >= "2.15.1") utils::globalVariables("x")



## TEST ZONE --------------------------------------------------------------

## "bouts3.mleFUN" <- function(x, p1, lambda1, p2, lambda2, lambda3)
## {
##     log(p1 * lambda1 * exp(-lambda1 * x) +
##         p2 * lambda2 * exp(-lambda2 * x) +
##         (1 - (p1 + p2)) * exp(-lambda3 * x))
## }

## "bouts3.ll" <- function(x) # 3-process Poisson
## {
##     function(p1, lambda1, p2, lambda2, lambda3) {
##         -sum(bouts3.mleFUN(x, p1, lambda1,
##                                      p2, lambda2, lambda3))
##     }
## }

## "bouts3.LL" <- function(x) # 3-process Poisson; transformed model
## {
##     function(p1, lambda1, p2, lambda2, lambda3) {
##         p1 <- unLogit(p1); p2 <- unLogit(p2)
##         lambda1 <- exp(lambda1)
##         lambda2 <- exp(lambda2)
##         lambda3 <- exp(lambda3)
##         -sum(bouts3.mleFUN(x, p1, lambda1,
##                                      p2, lambda2, lambda3))
##     }
## }

## ## Example code
## utils::example("diveStats", package="diveMove",
##                ask=FALSE, echo=FALSE)
## postdives <- tdrX.tab$postdive.dur
## postdives.diff <- abs(diff(postdives))

## ## Remove isolated dives
## postdives.diff <- postdives.diff[postdives.diff < 4000]
## lnfreq <- boutfreqs(postdives.diff, bw=0.1, plot=FALSE)
## startval <- boutinit(lnfreq, c(50, 400))
## p1 <- startval[[1]]["a"] / (startval[[1]]["a"] + startval[[2]]["a"] +
##                             startval[[3]]["a"])
## p2 <- startval[[2]]["a"] / (startval[[1]]["a"] + startval[[2]]["a"] +
##                             startval[[3]]["a"])

## ## Fit the reparameterized (transformed parameters) model
## ## Drop names by wrapping around as.vector()
## init.parms <- list(p1=as.vector(logit(p1)),
##                    lambda1=as.vector(log(startval[[1]]["lambda"])),
##                    p2=as.vector(logit(p2)),
##                    lambda2=as.vector(log(startval[[2]]["lambda"])),
##                    lambda3=as.vector(log(startval[[3]]["lambda"])))
## bout.fit1 <- bouts.mle(bouts3.LL, start=init.parms, x=postdives.diff,
##                        lower=c(0, -5, -20, -10, -10),
##                        upper=c(20, 0, 0, 0, 0),
##                        method="L-BFGS-B", control=list(trace=TRUE))
## coefs <- as.vector(coef(bout.fit1))

## ## Un-transform and fit the original parameterization
## init.parms <- list(p1=unLogit(coefs[1]), lambda1=exp(coefs[2]),
##                    p2=unLogit(coefs[3]), lambda2=exp(coefs[4]),
##                    lambda3=exp(coefs[5]))
## bout.fit2 <- bouts.mle(bouts3.ll, x=postdives.diff, start=init.parms,
##                        method="L-BFGS-B", lower=rep(1e-08, 5),
##                        control=list(parscale=c(1, 0.1, 1, 0.01, 0.001)))
## coefs <- as.vector(coef(bout.fit2))
## curve(log(coefs[1] * coefs[2] * exp(-coefs[2] * x) +
##           coefs[3] * coefs[4] * exp(-coefs[4] * x) +
##           (1 - (coefs[1] + coefs[3])) * exp(-coefs[5] * x)), 0, 10000)

## ## Simulations

## ## We choose these
## set.seed(10)
## n.sim <- 10000
## p1.sim <- 0.7
## lambda1.sim <- 0.05
## p2.sim <- 0.2
## lambda2.sim <- 0.005
## lambda3.sim <- 0.001
## pars.sim <- c(n.sim, p1.sim, lambda1.sim,
##               p2.sim, lambda2.sim, lambda3.sim)

## ## Try with boot for many simulations
## library(boot)
## nls.fun <- function(pars) {                  # pars=c(n, p, lambda1, lambda2)
##     x <- ifelse(runif(pars[1]) < pars[2], rexp(pars[1], pars[3]),
##                 rexp(pars[1], pars[4]))
##     x.hist <- boutfreqs(x, bw=5, plot=FALSE)
##     startval <- boutinit(x.hist, 80, plot=FALSE)
##     fit <- bouts2.nls(x.hist, start=startval, maxiter=500)
##     coefs.nls <- coef(fit)
##     c(coefs.nls[1]/(coefs.nls[1] + coefs.nls[3]),
##       coefs.nls[2], coefs.nls[4])
## }
## boot.nls <- boot(pars.sim, nls.fun, sim="parametric", R=1000)
## ## Bias with respect to known pars
## means.nls <- apply(boot.nls$t, 2, mean)
## se.nls <- apply(boot.nls$t, 2, function(x) sd(x)/sqrt(length(x)))
## bias.nls <- means.nls - pars.sim[-1]

## ## Fit with MLM
## mle.fun <- function(pars) {
##     x <- ifelse(runif(pars[1]) < pars[2], rexp(pars[1], pars[3]),
##                 rexp(pars[1], pars[4]))
##     x.hist <- boutfreqs(x, bw=5, plot=FALSE)
##     startval <- boutinit(x.hist, 80, plot=FALSE)
##     p <- startval$a1 / (startval$a1 + startval$a2)
##     init.parms <- list(p=logit(p), lambda1=log(startval$lambda1),
##                        lambda2=log(startval$lambda2))
##     fit1 <- bouts.mle(bouts2.LL, start=init.parms, x=x,
##                       method="L-BFGS-B", lower=c(-2, -5, -10))
##     coefs <- as.vector(coef(fit1))
##     init.parms <- list(p=unLogit(coefs[1]), lambda1=exp(coefs[2]),
##                        lambda2=exp(coefs[3]))
##     fit2 <- bouts.mle(bouts2.ll, x=x, start=init.parms,
##                       method="L-BFGS-B", lower=rep(1e-08, 3),
##                       control=list(parscale=c(1, 0.1, 0.01)))
##     coef(fit2)
## }
## boot.mle <- boot(pars.sim, mle.fun, sim="parametric", R=1000)
## ## Bias with respect to known pars
## means.mle <- apply(boot.mle$t, 2, mean)
## se.mle <- apply(boot.mle$t, 2, function(x) sd(x)/sqrt(length(x)))
## bias.mle <- means.mle - pars.sim[-1]

## ## Box plots of simulations
## boot.pars <- data.frame(method=factor(rep(c("MLM", "SDA"),
##                             each=c(nrow(boot.mle$t), nrow(boot.nls$t)))),
##                         p=c(boot.mle$t[, 1], boot.nls$t[, 1]),
##                         lambda1=c(boot.mle$t[, 2], boot.nls$t[, 2]),
##                         lambda2=c(boot.mle$t[, 3], boot.nls$t[, 3]))
## bplot.pars <- bwplot(p + lambda1 + lambda2 ~ method, data=boot.pars,
##                      layout=c(3, 1), as.table=TRUE, allow.multiple=TRUE,
##                      outer=TRUE, do.out=FALSE,
##                      scales=list(y="free", tck=c(0.8, 0), alternating=1,
##                          rot=c(0, 90), x=list(labels=c("MLM", "SDA"))),
##                      strip=strip.custom(bg="transparent",
##                          factor.levels=c(expression(italic(p)),
##                              expression(italic(lambda[f])),
##                              expression(italic(lambda[s])))),
##                      ylab="Estimated value", pch="|",
##                      panel=function(x, y, ...) {
##                          panel.bwplot(x, y, ...)
##                          if (panel.number() == 1)
##                              panel.abline(h=0.7, lty=2, ...)
##                          if (panel.number() == 2)
##                              panel.abline(h=0.05, lty=2, ...)
##                          if (panel.number() == 3)
##                              panel.abline(h=0.005, lty=2, ...)
##                      })

## trellis.device(pdf, file="par-bias.pdf", color=FALSE,
##                width=5, height=5)
## trellis.par.set(box.umbrella=list(lty=1))
## print(bplot.pars)
## dev.off()
