EL.statistic <- function(method, X, Y, d, t, bw = bw.nrd0)
  {
    EL.local(ELLR, X, Y, d, t, method=method, bw = bw)
  }


EL.plot <- function(method, X, Y, bw = bw.nrd0,
                    conf.level=NULL, simultaneous=FALSE,
                    bootstrap.samples=300, more.warnings=FALSE, ...)
  {
    numpoints <- 50
    call <- match.call(expand.dots = FALSE)
    optional <- call$...
    switch(match(method, c("fdiff", "qdiff", "qq", "pp", "roc")),
           {
             dxlab <- expression(x)
             dylab <- substitute(expression(F[a]-F[b]), list(a = call$Y, b = call$X))
           },
           {
             dxlab <- "P"
             dylab <- substitute(expression(F[a]^-1-F[b]^-1), list(a = call$Y, b = call$X))
           },
           {
             dxlab <- substitute(expression(F[a]^-1), list(a = call$Y))
             dylab <- substitute(expression(F[a]^-1), list(a = call$X))
           },
           {
             dxlab <- substitute(expression(F[a]), list(a = call$Y))
             dylab <- substitute(expression(F[a]), list(a = call$X))
           },
           {
             dxlab <- "False positive rate"
             dylab <- "True positive rate"
           }
           )
    tr <- EL.local(function(globals, X, Y) globals$trange(globals, X, Y), X, Y, method=method, bw=bw)
    if (simultaneous)
      {
        trlen <- tr[2] - tr[1]
        tr[1] <- tr[1] + 0.05 * trlen
        tr[2] <- tr[2] - 0.05 * trlen
      }
    t <- seq(tr[1], tr[2], length=numpoints)
    zz <- EL.smooth(method=method, X, Y, t, bw=bw, conf.level = conf.level,
                    simultaneous = simultaneous, bootstrap.samples=bootstrap.samples, more.warnings=more.warnings)
    if (!is.null(conf.level))
      dylim <- c(min(zz$conf.int[1,]), max(zz$conf.int[2,]))
    else
      dylim <- c(min(zz$estim), max(zz$estim))

    if (!("xlab" %in% names(optional)))
      optional <- c(optional, list(xlab=dxlab))
    if (!("ylab" %in% names(optional)))
      optional <- c(optional, list(ylab=dylab))
    if (!("ylim" %in% names(optional)))
      optional <- c(optional, list(ylim=dylim))
    eval(as.call(c(as.list(quote(plot(t, zz$estim, type='l'))), optional)))
    if (!is.null(conf.level))
      {
        eval(as.call(c(as.list(quote(lines(t, zz$conf.int[1,], lty="dashed"))), optional)))
        eval(as.call(c(as.list(quote(lines(t, zz$conf.int[2,], lty="dashed"))), optional)))
      }
  }

EL.smooth <- function(method, X, Y, t, bw = bw.nrd0,
                     conf.level=NULL, simultaneous=FALSE, bootstrap.samples=300,
                      more.warnings = FALSE)
  {
    warnings <- FALSE
    mm <- function(globals, tt)
      {
        d <- sapply(tt, function(t) deltasolve(globals, X, Y, t))
        bs <- NA
        conffun <- function(conf.level)
          {
            function(k) withCallingHandlers(ELconf(globals, X, Y, tt[k], conf.level, delta=d[k]),
                                            simpleWarning = function(e)
                                            {
                                              if (!more.warnings)
                                                {
                                                  warnings <<- TRUE
                                                  invokeRestart("muffleWarning")
                                                }
                                            })
          }        
        if (is.null(conf.level))
          ci <- NA
        else
          {
            if (!simultaneous)
              {
                ci <- sapply(1:length(tt), conffun(conf.level))
                attr(ci, "conf.level") <- conf.level
              }
            else
              {
                bootup <- function()
                  {
                    X1 <- sample(X, replace=TRUE)
                    Y1 <- sample(Y, replace=TRUE)
                    max(mapply(ELLR, delta=d, t=tt, MoreArgs=list(globals=globals, X=X1, Y=Y1)))
                  }
                bs <- sort(replicate(bootstrap.samples, bootup()))
                bs <- bs[ceiling(conf.level * bootstrap.samples)]
                ci <- sapply(1:length(tt), conffun(pchisq(bs, 1)))
              }
          }
        
        list(estimate=d, conf.int=ci, simultaneous.conf.int=simultaneous, bootstrap.crit = bs)
      }
    res <- EL.local(mm, t, method=method, bw = bw)
    if (warnings)
      warning("Estimates or confidence bands for some values could not be found. For a detailed report, run this command again with 'more.warnings = TRUE'.", call. = FALSE)
    res
  }

EL.means <- function(X, Y, mu = 0, conf.level=0.95)
  {
    call <- match.call()
    mm <- function(globals)
      {
        d <- deltasolve(globals, X, Y, 0)
        names(d) <- "Mean difference"
        ci <- ELconf(globals, X, Y, 0, conf.level, delta=d)
        attr(ci, "conf.level") <- conf.level
        ll <- ELLR(globals, X, Y, mu, 0)
        names(ll) <- "-2 * LogLikelihood"
        res <- list(estimate=d, conf.int=ci, p.value = 1 - pchisq(ll, 1), statistic = ll,
                    method="Empirical likelihood mean difference test",
                    null.value=mu, data.name = paste(deparse(call$X), " and ", deparse(call$Y)))
        class(res) <- "htest"
        res
      }
    EL.local(mm, method="mean")
  }

EL.Huber <- function(X, Y, mu=0, conf.level=0.95, scaleX=1, scaleY=1, VX=2.046, VY=2.046, k=1.35)
  {
    call <- match.call()
    mm <- function(globals)
      {
        globals$Hsigma1 <- scaleX
        globals$Hsigma2 <- scaleY
        globals$Hconst <- k
        globals$HV1 <- VX
        globals$HV2 <- VY
        d <- deltasolve(globals, X, Y, 0)
        names(d) <- "Huber estimator difference"
        ci <- ELconf(globals, X, Y, 0, conf.level, delta=d)
        ll <- ELLR(globals, X, Y, mu, 0)
        names(ll) <- "-2 * LogLikelihood"
        attr(ci, "conf.level") <- conf.level
        res <- list(estimate=d, conf.int=ci, p.value = 1 - pchisq(ll, 1), statistic = ll,
                    method="Empirical likelihood Huber estimator difference test",
                    null.value=mu, data.name = paste(deparse(call$X), " and ", deparse(call$Y)))
        class(res) <- "htest"
        res
      }
    EL.local(mm, method="huber")
  }

#### Globals setup

set.globals <- function(globals)
  {
    globals$lambda1 <- 0
    globals$lambda2 <- 0
    globals$wval1 <- NULL
    globals$wval2 <- NULL
    globals$lambdarange1 <- 0
    globals$lambdarange2 <- 0
    globals$degenerate <- FALSE
    globals$precision <- 1e-4
    globals
  }

#### Bandwidth setup

set.bw <- function(globals, bw)
  {
    globals$bw.X <- 0
    globals$bw.Y <- 0
    globals$bw <- bw
    globals$initbw <- bw.init
    globals
  }

bw.init <- function(globals, X, Y)
  {
    if ((globals$bw.X != 0) && (globals$bw.Y != 0))
      return(globals)
    
    if (is.function(globals$bw))
      {
        globals$bw.X <- globals$bw(X)
        globals$bw.Y <- globals$bw(Y)
      }
    else
      {
        globals$bw.X <- globals$bw[1]
        globals$bw.Y <- globals$bw[2]
      }
    globals
  }

#### Kernel setup

set.kernel <- function(globals)
  {
    globals$krnl <- kernel.norm
    globals$dkrnl <- kernel.dnorm
    globals$ikrnl <- kernel.inorm
    globals
  }

kernel.norm <- function(t)
	pnorm(t)

kernel.dnorm <- function(t)
	dnorm(t)

kernel.inorm <- function(t)
	qnorm(t)

#### Method setup

set.method <- function(globals, type)
  {
    if (!(type %in% c("pp", "qq", "roc", "mean", "fdiff", "qdiff", "huber")))
      stop(paste("Unknown EL method: ", type, "\n",
                 "Recognized EL methods are: pp, qq, roc, mean, fdiff, qdiff."))
    

    globals$w1 <- get(paste("w1.", type, sep=""))
    globals$w2 <- get(paste("w2.", type, sep=""))
    globals$alpha1 <- get(paste("alpha1.", type, sep=""))
    globals$alpha2 <- get(paste("alpha2.", type, sep=""))
    globals$trange <- get(paste("trange.", type, sep=""))
    globals$deltarange <- get(paste("deltarange.", type, sep=""))
    globals$thetarange <- get(paste("thetarange.", type, sep=""))
    globals$deltarange.opt <- get(paste("deltarange.opt.", type, sep=""))
    globals$thetarange.opt <- get(paste("thetarange.opt.", type, sep=""))
    globals$use.smoothing <- get(paste("use.smoothing.", type, sep=""))
    globals
  }

### Mean
w1.mean <- function(globals, X, theta, delta, t)
  X - theta

w2.mean <- function(globals, X, theta, delta, t)
  X - theta + delta

alpha1.mean <- function(globals, X, theta, delta, t)
  -1

alpha2.mean <- function(globals, X, theta, delta, t)
  -1

trange.mean <- function(globals, X, Y)
  stop("Mean differences are not dependent on the parameter t.")

deltarange.mean <- function(globals, X, Y, t)
  {
    tr <- c(min(X), max(X))
    tol <- 0.01 * (tr[2] - tr[1])
    c(tr[1] -max(Y) + tol, tr[2]-min(Y)- tol)
  }

deltarange.opt.mean <- deltarange.mean

thetarange.mean <- function(globals, X, Y, delta, t)
  {
    tr <- c(min(X), max(X))
    dr <- c(min(Y) + delta, max(Y) + delta)
    tol <- 0.01 * (tr[2] - tr[1])
    c(max(tr[1], dr[1]) + tol, min(tr[2], dr[2]) - tol)
  }

thetarange.opt.mean <- thetarange.mean

use.smoothing.mean <- FALSE

### Huber estimator

w1.huber <- function(globals, X, theta, delta, t)
  {
    globals$Hsigma <- globals$Hsigma1
    globals$HV <- globals$HV1
    Hsmooth(globals, (X - theta) / globals$Hsigma)
  }

w2.huber <- function(globals, X, theta, delta, t)
  {
    globals$Hsigma <- globals$Hsigma2
    globals$HV <- globals$HV2
    Hsmooth(globals, (X - theta + delta) / globals$Hsigma)
  }

alpha1.huber <- function(globals, X, theta, delta, t)
  {
    globals$Hsigma <- globals$Hsigma1
    globals$HV <- globals$HV1
    - DHsmooth(globals, (X - theta) / globals$Hsigma) / globals$Hsigma
  }

alpha2.huber <- function(globals, X, theta, delta, t)
  {
    globals$Hsigma <- globals$Hsigma2
    globals$HV <- globals$HV2
    - DHsmooth(globals, (X - theta + delta) / globals$Hsigma) / globals$Hsigma
  }

trange.huber <- function(globals, X, Y)
  stop("Huber estimator differences are not dependent on the parameter t.")

deltarange.huber <- function(globals, X, Y, t)
  {
    tr <- c(min(X), max(X))
    c(tr[1] -max(Y) + globals$precision, tr[2]-min(Y)- globals$precision)
  }

deltarange.huber <- deltarange.mean

deltarange.opt.huber <- deltarange.huber

thetarange.huber <- function(globals, X, Y, delta, t)
  {
    tr <- c(min(X) + globals$precision, max(X) - globals$precision)
    dr <- c(min(Y) + delta, max(Y) + delta)
    c(max(tr[1], dr[1]-globals$precision), min(tr[2], dr[2]+globals$precision))
  }

thetarange.huber <- thetarange.mean

thetarange.opt.huber <- thetarange.huber

use.smoothing.huber <- FALSE

Hsmooth <- function(globals, x)
  {
    x.l <- (x - globals$Hconst) / globals$HV
    x.r <- (x + globals$Hconst) / globals$HV
    pnorm(x.l) * (globals$Hconst - x) + pnorm(x.r) * (globals$Hconst + x) -
      globals$Hconst + globals$HV * (dnorm(x.r) - dnorm(x.l))
  }

DHsmooth <- function(globals, x)
  {
    x.l <- (x - globals$Hconst) / globals$HV
    x.r <- (x + globals$Hconst) / globals$HV
    pnorm(x.r) - pnorm(x.l)
  }

### ROC

w1.roc <- function(globals, X, theta, delta, t)
  globals$krnl((theta - X) / globals$bw.X) - (1 - delta)

w2.roc <- function(globals, X, theta, delta, t)
  globals$krnl((theta - X) / globals$bw.Y) - (1 - t)

alpha1.roc <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X

alpha2.roc <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta - X) / globals$bw.Y) / globals$bw.Y

trange.roc <- function(globals, X, Y)
  {
    c(0.01, 0.99)
  }

deltarange.roc <- function(globals, X, Y, t)
  {
    dY <- globals$bw.Y * globals$ikrnl(t)
    diffs <- c(max(X) - min(Y) + dY, min(X) - max(Y) + dY)
    range <- globals$krnl(diffs / globals$bw.X)
    c(min(range)+globals$precision, max(range)-globals$precision)
  }

deltarange.opt.roc <- function(globals, X, Y, t)
  {
    kk <- 1 - ecdf(X)(c(quantile(Y, 1-t) - 3 * globals$bw.X, quantile(Y, 1-t) + 3 * globals$bw.X))
    c(max(0.001, kk[2] - 0.001), min(0.999, kk[1] + 0.001))
  }

thetarange.roc <- function(globals, X, Y, delta, t)
  {
    dX <- globals$bw.X * globals$ikrnl(delta)
    dY <- globals$bw.Y * globals$ikrnl(t)
    dd <- c(max(min(X) - dX, min(Y) - dY), min(max(X) - dX, max(Y) - dY))
    c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))
  }

thetarange.opt.roc <- function(globals, X, Y, delta, t)
  {
    extreme <- thetarange.roc(globals, X, Y, delta, t)
    
    tol <- 2 * globals$bw.Y
    qY <- quantile(Y, 1-t)
    c(max(extreme[1], qY - tol), min(extreme[2], qY + tol))
  }

use.smoothing.roc <- TRUE


### P-P plots


w1.pp <- function(globals, X, theta, delta, t)
  globals$krnl((theta - X) / globals$bw.X) - delta

w2.pp <- function(globals, X, theta, delta, t)
  globals$krnl((theta - X) / globals$bw.Y) - t

alpha1.pp <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X

alpha2.pp <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta - X) / globals$bw.Y) / globals$bw.Y

trange.pp <- function(globals, X, Y)
  {
    c(0.01, 0.99)
  }

deltarange.pp <- function(globals, X, Y, t)
  {
    c(0.001, 0.999)
  }

deltarange.opt.pp <- function(globals, X, Y, t)
  {
    kk <- ecdf(X)(c(quantile(Y, t) - 3 * globals$bw.X, quantile(Y, t) + 3 * globals$bw.X))
    c(max(0.001, kk[1] - 0.001), min(0.999, kk[2] + 0.001))
  }

thetarange.pp <- function(globals, X, Y, delta, t)
  {
    dX <- globals$bw.X * globals$ikrnl(delta)
    dY <- globals$bw.Y * globals$ikrnl(t)
    dd <- c(max(min(X) + dX, min(Y) + dY), min(max(X) + dX, max(Y) + dY))
    c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))
  }

thetarange.opt.pp <- function(globals, X, Y, delta, t)
  {
    dX <- globals$bw.X * globals$ikrnl(delta)
    dY <- globals$bw.Y * globals$ikrnl(t)
    dd <- c(max(min(X) + dX, min(Y) + dY), min(max(X) + dX, max(Y) + dY))
    extreme <- c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))

    tol <- 2 * globals$bw.Y
    qY <- quantile(Y, t) 
    c(max(extreme[1],  qY - tol), min(extreme[2], qY + tol))
  }

use.smoothing.pp <- TRUE

### Q-Q plots
w1.qq <- function(globals, X, theta, delta, t)
  globals$krnl((delta - X) / globals$bw.X) - theta

w2.qq <- function(globals, X, theta, delta, t)
  globals$krnl((t - X) / globals$bw.Y) - theta

alpha1.qq <- function(globals, X, theta, delta, t)
  -1

alpha2.qq <- function(globals, X, theta, delta, t)
  -1

trange.qq <- function(globals, X, Y)
  c(min(Y), max(Y))

deltarange.qq <- function(globals, X, Y, t)
  {
    c(min(X) + globals$precision, max(X) - globals$precision)
  }

deltarange.opt.qq <- function(globals, X, Y, t)
  {
    extreme <- deltarange.qq(globals, X, Y, t)

    tol <- globals$bw.Y
    tr <- ecdf(Y)(c(t - tol, t + tol))
    upper <- quantile(X, tr[2]) +  globals$bw.X
    lower <- quantile(X, tr[1]) -  globals$bw.X
    c(max(lower, extreme[1]), min(upper, extreme[2]))
  }

thetarange.qq <- function(globals, X, Y, delta, t)
  {
    ww <- w2.qq(globals, Y, 0, 0, t)
    tr <- c(min(ww), max(ww))
    ww <- w1.qq(globals, X, 0, delta, 0)
    dr <- c(min(ww), max(ww))
    c(max(tr[1], dr[1]) + globals$precision, min(tr[2], dr[2]) - globals$precision)
  }

thetarange.opt.qq <- thetarange.qq

use.smoothing.qq <- TRUE


### F difference
w1.fdiff <- function(globals, X, theta, delta, t)
  globals$krnl((t - X) / globals$bw.X) - theta

w2.fdiff <- function(globals, X, theta, delta, t)
  globals$krnl((t - X) / globals$bw.Y) - theta - delta

alpha1.fdiff <- function(globals, X, theta, delta, t)
  -1

alpha2.fdiff <- function(globals, X, theta, delta, t)
  -1

trange.fdiff <- function(globals, X, Y)
  c(max(c(min(X), min(Y))), min(c(max(X), max(Y))))

deltarange.fdiff.prim <- function(globals, X, Y, t, prec)
{
  ww <- w1.fdiff(globals, X, 0, 0, t)
  tr <- c(min(ww), max(ww))
  ww <- w2.fdiff(globals, Y, 0, 0, t)
  c(min(ww) - tr[2] + prec, max(ww) - tr[1] - prec)
}

deltarange.fdiff <- function(globals, X, Y, t)
  deltarange.fdiff.prim(globals, X, Y, t, globals$precision)

deltarange.opt.fdiff <- function(globals, X, Y, t)
  deltarange.fdiff.prim(globals, X, Y, t, 1e-7)

thetarange.fdiff <- function(globals, X, Y, delta, t)
{
  ww <- w1.fdiff(globals, X, 0, 0, t)
  tr <- c(min(ww), max(ww))
  ww <- w2.fdiff(globals, Y, 0, delta, t)
  tr2 <- c(min(ww), max(ww))
  c(max(tr2[1], tr[1]) + 1e-7, min(tr2[2], tr[2]) - 1e-7)
}

thetarange.opt.fdiff <- thetarange.fdiff

use.smoothing.fdiff <- TRUE


### Quantile difference

w1.qdiff <- function(globals, X, theta, delta, t)
  globals$krnl((theta - X) / globals$bw.X) - t

w2.qdiff <- function(globals, X, theta, delta, t)
  globals$krnl((theta + delta - X) / globals$bw.Y) - t

alpha1.qdiff <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X

alpha2.qdiff <- function(globals, X, theta, delta, t)
  globals$dkrnl((theta + delta - X) / globals$bw.Y) / globals$bw.Y

trange.qdiff <- function(globals, X, Y)
  {
    minlen <- min(length(X), length(Y))
    c(1/minlen, 1 - 1/minlen)
  }

deltarange.qdiff <- function(globals, X, Y, t)
  {
    c(min(Y) - max(X) + globals$precision, max(Y) - min(X) - globals$precision)
  }

deltarange.opt.qdiff <- function(globals, X, Y, t)
  {
    extreme <- deltarange.qdiff(globals, X, Y, t)

    tol <- 2 * globals$bw.Y
    qdiff <- quantile(Y, t) - quantile(X, t)
    c(max(extreme[1], qdiff - tol), min(extreme[2], qdiff + tol))
  }

thetarange.qdiff <- function(globals, X, Y, delta, t)
  {
    dX <- globals$bw.X * globals$ikrnl(1-t)
    dY <- globals$bw.Y * globals$ikrnl(1-t)
    c(max(min(X) - dX, min(Y) - delta - dY) + globals$precision, min(max(X) - dX, max(Y) - delta - dY) - globals$precision)
  }

thetarange.opt.qdiff <- function(globals, X, Y, delta, t)
  {
    extreme <- thetarange.qdiff(globals, X, Y, delta, t)
    
    tol <- 2 * globals$bw.X
    qX <- quantile(X, t)
    c(max(extreme[1], qX - tol), min(extreme[2], qX + tol))    
  }

use.smoothing.qdiff <- TRUE


#### Generalized calling

EL.local <- function(fun, ..., bw = bw.nrd0, method)
  {
    globals <- list()
    globals <- set.method(globals, method)
    globals <- set.globals(globals)
    globals <- set.kernel(globals)
    globals <- set.bw(globals, bw)
    fun(globals, ...)
  }

#### Calculations

ELconf <- function(globals, X, Y, t, p.level=0.95, delta=NULL)
  {
    globals <- globals$initbw(globals, X, Y)
    if(is.null(delta))
      delta <- deltasolve(globals, X, Y, t)

    if (is.na(delta))
      return(c(NA, NA))
    
    critval <- qchisq(p.level, 1)
    ELlim <- function(delt)
      {
        ELLR(globals, X, Y, delt, t) - critval
      }

    drange <- globals$deltarange(globals, X, Y, t)
    if (drange[2] < drange[1])
      return(NA)

    if (ELlim(delta) < 0)
      {
        if (ELlim(drange[1]) > 0)
          {
            lo <- uniroot(ELlim, c(drange[1], delta))$root
          }
        else
          {
            warning(paste("Could not find lower confidence bound for t =", t, "."), call. = FALSE)
            lo <- drange[1]
          }
        
        if (ELlim(drange[2]) > 0)
          {
            hi <- uniroot(ELlim, c(delta, drange[2]))$root
          }
        else
          {
            warning(paste("Could not find upper confidence bound for t =", t, "."), call. = FALSE)
            hi <- drange[2]
          }
      }
    else
      {
        warning(paste("Could not find estimate at t =", t, "."), call. = FALSE)
        hi <- delta
        lo <- delta
      }
    c(lo, hi)
  }

deltasolve <- function(globals, X, Y, t)
  {
    globals <- globals$initbw(globals, X, Y)
    dr <- globals$deltarange.opt(globals, X, Y, t)
    if (dr[2] < dr[1])
      return(NA)
    if (dr[2] - dr[1] < 1e-7)
      return(dr[1])
    globals$tr <- globals$thetarange
    globals$thetarange <- globals$thetarange.opt
    oo <- suppressWarnings(optimize(function(d) ELLR(globals, X, Y, d, t), dr))
    globals$thetarange <- globals$tr
    if (is.finite(oo$objective))
      oo$minimum
    else
      NA
  }

ELLR <- function(globals, X, Y, delta, t)
  {
    globals <- globals$initbw(globals, X, Y)
    ts <- thetasolve(globals, X, Y, delta, t)
    theta <- ts$theta
    globals <- ts$globals
    if (!globals$degenerate)
      {
        lr <- 2* (sum(log(1 + globals$lambda1 * globals$wval1)) +
                  sum(log(1 + globals$lambda2 * globals$wval2)))
        if (lr < 1e-7)
          0
        else
          lr
      }
    else
      Inf
  }

thetasolve <- function(globals, X, Y, delta, t)
  {
    globals <- globals$initbw(globals, X, Y)
    len1 <- length(X)
    len2 <- length(Y)
    tgrad <- function(theta)
      {
        globals$degenerate <<- FALSE
        globals$wval1 <<- globals$w1(globals, X, theta, delta, t)
        globals$wval2 <<- globals$w2(globals, Y, theta, delta, t)
        alphaval1 <- globals$alpha1(globals, X, theta, delta, t)
        if (is.na(alphaval1[2]))
          alphaval1 <- rep.int(alphaval1, len1)
        alphaval2 <- globals$alpha2(globals, Y, theta, delta, t)
        if (is.na(alphaval2[2]))
          alphaval2 <- rep.int(alphaval2, len2)
        res <- .C("theta_equation",
                  as.integer(len1),
                  as.double(globals$wval1),
                  as.double(alphaval1),
                  as.integer(len2),
                  as.double(globals$wval2),
                  as.double(alphaval2),
                  l1 = double(1),
                  l2 = double(1),
                  res = double(1))
        globals$lambda1 <<- res$l1
        globals$lambda2 <<- res$l2
        res$res
      }

    tr <- globals$thetarange(globals, X, Y, delta, t)
    if (tr[1] > tr[2])
      {
        globals$degenerate <- TRUE
        return(list(theta = NA, globals=globals))
      }
    theta <- suppressWarnings(try(uniroot(tgrad, tr)$root, silent=TRUE))
    list(theta = theta, globals = globals)
  }
