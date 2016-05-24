mkCPSP <- function(st,
                   from=floor(min(st)),
                   to=ceiling(max(st))
                   ) {

  ## mkCPSPFct returns a CountingProcessSamplePath object
  ## Arguments:
  ##  st: A spike train object or a vector with strictly increasing elements
  ##      containing the spike times (measured in s)
  ##  from: The time (in s) at which the observations started
  ##  to: The time (in s) at which the observations ended
  ## Return
  ##  a list with the following components:
  ##    cpspFct: a function of a single time argument, t, instance of
  ##    the class "CountingProcessSamplePath". 
  ##    The function is cadlag (right continous and admits a limit on the left) between
  ##    from and to. The function returns NA outside of this interval. For a vector 
  ##    argument, a vector is returned. The function contains the original data in its
  ##    closure. The matched call is contained in attribute "call".
  ##    from: argument from
  ##    to: argument to
  ##    call: the matched call

  ## check that st is suitable
  st <- unclass(st)
  if (any(diff(st) <= 0)) stop("Wrong st.")
  
  cpspFct <- function(t) {
    if (missing(t)) t <- st
    result <- numeric(length(t))
    good <- from <= t & t <= to
    result[!good] <- NA
    result[good] <- sapply(t[good], function(x) sum(st <= x))
    result
  }

  ppspFct <- function() st

  spikeTrainFct <- function() as.spikeTrain(st)
  
  result <- list(cpspFct=cpspFct,
                 ppspFct=ppspFct,
                 spikeTrainFct=spikeTrainFct,
                 from=from,
                 to=to,
                 call=match.call()
                 )
  class(result) <- "CountingProcessSamplePath"
  result
  
}

print.CountingProcessSamplePath <- function(x,
                                            digits=5,
                                            ...) {

  cat(paste("A CountingProcessSamplePath object\n"))
  xx <- x$ppspFct()
  isi <- diff(xx)
  isi.m <- min(isi)
  isi.M <- max(isi)
  isi.mu <- mean(isi)
  isi.sd <- sd(isi)
  n <- length(xx)
  from <- x$from
  to <- x$to
  cat(paste("with",n,"events, from:",from,"to:",to,"\n"))
  cat(paste("Smallest inter event interval:",round(isi.m,digits=digits),
            ", largest inter event interval:",round(isi.M,digits=digits),
            "\n"))
  cat(paste("Mean inter event interval:",round(isi.mu,digits=digits),
      ", sd:", round(isi.sd,digits=digits),"\n"))
  
}

plot.CountingProcessSamplePath <- function(x,y,
                                           col,lwd,
                                           xlim,ylim,
                                           xlab,ylab,
                                           xaxs,yaxs,
                                           main,
                                           ...) {

  ## plot.CountingProcessSamplePath plot method for CountingProcessSamplePath objects
  ## Arguments:
  ##  x: A CountingProcessSamplePath object
  ##  y: Not used but required for plot methods.
  ##  xlim, ylim, xlab, ylab, main, lwd, col: Same as in plot
  ##  ...: additional arguments passed to plot

  xx <- x$ppspFct()
  n <- length(xx)
  if (missing(xlim)) xlim <- c(x$from,
                               x$to
                               )
  xx <- c(xx,xlim[2])
  if (missing(ylim)) ylim <- c(sum(xx<xlim[1]),sum(xx<=xlim[2]))
  if (missing(xlab)) xlab <- ""
  if (missing(ylab)) ylab <- "Cumulative number of evts."
  if (missing(main)) main <- "Counting Process Sample Path"
  if (missing(xaxs)) xaxs <- "i"
  if (missing(yaxs)) yaxs <- "i"
  plot(xlim,ylim,type="n",
       xlab=xlab,
       ylab=ylab,
       main=main,
       xaxs=xaxs,
       yaxs=yaxs,
       ...)
  yy <- 1:n
  if (missing(lwd)) lwd <- 1
  if (missing(col)) col <- 1
  segments(xx[-n],yy,xx[-1],yy,lwd=lwd,col=col)
  
}

lines.CountingProcessSamplePath <- function(x,
                                            ...) {

  ## lines.CountingProcessSamplePath lines method for CountingProcessSamplePath objects
  ## Arguments:
  ##  x: A CountingProcessSamplePath object
  ##  ...: additional arguments passed to segments
  
  xx <- evalq(st,envir=environment(x$cpspFct))
  n <- length(xx)
  xlim <- c(x$from,
            x$to
            )
  xx <- c(xx,xlim[2])
  yy <- 1:n
  segments(xx[-n],yy,xx[-1],yy,...)
  
}


as.CPSP <- function(x) {

  if (is.numeric(x)) return(mkCPSP(x))
  if (is.spikeTrain(x)) return(mkCPSP(unclass(x)))
  
}

summary.CountingProcessSamplePath <- function(object,
                                              exact = TRUE,
                                              lag.max = NULL,
                                              d = max(c(2, sqrt(length(object$ppspFct()))%/%5)),
                                              ...) {

  x <- object$ppspFct()
  n <- length(x)
  xn <- (x-object$from)/(x[n]-object$from)
  UniformGivenN <- ks.test(xn[-n],punif,exact=exact)$p.value

  b95 <- function(t) 0.299944595870772 + 2.34797018726827*sqrt(t)
  b99 <- function(t) 0.313071417065285 + 2.88963206734397*sqrt(t)

  ## Modification 2008 10 20 by C Pouzat
  ## y.w <- (seq(x) - x)/sqrt(object$to-object$from)
  y.w <- (x-x[1])[-1]
  ny <- length(y.w)
  y.w <- (y.w-1:ny)/sqrt(ny)
  x.w <- (1:ny)/ny

  Wiener95 <- all(-b95(x.w) < y.w & y.w < b95(x.w))
  Wiener99 <- all(-b99(x.w) < y.w & y.w < b99(x.w))

  xs <- sort(diff(x))
  BermanTest <- ks.test(xs,pexp,exact=exact)$p.value

  isi <- diff(x)
  isi.o <- rank(isi)/length(isi)
  isi.l <- length(isi)
  if (is.null(lag.max)) 
    lag.max <- round(10 * log10(isi.l))
  lag.max <- min(isi.l - 1, lag.max)
  grid <- seq(0, 1, length.out = d + 1)
  getChi2 <- function(lag) {
    isi.x <- isi.o[1:(isi.l - lag.max)]
    isi.y <- isi.o[(1 + lag):(isi.l - lag.max + lag)]
    isi.x <- as.integer(cut(isi.x, breaks = grid))
    isi.y <- as.integer(cut(isi.y, breaks = grid))
    counts <- matrix(0, nrow = d, ncol = d)
    for (i in seq(along.with = isi.x)) counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], 
                    isi.y[i]] + 1
    chisq.test(counts, ...)
  }
  chi2seq <- lapply(1:lag.max, getChi2)
  chi2.df <- chi2seq[[1]]$parameter
  ## minChi2 <- qchisq(0.025, df = chi2.df)
  ## maxChi2 <- qchisq(0.975, df = chi2.df)
  chi2V <- sapply(chi2seq, function(l) l$statistic)
  ## outOf95 <- chi2V < minChi2 | chi2V > maxChi2
  chi2.w <- cumsum((chi2V-chi2.df)/sqrt(2*chi2.df*lag.max))
  chi2.95 <- all(abs(chi2.w) < b95((1:lag.max)/lag.max))
  chi2.99 <- all(abs(chi2.w) < b99((1:lag.max)/lag.max))
                 
  vt <- varianceTime(x,CI=c(0.95,0.99))
  if (!is.null(dim(vt$ciLow))) {
    vt.out95 <- vt$s2 < vt$ciLow[1,] | vt$ciUp[1,] < vt$s2
    vt.out99 <- vt$s2 < vt$ciLow[2,] | vt$ciUp[2,] < vt$s2
  } else {
    vt.out95 <- vt$s2 < vt$ciLow[1] | vt$ciUp[1] < vt$s2
    vt.out99 <- vt$s2 < vt$ciLow[2] | vt$ciUp[2] < vt$s2
    warning("Only one window size for varianceTime")
  }
  
  result <- list(UniformGivenN=UniformGivenN,
                 Wiener95=Wiener95,
                 Wiener99=Wiener99,
                 BermanTest=BermanTest,
                 RenewalTest=list(chi2.95=chi2.95,chi2.99=chi2.99,total=lag.max),
                 varianceTime=vt,
                 varianceTimeSummary=c(total=length(vt$s2),out95=sum(vt.out95),out99=sum(vt.out99)),
                 n=n,
                 call=match.call()
                 )

  class(result) <- "CountingProcessSamplePath.summary"
  result
  
}

print.CountingProcessSamplePath.summary <- function(x,
                                                    digits=5,
                                                    ...) {

  cat(paste(" *** Test of uniformity on the time axis \n",
            "    Prob. of the Kolmogorov statistic under H0:",
            round(x$UniformGivenN,digits=digits),"\n"
            )
      )
  cat(paste(" *** Wiener process test \n",
            "    Inside 95% domain:",x$Wiener95,", inside 99% domain:",x$Wiener99,"\n"
            )
      )
  cat(paste(" *** Berman test  \n",
            "    Prob. of the Kolmogorov statistic under H0:",
            round(x$BermanTest,digits=digits),"\n"
            )
      )
  cat(paste(" *** Renewal test \n",
            "    Inside 95% domain:",
            x$RenewalTest[["chi2.95"]],
            ", inside 99% domain:",
            x$RenewalTest[["chi2.99"]],"\n",
            "    Maximum lag:", x$RenewalTest[["total"]],"\n"
            )
      )
  ## if (x$RenewalTest["out"]>1)
  ##   cat(paste("    ",x$RenewalTest["out"],"lags out on a total of:",x$RenewalTest["total"],"\n"))
  ## else
  ##  cat(paste("    ",x$RenewalTest["out"],"lag out on a total of:",x$RenewalTest["total"],"\n"))

  total <- x$varianceTimeSummary["total"]
  out95 <- x$varianceTimeSummary["out95"]
  out99 <- x$varianceTimeSummary["out99"]
  cat(paste(" *** Variance vs \"time\" with", total,"time windows:\n"))
  if (out95>1) cat(paste("    ",out95,"windows out at 95% level\n"))
  else cat(paste("    ",out95,"window out at 95% level\n"))
  if (out99>1) cat(paste("    ",out99,"windows out at 99% level\n"))
  else cat(paste("    ",out99,"window out at 99% level\n"))
  cat(paste(" *** The object contains", x$n,"events.\n"))

}


plot.CountingProcessSamplePath.summary <- function (x, y,
                                                    which = c(1,2,6,8),
                                                    main,
                                                    caption = c(
                                                      expression(paste("Uniform on ", Lambda," Test")), 
                                                      "Berman's Test",
                                                      "Log Survivor Function",
                                                      expression(paste(U[k+1]," vs ", U[k])), 
                                                      "Variance vs Mean Test",
                                                      "Wiener Process Test",
                                                      "Autocorrelation Fct.",
                                                      "Renewal Test"),
                                                    ask = FALSE,
                                                    lag.max = NULL,
                                                    d = max(c(2, sqrt(length(eval(x$call[[2]])$ppspFct()))%/%5)),
                                                    ...) {
  
    if (!inherits(x, "CountingProcessSamplePath.summary")) 
        stop("use only with \"CountingProcessSamplePath.summary\" objects")

    force(d)
    ## Get the spike times
    cpsp <- eval(x$call[[2]])
    xx <- cpsp$ppspFct()
    isi <- diff(xx)
    Y <- seq(xx)
    nbSpikes <- length(xx)
    b95 <- function(t) 0.299944595870772 + 2.34797018726827*sqrt(t)
    b99 <- function(t) 0.313071417065285 + 2.88963206734397*sqrt(t)
      
    show <- logical(8)
    show[which] <- TRUE
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    else {
        if (sum(show) == 2) 
            layout(matrix(1:2, nrow = 2))
        if (2 < sum(show) && sum(show) < 5) 
            layout(matrix(1:4, nrow = 2))
        if (sum(show) >= 5) 
            layout(matrix(1:9, nrow = 3))
    }
    mainGiven <- !missing(main)
    if (show[1]) {
      slopeKS <- length(xx)/max(xx)
      plot(as.numeric(xx), Y, type = "n", xlab = expression(Lambda), 
           ylab = expression(N(Lambda)),
           main = ifelse(mainGiven, main, caption[1]),
           sub = ifelse(mainGiven,caption[1],"")
           )
        abline(a = 0, b = slopeKS)
        abline(a = 1.36 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = -1.36 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = 1.63 * sqrt(nbSpikes), slopeKS, lty = 2)
        abline(a = -1.63 * sqrt(nbSpikes), slopeKS, lty = 2)
        lines(as.numeric(xx), Y, col = 2, lwd = 2)
    }
    lambda <- 1 - exp(-isi)
    if (show[2]) {
        plot(c(0, 1), c(0, 1), type = "n", xlab = expression(U[(k)]), 
            ylab = "ECDF", main = ifelse(mainGiven, 
                main, caption[2]), sub = ifelse(mainGiven, caption[2], 
                ""))
        abline(a = 0, b = 1)
        abline(a = 1.36/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = -1.36/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = 1.63/sqrt(nbSpikes - 1), 1, lty = 2)
        abline(a = -1.63/sqrt(nbSpikes - 1), 1, lty = 2)
        lines(sort(lambda), (1:(nbSpikes - 1))/(nbSpikes - 1), 
            col = 2, lwd = 2)
    }
    if (show[3]) {
        nI <- length(isi)
        Y <- (nI:1)/nI
        X <- sort(isi)
        Yth <- exp(-X)
        Y95p <- qbinom(0.975, nI, Yth)/nI
        Y95m <- qbinom(0.025, nI, Yth)/nI
        Y99p <- qbinom(0.995, nI, Yth)/nI
        Y99m <- qbinom(0.005, nI, Yth)/nI
        maxId <- max(which(Yth > 0.001))
        plot(c(0, X[maxId]), c(0.001, 1), type = "n", xlab = expression(Y[(k)]), 
            ylab = "Survivor Fct", main = ifelse(mainGiven, main, 
                caption[3]), sub = ifelse(mainGiven, caption[3], 
                ""), log = "y")
        lines(X, Y95p, lty = 2)
        lines(X, Y95m, lty = 2)
        lines(X, Y99p, lty = 2)
        lines(X, Y99m, lty = 2)
        lines(X, Y, col = 2, lwd = 2)
    }
    if (show[4]) {
        plot(lambda[-length(lambda)], lambda[-1], xlab = expression(U[k]), 
            ylab = expression(U[k + 1]), pch = 3, main = ifelse(mainGiven, 
                main, caption[4]), sub = ifelse(mainGiven, caption[4], 
                ""))
    }
    if (show[5]) {
        plot(x$varianceTime, style = "Ogata", xlab = "Window Length", 
            main = ifelse(mainGiven, main, caption[5]), sub = ifelse(mainGiven, 
                caption[5], ""))
    }
    if (show[6]) {
      y.w <- (xx-xx[1])[-1]
      n <- length(y.w)
      y.w <- c(0,(y.w-1:n)/sqrt(n))
      x.w <- (0:n)/n
      tt <- seq(0,1,0.01)
      ylim <- c(min(-b99(tt[length(tt)]),min(y.w)),
                max(b99(tt[length(tt)]),max(y.w))
                )
      plot(x.w, y.w, type = "n", xlab = expression(t),
           ylab = expression(X[t]^n), 
           ylim = ylim,
           main = ifelse(mainGiven, main, caption[6]),
           sub = ifelse(mainGiven, caption[6], ""),
           xaxs="i",yaxs="i"
           )
        abline(h = 0)
        lines(tt, b95(tt), lty = 2)
        lines(tt, -b95(tt), lty = 2)
        lines(tt, b99(tt), lty = 2)
        lines(tt, -b99(tt), lty = 2)
        lines(x.w, y.w, col = 2, lwd = 2)
    }
    if (show[7]) {
      acf(isi,lag.max=lag.max,
          main = ifelse(mainGiven, main, caption[7]),
           sub = ifelse(mainGiven, caption[7], ""),
          ...)
    }
    if (show[8]) {
      isi <- diff(xx)
      isi.o <- rank(isi)/length(isi)
      isi.l <- length(isi)
      if (is.null(lag.max)) 
        lag.max <- round(10 * log10(isi.l))
      lag.max <- min(isi.l - 1, lag.max)
      grid <- seq(0, 1, length.out = d + 1)
      getChi2 <- function(lag) {
        isi.x <- isi.o[1:(isi.l - lag.max)]
        isi.y <- isi.o[(1 + lag):(isi.l - lag.max + lag)]
        isi.x <- as.integer(cut(isi.x, breaks = grid))
        isi.y <- as.integer(cut(isi.y, breaks = grid))
        counts <- matrix(0, nrow = d, ncol = d)
        for (i in seq(along.with = isi.x)) counts[isi.x[i], isi.y[i]] <- counts[isi.x[i], 
                        isi.y[i]] + 1
        chisq.test(counts, ...)
      }
      chi2seq <- lapply(1:lag.max, getChi2)
      chi2.df <- chi2seq[[1]]$parameter
      chi2V <- sapply(chi2seq, function(l) l$statistic)
      y.w <- c(0,cumsum((chi2V-chi2.df)/sqrt(2*chi2.df*lag.max)))
      x.w <- (0:lag.max)/lag.max
      tt <- seq(0,1,0.01)
      ylim <- c(min(-b99(tt[length(tt)]),min(y.w)),
                max(b99(tt[length(tt)]),max(y.w))
                )
      sub <- ifelse(mainGiven,
                    paste(caption[8],", lag.max = ",lag.max,", d =",d,sep=""),
                    paste("lag.max = ",lag.max,", d = ",d,sep="")
                    )
      plot(x.w, y.w, type = "n", xlab = expression(t),
           ylab = expression(X[t]^n), 
           ylim = ylim,
           main = ifelse(mainGiven, main, caption[8]),
           sub = sub,
           xaxs="i",yaxs="i"
           )
        abline(h = 0)
        lines(tt, b95(tt), lty = 2)
        lines(tt, -b95(tt), lty = 2)
        lines(tt, b99(tt), lty = 2)
        lines(tt, -b99(tt), lty = 2)
        lines(x.w, y.w, col = 2, lwd = 2)
      
    }

}
