"chisq.setup" <-
  function(df=1,
           ncp=0,
           log.p=FALSE,
           xlim.in=c(0, qchisq.intermediate(p=1-.01, df=df, ncp=ncp, log.p=log.p)),
           ylim.in=range(c(0, 1.1*dchisq.intermediate(x=seq(max(0.5,df-2),df+2,.01), df=df, ncp=ncp, log=log.p))),
           main.in=main.calc, ylab.in="Chisq density",
           ...) {
    main.calc <- parse(text=
                         paste("Chisq~density: ~~ nu==", df,
                               sep=""))
    if (ncp > 0) paste(main.calc, " ncp=", ncp, sep="")
    if (ncp < 0) stop("non-centrality parameter must be non-negative.")

    plot(xlim=xlim.in, ylim=ylim.in,
         x=xlim.in, y=ylim.in,
         yaxt="n", type="n",
         las=1,
         xlab="",
         ylab=ylab.in,
         main=main.in,
         ...)
    abline(h=0, v=0)
    axis(4, las=1)
  }

"chisq.curve" <-
  function(df=1,
           ncp=0,
           log.p=FALSE,
           alpha=.05,
           critical.values=chisq.alpha,
           chisq=seq(0, par()$usr[2], length=109),
           shade="right", col=par("col"),
           axis.name="chisq",
           ...) {

    ## Valid values for shade are "right", "left", "inside", "outside",
    ## "none".  Default is "right" for one-sided critical.values and
    ## "outside" for two-sided critical values.  "none" is used to
    ## redraw an outline of the curve that would otherwise be obscured
    ## by a solid color from the shaded area of another curve.

    ## NA and missing alpha
    switch(1 + length(alpha),
           alpha <- c(0,0),
           if (is.na(alpha)) alpha <- c(0,0),
           alpha[is.na(alpha)] <- 0)

    chisq.alpha <- if (length(alpha)==1)
      qchisq.intermediate(p=1-alpha, df=df, ncp=0, log.p=log.p)
    else
      qchisq.intermediate(p=c(alpha[1], 1-alpha[2]), df=df, ncp=ncp, log.p=log.p)

    lines(y=dchisq.intermediate(x=chisq, df=df, ncp=ncp, log=log.p), x=chisq)

    if (missing(shade))
      shade <- switch(length(critical.values)+1,
                      "none",
                      "right",
                      "outside",
                      stop("Specify no more than 2 critical values."))

    critical.one <- TRUE
    if (length(critical.values)==1) {
      if (shade=="right") {
        x <- seq(critical.values, max(chisq), length=51)
        shaded.area <- 1-pchisq.intermediate(q=critical.values, df=df, ncp=ncp, log.p=log.p)
      }
      else {
        x <- seq(min(chisq), critical.values, length=51)
        shaded.area <- pchisq.intermediate(q=critical.values, df=df, ncp=ncp, log.p=log.p)
      }
    }
    if (length(critical.values)==2) {
      if (sum(alpha) > 1)
        stop(paste("left.alpha=", alpha[1],
                   " + right.alpha=", alpha[2],
                   " > 1", sep=""))
      if (shade=="outside") {
        critical.one <- FALSE
        x1 <- seq(min(chisq), critical.values[1], length=51)
        x2 <- seq(critical.values[2], max(chisq), length=51)
        shaded.area <- 1-diff(pchisq.intermediate(q=critical.values, df=df, ncp=ncp, log.p=log.p))
      }
      else { ## shade == "inside"
        x <- seq(critical.values[1], critical.values[2], length=51)
        shaded.area <- diff(pchisq.intermediate(q=critical.values, df=df, ncp=ncp, log.p=log.p))
      }
    }
    if (shade != "none") {
      if (critical.one)
        polygon(x=c(x[1], x, x[length(x)]),
                y=c(0, dchisq.intermediate(x=x, df=df, ncp=ncp,
                  log=log.p), 0),
                col=col)
      else {
        polygon(x=c(x1[1], x1, x1[length(x1)]),
                y=c(0, dchisq.intermediate(x=x1, df=df, ncp=ncp,
                  log=log.p), 0),
                col=col)
        polygon(x=c(x2[1], x2, x2[length(x2)]),
                y=c(0, dchisq.intermediate(x=x2, df=df, ncp=ncp,
                  log=log.p), 0),
                col=col)
      }
    }

    axis(1, at=critical.values, tck=-.09, labels=FALSE)
    left.margin <- .15*diff(par()$usr[1:2])
    mtext(side=1, at=par()$usr[1]-left.margin, line=1, text=axis.name)
    axis(1, at=critical.values, tick=FALSE, line=2,
         labels=round(critical.values, 3), col.axis=col)

    mtext(side=1, at=par()$usr[1]-left.margin, line=3, text=axis.name, col=col)
    if (shade != "none") {
      mtext(side=1, at=par()$usr[2]+left.margin, line=1,
            text="shaded area", cex=par()$cex)
      mtext(side=1, at=par()$usr[2]+left.margin, line=3,
            text=format(shaded.area, digits=3), cex=par()$cex, col=col)
    }
  invisible(NULL)
}

chisq.observed <- function(chisq.obs, col="green",
                       df=1,
                       ncp=0,
                       log.p=FALSE,
                       axis.name="chisq",
                       shade="right",
                       shaded.area=0,
                       display.obs=TRUE) {
  abline(v=chisq.obs, col=col, lty=5)
  chisq.obs2 <- c(chisq.obs, chisq.obs)
  arrows(chisq.obs2, par()$usr[3:4]+c(-.01,.01), chisq.obs2, par()$usr[3:4],
         xpd=TRUE, col=col, length=.1)
  axis(side=1, at=chisq.obs, labels=FALSE, col=col)
  axis(side=3, at=chisq.obs, labels=FALSE, col=col)
  mtext(side=3, text=round(chisq.obs,3), at=chisq.obs, line=.5, cex=par()$cex, col=col)

  ## shade=="right"  ## we outline the right region only
  if (shade=="right" || shade=="outside") {
  x <- seq(chisq.obs, par()$usr[2], length=51)
  shaded.area <- shaded.area +
    1-pchisq.intermediate(q=chisq.obs, df=df, ncp=ncp, log.p=log.p)
  left.margin <- .15*diff(par()$usr[1:2])
    if (display.obs) {
  mtext(side=1, at=par()$usr[2]+left.margin, line=4.5,
        text=format(shaded.area, digits=3), cex=par()$cex, col=col)
  mtext(side=1, text=round(chisq.obs,3), at=chisq.obs, line=4.5, cex=par()$cex, col=col)
  mtext(side=1, at=par()$usr[1]-left.margin, line=4.5, text=axis.name, col=col)
    }
}

  if (shade=="left" || shade=="outside") {
  x <- seq(par()$usr[1], chisq.obs, length=51)
  shaded.area <- shaded.area +
    pchisq.intermediate(q=chisq.obs, df=df, ncp=ncp, log.p=log.p)
  left.margin <- .15*diff(par()$usr[1:2])
    if (display.obs) {
  mtext(side=1, at=par()$usr[2]+left.margin, line=4.5,
        text=format(shaded.area, digits=3), cex=par()$cex, col=col)
  mtext(side=1, text=round(chisq.obs,3), at=chisq.obs, line=4.5, cex=par()$cex, col=col)
  mtext(side=1, at=par()$usr[1]-left.margin, line=4.5, text=axis.name, col=col)
}
}

  lines(x=c(x[1], x, x[length(x)], x[1]),
        y=c(0, dchisq.intermediate(x=x, df=df, ncp=ncp, log=log.p), 0, 0),
        col=col, lwd=3)
  invisible(shaded.area)
}

## source("~/HH-R.package/HH/R/chisq.curve.R")
