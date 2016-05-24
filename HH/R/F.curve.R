"F.setup" <-
  function(df1=1,
           df2=Inf,
           ncp=0,
           log.p=FALSE,
           xlim.in=c(0, 5),
           ylim.in=range(c(0, 1.1*df.intermediate(x=seq(.5,1.5,.01), df1=df1, df2=df2, ncp=ncp, log=log.p))),
           main.in=main.calc, ylab.in="F density",
           ...) {
    main.calc <- parse(text=
                         paste("F~density: ~~ nu[1]==", df1,
                               "~~~~~nu[2]==", df2, sep=""))
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

"F.curve" <-
  function(df1=1,
           df2=Inf,
           ncp=0,
           log.p=FALSE,
           alpha=.05,
           critical.values=f.alpha,
           f=seq(0, par()$usr[2], length=109),
           shade="right", col=par("col"),
           axis.name="f",
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

    f.alpha <- if (length(alpha)==1)
      qf.intermediate(p=1-alpha, df1=df1, df2=df2, ncp=ncp, log.p=log.p)
    else
      qf.intermediate(p=c(alpha[1], 1-alpha[2]), df1=df1, df2=df2, ncp=ncp, log.p=log.p)

    lines(y=df.intermediate(x=f, df1=df1, df2=df2, ncp=ncp, log=log.p), x=f)

    if (missing(shade))
      shade <- switch(length(critical.values)+1,
                      "none",
                      "right",
                      "outside",
                      stop("Specify no more than 2 critical values."))

    critical.one <- TRUE
    if (length(critical.values)==1) {
      if (shade=="right") {
        x <- seq(critical.values, max(f), length=51)
        shaded.area <- 1-pf.intermediate(q=critical.values, df1=df1, df2=df2, ncp=ncp, log.p=log.p)
      }
      else {
        x <- seq(min(f), critical.values, length=51)
        shaded.area <- pf.intermediate(q=critical.values, df1=df1, df2=df2, ncp=ncp, log.p=log.p)
      }
    }
    if (length(critical.values)==2) {
      if (sum(alpha) > 1)
        stop(paste("left.alpha=", alpha[1],
                   " + right.alpha=", alpha[2],
                   " > 1", sep=""))
      if (shade=="outside") {
        critical.one <- FALSE
        x1 <- seq(min(f), critical.values[1], length=51)
        x2 <- seq(critical.values[2], max(f), length=51)
        shaded.area <- 1-diff(pf.intermediate(q=critical.values, df1=df1, df2=df2, ncp=ncp, log.p=log.p))
      }
      else { ## shade == "inside"
        x <- seq(critical.values[1], critical.values[2], length=51)
        shaded.area <- diff(pf.intermediate(q=critical.values, df1=df1, df2=df2, ncp=ncp, log.p=log.p))
      }
    }
    if (shade != "none") {
      if (critical.one)
        polygon(x=c(x[1], x, x[length(x)]),
                y=c(0, df.intermediate(x=x, df1=df1, df2=df2, ncp=ncp,
                  log=log.p), 0),
                col=col)
      else {
        polygon(x=c(x1[1], x1, x1[length(x1)]),
                y=c(0, df.intermediate(x=x1, df1=df1, df2=df2, ncp=ncp,
                  log=log.p), 0),
                col=col)
        polygon(x=c(x2[1], x2, x2[length(x2)]),
                y=c(0, df.intermediate(x=x2, df1=df1, df2=df2, ncp=ncp,
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

F.observed <- function(f.obs, col="green",
                       df1=1,
                       df2=Inf,
                       ncp=0,
                       log.p=FALSE,
                       axis.name="f",
                       shade="right",
                       shaded.area=0,
                       display.obs=TRUE) {
  f.obs2 <- c(f.obs, f.obs)
  if (display.obs) {
    abline(v=f.obs, col=col, lty=5)
    arrows(f.obs2, par()$usr[3:4]+c(-.01,.01), f.obs2, par()$usr[3:4],
           xpd=TRUE, col=col, length=.1)
    axis(side=1, at=f.obs, labels=FALSE, col=col)
    axis(side=3, at=f.obs, labels=FALSE, col=col)
    mtext(side=3, text=round(f.obs,3), at=f.obs, line=.5, cex=par()$cex, col=col)
  }

  ## shade=="right"  ## we outline the right region only
  if (shade=="right" || shade=="outside") {
    x <- seq(f.obs, par()$usr[2], length=51)
    shaded.area <- shaded.area +
      1-pf.intermediate(q=f.obs, df1=df1, df2=df2, ncp=ncp, log.p=log.p)
    left.margin <- .15*diff(par()$usr[1:2])
    if (display.obs) {
      mtext(side=1, text=round(f.obs,3), at=f.obs, line=4.5, cex=par()$cex, col=col)
      mtext(side=1, at=par()$usr[1]-left.margin, line=4.5, text=axis.name, col=col)
    }
  }

  if (shade=="left" || shade=="outside") {
    x <- seq(par()$usr[1], f.obs, length=51)
    shaded.area <- shaded.area +
      pf.intermediate(q=f.obs, df1=df1, df2=df2, ncp=ncp, log.p=log.p)
    left.margin <- .15*diff(par()$usr[1:2])
    if (display.obs) {
      mtext(side=1, text=round(f.obs,3), at=f.obs, line=4.5, cex=par()$cex, col=col)
      mtext(side=1, at=par()$usr[1]-left.margin, line=4.5, text=axis.name, col=col)
    }
  }

  ## mtext(side=1, at=par()$usr[2]+left.margin, line=4.5,
  ##      text=format(shaded.area, digits=3), cex=par()$cex, col=col)
  lines(x=c(x[1], x, x[length(x)], x[1]),
        y=c(0, df.intermediate(x=x, df1=df1, df2=df2, ncp=ncp, log=log.p), 0, 0),
        col=col, lwd=3)
  invisible(shaded.area)
}

## source("~/HH-R.package/HH/R/F.curve.R")
