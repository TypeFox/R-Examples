NTplot <- function(mean0, ...) {
  UseMethod("NTplot")
}

NTplot.htest <- function(mean0, ..., shiny=FALSE, NTmethod="htest") {
  if (shiny)
    shiny.NormalAndTplot(mean0, ..., NTmethod=NTmethod)
  else
    NormalAndTplot.htest(mean0, ...)
}

NTplot.NormalAndTplot <- function(mean0, ..., shiny=FALSE) {
  if (shiny)
    shiny.NormalAndTplot(mean0, ...)
  else {
    calllist <- attr(mean0, "call.list")
    dotdotdot <- list(...)
    calllist[names(dotdotdot)] <- dotdotdot
    do.call("NormalAndTplot", calllist)
  }
}

NTplot.default <- function(mean0=0, ..., shiny=FALSE, distribution.name=c("normal","z","t","binomial")) {
  distribution.name <- match.arg(distribution.name)

  xbar <- list(...)$xbar
  if (distribution.name == "binomial" &&
      missing(mean0) &&
      (is.null(xbar) || is.na(xbar)) &&
      !shiny)
    return(normalApproxBinomial(...))

  if (shiny)
    do.call("shiny.NormalAndTplot",
            c(list(mean0=mean0, distribution.name=distribution.name), list(...)))
  else
    NormalAndTplot(mean0, ..., distribution.name=distribution.name)
}

NormalAndTplot <- function(mean0, ...)
  UseMethod("NormalAndTplot")

NormalAndTplot.NormalAndTplot <- function(mean0, ...)
  do.call("NormalAndTplot", c(attr(mean0, "call.list"), list(...)))



NTplot.power.htest <-
  function(mean0, ..., shiny=FALSE, xbar=NA, ## these input values are used
           mean1, n, df, sd, distribution.name, sub, ## these input values ignored
           alpha.left, alpha.right, number.vars) { ## these input values ignored
    NT <- mean0
    sides <- ifelse(NT$alternative=="two.sided", 2, 1) ## "one.sided"

    if (is.null(NT$delta)) stop("This function is not yet working for the two-sample case\n",
                                NT$method,
                                call.=FALSE)

    NT.call.list <- list(mean1=NT$delta,
                         n=NT$n,
                         df=if (NT$method=="Two-sample t test power calculation")
                              2*(NT$n-1)
                            else
                              NT$n-1,
                         ## sd=NT$sd,
                         ## stderr=if (NT$method=="Two-sample t test power calculation")
                         ##          NT$sd*sqrt(2)/sqrt(NT$n)
                         ##        else
                         ##          NT$sd/sqrt(NT$n),
                         sd=if (NT$method=="Two-sample t test power calculation")
                              NT$sd*sqrt(2)
                            else
                              NT$sd,
                         distribution.name="t",
                         sub=if (is.null(NT$note)) NT$method else paste(NT$method, NT$note, sep="\n"),
                         alpha.left=if (sides==2)
                                      NT$sig.level/2
                                    else {if (NT$delta < 0)
                                            NT$sig.level
                                          else
                                            0},
                         alpha.right=if (sides==2)
                                       NT$sig.level/2
                                     else {if (NT$delta > 0)
                                             NT$sig.level
                                           else
                                             0},
                         number.vars=if (NT$method=="One-sample t test power calculation")
                                       1
                                     else
                                       2
                         )
    stderr <- NT.call.list$sd/sqrt(NT$n)
    NT.call.list$xlim <- range(0, NT$delta, xbar, na.rm=TRUE) + c(-3,3)*stderr
    NT.call.list$xbar <- xbar
    NT.call.list$NTmethod <- "power.htest"

    result <- do.call("NormalAndTplot.default", c(NT.call.list, list(...)))

    if (shiny) {
      ## attr(result, "call.list")$stderr <-
      ##   attr(result, "call.list")$stderr * sqrt(attr(result, "call.list")$n)
      shiny.NormalAndTplot(result, NTmethod="power.htest")
    }
    else
      result
  }

globalVariables(c('w','obs'))

NormalAndTplot.default <- function(mean0=0,
                                   mean1=NA,
                                   xbar=NA,
                                   df=Inf, n=1,
                                   sd=1,
                                   xlim=c(-3, 3)*sd/sqrt(n) + range(c(mean0, mean1, xbar), na.rm=TRUE),
                                   ylim,
                                   alpha.right=.05, alpha.left=0,
                                   float=TRUE, ntcolors="original",
                                   digits=4, digits.axis=digits,
                                   digits.float=digits,
                                   distribution.name=c("normal","z","t","binomial"),
                                   type=c("hypothesis", "confidence"),
                                   zaxis=FALSE, z1axis=FALSE,
                                   cex.z=.5, cex.prob=.6, cex.top.axis=1,
                                   main=NA, xlab, ylab,
                                   prob.labels=(type=="hypothesis"),
                                   xhalf.multiplier=1,
                                   yhalf.multiplier=1,
                                   cex.main=1,
                                   key.axis.padding=4.5,
                                   number.vars=1,
                                   sub=NULL,
                                   NTmethod="default",
                                   power=FALSE,
                                   beta=FALSE,
                                   ...) {

  type <- match.arg(type)
  if (type == "confidence") {
    if (!is.na(xbar)) mean0 <- xbar
    if (is.na(xbar)) xbar <- mean0
  }

  stderr <- sd/sqrt(n)

  if (is.list(zaxis)) {
    zaxis.list <- zaxis
    if (!all(c("at", "labels") %in% names(zaxis)))
      stop("The list 'zaxis' must contain 'at' and 'labels' components.", call.=FALSE)
    zaxis <- TRUE
  } else zaxis.list <- list()

  if (is.list(z1axis)) {
    z1axis.list <- z1axis
    if (!all(c("at", "labels") %in% names(z1axis)))
      stop("The list 'z1axis' must contain 'at' and 'labels' components.", call.=FALSE)
    z1axis <- TRUE
  } else z1axis.list <- list()

  distribution.name <- match.arg(distribution.name)
  if (distribution.name=="t" && is.infinite(df)) distribution.name <- "normal"
  if ((distribution.name=="z" || distribution.name=="normal" || distribution.name=="binomial") &&
      !is.infinite(df)) df <- Inf

  sided <- "both"
  if (alpha.left > 0 && alpha.right == 0) sided <- "left"
  if (alpha.left == 0 && alpha.right > 0) sided <- "right"
  ## if (alpha.right > 0 && alpha.left  > 0) sided <- "both"

  ncp <- (mean1-mean0)/stderr
  if (sided == "both") ncp <- abs(ncp)

  switch(distribution.name,
         normal=,
         z={
           dfunction <- function(x, mean=mean, sd=stderr, ...) dnorm(x, mean=mean, sd=sd)
           qfunction <- function (..., df, ncp) qnorm(...)
           ## pfunction <- function (..., df, ncp) pnorm(...)
           sigma.p1 <- stderr
           pfunction <- function (q, ..., df, ncp=0, sigma.p1) pnorm(q-ncp, ...)
           ## pnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
         },
         t={
           dfunction <- function(x, mean=mean, sd=stderr, ...) dt(x=(x-mean)/sd, ...) / sd
           qfunction <- qt
           sigma.p1 <- stderr
           pfunction <- function(q, ..., df, ncp=0, sigma.p1) pt(q=q, df=df, ncp=ncp, ...)
           ## pt <- function (q, df, ncp, lower.tail = TRUE, log.p = FALSE)
         },
         binomial={
           dfunction <- function(x, mean=mean, sd=stderr, ...) dnorm(x, mean=mean, sd=sd)
           qfunction <- function (..., df, ncp) qnorm(...)
           if (number.vars==1) {  ## This matches power.t.test
             p1 <- mean1
             sigma.p1 <- sqrt(p1*(1-p1)/n)
             pfunction <- function (q, ..., df, ncp=0, sigma.p1) pnorm((q-ncp)*stderr/sigma.p1, ...)
           }
           else { ## number.vars==2  ## This is not right.  I don't understand power.prop.test yet.
             p1 <- mean1
             sigma.p1 <- sqrt(p1*(1-p1)/n)
             pfunction <- function (q, ..., df, ncp=0, sigma.p1) pnorm((q-ncp)*stderr/sigma.p1, ...)
           }
         }
         )

  green127 <- ColorWithAlpha("green")
  blue127 <- ColorWithAlpha("blue")
  black127 <- ColorWithAlpha("black")

  if (length(ntcolors)==1) {

  if (ntcolors == "original") {
    col.alpha             <- "blue"
    col.notalpha          <- "lightblue"
    col.beta              <- "red"
    col.power             <- "pink"
    col.pvalue            <- "green"
    col.pvaluetranslucent <- green127
    col.critical          <- "gray50"
    col.border            <- black127
    col.text              <- "black"
    col.conf              <- "lightgreen"
  }
  if (ntcolors == "stoplight") {
    col.alpha             <- "red"
    col.notalpha          <- "honeydew2"
    col.beta              <- "orange"
    col.power             <- "pink"
    col.pvalue            <- "blue"
    col.pvaluetranslucent <- blue127
    col.critical          <- "gray50"
    col.border            <- black127
    col.text              <- "black"
    col.conf              <- "lightgreen"
  }
} else {

  ntcolors.names <- c("col.alpha", "col.notalpha", "col.beta", "col.power",
                      "col.pvalue", "col.pvaluetranslucent", "col.critical", "col.border",
                      "col.text", "col.conf")
  if (!all(match(ntcolors.names, names(ntcolors), nomatch=FALSE)))
    stop("The 'ntcolors' argument must be a named character vector with names: \n",
         paste(ntcolors.names, collapse=", "), call.=FALSE)

  for (i in names(ntcolors)) assign(i, ntcolors[i])
}


  if (type == "confidence") {
    col.alpha <- "white"
    col.notalpha <- col.conf
  }


  if (missing(ylim)) ylim <- c(0, dfunction(x=mean0, mean=mean0, sd=stderr, df=df) * 1.04)


  if (missing(ylab))
    ylab <- switch(distribution.name,
                   binomial=if (type=="hypothesis")
                              list(expression(phi(z)/sigma[p[0]]), cex=1.5, rot=0)
                            else
                              list(expression(phi(z)/s[hat(p)]), cex=1.5, rot=0),
                   normal=,
                   z=list(c(
                     expression(phi(z)/sigma[bar(x)]),
                     expression(phi(z)/sigma[bar(x)[1]-bar(x)[2]]))[number.vars],
                     cex=1.5, rot=0),
                   t=list(c(
                     expression(f[nu](t)/s[bar(x)]),
                     expression(f[nu](t)/s[bar(x)[1]-bar(x)[2]]))[number.vars],
                    cex=1.5, rot=0)
                   )

  if (distribution.name != "binomial" && missing(xlab))
    xlab <- if (number.vars==1)
              expression(w == bar(x))
            else
              expression(w == bar(x)[1] - bar(x)[2])

  ##  main <- Main(mean0, mean1, xbar, sd, n, df)
  if (missing(main) || is.null(main)) ##|| is.na(main))
    main <- list(MainSimpler(mean0, mean1, xbar, stderr, n, df, distribution.name,
                             digits=digits.axis, number.vars=number.vars, type=type),
                 cex=cex.main)

  Setup <- Base(dfunction, xlim, ylim,
                ylab=ylab,
                xlab=xlab,
                main=main, ...,
                axis.bottom=1+.13*(zaxis+z1axis)*cex.z,
                key.axis.padding=key.axis.padding,
                number.vars=number.vars,
                sub=sub)
  Setup.xlim <- Setup$x.limits
  Setup.ylim <- Setup$y.limits


  zc.right <- qfunction(p=alpha.right, lower=FALSE, df=df)
  xbarc.right <- zc.right * stderr + mean0

  if (is.infinite(xbarc.right) || is.na(xbarc.right)) xbarc.right <- Setup.xlim[2]
  zc.left <- qfunction(p=alpha.left, lower=TRUE, df=df)
  xbarc.left <- zc.left * stderr + mean0
  if (is.infinite(xbarc.left) || is.na(xbarc.left)) xbarc.left <- Setup.xlim[1]
## recover()
  Border0 <- Border(dfunction, Setup.xlim[1], Setup.xlim[2], base=TRUE, border=col.border, mean=mean0, stderr=stderr, df=df)
  if (!is.na(mean1))
    Border1 <- switch(distribution.name,
                      t=Border(dfunction, Setup.xlim[1], Setup.xlim[2], base=TRUE, border=col.border,
                        mean=mean0,  ## mean0 is correct, ncp makes the adjustment
                        stderr=stderr, df=df, ncp=ncp),
                      normal=,
                      z=Border(dfunction, Setup.xlim[1], Setup.xlim[2], base=TRUE, border=col.border, mean=mean1, stderr=stderr, df=df, ncp=ncp),
                      binomial=Border(dfunction, Setup.xlim[1], Setup.xlim[2], base=TRUE, border=col.border, mean=mean1, stderr=sigma.p1, df=df, ncp=ncp)
                      )
  Area0.left   <- Area(dfunction, Setup.xlim[1], xbarc.left, base=TRUE, col=col.alpha, mean=mean0, stderr=stderr, df=df)
  Area0.middle <- Area(dfunction, xbarc.left, xbarc.right, base=TRUE, col=col.notalpha, mean=mean0, stderr=stderr, df=df)
  Area0.right  <- Area(dfunction, xbarc.right, Setup.xlim[2], base=TRUE, col=col.alpha, mean=mean0, stderr=stderr, df=df)
  Middle0 <- Vertical(mean0, col=col.notalpha, lwd=4)

if (!is.na(mean1)) {
  ## Area1.left   <- Area(dfunction, Setup.xlim[1], xbarc.left, base=TRUE, col=col.power, mean=mean1, stderr=stderr, df=df, ncp=ncp)
  ## Area1.middle <- Area(dfunction, xbarc.left, xbarc.right, base=TRUE, col=col.beta, mean=mean1, stderr=stderr, df=df, ncp=ncp)
  ## Area1.right  <- Area(dfunction, xbarc.right, Setup.xlim[2], base=TRUE, col=col.power, mean=mean1, stderr=stderr, df=df, ncp=ncp)
  switch(distribution.name,
         t={ ## mean0 is correct, ncp makes the adjustment
           Area1.left   <- Area(dfunction, Setup.xlim[1], xbarc.left, base=TRUE, col=col.power, mean=mean0, stderr=stderr, df=df, ncp=ncp)
           Area1.middle <- Area(dfunction, xbarc.left, xbarc.right, base=TRUE, col=col.beta, mean=mean0, stderr=stderr, df=df, ncp=ncp)
           Area1.right  <- Area(dfunction, xbarc.right, Setup.xlim[2], base=TRUE, col=col.power, mean=mean0, stderr=stderr, df=df, ncp=ncp)
         },
         normal=,
         z={
           Area1.left   <- Area(dfunction, Setup.xlim[1], xbarc.left, base=TRUE, col=col.power, mean=mean1, stderr=stderr, df=df, ncp=ncp)
           Area1.middle <- Area(dfunction, xbarc.left, xbarc.right, base=TRUE, col=col.beta, mean=mean1, stderr=stderr, df=df, ncp=ncp)
           Area1.right  <- Area(dfunction, xbarc.right, Setup.xlim[2], base=TRUE, col=col.power, mean=mean1, stderr=stderr, df=df, ncp=ncp)
         },
         binomial={
           Area1.left   <- Area(dfunction, Setup.xlim[1], xbarc.left, base=TRUE, col=col.power, mean=mean1, stderr=sigma.p1, df=df, ncp=ncp)
           Area1.middle <- Area(dfunction, xbarc.left, xbarc.right, base=TRUE, col=col.beta, mean=mean1, stderr=sigma.p1, df=df, ncp=ncp)
           Area1.right  <- Area(dfunction, xbarc.right, Setup.xlim[2], base=TRUE, col=col.power, mean=mean1, stderr=sigma.p1, df=df, ncp=ncp)
         })
  Middle1 <- Vertical(mean1, col=col.power, lwd=2)
}

  if (!is.na(xbar)) {
    if (sided == "left") {
      xbar.left <- xbar
      xbar.right <- Inf ## Setup.xlim[2]
      xbar.otherside <- Inf ## xbar  ## place holder
    }
    if (sided == "right") {
      xbar.left <- -Inf ## Setup.xlim[1]
      xbar.right <- xbar
      xbar.otherside <- -Inf ##xbar  ## place holder
    }
    if (sided == "both") {
      if (xbar >= mean0) {
        xbar.left <- mean0 - (xbar - mean0)
        xbar.right <- xbar
        xbar.otherside <- xbar.left
      }
      else {
        xbar.left <- xbar
        xbar.right <- mean0 + (mean0 - xbar)
        xbar.otherside <- xbar.right
      }
    }

      top.axis.vector <- c(xbar=as.numeric(xbar), xbar.otherside=as.numeric(xbar.otherside),
                           mean0=as.numeric(mean0), mean1=as.numeric(mean1),
                           xbarc.left=as.numeric(xbarc.left), xbarc.right=as.numeric(xbarc.right))
    NotInf <- !is.infinite(top.axis.vector)
    top.axis.vector[NotInf] <- zapsmall(top.axis.vector[NotInf], digits=digits)

    Empty <- layer(panel.points(x=vector(), y=vector()))

    Borderxbar.left <- if (is.infinite(xbar.left))
                         Empty
                       else
                         Border(dfunction, Setup.xlim[1], xbar.left, base=TRUE, border=col.pvalue, lwd=2, mean=mean0, stderr=stderr, df=df)
    Borderxbar.right <- if (is.infinite(xbar.right))
                          Empty
                        else
                          Border(dfunction, xbar.right, Setup.xlim[2], base=TRUE, border=col.pvalue, lwd=2, mean=mean0, stderr=stderr, df=df)
    Vertical.xbar <- Vertical(xbar, lty=2, lwd=2, col.pvalue)
    Vertical.xbar.otherside <- Vertical(xbar.otherside, lty=2, lwd=1, col.pvalue)
    Areaxbar.left <-  if (is.infinite(xbar.left))
                        Empty
                      else
                        Area(dfunction, Setup.xlim[1], xbar.left, base=TRUE, col=col.pvaluetranslucent, mean=mean0, stderr=stderr, df=df)
    Areaxbar.right <- if (is.infinite(xbar.right))
                        Empty
                      else
                        Area(dfunction, xbar.right, Setup.xlim[2], base=TRUE, col=col.pvaluetranslucent, mean=mean0, stderr=stderr, df=df)

    ## phantom <- function() {} ## placeholder to prevent R CMD check from complaining
    xbar.expr <- as.expression(substitute(xbarsymbol==xbar, c(alist(xbarsymbol=w["obs"]), ##bar(x)["obs"]),
        list(xbar=format(top.axis.vector["xbar"], digits=digits.axis)))))
    xbar.otherside.expr <- as.expression(substitute(xbarsymbol==xbar.otherside, c(alist(xbarsymbol=w["otherside"]), ##bar(x)["otherside"]),
        list(xbar.otherside=format(top.axis.vector["xbar.otherside"], digits=digits.axis)))))
    Axis.xbar <- AxisNormal(at=xbar, labels=xbar.expr, line.col=col.pvalue, text.col=col.text, line.lwd=1, tck=5*cex.top.axis, text.cex=cex.top.axis) +
      AxisNormal(side="bottom", at=xbar, labels=format(xbar, digits=digits.axis), line.col=col.pvalue, text.col="transparent", line.lwd=1)
    Axis.xbar.otherside <- AxisNormal(at=xbar.otherside, labels=xbar.otherside.expr, line.col=col.pvalue, text.col=col.text, line.lwd=1, tck=5*cex.top.axis, text.cex=cex.top.axis) +
      AxisNormal(side="bottom", at=xbar.otherside, labels=format(xbar.otherside, digits=digits.axis), line.col=col.pvalue, text.col="transparent", line.lwd=1)
  }
  else
    {
      xbar.left <- NA
      xbar.right <- NA
      xbar.otherside <- NA
      top.axis.vector <- c(xbar=as.numeric(xbar), xbar.otherside=as.numeric(xbar.otherside),
                           mean0=as.numeric(mean0), mean1=as.numeric(mean1),
                           xbarc.left=as.numeric(xbarc.left), xbarc.right=as.numeric(xbarc.right))
      NotInf <- !is.infinite(top.axis.vector)
      top.axis.vector[NotInf] <- zapsmall(top.axis.vector[NotInf], digits=digits)
    }

  if (type == "hypothesis") {
    mean0.alist <- alist(mu0=mu[0])
    mean1.alist <- alist(mu1=mu[a])
    if (number.vars == 2) {
      mean0.alist <- alist(mu0=(mu[1]-mu[2])[0])
      mean1.alist <- alist(mu1=(mu[1]-mu[2])[a])
    }
    xbarc.left.alist <- alist(xbarsymbol=w[c]) ##bar(x)[c])
    xbarc.right.alist <- alist(xbarsymbol=w[c]) ##bar(x)[c])
  }
  else {
    mean0.alist <-  alist(mu0=w[obs]) ##bar(x))
    mean1.alist <-  alist(mu1=bar(x)) ## placeholder
    if (number.vars == 2) {
      mean0.alist <- alist(mu0=w[obs]) ##(bar(x)[1]-bar(x)[2]))
      mean1.alist <- alist(mu1=(bar(x)[1]-bar(x)[2])) ## placeholder
    }
    xbarc.left.alist <- alist(xbarsymbol=w["LCL"]) ##mu["LCL"])
    xbarc.right.alist <- alist(xbarsymbol=w["UCL"]) ##mu["UCL"])
  }
  mean0.expr <- as.expression(substitute(mu0==mean0, c(mean0.alist, list(mean0=format(top.axis.vector["mean0"], digits=digits.axis)))))
  mean1.expr <- as.expression(substitute(mu1==mean1, c(mean1.alist, list(mean1=format(top.axis.vector["mean1"], digits=digits.axis)))))
  Axis.0 <- AxisNormal(at=mean0, labels=mean0.expr, line.col=col.notalpha, line.lwd=2, tck=1*cex.top.axis, text.cex=cex.top.axis) +
    AxisNormal(side="bottom", at=mean0, labels=format(mean0, digits=digits.axis), line.col=col.notalpha, text.col="transparent", line.lwd=2)
  Axis.1 <- AxisNormal(at=mean1, labels=mean1.expr, line.col=col.power, line.lwd=2, tck=1*cex.top.axis, text.cex=cex.top.axis) +
    AxisNormal(side="bottom", at=mean1, labels=format(mean1, digits=digits.axis), line.col=col.power, text.col="transparent", line.lwd=2)
  xbarc.left.expr <-  as.expression(substitute(xbarsymbol==xbar, c(xbarc.left.alist, list(xbar=format(top.axis.vector["xbarc.left"], digits=digits.axis)))))
  xbarc.right.expr <- as.expression(substitute(xbarsymbol==xbar, c(xbarc.right.alist, list(xbar=format(top.axis.vector["xbarc.right"], digits=digits.axis)))))
  Axis.xbarc.right <- AxisNormal(at=xbarc.right, labels=xbarc.right.expr, line.col=col.critical, text.col=col.text, tck=3.00*cex.top.axis, text.cex=cex.top.axis) +
    AxisNormal(side="bottom", at=xbarc.right, labels=format(xbarc.right, digits=digits.axis), line.col=col.critical, text.col="transparent", line.lwd=1.25)
  Axis.xbarc.left <- AxisNormal(at=xbarc.left, labels=xbarc.left.expr, line.col=col.critical, text.col=col.text, tck=3.00*cex.top.axis, text.cex=cex.top.axis) +
    AxisNormal(side="bottom", at=xbarc.left, labels=format(xbarc.left, digits=digits.axis), line.col=col.critical, text.col="transparent", line.lwd=1.25)
  Vertical.xbarc.left <- Vertical(xbarc.left, col=col.critical, lty=2)
  Vertical.xbarc.right <- Vertical(xbarc.right, col=col.critical, lty=2)


  if (!is.na(xbar)) {
    pvalue.right <- if (sided != "left") pfunction((xbar.right-mean0)/stderr, df=df, lower=FALSE, sigma.p1=stderr) else NA
    pvalue.left <- if (sided != "right") pfunction((xbar.left-mean0)/stderr, df=df, lower=TRUE, sigma.p1=stderr) else NA
    pvalue <- switch(sided,
                     right=pvalue.right,
                     left= pvalue.left,
                     both= pvalue.right + pvalue.left)
  } else {
    pvalue.right <- NA
    pvalue.left <- NA
    pvalue <- NA
  }

  if (!is.na(mean1)) {
    power.right <- if (sided != "left")
                     ## pfunction((xbarc.right-mean1)/stderr, df=df, lower=FALSE, ncp=ncp)
                     pfunction((xbarc.right-mean0)/stderr, df=df, lower=FALSE, ncp=ncp, sigma.p1=sigma.p1)
                   else
                     NA
    power.left <- if (sided != "right")
                    ## pfunction((xbarc.left-mean1)/stderr, df=df, lower=TRUE, ncp=ncp)
                    pfunction((xbarc.left-mean0)/stderr, df=df, lower=TRUE, ncp=ncp, sigma.p1=sigma.p1)
                  else
                    NA
    power.total <- switch(sided,
                          right=power.right,
                          left= power.left,
                          both= power.right + power.left)
  } else {
    power.right <- NA
    power.left <- NA
    power.total <- NA
  }

  Floats <- Float(sided, type,
                  xbarc.left, xbarc.right, xbar,
                  alpha=alpha.left + alpha.right,
                  beta=1-power.total,
                  power=power.total,
                  pvalue=pvalue,
                  conf=1-(alpha.left + alpha.right),
                  Setup.xlim, Setup.ylim,
                  mean0, mean1,
                  col.alpha, col.beta, col.power, col.pvalue, col.conf, cex.prob,
                  prob.labels,
                  xhalf.multiplier, yhalf.multiplier, digits=digits.float)



  if (type == "hypothesis") {
    result <- Setup
    if (sided != "right") result <- result + Vertical.xbarc.left
    if (sided != "left")  result <- result + Vertical.xbarc.right
    result <- result + Area0.middle + Middle0
    if (!is.na(mean1)) {
      result <- result + Area1.middle
      if (sided != "right") {
        if (power.left > alpha.left)
          result <- result + Area1.left + Area0.left
        else
          result <- result + Area0.left + Area1.left
      }
      if (sided != "left") {
        if (power.right > alpha.right)
          result <- result + Area1.right + Area0.right
        else
          result <- result + Area0.right + Area1.right
      }
      if (!is.na(mean1)) result <- result + Middle1 + Border0
    }
    else
      result <- result + Area0.left + Area0.right + Border0
    ## if (!is.na(mean1)) result <- result + Area1.left + Area1.middle + Area1.right + Middle1
    ## result <- result + Area0.left + Area0.right + Border0
    if (!is.na(mean1)) result <- result + Border1
    if (!is.na(xbar)) {
      result <- result + Borderxbar.left + Borderxbar.right + Axis.xbar + Vertical.xbar + Areaxbar.left + Areaxbar.right
      if (!is.na(mean1)) result <- result + Middle1
    }
    result <- result + Zeroline
    if (!is.na(xbar) && sided == "both") result <- result + Axis.xbar.otherside + Vertical.xbar.otherside
    if (sided != "right") result <- result + Axis.xbarc.left
    if (sided != "left")  result <- result + Axis.xbarc.right
    if (!is.na(mean1)) result <- result + Axis.1
    result <- result + Axis.0
    if (float) result <- result + Floats
    if (zaxis) {
      if (length(zaxis.list) == 0) {
        z.pretty <- pretty((Setup.xlim - mean0)/stderr)
        z.at <- stderr * z.pretty + mean0
        z.labels <- signif(z.pretty, digits)
      }
      else {
        z.at <- zaxis.list$at
        z.labels <- zaxis.list$labels
      }
      result <- result + layer(
        {
          panel.axis("bottom", outside=TRUE, tck=2+cex.z, at=z.at, labels=z.labels,
                     text.cex=cex.z, rot=0, line.col="transparent")
          panel.text(x=convertX(unit(-3,   "strwidth",  data="z"), unitTo="native", valueOnly=TRUE),
                     y=convertY(unit(-2-cex.z, "strheight", data="z"), unitTo="native", valueOnly=TRUE),
                     if (distribution.name=="t") "t" else "z",
                     cex=cex.z)
        },
        data=list(z.at=z.at, z.labels=z.labels, distribution.name=distribution.name, cex.z=cex.z))

      if (z1axis && !is.na(mean1)) {
        stderrz1 <- if (distribution.name=="binomial") sigma.p1 else stderr
        if (length(z1axis.list) == 0) {
          z1.pretty <- pretty((Setup.xlim - mean1)/stderrz1)
          z1.at <- stderrz1 * z1.pretty + mean1
          z1.labels <- signif(z1.pretty, digits)
        }
        else {
          z1.at <- z1axis.list$at
          z1.labels <- z1axis.list$labels
        }
        result <- result + layer(
          {
            panel.axis("bottom", outside=TRUE, tck=3+2*cex.z, at=z1.at, labels=z1.labels,
                       text.cex=cex.z, rot=0, line.col="transparent")
            panel.text(x=convertX(unit(-3,   "strwidth",  data="z"), unitTo="native", valueOnly=TRUE),
                       y=convertY(unit(-2.75-2*cex.z, "strheight", data="z"), unitTo="native", valueOnly=TRUE),
                       if (distribution.name=="t") expression(t[1]) else expression(z[1]), cex=cex.z)
          },
          data=list(z1.at=z1.at, z1.labels=z1.labels, distribution.name=distribution.name, cex.z=cex.z))
      }
    }
  }
  else { ## confidence interval
    result <- Setup
    if (sided != "right") result <- result + Vertical.xbarc.left
    if (sided != "left")  result <- result + Vertical.xbarc.right
    result <- result + Area0.middle + Middle0
    result <- result + Area0.middle + Border0
    result <- result + Zeroline
    if (sided != "right") result <- result + Axis.xbarc.left
    if (sided != "left")  result <- result + Axis.xbarc.right
    result <- result + Axis.0
    if (float) result <- result + Floats
    if (zaxis) {
      if (length(zaxis.list) == 0) {
        z.pretty <- pretty((Setup.xlim - mean0)/stderr)
        z.at <- stderr * z.pretty + mean0
        z.labels <- signif(z.pretty, digits)
      }
      else {
        z.at <- zaxis.list$at
        z.labels <- zaxis.list$labels
      }
      result <- result + layer(
        {
          panel.axis("bottom", outside=TRUE, tck=2+cex.z, at=z.at, labels=z.labels,
                     text.cex=cex.z, rot=0, line.col="transparent")
          panel.text(x=convertX(unit(-3,   "strwidth",  data="z"), unitTo="native", valueOnly=TRUE),
                     y=convertY(unit(-2-cex.z, "strheight", data="z"), unitTo="native", valueOnly=TRUE),
                     ifelse(distribution.name=="t", "t", "z"), cex=cex.z)
        },
        data=list(z.at=z.at, z.labels=z.labels, distribution.name=distribution.name, cex.z=cex.z))
    }
  }

  result.table <- NTplotTable(distribution.name= distribution.name,
                                   type=              type,
                                   mean0=             top.axis.vector["mean0"],
                                   mean1=             top.axis.vector["mean1"],
                                   xbar=              top.axis.vector["xbar"],
                                   ## sd=                sd,
                                   df=                df,
                                   n=                 n,
                                   alpha.right=       alpha.right,
                                   alpha.left=        alpha.left,
                                   stderr=            stderr,
                                   sigma.p1=          sigma.p1,
                                   zc.right=          zc.right,
                                   xbarc.right=       top.axis.vector["xbarc.right"],
                                   zc.left=           zc.left,
                                   xbarc.left=        top.axis.vector["xbarc.left"],
                                   sided=             sided,
                                   xbar.left=         xbar.left,
                                   xbar.right=        xbar.right,
                                   xbar.otherside=    top.axis.vector["xbar.otherside"],
                                   pvalue.left=       pvalue.left,
                                   pvalue.right=      pvalue.right,
                                   pvalue=            pvalue,
                                   power.left=        power.left,
                                   power.right=       power.right,
                                   power.total=       power.total,
                                   conf.left=         alpha.left, ## ifelse(sided != "right", 1 - alpha.right, alpha.right),
                                   conf.right=        alpha.right, ## ifelse(sided != "left", 1 - alpha.left, alpha.left),
                                   conf=              1 - (alpha.left + alpha.right),
                                   mean0.alist=       mean0.alist,
                                   mean1.alist=       mean1.alist
                                     )
  attr(result,"table") <- result.table$normalTable
  attr(result,"prob") <- result.table$prob
  attr(result,"scales") <- result.table$scales

  ExpressionOrText <- function(x) {
    if (is.character(x)) return(x)
    xdp <-
      if (length(x)>1)
        deparse(x[[1]], width.cutoff=500)
      else
        deparse(x, width.cutoff=500)
    xdp
  }

  ## attr(result,"call") <- deparse(match.call())  ## works for commandline call, not from shiny
  attr(result,"call") <- paste("NTplot(mean0=", ifelse(type=="hypothesis", mean0, NA),
                               ", mean1=",mean1,
                               ", xbar=", ifelse(type=="confidence", mean0, xbar),
                               ## ", sd=",sd,
                               ", df=",df,
                               ", n=",n,
                               ", sd=",sd,
                               ", xlim=c(",xlim[1],",",xlim[2],")",
                               ", ylim=c(",ylim[1],",",ylim[2],")",
                               ", alpha.right=",alpha.right,
                               ", alpha.left=",alpha.left,
                               ", float=",float,
                               ", ntcolors=\"",ntcolors,"\"",
                               ", digits=",digits,
                               ", distribution.name=\"",distribution.name,"\"",
                               ", type=\"",type,"\"",
                               ", zaxis=",zaxis,
                               ", z1axis=",z1axis,
                               ", cex.z=",cex.z,
                               ", cex.prob=",cex.prob,
                               ", main=",ExpressionOrText(main),
                               #if (length(main)>1)
                               #             deparse(main[[1]], width.cutoff=500)
                               #           else
                               #             deparse(main, width.cutoff=500),
                               ", xlab=",ExpressionOrText(xlab),
                               #if (length(xlab)>1)
                               #             deparse(xlab[[1]], width.cutoff=500)
                               #           else
                               #             deparse(xlab, width.cutoff=500),
                               ##", ylab=",ylab,
                               ", prob.labels=",prob.labels,
                               ##", xhalf.multiplier=",xhalf.multiplier,
                               ", number.vars=",number.vars,
                               ", sub=",ifelse(is.null(sub),"NULL",paste("\"", sub, "\"", sep="")),
                               ", NTmethod=\"",NTmethod, "\"",
                               ", power=",power,
                               ", beta=",beta,
                               ")",
                               sep="")

  attr(result,"call.list") <- list(mean0=ifelse(type=="hypothesis", mean0, NA),
                                   mean1=mean1,
                                   xbar=ifelse(type=="confidence", mean0, xbar),
                                   ## sd=sd,
                                   df=df,
                                   n=n,
                                   sd=sd,
                                   xlim=xlim,
                                   ylim=ylim,
                                   alpha.right=alpha.right,
                                   alpha.left=alpha.left,
                                   float=float,
                                   ntcolors=ntcolors,
                                   digits=digits,
                                   distribution.name=distribution.name,
                                   type=type,
                                   zaxis=zaxis,
                                   z1axis=z1axis,
                                   cex.z=cex.z,
                                   cex.prob=cex.prob,
                                   main=main,
                                   xlab=xlab,
                                   ## ylab=ylab,
                                   prob.labels=prob.labels,
                                   ## xhalf.multiplier=xhalf.multiplier,
                                   number.vars=number.vars,
                                   sub=sub,
                                   NTmethod=NTmethod,
                                   power=power,
                                   beta=beta
                                   )

  attr(result,"colors")=c(
    col.alpha            =col.alpha,
    col.notalpha         =col.notalpha,
    col.beta             =col.beta,
    col.power            =col.power,
    col.pvalue           =col.pvalue,
    col.pvaluetranslucent=col.pvaluetranslucent,
    col.critical         =col.critical,
    col.border           =col.border,
    col.text             =col.text,
    col.conf             =col.conf)

  class(result) <- c("NormalAndTplot", class(result))

  if (type=="hypothesis" && (power || beta))
    result <- powerplot(result, power=power, beta=beta, ...)

  result
}


print.NormalAndTplot <- function(x, tablesOnPlot=TRUE, plot=TRUE,
                                 scales=FALSE, prob=FALSE, call=FALSE,
                                 ..., cex.table=.7, digits=getOption("digits"),
                                 position.2=.17) {

  if (scales) {
    cat("\nscales\n")
    print(attr(x, "scales"))
  }
  if (prob) {
    cat("\nprobabilities\n")
    print(attr(x, "prob"))
  }
  if (call)
    cat("\ncall\n", attr(x, "call"), "\n")

  if (plot) {

    if (!tablesOnPlot) {
      return(NextMethod("print"))
    }

    if (tablesOnPlot && !is.null(list(...)$position))
      stop("position= argument is incompatible with tablesOnPlot=TRUE")

    NextMethod("print", position=c(0, position.2, 1, 1))
    ## lattice:::print.trellis(x, position=c(0, position.2, 1, 1))

## new
    NTplotTheme <- ttheme_default(
      core   =list(fg_params=list(parse=FALSE, cex = cex.table, hjust=1, x=0.9 ), bg_params = list(fill="grey98", lwd=1.5, col="white")),
##    core   =list(fg_params=list(parse=FALSE, cex = cex.table, hjust=0.1, x=.2), bg_params = list(fill="grey98", lwd=1.5, col="white")),
      colhead=list(fg_params=list(parse=TRUE,  cex = cex.table, fontface="bold"), bg_params = list(fill="grey95", lwd=1.5, col="white")),
      rowhead=list(fg_params=list(parse=TRUE,  cex = cex.table, fontface="bold"), bg_params = list(fill="grey95", lwd=1.5, col="white"))
    )
## new end

    old.digits <- options(digits=digits)

    pushViewport(viewport(x = 0, y = 0,
                          width = .6,
                          height = .2,
                          just = c("left", "bottom")))

    axs <- attr(x, "scales")
    axsf <- axs
    ## axsf[ 1,] <- format(zapsmall(axs[ 1,,drop=FALSE], digits=digits),
    ##                     digits=digits)
    ## axsf[-1,] <- format(zapsmall(axs[-1,,drop=FALSE], digits=digits),
    #                      digits=digits)
    axsf <- format(zapsmall(axs, digits=digits),
                   digits=digits)
## ## old
##     gridExtra::grid.table(axsf,
##                           parse=TRUE,
##                           core.just="right", row.just="center", col.just="center",
##                           gpar.corefill = gpar(fill = "grey98", col = "white"),
##                           gpar.rowfill = gpar(fill = "grey95", col = "white"),
##                           gpar.colfill = gpar(fill = "grey95", col = "white"),
##                           gpar.rowtext = gpar(cex = cex.table, fontface = "bold"),
##                           gpar.coltext = gpar(cex = cex.table, fontface = "bold"),
##                           gpar.coretext = gpar(cex = cex.table))
## ## old end
## new
    gridExtra::grid.table(axsf,
                          theme=NTplotTheme)
## new end
    popViewport()

    pushViewport(viewport(x = .55, y = 0,
                          width = .4,
                          height = .2,
                          just = c("left", "bottom")))
## ## old
##     gridExtra::grid.table(format(round(attr(x, "prob"), digits=digits), nsmall=4),
##                           parse=TRUE,
##                           core.just="right", row.just="center", col.just="center",
##                           gpar.corefill = gpar(fill = "grey98", col = "white"),
##                           gpar.rowfill = gpar(fill = "grey95", col = "white"),
##                           gpar.colfill = gpar(fill = "grey95", col = "white"),
##                           gpar.rowtext = gpar(cex = cex.table, fontface = "bold"),
##                           gpar.coltext = gpar(cex = cex.table, fontface = "bold"),
##                           gpar.coretext = gpar(cex = cex.table))
## ## old end
## new
    gridExtra::grid.table(format(round(attr(x, "prob"), digits=digits), nsmall=4),
                          theme=NTplotTheme)
## new end
    popViewport()

    options(old.digits)
  }

  invisible(x)
}

