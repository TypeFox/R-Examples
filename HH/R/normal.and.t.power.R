powerplot <- function (nt, ...) {
  UseMethod("powerplot")
}

## NTpbplot
powerplot.NormalAndTplot <- function(nt, power=TRUE, beta=FALSE, ...,
                                     hh=if (power && beta) c(6,2,2) else c(6,2)) {

  if (length(nt$condlevels[[1]]) > 1) {
    warning("The input object ", substitute(nt), " has more than one panel.  Only the first panel is used.", call.=FALSE)
    nt <- nt[1]
  }

  ntp <- if (power) NormalAndTPower(nt, "power", ...)
  ntb <- if (beta)  NormalAndTPower(nt, "beta",  ...)

  nt.xlab <- nt$xlab
  nt <- update(nt, xlab=" ") +
    layer(panel.axis(at=mean(current.panel.limits()$x),
                     labels=nt.xlab,
                     outside=TRUE, rot=0, line.col="transparent"),
          data=list(nt.xlab=nt.xlab))

  if (power && beta)
    cdpb <- c(density=nt, power=ntp, beta=ntb, layout=c(1,3), y.same=FALSE, x.same=TRUE)
  if ( power && !beta)
    cdpb <- c(density=nt, power=ntp,           layout=c(1,2), y.same=FALSE, x.same=TRUE)
  if (!power &&  beta)
    cdpb <- c(density=nt,            beta=ntb, layout=c(1,2), y.same=FALSE, x.same=TRUE)

  result <-
  resizePanels(h=hh,
               update(cdpb,
                      as.table=TRUE,
                      between=list(y=c(
                                     3+1.25*attr(nt, "call.list")$zaxis +
                                       1.25*attr(nt, "call.list")$z1axis,
                                     2)),
                      strip=FALSE, strip.left=TRUE,
                      scales=list(y=list(rot=0), x=list(limits=nt$x.limits)),
                      ylab=NULL)
               )

  result$y.scales$at <- list(FALSE,NULL,NULL)
  result$y.scales$labels <- list(FALSE,NULL,NULL)
  for (pp in seq(along=result$condlevels[[1]])[-1]) {
    result$y.limits[[pp]] <- c(0,1)
    result$y.scales$at[[pp]] <- ntp$y.scales$at
    result$y.scales$labels[[pp]] <- ntp$y.scales$labels
  }


  attributes(result)[c("table", "prob", "scales", "call" )] <-
    attributes(nt)[c("table", "prob", "scales", "call" )]
  result
}


NormalAndTPower <- function(nt,
                            which=c("power","beta"),
                            digits=4, digits.top.axis=digits,
                            digits.left=digits,
                            col.power=attr(nt, "color")["col.power"],
                            col.beta=attr(nt, "color")["col.beta"],
                            cex.top.axis=1, cex.left.axis=1,
                            lwd.reference=4, lwd.line=2,
                            main=which, ...) {
  which <- match.arg(which)

  number.vars <- attr(nt, "call.list")$number.vars
  mean0.alist <- alist(mu0=mu[0])
  mean1.alist <- alist(mu1=mu[a])
  if (number.vars == 2) {
    mean0.alist <- alist(mu0=(mu[1]-mu[2])[0])
    mean1.alist <- alist(mu1=(mu[1]-mu[2])[a])
  }

  mean0.char <- as.character(mean0.alist)
  mean1.char <- as.character(mean1.alist)

  tables <- attr(nt, "table")
  ## scales <- attr(nt, "scales")
  ## prob <- attr(nt, "prob")

 ## nt.call <- as.list(nt$call[[2]])[-1]

  nt.type <- ifelse(is.na(tables["Confidence", "Probability"]), "hypothesis", "confidence")
  if (nt.type == "confidence")
    stop("power isn't meaningful for confidence graphs", call.=FALSE)

  nt.mean0 <- tables["bar(x)", mean0.char] ## "mu[0]"]
  nt.mean1 <- tables["bar(x)", mean1.char] ## "mu[1]"]
  nt.xbar  <- tables["bar(x)", 'w["obs"]']
  nt.stderr  <- tables["bar(x)", "sigma"]
##  nt.sd <- ifelse(is.null(nt.sd), 1, nt.sd)

  nt.n  <- tables["bar(x)", "n"]
##  nt.n <- ifelse(is.null(nt.n), 1, nt.n)

  nt.df  <- tables["bar(x)", "df"]
##  nt.df <- ifelse(is.null(nt.df), Inf, nt.df)

  nt.distribution <- ifelse(is.infinite(nt.df), "z", "t")

  nt.x <- seq(nt$x.limits[1], nt$x.limits[2], length=101)

  nt.crit <- attr(nt,"table")["bar(x)", c("w[crit.L]", "w[crit.R]")]
##  nt.crit <- nt.crit[nt.crit != nt$x.limits]

  nt.mean0 <- attr(nt,"table")["bar(x)", mean0.char]

  alphaNotMissing <- !is.na(tables["alpha", c("w[crit.L]", "w[crit.R]")])
  nt.sided <- if (sum(alphaNotMissing) == 2)
                "both"
              else
                ifelse(alphaNotMissing[1], "left", "right")

  switch(nt.distribution,
         normal=,
         z=pfunction <- function(q, ..., df, ncp=0) pnorm(q-ncp, ...),
         t=pfunction <- pt
         )

  warn.old <- options(warn=-1) ## pt gives precision warnings at the extreme tails.
  nt.ypower <- switch(nt.sided,
                      left=pfunction(
                        q=(nt.crit[1]-nt.mean0)/nt.stderr,
                        df=nt.df,
                        lower.tail=TRUE, ncp=(nt.x-nt.mean0)/nt.stderr),
                      right= 1 - pfunction(
                        q=(nt.crit[2]-nt.mean0)/nt.stderr,
                        df=nt.df,
                        lower.tail=TRUE, ncp=(nt.x-nt.mean0)/nt.stderr),
                      both=pfunction(
                        q=(nt.crit[1]-nt.mean0)/nt.stderr,
                        df=nt.df,
                        lower.tail=TRUE, ncp=(nt.x-nt.mean0)/nt.stderr)
                      + 1 - pfunction(
                        q=(nt.crit[2]-nt.mean0)/nt.stderr,
                        df=nt.df,
                        lower.tail=TRUE, ncp=(nt.x-nt.mean0)/nt.stderr)
                      )
  options(warn.old)

  ## nt.alpha <- attr(nt,"table")["alpha",]
  nt.beta  <- attr(nt,"table")["beta","Probability"]
  nt.power <- attr(nt,"table")["power","Probability"]
  ## nt.p     <- attr(nt,"table")["p",]

  if (which == "power") {
    ylab=list(expression(atop(scriptstyle("power ="), 1-beta)), rot=0)
    nt.powerbeta <- nt.power
    col.powerbeta <- col.power
  }
  else {
    nt.ypower <- 1 - nt.ypower
    ylab=list(expression(beta), rot=0)
    nt.powerbeta <- nt.beta
    col.powerbeta <- col.beta
  }

  mean1.expr <- as.expression(substitute(mu1==mean1, c(mean1.alist, ## mu[1]),
      list(mean1=format(nt.mean1, digits=digits.top.axis)))))

  xyplot(xlim=nt$x.limits, ylim=c(0,1), type="l", col="gray60", lwd=lwd.line,
         scales=list(y=list(at=seq(0,1,.25), labels=seq(0,1,.25))),
         nt.ypower ~ nt.x, main=main, xlab=as.expression(mean1.alist), ##  expression(mu[1]),
         ylab=ylab,
         par.settings=list(clip=list(panel=FALSE))) +
           layer({
             panel.abline(h=nt.powerbeta, v=nt.mean1,
                          lty=3, col=col.powerbeta, lwd=lwd.reference)
             panel.axis("left", at=nt.powerbeta,
                       labels=format(nt.powerbeta, digits=digits.left),
                        text.cex=cex.left.axis,
                        tck=6, outside=TRUE,
                        line.col=col.powerbeta, line.lwd=lwd.reference, line.lty=3)
             panel.axis("top", at=nt.mean1,
                        labels=mean1.expr,
                        rot=0, outside=TRUE,
                        text.cex=cex.top.axis,
                        line.col=col.powerbeta, line.lwd=lwd.reference, line.lty=3)
           },
                 data=list(nt.mean1=nt.mean1, nt.powerbeta=nt.powerbeta,
                   col.powerbeta=col.powerbeta, digits.top.axis=digits.top.axis,
                   digits.left=digits.left, cex.top.axis=cex.top.axis,
                   cex.left.axis=cex.left.axis,
                   lwd.reference=lwd.reference,
                           mean1.expr=mean1.expr))
}
