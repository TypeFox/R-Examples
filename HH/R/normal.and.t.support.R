## library(lattice)
## library(latticeExtra)

Base <- function(dfunction,
                 xlim,
                 ylim,
                 ylab=deparse(substitute(dfunction)),
                 xlab="x",
                 main=NULL,
                 sub=NULL,
                 number.vars=1,
                 key.axis.padding=4.5,
                 axis.bottom=1,
                 ...,
                 par.settings=list(
                   clip=list(panel=FALSE),
                   layout.heights=list(key.axis.padding=key.axis.padding,
                                       axis.bottom=axis.bottom,
                                       axis.xlab.padding=2,
                                       bottom.padding=ifelse(number.vars==1, 1, 2.5)),
                   layout.widths=list(left.padding=7)
                 )) {
  xyplot(ylim ~ xlim, type="n", ylab=ylab, xlab=xlab, main=main, sub=sub, ...,
         par.settings=par.settings)
}

Curve <- function(dfunction, xlow, xhigh,
                  col=NA, border="black", lwd=1,
                  base=FALSE, closed=TRUE, ..., mean=0, stderr=1, df=Inf, ncp=0) {
  xx <- seq(xlow, xhigh, length=201)
  yy <- dfunction(xx, mean=mean, sd=stderr, df=df, ncp=ncp)
  if (closed) {
    if (base) {
      xx <- c(xx[1], xx, xx[length(xx)])
      yy <- c(0,     yy, 0             )
    }
    layer(panel.polygon(x=xx,
                        y=yy,
                        col=col,
                        border=border, lwd=lwd),
          data=list(xx=xx, yy=yy, col=col, border=border, lwd=lwd))
  }
  else
        layer(panel.lines(x=xx,
                        y=yy,
                        col=border,
                        lwd=lwd),
          data=list(xx=xx, yy=yy, col=col, border=border, lwd=lwd))
}


Area <- function(..., border=NA) {
  Curve(..., border=border)
}

Border <- function(..., col=NA) {
  Curve(..., col=col, closed=FALSE)
}

Vertical <- function(v, lty=4, col="black", lwd=1) {
  layer(panel.abline(v=v, lty=lty, col=col, lwd=lwd),
        data=list(v=v, lty=lty, col=col, lwd=lwd))
}

Zeroline <- layer(panel.refline(h=0, col="gray70"))

AxisNormal <- function(side="top", at, labels=TRUE, outside=TRUE, rot=0,
                       line.col="black", text.col="black", tck=1, text.cex=1, ...) {
  layer(panel.axis(side=side, at=at, labels=labels, outside=outside, rot=rot,
                   line.col=line.col, text.col=text.col, tck=tck, text.cex=text.cex),
        data=list(side=side, at=at, labels=labels, outside=outside, rot=rot,
          line.col=line.col, text.col=text.col, tck=tck, text.cex=text.cex))
}



Float <- function(sided, type,
                  xbarc.left, xbarc.right, xbar,
                  alpha, beta, power, pvalue,
                  conf,
                  xlim, ylim,
                  mean0, mean1,
                  col.alpha, col.beta, col.power, col.pvalue, col.conf, cex.prob=.6,
                  prob.labels,
                  xhalf.multiplier, yhalf.multiplier, digits) {
  labels <- if (prob.labels)
              as.expression(
                c(substitute(symbol==value, c(alist(symbol=alpha),  list(value=round(alpha,  digits)))),
                  substitute(symbol==value, c(alist(symbol=beta),   list(value=round(beta,   digits)))),
                  substitute(symbol==value, c(alist(symbol=1-beta), list(value=round(power,  digits)))),
                  substitute(symbol==value, c(alist(symbol="p"),      list(value=round(pvalue, digits)))),
                  substitute(symbol==value, c(alist(symbol="Conf"),   list(value=round(conf,   digits)))))
              )
            else
              c(alp=round(alpha, digits),
                bet=round(beta, digits),
                pow=round(power, digits),
                pva=round(pvalue, digits),
                cnf=round(conf, digits))

  if (sided=="left" ||
      (sided=="both" && (!is.na(mean1) && mean1 < mean0) || (!is.na(xbar) && xbar < mean0))
      ) {
    AddSubtract  <- c(-1, 1, -1, -1.5, 1.5)*2.5
    xbarc <- xbarc.left
  } else {
    AddSubtract  <- c(1, -1, 1, 1.5, -1.5)*2.5
    xbarc <- xbarc.right
  }


  xvalues  <- c(alpha=xbarc, beta=xbarc, power=xbarc, pvalue=xbar, conf=xbarc) +
    AddSubtract * diff(xlim)/40
  yvalues <- c(alpha=3, beta=5, power=7, pvalue=1, conf=3) * diff(ylim)/30

  border <- c(col.alpha=col.alpha,
              col.beta=col.beta,
              col.power=col.power,
              col.pvalue=col.pvalue,
              col.conf=col.conf)
  xhalf <- if (prob.labels)
             cex.prob * diff(xlim)/9 * xhalf.multiplier
           else
             cex.prob * diff(xlim)/15 *xhalf.multiplier
  yhalf <- cex.prob * diff(ylim)/30 *yhalf.multiplier

  subscripts <- c(alpha=(type=="hypothesis"),
                  beta=(type=="hypothesis") && (!is.na(mean1)),
                  power=(type=="hypothesis") && (!is.na(mean1)),
                  pvalue=(type=="hypothesis") && (!is.na(xbar)),
                  conf=(type=="confidence"))
  xvalues <- xvalues [subscripts]
  yvalues <- yvalues [subscripts]
  labels  <- labels  [subscripts]
  border  <- border  [subscripts]

  layer({
    panel.rect(xvalues-xhalf, yvalues-yhalf, xvalues+xhalf, yvalues+yhalf,
               col="white", border=border, lwd=2)
    old.digits <- options(digits=digits)
    panel.text(xvalues, yvalues, labels=labels, cex=cex.prob)
    options(old.digits)
  }, data=list(
       xvalues=xvalues, yvalues=yvalues, labels=labels, border=border,
       xhalf=xhalf, yhalf=yhalf, cex.prob=cex.prob, digits=digits)
        )
}



globalVariables(c('mu', 'bar', 'x', 'sigma', 'nu', 's'))

MainSimpler <- function(mean0, mean1, xbar, stderr, n, df, distribution.name, digits, number.vars, type) {

  argnames <- alist(mu0=mu[0], mu1=mu[1], ## xbarsymbol=bar(x),
                    sesymbol=sigma[bar(x)],
                    sesymbolt=s[bar(x)],
                    sesymbolph=sigma[p[0]],
                    sesymbolpc=s[hat(p)],
                    nn=n, dfdf=nu)
  if (number.vars==2)
    argnames[c("sesymbol","sesymbolt")] <-
      alist(sesymbol=sigma[bar(x)[1]-bar(x)[2]],
            sesymbolt=s[bar(x)[1]-bar(x)[2]])

  argvals <- list(mean0=mean0, mean1=mean1, xbar=xbar,
                  n=n, df=df, se=format(stderr, digits=digits))
  argsboth <- c(argnames, argvals)

  main <- switch(distribution.name,
                 t=substitute("t: " * sesymbolt==se * ", " ~ nn==n * ", " ~ dfdf==df, argsboth),
                 normal=,
                 z=substitute("normal: " * sesymbol==se * ", " ~ nn==n, argsboth),
                 binomial=if (type=="hypothesis")
                            substitute("normal approximation to the binomial: " * sesymbolph==se * ", " ~ nn==n, argsboth)
                          else
                            substitute("normal approximation to the binomial: " * sesymbolpc==se * ", " ~ nn==n, argsboth)
                 )

  as.expression(main)
}






ColorWithAlpha <- function(colorname, alpha=127) { ## designed for scalar color name from colors()
  rgb(t(col2rgb(colorname)), maxColorValue=255, alpha=alpha)
}
