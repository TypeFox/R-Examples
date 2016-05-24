## Time Series Plotting Functions
## Displays for Direct Comparison of ARIMA Models
## The American Statistician, May 2002, Vol. 56, No. 2, pp. 131-138
## Richard M. Heiberger, Temple University
## Paulo Teles, Faculdade de Economia do Porto

## The contents are copyrighted under the Gnu Public License.
##
## July 2007 revision.  Works with R-2.5.1 and S-Plus 8.
##

## November 2014 revision:
## Change default arrangement for printing panels in plots.
## Remove S-Plus support.
## Redo with latticeExtra::useOuterStrips

## x is an "arima" object in S-Plus and an "Arima" object in R
npar.arma <- function(x, arima=FALSE) {
  if.R(r=npar.rarma(x$arma, arima=arima),
       s=npar.sarma(x$model, arima=arima))
}


## S-Plus (really an input style model)
npar.sarma <- function(model, arima=FALSE) {
  if (!is.null(names(model))) model <- list(model)
  subscr <- if (arima) 1:3 else -2
  sum(sapply(model, function(x, subscr) x$order[subscr], subscr=subscr))
}


## R (the arma result)
npar.rarma <- function(arma, arima=FALSE) {
  if (arima) sum(arma[-5]) else sum(arma[1:4])
}



acfplot <- function(rdal, type="acf",
                    main=paste("ACF of std.resid:", rdal$series,
                               "   model:",         rdal$model),
                    lag.units=rdal$tspar["frequency"],
                    lag.lim=range(pretty(rdal[[type]]$lag)*lag.units),
                    lag.x.at=pretty(rdal[[type]]$lag)*lag.units,
                    lag.x.labels={tmp <- lag.x.at
                                  tmp[as.integer(tmp)!=tmp] <- ""
                                  tmp},
                    lag.0=TRUE,
                    xlim=xlim.function(lag.lim/lag.units),
                    ...) {
  if (!is.R() && is.null(lag.units)) stop('time series must be an "rts".')
  xlim.function <- function(lu) {incr <- .02*lu[2]
                                 lu + c(-1,1)*incr}
  rdal[[type]]$tmp <-
    ordered(paste("ar:",rdal[[type]]$ar," ma:",rdal[[type]]$ma, sep=""),
            levels=outer(unique(rdal[[type]]$ma),
                         unique(rdal[[type]]$ar),
                         FUN=function(ma,ar) paste("ar:",ar," ma:",ma, sep="")))
  rdal[[type]]$subset <- (lag.0 | rdal[[type]]$lag != 0)

  xyplot(acf ~ lag | tmp, data=rdal[[type]], subset=subset,
         n.used=rdal$n.used,
         layout=c(length(unique(rdal[[type]]$ma)),
                  length(unique(rdal[[type]]$ar))),
         panel=panel.acf,
         strip=function(...) strip.default(..., style = 1),
         as.table=TRUE,
         scales=list(alternating=FALSE,
                     x=list(at=lag.x.at/lag.units, labels=lag.x.labels)),
         main=main, xlim=xlim, ...)
}

residplot <- function(rdal,
                      main=paste("std.resid:", rdal$series,
                                 "   model:",  rdal$model),
                      ...) {
  rdal$std.resid$tmp <-
    ordered(paste("ar:",rdal$std.resid$ar," ma:",rdal$std.resid$ma, sep=""),
            levels=outer(unique(rdal$std.resid$ma),
                         unique(rdal$std.resid$ar),
                         FUN=function(ma,ar) paste("ar:",ar," ma:",ma, sep="")))

  xyplot(resid ~ time | tmp, data=rdal$std.resid,
         layout=c(length(unique(rdal$std.resid$ma)),
                  length(unique(rdal$std.resid$ar))),
         panel=panel.std.resid,
         strip=function(...) strip.default(..., style = 1),
         as.table=TRUE,
         scales=list(alternating=FALSE),
         main=main, ...)
}

gofplot <- function(rdal,
                    main=paste("P-value for gof:", rdal$series,
                               "   model:",         rdal$model),
                    lag.units=rdal$tspar["frequency"],
                    lag.lim=range(pretty(rdal$gof$lag)*lag.units),
                    lag.x.at=pretty(rdal$gof$lag)*lag.units,
                    lag.x.labels={tmp <- lag.x.at
                                  tmp[as.integer(tmp)!=tmp] <- ""
                                  tmp},
                    xlim=xlim.function(lag.lim/lag.units),
                    pch=16, ...) {
xlim.function <- function(lu) {incr <- .02*lu[2]
                                 lu + c(-1,1)*incr}
  rdal$gof$tmp <-
    ordered(paste("ar:",rdal$gof$ar," ma:",rdal$gof$ma, sep=""),
            levels=outer(unique(rdal$gof$ma),
                         unique(rdal$gof$ar),
                         FUN=function(ma,ar) paste("ar:",ar," ma:",ma, sep="")))
  xyplot(p ~ lag | tmp, data=rdal$gof,
         pch=pch,
         ylim=range(-.05, .05, rdal$gof$p, na.rm=TRUE),
         layout=c(length(unique(rdal$gof$ma)),
                  length(unique(rdal$gof$ar))),
         panel=panel.gof,
         strip=function(...) strip.default(..., style = 1),
         as.table=TRUE,
         scales=list(alternating=FALSE,
                     x=list(at=lag.x.at/lag.units, labels=lag.x.labels)),
         main=main, xlim=xlim, ...)
}

panel.acf <- function(..., n.used) {
# based on acf.list piece of arima.diag.plot
  panel.abline(h=0, col="gray50")
  if (length(n.used)==1) this.n <- n.used
  else {
    this.frame <- sys.parent()
    this.cell <-
           packet.number()
    this.n <- t(n.used)[this.cell]
  }
  c1 <- 2/sqrt(this.n)
  panel.abline(h= c(-c1,c1), lty=2, lwd=.75, col="gray50")
  panel.xyplot(..., type="h")
}


tsdiagplot <- function(x,
                       p.max=2, q.max=p.max,
                       model=c(p.max, 0, q.max), ## S-Plus
                       order=c(p.max, 0, q.max), ## R
                       lag.max=36, gof.lag=lag.max,
                       armas=arma.loop(x, order=order,
                           series=deparse(substitute(x)), ...),
                       diags=diag.arma.loop(armas, x,
                                            lag.max=lag.max,
                                            gof.lag=gof.lag),
                       ts.diag=rearrange.diag.arma.loop(diags),
                       lag.units=ts.diag$tspar["frequency"],
                       lag.lim=range(pretty(ts.diag$acf$lag))*lag.units,
                       lag.x.at=pretty(ts.diag$acf$lag)*lag.units,
                       lag.x.labels={tmp <- lag.x.at
                                  tmp[as.integer(tmp)!=tmp] <- ""
                                  tmp},
                       lag.0=TRUE,
                       main, lwd=0,
                       ...) { ## ... is seasonal in R (if used)
  sum.armas <- summary(armas)

  sigma2 <- aicsigplot(sum.armas$sigma2, "sigma2", sum.armas$series,
                             model=sum.armas$model,
                             xlab=NULL, ylab=NULL,
                             main="sigma2")
  resid  <- residplot(ts.diag,
                            xlab=NULL, ylab=NULL, lwd=lwd,
                            main="standardized residuals")
  aic    <- aicsigplot(sum.armas$aic,    "aic",    sum.armas$series,
                             model=sum.armas$model,
                             xlab=NULL, ylab=NULL,
                             main="aic")
  acf    <-       acfplot(ts.diag, xlab=NULL, ylab=NULL,
                          lag.units=lag.units,
                          lag.lim=lag.lim,
                          lag.x.at=lag.x.at,
                          lag.x.labels=lag.x.labels,
                          lag.0=lag.0,
                          main="ACF of standardized residuals")
  pacf   <-       acfplot(ts.diag, type="pacf", xlab=NULL, ylab=NULL,
                          lag.units=lag.units,
                          lag.lim=lag.lim,
                          lag.x.at=lag.x.at,
                          lag.x.labels=lag.x.labels,
                          lag.0=lag.0,
                          main="PACF of standardized residuals")
    y.range <- range(c(acf$y.limits, pacf$y.limits))
    acf$y.limits  <- y.range
    pacf$y.limits <- y.range

  gof    <- gofplot(ts.diag, xlab=NULL, ylab=NULL,
                          lag.units=lag.units,
                          lag.lim=lag.lim,
                          lag.x.at=lag.x.at,
                          lag.x.labels=lag.x.labels,
                          main="P-value for gof", cex=.25)

  method <- armas[[nrow(armas), ncol(armas)]]$method
  mm <- arima.model(armas[[nrow(armas), ncol(armas)]])
  mm[[1]]$order[c(1,3)] <- c("p","q")
  mmc <- as.character(mm)
  series <- sum.armas$series
  main.diag <- paste("series:", series, "  model:", mmc, "  by", method)
  result <- list(
                 sigma2=sigma2,
                 aic=aic, resid=resid, acf=acf, pacf=pacf, gof=gof,
                 main.diag=main.diag, series=series, model=mmc, method=method)
  if (!missing(main)) result$main.diag <- main
  class(result) <- "tsdiagplot"
  result
}


print.tsdiagplot <- function(x, ..., portrait=FALSE) {
  pages <- list(...)$pages
  if (!portrait) {
    if(is.null(pages) || pages==1)
      print1.tsdiagplot(x)
    else
      print2.tsdiagplot(x)
  }
  else {
    print3.tsdiagplot(x)
  }
  invisible(x)
}


print1.tsdiagplot <- function(x) {
##  is.main.title <- as.logical(match("main",names(x),FALSE))
##  if (is.main.title)
##  oldpar.oma <- par(oma=c(0,0,3,0))

xrg.left <- .44
xa.left <- .73
y.top.bottom <- .48
y.bottom.top <- .47

  print(x$resid, position=c( xrg.left, .00,          1.00,  y.bottom.top), more=TRUE)
  print(x$aic,   position=c( xa.left,  y.top.bottom, 1.00,  .95),          more=TRUE)
  print(x$acf,   position=c( .00,      y.top.bottom,  .45,  .95),          more=TRUE)
  print(x$pacf,  position=c( .00,      .00,           .45,  y.bottom.top), more=TRUE)
  print(x$gof,   position=c( xrg.left, y.top.bottom,  .75,  .95),         more=FALSE)

  title.trellis(x$main.diag)

##   if (is.main.title) {
##     par(oma=c(0,0,0,0))
##     title.trellis(x$main)
##     par(oldpar.oma)
##   }

  invisible(x)
}

print2.tsdiagplot <- function(x) {
  is.main.title <- as.logical(match("main", names(x), FALSE))
## page 1
  if (is.main.title)
    oldpar.oma <- par(oma=c(0,0,3,0))
  else
    oldpar.oma <- par(oma=par()$oma)

acf.bottom <- .48
pacf.top <- .47
aic.sigma2.top <- .54


  print(x$acf,   position=c( .00, acf.bottom,  .45,  .95), more=TRUE)
  print(x$pacf,  position=c( .00, .00,  .45,  pacf.top), more=TRUE)

  print(x$gof,   position=c( .50, .55,  .95,  .95), more=TRUE)

  print(x$aic,   position=c( .45, .00,  .75,  aic.sigma2.top), more=TRUE)
  print(x$sigma2,position=c( .70, .00, 1.00,  aic.sigma2.top), more=FALSE)

  title.trellis(x$main.diag)

  if (is.main.title) {
    par(oma=c(0,0,0,0))
    title.trellis(x$main)
  }

## page 2
  if (is.main.title)
    par(oma=c(0,0,6,0))
  else
    par(oma=c(0,0,3,0))

  print(x$resid)

  if (is.main.title)
    par(oma=c(0,0,3,0))
  else
    par(oma=c(0,0,0,0))
  title.trellis(x$main.diag)

  if (is.main.title) {
    par(oma=c(0,0,0,0))
    title.trellis(x$main)
  }

  par(oldpar.oma)

## return
  invisible(x)
}


print3.tsdiagplot <- function(x) {
## bottom
  print(x$resid, position=c( .00, .00, 1.00, .35), more=TRUE)
## right middle
  print(rearrangeAICplot(x$aic),
                 position=c( .60, .35, 1.00, .65), more=TRUE)
## left top
  print(update(x$acf, lwd=2.5),   position=c( .00, .65,  .60, .95), more=TRUE)
## left middle
  print(update(x$pacf, lwd=2.5),  position=c( .00, .35,  .60, .65), more=TRUE)
## right top
  print(update(x$gof, cex=.7),   position=c( .60, .65, 1.00, .95), more=FALSE)
## title
  title.trellis(x$main.diag)

  invisible(x)
}

rearrangeAICplot <- function(aic) {
  newlegend <- aic$legend
  aic$legend <- NULL
  names(newlegend) <- "bottom"
  newlegend$bottom$args$key$space <- "bottom"
  newlegend$bottom$args$key$title <- "ma               ar"
  newlegend$bottom$args$key$cex.title <- .8
  newlegend$bottom$args$key$text$col    <- newlegend$bottom$args$key$lines$col
  newlegend$bottom$args$key$text$cex    <- .8
  newlegend$bottom$args$key$column <- length(newlegend$bottom$args$key$text[[1]])
  newlegend$bottom$args$key$lines$lwd <- 2
##recover()
  update(aic,
         layout=c(2,1),
         between=list(x=1),
         legend=newlegend,
         xlab=c("ar\n","ma\n"),
         lwd=2)
}

if (FALSE) {
print3.tsdiagplot(elnino.diagplot)

ddco2.diagplot <-
tsdiagplot(armas=ddco2.loop, diags=ddco2.diags,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)

print3.tsdiagplot(ddco2.diagplot)

pdf("test.pdf", height=14, width=9)
print3.tsdiagplot(ddco2.diagplot)
dev.off()

}



tsacfplots <- function(x,
                       ylab=deparse(substitute(x)),
                       x.name=ylab[[1]],
                       main=paste("Series:", x.name),
                       lag.at=NULL,
                       lag.max=NULL,
                       lag.units=NULL,
                       lag.0=TRUE,
                       ...) {
  xplot <- seqplot(x, ylab=ylab,
                   main=NULL,
                   h=0, ...)

  if (missing(lag.at))
    acf.plots <- acf.pacf.plot(x, series=x.name, lag.max=lag.max, lag.0=lag.0,
                               main=NULL, ...)
  else
    acf.plots <- acf.pacf.plot(x, series=x.name, lag.max=lag.max, lag.0=lag.0,
                               main=NULL,
                               lag.at=lag.at, lag.units=lag.units, ...)

  result <- list(xplot=xplot, acf.plots=acf.plots, main=main)
  class(result) <- "tsacfplots"
  result
}

print.tsacfplots <- function(x,
                             ts.pos=c(.00, .00,  .70, 1.00),
                             acf.pos=c(.65, .10, 1.00,  .90),
                             ...,
                             portrait=FALSE,
                             ts.pos.portrait=c(0, .3, 1, 1),
                             acf.pos.portrait=c(.1, 0, .9, .35)) {
  if (!portrait) {
    print(x$xplot,     position=ts.pos,  more=TRUE)
    print(x$acf.plots, position=acf.pos, more=FALSE)
    title.trellis(x$main)
  }
  else {
    print(x$xplot,                   position=ts.pos.portrait,  more=TRUE)
    print(update(x$acf.plots, layout=c(2,1), between=list(x=2), lwd=2),
          position=acf.pos.portrait, more=FALSE)
    title.trellis(x$main)
  }
  invisible(x)
}

acf.pacf.plot <- function(x,
                          ylab=NULL,
                          series=deparse(substitute(x)),
                          main=paste("ACF and PACF:", series),
                          lag.max,
                          lag.units=frequency(x),
                          lag.at=pretty(apacf$lag),
                          lag.labels=lag.at*lag.units,
                          lag.0=TRUE,
                          strip=TRUE, strip.left=FALSE,
                          ...) {
  if (missing(lag.max) || is.null(lag.max))
    acf.x <- acf(x, plot=FALSE)
  else
    acf.x <- acf(x, plot=FALSE, lag.max=lag.max)
  acf.x$series <- series

  if (missing(lag.max) || is.null(lag.max))
    pacf.x <- acf(x, type="partial", plot=FALSE)
  else
    pacf.x <- acf(x, type="partial", plot=FALSE, lag.max=lag.max)
  pacf.x$series <- series

  apacf <- data.frame(acf=c(acf.x$acf, 1, pacf.x$acf),
                      lag=c(acf.x$lag, 0, pacf.x$lag),
                      type=c(rep("acf",    length( acf.x$acf)),
                             rep("pacf", 1+length(pacf.x$acf))))
  apacf$subset <- (lag.0 | apacf$lag != 0)

  xyplot(acf ~ lag | type, data=apacf, subset=subset,
         n.used=acf.x$n.used,
         panel=panel.acf,
         as.table=TRUE,
         strip=strip, strip.left=strip.left,
         between=list(y=1),
         scales=list(alternating=FALSE, x=list(at=lag.at, labels=lag.labels)),
         ylab=ylab,
         main=main,
         layout=c(1,2), ...)
}




arma.loop <- function(x,
                      model,             ## S-Plus
                      order, seasonal,   ## R
                      series=deparse(substitute(x)), ...)
{
if( !missing(model) || missing(order) )  ## missing(seasonal) is OK
          stop("Please use valid arguments for arima.")
model <-
            if (missing(seasonal))
              list(list(order=order))
            else
              list(list(order=order), seasonal)
  ## The ... is used to pass arguments to arima in R
  if (!is.null(names(model))) model <- list(model)
  p.max <- model[[1]]$order[1]
  q.max <- model[[1]]$order[3]
  z <- matrix(vector("list", (p.max+1)*(q.max+1)), p.max+1, q.max+1,
              dimnames=list(as.character(0:p.max), as.character(0:q.max)))
  dz <- dimnames(z)
  for (q in dz[[2]]) for (p in dz[[1]]) {
    model[[1]]$order[1] <- as.numeric(p)
    model[[1]]$order[3] <- as.numeric(q)

      dps <- NULL
      if (!is.list(model))
        dpm <- deparse(model, width.cutoff=100, control=NULL)
      else {
        if (is.list(model) && length(model)==1) {
          if (is.list(model[[1]]))
            dpm <- deparse(model[[1]][[1]], width.cutoff=100, control=NULL)
          else
            dpm <- deparse(model[[1]], width.cutoff=100, control=NULL)
        }
        else
          {dpm <- deparse(model[[1]][[1]], width.cutoff=100, control=NULL)
           dps <- deparse(model[[2]], width.cutoff=100, control=NULL)}
      }
      arima.text <- paste("arima(x, order=",
                          dpm,
                          ## deparse(model[[1]][[1]],
                          ##         width.cutoff=100, control=NULL),
                          if (!is.null(dps)) paste(", seasonal=", dps),
                          ", ...)")
      tmp <-  try(eval(parse(text=arima.text)))
      input.method <- list(...)$method
      tmp$method <- ifelse(is.null(input.method), "CSS-ML", input.method)


    tmp$series <- series
    z[[p, q]] <- tmp
  }
  class(z) <- "arma.loop"
  z
}

"[.arma.loop" <- function(x, ..., drop = TRUE)
{
  class.x <- class(x)
  x <- unclass(x)
  result <- NextMethod("[", drop=drop)
  class(result) <- class.x
  result
}

seqplot.forecast <- function(...)
  .Defunct("seqplotForecast", package="HH")

seqplotForecast <- function(xts, forecast, multiplier=1.96,
                            series=deparse(substitute(observed)), ylim,
                            CI.percent=round((1-2*(1-pnorm(multiplier)))*100,2),
                            main = paste(
                              series, " with forecast + ",
                              CI.percent, "% CI", sep=""),
                            xlab=NULL, ylab=NULL,
                            ...) ## x.at, xlim
{
  ## on input, xts is observed data
    fm <- forecast$pred
    fs <- forecast$se
    xts <- ts(c(xts, fm), time(xts)[1], frequency=frequency(xts))

  upperbound <- fm + multiplier * fs
  lowerbound <- fm - multiplier * fs
  if (missing(ylim))
    ylim <- range(xts, upperbound, lowerbound)
  seqplot(xts=xts,
          ylim=ylim,
          xlab=xlab,
          ylab=ylab,
          main=main,
          x.pred=time(fm),
          forecast=fm,
          upperbound=upperbound,
          lowerbound=lowerbound,
          ...)
}




aicsigplot <- function(z, z.name=deparse(substitute(z)), series.name="ts",
                       model=NULL,
                       xlab="", ylab=z.name,
                       main=paste(z.name,  series.name, model),
                       layout=c(1,2), between=list(x=1,y=1), ...)
{
  ar <- as.numeric(dimnames(z)[[1]][row(z)])
  ma <- as.numeric(dimnames(z)[[2]][col(z)])
  rc <- max(dim(z))

  zz <- c(z,z)
  f1 <- c(ar,ma)
  f2 <- c(ma,ar)
  f3 <- factor(rep(1:2, rep(length(z),2)),
               labels=c(paste(z.name, "~ ar | ma"),
                        paste(z.name, "~ ma | ar")))
  xyplot(zz ~ f1 | f3,
         panel=function(x,y,subscripts,groups, ...) {
               panel.superpose(x,y,subscripts,groups, ...)
             },
         groups=f2, type="b", pch=16, cex=.5,
         layout=layout,
         ylab=ylab,
         xlab=xlab,
#        xlab=" ",
         main=main,
         strip=function(...) strip.default(..., style=1),
         scales=list(x=list(at=unique(f1), alternating=FALSE, relation="free"),
                     y=list(alternating=FALSE)),
         key=list(text=list(format(sort(unique(c(ar,ma))))),
           space="right",
           size=2, cex=.75, between=1,
           border=1,
           lines=Rows(trellis.par.get("superpose.line"), 1:rc)),
         ...
         )
}



summary.arma.loop <- function(object, ...) {
  dimnames.object <- dimnames(object)
  sigma2 <- array(NA, dim=dim(object), dimnames=dimnames.object)
  aic <- sigma2
  coef <- array(vector("list", prod(dim(object))),
                dim=dim(object), dimnames=dimnames.object)
  t.coef <- coef
  model.object <- array(list(), dim=dim(object), dimnames=dimnames.object)
  model <- array("", dim=dim(object), dimnames=dimnames.object)
  for (q in dimnames.object[[2]]) for (p in dimnames.object[[1]]) {
    model.object[[p,q]] <- arima.model(object[[p,q]])
    model[[p,q]]   <- as.character(model.object[[p,q]])
    if (is.null(object[[p,q]]$sigma2)) next
    sigma2[p,q] <- object[[p,q]]$sigma2
    if (is.null(object[[p,q]]$aic)) next
    aic[p,q]    <- object[[p,q]]$aic

    if (npar.sarma(model.object[[p,q]], arima=TRUE) == 0) next
    mi.coef <- coef(object[[p,q]])
    coef[[p,q]]  <- mi.coef
    rec.sd <- diag.maybe.null(object[[p,q]]$var.coef)^-.5
    t.coef[[p,q]] <- mi.coef * rec.sd
  }
  class(coef) <- "arma.loop.list"
  names(coef) <- model
  class(t.coef) <- "arma.loop.list"
  names(t.coef) <- model

  mm <- model.object[[length(object)]]
  mm[[1]]$order[c(1,3)] <- c("p","q")
  model <- as.character(mm)

  list(series=object[[2]]$series, model=model, sigma2=sigma2, aic=aic,
       coef=coef, t.coef=t.coef)
}

print.arma.loop <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.arma.loop.list <- function(object, ...) {
  rr <- length(object)
  cc <- length(object[[rr]])
  result <- matrix(nrow=rr, ncol=cc, dimnames=list(names(object), names(object[[rr]])))
  for (i in 1:rr) {
    if (length(object[[i]]) == 0) next
    result[i,names(object[[i]])] <- object[[i]]
  }
  result
}

print.arma.loop.list <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}



diag.arma.loop <- function(z, x=stop("The time series x is needed in S-Plus when p=q=0."),
                           lag.max=36, gof.lag=lag.max) {
  dz <- dimnames(z)
  d <- array(vector("list", length(z)), dim=dim(z), dimnames=dz)
  for (q in dz[[2]]) for (p in dz[[1]]) {

    ## based on stats:::tsdiag.Arima
    arima.diag <- NA ## placeholder to make R-2.6.0dev happy
    object <- z[[p,q]]
    tmp <- list()
    tmp$series <- object$series
    rs <- object$resid
    if (is.null(rs)) next
    tmp$acf.list <- acf(rs, lag.max=lag.max, plot = FALSE, na.action = na.pass)
    tmp$pacf.list <- pacf(rs, lag.max=lag.max, plot = FALSE, na.action = na.pass)
    tmp$std.resid <- rs/sqrt(object$sigma2)
    nlag <- gof.lag
    pval <- numeric(nlag)
    gof <- t(sapply(1:nlag, Box.test, x=rs,
                    type = "Ljung-Box"))[,c("statistic","parameter","p.value","parameter")]
    dimnames(gof)[[2]][c(2,4)] <- c("df","lag")
    tmp$gof <- data.frame(apply(gof, 2, unlist))
         tmp$gof$lag <- tmp$gof$lag/object$arma[5] ## frequency
    d[[p,q]] <- tmp
  }
  mm <- arima.model(z[[length(z)]])
  mm[[1]]$order[c(1,3)] <- c("p","q")
  attr(d,"model") <- as.character(mm)
  class(d) <- "diag.arma.loop"
  d
}

"[.diag.arma.loop" <- function(x, ..., drop = TRUE)
{
  class.x <- class(x)
  x <- unclass(x)
  result <- NextMethod("[", drop=drop)
  class(result) <- class.x
  attr(result,"model") <- attr(x,"model")
  result
}

rearrange.diag.arma.loop <- function(z) {
  dz <- dimnames(z)
  acf     <- list(acf=vector(), lag=vector(), ar=vector(), ma=vector())
  pacf    <- acf
  n.used  <- array(0, dim=dim(z), dimnames=dz)
  gof     <- list(gof=vector(), df=vector(), p=vector(), lag=vector(),
                  ar=vector(), ma=vector())
  resid   <- list(resid=vector(), time=vector(), ar=vector(), ma=vector())

  for (q in dz[[2]]) for (p in dz[[1]]) {
    if (length(z[[p,q]]) == 0) {

      acf$acf <- c(acf$acf, NA)
      acf$lag <- c(acf$lag, NA)
      acf$ar  <- c(acf$ar, p)
      acf$ma  <- c(acf$ma, q)

      pacf$acf <- c(pacf$acf, NA)
      pacf$lag <- c(pacf$lag, NA)
      pacf$ar  <- c(pacf$ar, p)
      pacf$ma  <- c(pacf$ma, q)

      n.used[p,q] <- 0

      gof$gof <- c(gof$gof, NA)
      gof$df  <- c(gof$df,  NA)
      gof$p   <- c(gof$p,   NA)
      gof$lag <- c(gof$lag, NA)
      gof$ar  <- c(gof$ar, p)
      gof$ma  <- c(gof$ma, q)

      resid$resid <- c(resid$resid, NA)
      resid$time  <- c(resid$time, NA)
      resid$ar    <- c(resid$ar, p)
      resid$ma    <- c(resid$ma, q)

      next
    }
    acf$acf <- c(acf$acf, z[[p,q]]$acf.list$acf)
    acf$lag <- c(acf$lag, z[[p,q]]$acf.list$lag)
    acf$ar  <- c(acf$ar, rep(p, length(z[[p,q]]$acf.list$acf)))
    acf$ma  <- c(acf$ma, rep(q, length(z[[p,q]]$acf.list$acf)))

    pacf$acf <- c(pacf$acf, 1, z[[p,q]]$pacf.list$acf)
    pacf$lag <- c(pacf$lag, 0, z[[p,q]]$pacf.list$lag)
    pacf$ar  <- c(pacf$ar, rep(p, 1+length(z[[p,q]]$pacf.list$acf)))
    pacf$ma  <- c(pacf$ma, rep(q, 1+length(z[[p,q]]$pacf.list$acf)))

    n.used[p,q] <- z[[p,q]]$acf.list$n.used

    gof$gof <- c(gof$gof, z[[p,q]]$gof$statistic)
    gof$df  <- c(gof$df,  z[[p,q]]$gof$df)
    gof$p   <- c(gof$p,   z[[p,q]]$gof$p.value)
    gof$lag <- c(gof$lag, z[[p,q]]$gof$lag)
    gof$ar  <- c(gof$ar, rep(p, length(z[[p,q]]$gof$df)))
    gof$ma  <- c(gof$ma, rep(q, length(z[[p,q]]$gof$df)))

    resid$resid <- c(resid$resid, as.numeric(z[[p,q]]$std.resid))
    resid$time  <- c(resid$time, as.numeric(time(z[[p,q]]$std.resid)))
    resid$ar    <- c(resid$ar, rep(p,length(z[[p,q]]$std.resid)))
    resid$ma    <- c(resid$ma, rep(q,length(z[[p,q]]$std.resid)))
  }
  lz <- length(z)
        tspar <- function(x) {
           result <- tsp(x)
           result[2] <- 1/result[3]
           names(result) <-c("start", "deltat", "frequency")
           result
         }

  list(series=z[[lz]]$series, tspar=tspar(z[[lz]]$std.resid),
       model=attr(z,"model"),
       acf=acf, pacf=pacf, n.used=n.used, gof=gof, std.resid=resid)
}


diag.maybe.null <- function(x, ...) {
  if (length(x) == 0) return(x)
  if ((length(x) == 1) && is.na(x)) return(x)
  else diag(x, ...)
}



panel.std.resid <- function(...) {
# based on std.resid piece of arima.diag.plot
  panel.xyplot(..., type="h")
  panel.abline(h=0)
}

panel.gof <- function(...) {
# based on gof piece of arima.diag.plot
  panel.xyplot(...)
  panel.abline(h=0.05, lty=2)
}


# seqplot permits abline parameters

seqplot <- function(xts, ...)
UseMethod("seqplot")


seqplot.ts <- function(xts, pch.seq=letters, groups=as.numeric(cycle(xts)),
                        x.at=pretty(time(xts)[groups==min(groups)]),
                        x.labels,
                        ylab=deparse(substitute(xts)),
                        ...) {
  groups <- rep(groups, length=length(xts))
  if(missing(x.at) && length(unique(groups))==1 && missing(pch.seq)) {
    x.at <- pretty(time(xts))
    pch.seq <- 16
  }
  if(missing(x.at) && (length(xts) <= frequency(xts)*5))
    x.at <- pretty(time(xts))
  if(missing(pch.seq) && (frequency(xts) == 12))
    pch.seq <- substring(month.name, 1, 1)
  if(missing(pch.seq) && (frequency(xts) == 7)) {
    day.name <- weekdays(seq(as.Date("2007-07-29"), by="day", length=7))  ## "Sunday":"Saturday"
    pch.seq <- substring(day.name, 1, 1)
  }
  if (missing(x.labels)) x.labels <- format(x.at)
  seqplot.default(xts, pch.seq=pch.seq, groups=groups,
                  scales=list(x=list(at=x.at, labels=x.labels)),
                  ylab=ylab,
                  ...)
}




seqplot.default <- function(xts,
                            pch.seq=letters,
                            groups=as.numeric(cycle(xts)),
                            a=NULL, b=NULL, h=NULL, v=NULL,
                            ylab=deparse(substitute(xts)),
                            xlab="Time",
                            lwd=1, lty=c(1,3),
                            type="b",
                            col=trellis.par.get("superpose.symbol")$col,
                            col.line="gray60",
                            ...) {
  xyplot(xts ~ as.numeric(time(xts)), groups=rep(groups, length=length(xts)),
         xlab=xlab,
         ylab=ylab,
         lwd=lwd, lty.in=lty,
         pch=rep(pch.seq, length=max(groups)),
         type=type,
         ab.params=list(a=a,b=b,v=v,h=h),
         col=col, col.line=col.line,
         panel=function(x, y, subscripts, groups, pch, ab.params,
                        col=col, col.line=col.line, ...,
                        type, x.pred, forecast, upperbound, lowerbound, lty.in) {
           ab.params <- ab.params[!sapply(ab.params, is.null)]
           if (length(ab.params)) do.call("panel.abline", c(ab.params, list(col="gray70")))
           panel.xyplot(x, y, type="l", col=col.line, lty=lty.in[1], ...)
           if (type!="l")
             panel.superpose(x, y, subscripts, groups, pch=pch, col=col, ...)
           if (!missing(upperbound))
             panel.xyplot(x=x.pred, y=upperbound, type="l", lty=lty.in[2])
           if (!missing(lowerbound))
             panel.xyplot(x=x.pred, y=lowerbound, type="l", lty=lty.in[2])
         },
         ...)
}


gof.calculation <- function(acf.list, gof.lag, n, n.parms) {
### copied from arima.diag
  gof <- numeric(0)
  if(gof.lag > 0) {
    gof <- n * cumsum(acf.list$acf[2:(n.parms + gof.lag + 1)]^2)
    gof.df <- 1:gof.lag
    gof.lags <- acf.list$lag[(n.parms + 1):(n.parms + gof.lag)]
    gof.p.value <- numeric(gof.lag)
    gof <- gof[(n.parms + 1):(n.parms + gof.lag)]
    for(i in 1:gof.lag)
      gof.p.value[[i]] <- (1 - pchisq(gof[[i]], gof.df[[i]]))
    gof <- list(statistic = gof, df = gof.df,
                p.value = gof.p.value, lag = gof.lags)
  }
}


