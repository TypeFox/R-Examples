plotCL <- function(model, what="c", fit=TRUE, swap=FALSE, series=NULL, sex=NULL, years=NULL, lengths=NULL, axes=TRUE,
                   same.limits=TRUE, log=FALSE, base=10, eps.log=1e-5, main="", xlab="", ylab="", cex.main=1.2,
                   cex.lab=1, cex.axis=0.8, cex.strip=0.8, col.strip="gray95", strip=strip.custom(bg=col.strip),
                   las=!fit, tck=c(1,fit)/2, tick.number=5, lty.grid=3, col.grid="gray", pch=16, cex.points=0.5,
                   col.points="black", lty.lines=1, lwd.lines=2, col.lines=c("red","blue"), plot=TRUE, ...)
{
  ## 1  Define functions
  panel.bubble <- function(x, y, ...)  # bubble plot obs in one single-sex panel
  {
    panel.abline(v=pretty(x,tick.number), h=pretty(y,tick.number), lty=lty.grid, col=col.grid)
    panel.xyplot(x, y, ...)
  }
  panel.obs <- function(x, y, ...)  # overlay male and female lines in year panels
  {
    panel.abline(v=pretty(x,tick.number), lty=lty.grid, col=col.grid)
    panel.superpose(x, y, ...)
  }
  panel.fit <- function(x, y, subscripts, col.points, col.lines, ...)  # overlay obs and fit in sex:year panels
  {
    panel.abline(v=pretty(x,tick.number), lty=lty.grid, col=col.grid)
    panel.superpose.2(x, y, subscripts=subscripts, col.symbol=col.points[subscripts], col.line=col.lines[subscripts],
                      ...)
  }

  ## 2  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  what <- match.arg(what, c("c","s"))
  relation <- if(same.limits) "same" else "free"
  las <- as.numeric(las)

  ## 3  Prepare data (extract, rearrange, filter, transform)
  if(what == "c")
  {
    if(!any(names(model)=="CLc"))
    {
      what <- "s"
      warning("Element 'CLc' (commercial C@L) not found; assuming what=\"s\" was intended")
    }
    else
      x <- model$CLc
  }
  if(what == "s")  # value of 'what' may have changed
  {
    if(!any(names(model)=="CLs"))
      stop("Element 'CLs' (survey C@L) not found; try plotCL(x, what=\"c\") to plot commercial C@L")
    x <- model$CLs
  }
  x <- data.frame(Series=rep(x$Series,2), Year=rep(x$Year,2), SS=rep(x$SS,2), Sex=rep(x$Sex,2), Length=rep(x$Length,2),
                  ObsFit=c(rep("Obs",nrow(x)),rep("Fit",nrow(x))), P=c(x$Obs,x$Fit))
  if(is.null(series))
    series <- unique(x$Series)
  if(is.null(sex))
    sex <- unique(x$Sex)
  if(is.null(years))
    years <- unique(x$Year)
  if(is.null(lengths))
    lengths <- unique(x$Length)
  if(length(series) > 1)
  {
    series <- series[1]
    warning("More than one C@L series found; assuming series=\"", series, "\" was intended", sep="")
  }
  ok.series  <- x$Series %in% series;  if(!any(ok.series))  stop("Please check if the 'series' argument is correct")
  ok.sex     <- x$Sex    %in% sex;     if(!any(ok.sex))     stop("Please check if the 'sex' argument is correct")
  ok.years   <- x$Year   %in% years;   if(!any(ok.years))   stop("Please check if the 'years' argument is correct")
  ok.lengths <- x$Length %in% lengths; if(!any(ok.lengths)) stop("Please check if the 'lengths' argument is correct")
  x <- x[ok.series & ok.sex & ok.years & ok.lengths,]
  if(nrow(x) == 0)
    stop("Empty data frame; please check if subsetting args (series, sex, years, lengths) had mistakes")
  nsexes <- length(unique(x$Sex))
  if(log)
    x$P <- log(x$P+eps.log, base=base)

  ## 4  Prepare plot (set pars, vectorize args, create list args)
  col.points <- rep(col.points, length.out=2)
  col.lines <- rep(col.lines, length.out=2)
  mymain <- list(label=main, cex=cex.main)
  myxlab <- list(label=xlab, cex=cex.lab)
  myylab <- list(label=ylab, cex=cex.lab)
  myrot <- switch(as.character(las), "0"=list(x=list(rot=0),y=list(rot=90)), "1"=list(x=list(rot=0),y=list(rot=0)),
                  "2"=list(x=list(rot=90),y=list(rot=0)), "3"=list(x=list(rot=90),y=list(rot=90)))
  myscales <- c(list(draw=axes,relation=relation,cex=cex.axis,tck=tck,tick.number=tick.number), myrot)
  mystrip <- strip.custom(bg=col.strip)
  mytext <- list(cex=cex.strip)

  ## 5  Create trellis object
  fixed.ylim <- FALSE
  if(nsexes==1 && !fit)
  {
    x <- x[x$ObsFit=="Obs",]
    col.points <- ifelse(x$P==0|is.na(x$P), "transparent", col.points)
    mycex <- cex.points*sqrt(x$P/mean(x$P,na.rm=TRUE)) + 1/1000  # cex=0 is illegal on PDF device
    myformula <- if(!swap) Year~Length|switch(what,c="Commercial C@L","Survey C@L")
    else Length~Year|switch(what,c="Commercial C@L","Survey C@L")
    graph <- xyplot(myformula, data=x, panel=panel.bubble,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    pch=pch, cex=mycex, col=col.points, ...)
    graph$y.limits <- rev(graph$y.limits)
    fixed.ylim <- TRUE
  }
  if(nsexes==1 && fit)
  {
    myformula <- if(!swap) P~Length|factor(Year) else P~Year|factor(Length)
    graph <- xyplot(myformula, data=x, groups=x$ObsFit, panel=panel.fit, type=c("l","p"), as.table=TRUE,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    pch=pch, cex=cex.points, col.points=col.points[factor(x$Sex)], lty=lty.lines, lwd=lwd.lines,
                    col.lines=col.lines[factor(x$Sex)], ...)
  }
  if(nsexes==2 && !fit)
  {
    x <- x[x$ObsFit=="Obs",]
    myformula <- if(!swap) P~Length|factor(Year) else P~Year|factor(Length)
    graph <- xyplot(myformula, data=x, groups=x$Sex, panel=panel.obs, type="l", as.table=TRUE,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    lty=lty.lines, lwd=lwd.lines, col=col.lines, ...)
  }
  if(nsexes==2 && fit)
  {
    myformula <- if(!swap) P~Length|factor(Year)*factor(Sex) else P~Year|factor(Length)*factor(Sex)
    graph <- xyplot(myformula, data=x, groups=x$ObsFit, panel=panel.fit, type=c("l","p"), as.table=TRUE,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    pch=pch, cex=cex.points, col.points=col.points[factor(x$Sex)], lty=lty.lines, lwd=lwd.lines,
                    col.lines=col.lines[factor(x$Sex)], ...)
  }
  if(!log && !fixed.ylim)  # leave ylim alone if log-transformed or bubble plot
  {
    if(is.list(graph$y.limits))                                                 # set lower ylim to 0
      graph$y.limits <- lapply(graph$y.limits, function(y){y[1]<-0;return(y)})  # multi-panel plot
    else
      graph$y.limits[1] <- 0                                                    # single-panel plot
  }

  ## 6  Finish
  if(plot)
  {
    print(graph)
    invisible(x)
  }
  else
  {
    invisible(graph)
  }
}
