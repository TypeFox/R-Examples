
#  tfplot provides generic methods for plotting time series objects. 
# Plot methods will probably do some processing
#  and eventually call tfOnePlot for each panel.

tfplot <- function(x, ...)  UseMethod("tfplot")

tfplot.default <- function(x, ..., tf=tfspan(x , ...), start=tfstart(tf), end=tfend(tf),
	series=seq(nseries(x)), 
	Title=NULL, title=Title, subtitle=NULL,
	lty = 1:5, lwd = 1, pch = 1, col = 1:6, cex = NULL,
	xlab=NULL, ylab=seriesNames(x), xlim = NULL, ylim = NULL,
	graphs.per.page=5, par=NULL, reset.screen=TRUE,
	Xaxis="auto", L1=NULL,
	YaxisL=TRUE, YaxisR=FALSE, Yaxis.lab.rot = "vertical",
	splitPane=NULL,
	lastObs=FALSE, source=NULL,
	footnote=NULL, footnoteLeft=footnote, footnoteRight=NULL,
	legend=NULL, legend.loc="topleft")
    {#  ... before other args means abbreviations do not work, but otherwise
     # positional matching seems to kick in and an object to be plotted gets used
     #  for start.
     if (!is.tframed(x)) {
       if (is.matrix(x) || is.vector(x)) x <- ts(x)
       else return(plot(x))
       }
     if(inherits(x, "TSmodel"))
        stop("tfplot does not know how to plot a model. ",
             "Consider simulating the model: tfplot(simulate(model)) ",
             "or evaluating the model with data: tfplot(l(model, data)).")
     if( !is.numeric(x) )
        stop("tfplot.default does not know how to plot this object.")
     old.par <- par(par)
     on.exit(par(old.par)) 
     mr <- par()$mar
     names <- seriesNames(x)
     Ngraphs <- min(length(series), graphs.per.page)
     if( (!is.list(xlim)) && (2 == length(xlim)))
              xlim <- rep(list(xlim), length(series))
     if( (!is.list(ylim)) && (2 == length(ylim)))
              ylim <- rep(list(ylim), length(series))
     if(reset.screen)  {
        par(mfcol = c(Ngraphs, 1), mar=mr, no.readonly=TRUE)
	}  
#     tf <- tframe(tfwindow(x, start=start, end=end))
# would be nice if this could expand tf (tfwindow only truncates - need a
# replacement that expands too.)
     N <- nseries(x)
     if(length(xlab)          < N) xlab          <- rep(xlab, N)
     if(length(subtitle)      < N) subtitle      <- rep(subtitle, N)
     if(length(footnoteRight) < N) footnoteRight <- rep(footnoteRight, N)
     if(length(footnoteLeft)  < N) footnoteLeft  <- rep(footnoteLeft, N)
     if(length(source)        < N) source        <- rep(source, N)
     if(length(legend)        < N) legend        <- rep(legend, N)
     if(length(legend.loc)    < N) legend.loc    <- rep(legend.loc, N)
     if(length(YaxisR)        < N) YaxisR        <- rep(YaxisR, N)
 
     for (i in series)
       {if(mode(i)=="character") i <- match(i, names)
	z <-  selectSeries(x, series=i)
        for (d in list(...))
    	   z <- tbind(z, selectSeries(d, series=i)) 
	lgd <- if (is.matrix(legend))legend[,i] else legend
	tfOnePlot(z, tf=tf, start=start, end=end, subtitle=subtitle[i],
	          lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
		  xlab=xlab[i], ylab=ylab[i], xlim=xlim[[i]], ylim=ylim[[i]],
		  Xaxis=Xaxis, L1=L1,
		  YaxisL=YaxisL, YaxisR=YaxisR[i], Yaxis.lab.rot=Yaxis.lab.rot,
		  splitPane=splitPane,lastObs=lastObs, source=source[i],
		  footnoteLeft=footnoteLeft[i], footnoteRight=footnoteRight[i],
		  legend=lgd, legend.loc=legend.loc[i])
        if(!is.null(title) && (i==1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = title)	
    	}
    
  invisible()
 }

tfOnePlot <- function(x, tf=tframe(x), start=tfstart(tf), end=tfend(tf), 
         Title=NULL, title=Title, subtitle=NULL, lty=1:5, lwd=1, pch=1, col=1:6, cex=NULL,
        xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, par=NULL,
	Xaxis="auto", L1=NULL,
	YaxisL=TRUE, YaxisR=FALSE, Yaxis.lab.rot = "vertical",
	splitPane=NULL,
	lastObs=FALSE,  
	source=NULL,
	footnote=NULL, footnoteLeft=footnote, footnoteRight=NULL,
	legend=NULL, legend.loc="topleft"){
  if (inherits(x, "zooreg")) x <- as.ts(x)
  old.par <- par(par)
  on.exit(par(old.par)) 
  mr <- par()$mar

  if (!is.tframed(x)) plot(x, start=start, end=end, 
              lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
              xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  else {
      #if(is.null(source)) source <- 
     #        if(is.null(options("TSsource"))) NULL else(options()$TSsource)(x)
     if (!is.null(start)) x <- tfwindow(x, start=start, warn=FALSE)
     if (!is.null(end))   x <- tfwindow(x, end=end, warn=FALSE)
     if(is.null(xlab)) xlab <- ""
     if(is.null(ylab)) ylab <- paste(seriesNames(x), collapse="  ")
     if(is.null(ylim)) ylim <- range(x, na.rm=TRUE)
     tline <- time(x)
     if( inherits(tline, "ts")) tline <- unclass(tline)
     # formerly matplot with tline not a matrix was used, but this does
     # not plot (non-ts) dates as well as plot().
     if (lastObs) {
 	dt <- end(x)
	if(!is.null(frequency(x))) {
	  if(frequency(x) == 12)
	      dt <- paste(c("Jan", "Feb","Mar","Apr","May","Jun","Jul",
 	          "Aug","Sep","Oct","Nov","Dec")[dt[2]],dt[1], collapse=" ")
 	  else if(frequency(x) == 4)
	      dt <- paste(c("Q1","Q2","Q3","Q4")[dt[2]],dt[1], collapse=" ")
	}
 	last <- paste("Last observation:", dt)
       }
     # zoo freq==1 could be anything, so let zoo handle it
     noAuto <- inherits(x, "zoo") ||is.null(Xaxis) || !(frequency(x) %in% c(1,4,12))
     N <- nseries(x)
     if (1 == N) dim(x) <- c(length(x),1)
     else {
        if (length(lty) < N) lty <- rep(lty,length.out=N)
        if (length(lwd) < N) lwd <- rep(lwd,length.out=N)
        if (length(pch) < N) pch <- rep(pch,length.out=N)
        if (length(col) < N) col <- rep(col,length.out=N)
	}

     # add extra space for titles with a new line character
     m3 <- mr[3]
     if (is.character(title)    && grepl("\n", title))    m3 <- m3 + 1 
     if (is.character(subtitle) && grepl("\n", subtitle)) m3 <- m3 + 1 
     
     mr[3] <- m3
     
     if(is.null(splitPane)){
        par(mar=mr)
        plot(tline, x[,1], type="l", lty=lty, lwd=lwd, pch=pch, col=col, 
	    cex=cex, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
   	    xaxt = if(noAuto) "s" else "n", yaxt = "n")
   	
   	if (2 <= N) for (i in 2:N) lines(tline, x[,i],
   	  type="l", lty=lty[i], lwd=lwd[i], pch=pch[i], col=col[i])

   	#if(noAuto) axis(side=1) in some cases this does not seem to work
	#as well as specifying plot(xaxt = "s" )
	
   	if(!noAuto) 
	  if("auto" == Xaxis) tfXaxis(tline, L1=L1)
   	  else stop("Xaxis specification invalid.")

   	tfYaxis(YaxisL=YaxisL, YaxisR=YaxisR, Yaxis.lab.rot=Yaxis.lab.rot)
	}
     else { # splitPane
	mrL <- mrR <- mr
	mrL[4] <- 0
	mrR[2] <- 0
	if (YaxisR) mrR[4] <- mrR[4] + 1  # add for right labels (default is 2.1 )
	#mn <- min(x)
    	#mn <- mn - 0.01 * abs(mn)
    	#mx <- max(x)
    	#mx <- mx  + 0.01 * abs(mx)
    	#left side, screen(1)
    	par(fig=c(0, .65, 0,1), mar=mrL)
   	plot(tline, x[,1], type="l", lty=lty, lwd=lwd, pch=pch, col=col, 
	    cex=cex, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, 
   	    xaxt = if(noAuto) "s" else "n", yaxt = "n")
   	
   	if (2 <= N) for (i in 2:N) lines(tline, x[,i],
   	  type="l", lty=lty[i], lwd=lwd[i], pch=pch[i], col=col[i])

   	if(!noAuto) 
	  if("auto" == Xaxis) tfXaxis(tline, L1=L1)
   	  else stop("Xaxis specification invalid.")

   	tfYaxis(YaxisL=YaxisL, YaxisR=FALSE, Yaxis.lab.rot=Yaxis.lab.rot)
 
   	#right side, screen(2)
   	b  <-  tfwindow(x, start=tline[length(tline)] -(splitPane-1)/frequency(x))
   	bt <- time(b)
   	if( inherits(bt, "ts")) bt <- unclass(bt)
	par(fig=c(.7, 1, 0,1), new=TRUE, mar=mrR)
	plot(bt, b[,1], type="l", lty=lty, lwd=lwd, pch=pch, col=col, 
	    cex=cex, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, 
   	    xaxt = if(noAuto) "s" else "n", yaxt = "n")
   	
   	if (2 <= N) for (i in 2:N) lines(bt, b[,i],
   	  type="l", lty=lty[i], lwd=lwd[i], pch=pch[i], col=col[i])

   	if(!noAuto) 
	  if("auto" == Xaxis) tfXaxis(bt, L1=L1)
   	  else stop("Xaxis specification invalid.")

   	tfYaxis(YaxisL=FALSE, YaxisR=YaxisR, Yaxis.lab.rot=Yaxis.lab.rot)
	# Now set back to full device for title and footnotes.
	# setting usr works around what appears to be  bug. (The footnotes
	# do not get set properly relative to the center.)
	par(fig=c(0, 1, 0,1), new=FALSE, mar=mr, usr=c(0,1,0,1))
	}
     }
  if (!is.null(title) && (is.null(options()$PlotTitles) ||
      options()$PlotTitles)) title(main = title)      
  if (!is.null(subtitle) && (is.null(options()$PlotSubtitles) ||
      options()$PlotSubtitles)) title(main = subtitle, line=0.5, 
        cex.main=0.8 *par("cex.main"), font.main=0.5 *par("font.main"))       
  if (!is.null(source) && (is.null(options()$PlotSources) ||
      options()$PlotSourcse)) 
               mtext(source, side=1, line = 2, adj=0, cex=0.7)    
  if (lastObs) mtext(last,   side=1, line = 2, adj=1, cex=0.7)
   # footnote will go on another line with \n
  if (!is.null(footnoteLeft) && (is.null(options()$PlotFootnotes) ||
      options()$PlotFootnotes)) 
           mtext(footnoteLeft, side=1, line = 3, adj=0, cex=0.7)      
  if (!is.null(footnoteRight) && (is.null(options()$PlotFootnotes) ||
      options()$PlotFootnotes)) 
           mtext(footnoteRight, side=1, line = 3, adj=1, cex=0.7)     
  if (!is.null(legend)) legend(legend.loc, inset = c(0.05, .05), 
     col=col, lty=lty, cex=0.7, legend=legend)
  invisible(x)
 }


tfYaxis <- function (YaxisL=TRUE, YaxisR=FALSE, Yaxis.lab.rot = "vertical" ) {
   if (Yaxis.lab.rot == "vertical" )   las <- 3
   if (Yaxis.lab.rot == "horizontal" ) las <- 1
   ysc <- par()$yaxp
   if(ysc[3] < 10) blank <- FALSE else TRUE #unlabelled intermediate tick marks
   if (! blank) {
      # calculate tick marks in between labeled ones. If the labels are at intervals of 
      #   5 ..9 this should give 1. At 0.5 ... 0.9 this should give 0.1
      by <- 10^(ceiling(log10(((ysc[2] - ysc[1]) / ysc[3]))-1))
      at1 <- seq(ysc[1], ysc[2], by = by)
      }
   if (YaxisL) {
      axis(side = 2, las=las)
      if (! blank) axis(side = 2, at= at1, tcl=-0.3, labels=FALSE)
      }
   
   if ((is.numeric(YaxisR)) || YaxisR) { # numeric or T
      scR <- axis(side = 4, las=las, labels=!is.numeric(YaxisR))
      if (is.numeric(YaxisR))
        axis(side = 4, las=las, at=scR, labels=as.graphicsAnnot(YaxisR * scR))
      if (! blank) axis(side = 4, at= at1, tcl=-0.3, labels=FALSE)
      }
   invisible()
   }


tfXaxis <- function (x, L1 = NULL) {
   at1 <- time(x)
   fr <- frequency(x)
   mgp1 <- 0.5  # ignored with lab1 FALSE
   mgp2 <- 1.2  # ignored with lab1 FALSE
   
   # note that fr not in these cases is usually given to UseMethod in tfOnePlot
   # and should never call tfXaxis.
   if (is.null(fr))  stop("frequency must not be null for auto axis calculation.")
   else if (fr == 1 ) {
     blank <- TRUE 
     at2 <- at1
     }
   else if (fr == 4 ) {
     if(is.null(L1)) L1 <- c("Q1","Q2","Q3","Q4")
     at2 <- unique(floor(at1))
     blank <- length(at2) > 8 
     }
   else if (fr == 12) {
     if(is.null(L1)) L1 <- c("J","F","M","A","M","J","J","A","S","O","N","D")
     at2 <- unique(floor(at1))
     blank <- length(at2) > 4 
     }
   else stop("frequency ", fr, "is not handled by auto axis calculation.")

   # omit period labels when it gets too crowded    
   if (blank) lab1 <- FALSE # omit period labels
   else {
     lab1 <- rep(L1, length.out=fr * length(at2))
     # arrange for not starting in first period of year
     s <- fr - sum(at1 < (1 + at2[1] - 1e-5)) # less eps because date may fail < for Jan
     if (0 < s)  lab1 <- lab1[-seqN(s)]
     lab1 <- lab1[seqN(Tobs(x))]
     }
   axis(side = 1, at = at1,   cex.axis=0.6, labels = lab1, mgp = c(3, mgp1, 0))
   axis(side = 1, at = at2,	  tcl=-0.8, labels = FALSE)
   axis(side = 1, at = 0.5 + at2, tcl=0,    labels=at2, mgp = c(3, mgp2, 0))
   invisible(x)
}


diffLog <- function(obj, lag = 1, base = exp(1),
              names=paste("diff of log of ", seriesNames(obj))) 
   UseMethod("diffLog")
 
diffLog.default <- function(obj, lag = 1, base = exp(1),
              names=paste("diff of log of ", seriesNames(obj)))
{#Calculate the difference from lag periods prior for log of data.
 obj <- diff(log(obj, base = base), lag = lag)
 if(is.null(options()$ModSeriesNames) || options()$ModSeriesNames)
        seriesNames(obj) <- names
 obj
}


# Note 1. This is generic so methods can be defined on series within an object,
# as in TSdata and TSestModel. For different types of series, eg zoo, it
# should not be necessary to define methods for this, that should be done with
# lower level utilities like tfL and diff, which should be used here.
ytoypc <- function(obj, names=paste("y to y %ch", seriesNames(obj))) 
   UseMethod("ytoypc")
 
ytoypc.default <- function (obj, names=paste("y to y %ch", seriesNames(obj)) ){
   obj <- percentChange(obj, lag = tffrequency(obj))
   if(is.null(options()$ModSeriesNames) || options()$ModSeriesNames)
        seriesNames(obj) <- names
   obj
}


# See Note 1 above
percentChange <- function(obj, ...) UseMethod("percentChange")

percentChange.default <- function(obj, base=NULL, lag=1, 
      cumulate=FALSE, e=FALSE, ...)
{#  (... further arguments, currently disregarded)
   cls <- class(obj)
   # note next has to be applied to a shorter object in the end
   #if (is.tframed(obj)) tf <- list(end=tfend(obj), frequency=tffrequency(obj))
   #else tf <- NULL
   if (!is.tframed(obj)) stop("percentChange only works on tframed objects.")
   tf <- tframe(diff(obj, lag=lag)) # not very efficient
   if (is.null(dim(obj)))
     {vec <- TRUE
      obj <- matrix(obj, length(obj),1)
     }
   else vec <- FALSE
   mm <- unclass(rbind(base,obj))
   if (any(cumulate))
          mm[,cumulate] <-apply(mm[,cumulate,drop=FALSE],2,cumsum)
   if (any(e)) mm[,e] <- exp(mm[,e,drop=FALSE])
   N <- NROW(mm)
   pchange <-100*(mm[(lag+1):N,,drop=FALSE] - 
                    mm[1:(N-lag),,drop=FALSE])/mm[1:(N-lag),,drop=FALSE]
   if (vec) pchange <- pchange[,1]
   tframed(pchange, tf) 
}

# See Note 1 above
annualizedGrowth <- function(obj, ...) UseMethod("annualizedGrowth")

annualizedGrowth.default <- function(obj, lag=1, freqLagRatio=frequency(obj)/lag,
        names=paste("Annual Growth of", seriesNames(obj)), ...) {
  if (!is.tframed(obj)) stop("annualizedGrowth only works on tframed objects.")
  d <- tfL(obj, p= lag)
  r <- 100*((obj / d )^freqLagRatio - 1)
  if(is.null(options()$ModSeriesNames) || options()$ModSeriesNames)
        seriesNames(r) <- names
  tframed(r, tframe(diff(obj, lag=lag))) # not very efficient
  }

addDate <- function(date, periods, freq)
  {if (is.null(periods)) periods <- 0
   c(date[1]+(date[2]+periods-1)%/%freq, 1+(date[2]+periods-1)%%freq)
  }


tsScan <- function(file="", skip=1, nseries=1, sep=",", 
           na.strings=c("NA", "NC", "ND"), ...)
   {# all args passed to scan. Expects a file with (default one) title line 
    # to skip and data in three columns (default separated with commas):
    #   year, period, data[;1], data[;2], ..., data[;nseries]
    # and builds a ts with freq set to max(period)
    z <- scan(file=file, skip=skip, what=as.list(c(seq(2), double(nseries))),
              sep=sep, na.strings=na.strings, ...)
    zz <- NULL
     for  (i in 1:nseries) zz <- cbind(zz, z[[2+i]])
    ts(zz, start=c(z[[1]][1], z[[2]][1]), frequency=max(z[[2]]))
    }

tsWrite <- function(x, file="data", header=TRUE, sep=",", digits=16)
   {# all args passed to scan. Expects a ts or mts.
    # write file with (default one) title line. 
    # then data in three columns (efault separated with commas):
    #   year, period, data[;1], data[;2], ...
    if (header) write(paste("year", "period", 
                   paste(seriesNames(x), collapse=sep),sep=sep), file=file)
    yr  <- floor(time(x))
    pr  <- 1+ (time(x) %% 1) * frequency(x)
    dg <- options(digits=digits)
    on.exit(options(dg))
    x <- as.matrix(x)
    write(t(cbind(yr, pr, x)), file=file, ncolumns = 2 + ncol(x), sep=sep, append=header)
    }


tfVisPlot <- function (x, tf = tframe(x), start = tfstart(tf), end = tfend(tf), 
                  options=list(title=NULL), ...){   
    # x a multivariate series. ... passed to gvisLineChart
    if (!is.null(start)) x <- tfwindow(x, start = start, warn = FALSE)
    if (!is.null(end))   x <- tfwindow(x, end = end, warn = FALSE)
    tm <- time(x)
    y <- floor(time(tm))
    if (frequency(x) == 4) {
    	p <- cycle(tm)
    	pp <- c("Q1", "Q2", "Q3", "Q4")[p]
    	dt <- paste(pp, y)
    	seriesData <- data.frame(date = dt, x)
    }
    else if (frequency(x) == 12) {
    	p <- cycle(tm)
    	pp <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
    	    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[p]
    	dt <- paste(pp, y)
    	seriesData <- data.frame(date = dt, x)
    }
    else seriesData <- data.frame(date = tm, x)
    
    nm <- seriesNames(x)
    names(seriesData) <- c("date", nm)
    
    if (!requireNamespace("googleVis")) stop("tfVisPlot requires googleVis")
    else  plot(googleVis::gvisLineChart(seriesData, xvar="date", yvar=nm,
         options=options, ...))
    cat("look for chart in web browser.\n")
    invisible(x)
}

