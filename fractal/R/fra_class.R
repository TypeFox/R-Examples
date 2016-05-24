################################################
## FRACTAL nonlinear dynamics constructor
## functions and corresponding methods
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class chaoticInvariant
## Constructor function: chaoticInvariant
## Methods:
##
##   eda.plot.chaoticInvariant
##   plot.chaoticInvariant
##   print.chaoticInvariant
##   print.summary.chaoticInvariant
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class FDSimulate
## Constructor function: FDSimulate
## Methods:
##
##   plot.FDSimulate
##   print.FDSimulate
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class fractalBlock
## Constructor function: fractalBlock
## Methods:
##
##   eda.plot.fractalBlock
##   plot.fractalBlock
##   print.fractalBlock
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class FNN
## Constructor function: FNN
## Methods:
##
##   plot.FNN
##   print.FNN
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class embedSeries
## Constructor function: embedSeries
## Methods:
##
##   [.embedSeries
##   as.matrix.embedSeries
##   plot.embedSeries
##   print.embedSeries
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class lyapunov
## Constructor function: lyapunov
## Methods:
##
##   plot.lyapunov
##   print.lyapunov
##   print.summary.lyapunov
##   summary.lyapunov
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class spaceTime
## Constructor function: spaceTime
## Methods:
##
##   as.matrix.spaceTime
##   eda.plot.spaceTime
##   plot.spaceTime
##   print.spaceTime
##
################################################

################################################
##
## Class chaoticInvariant
## Constructor function: chaoticInvariant
## Methods:
##
##   eda.plot.chaoticInvariant
##   plot.chaoticInvariant
##   print.chaoticInvariant
##   print.summary.chaoticInvariant
##
################################################

###
# Constructor: chaoticInvariant
###

"chaoticInvariant" <- function(x, dimension=NA, n.embed=NA, epsilon=NA,
  tlag=NA, olag=NA, resolution=NA, n.reference=NA, metric=NA,
  series.name="", xlab ="", ylab= "", invariant="", n.neighbor=NA,
  fit=lmsreg, series=NA)
{
  # check inputs
  if (!is.list(x) || (length(x) != 2))
    stop("Required input must be a list containing two matrices")
  if (!is.matrix(x[[1]]))
    stop("First object of required input must be a matrix")
  if (!is.integer(dimension))
    stop("The dimension argument must be a integer scalar or a vector of integers")
  checkScalarType(n.embed,"integer")
  checkScalarType(n.reference,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(olag,"integer")
  checkScalarType(series.name,"character")
  if (!is.vector(series))
    stop("The series argument must be a vector")
  checkScalarType(metric,"numeric")
  checkScalarType(xlab,"character")
  checkScalarType(ylab,"character")
  checkScalarType(invariant,"character")
  if (!is.na(resolution))
    checkScalarType(resolution,"numeric")
  if (!is.na(n.neighbor))
    checkScalarType(n.neighbor,"integer")
  if (!any(is.na(epsilon)) && !is.vector(epsilon))
    stop("The epsilon argument must be an vector without any NAs")
  if (ncol(x[[1]]) != length(dimension))
    stop("The dimension argument does not match the number of columns ",
      "in the supplied data matrix")
  if (nrow(x[[1]]) != length(x[[2]]))
    stop("The length of the scale vector must match the number of rows ","
      in the supplied data matrix")

  # initialize variables
  stat  <- x[[1]]
  scale <- as.vector(x[[2]])

  # calculate various statistics on the data
  fits <- apply(stat, MARGIN=2, function(x, scale, fit){
    linearFit(scale, x, fit=fit, angle.tolerance=5, aspect=FALSE, method="widest")},
	scale=scale, fit=fit)

  # extract slope coefficients
  # and form difference statistics
  slopes <- unlist(lapply(fits,
    function(x){as.numeric(coefficients(x)[2])}))

  # test for information dimension analysis
  is.d1 <- is.element(invariant, "information dimension")

  dlag   <- 2
  dstat  <- apply(stat, MARGIN=2,
    function(x, scale, dlag, is.d1){

      z <- rep(NA, length(x) - dlag)
      good <- which(!is.na(x))

      fit <- ifelse1(is.d1,
      	loess.smooth(x[good], scale[good], evaluation=length(x[good]), degree=2),
        loess.smooth(scale[good], x[good], evaluation=length(x[good]), degree=2))

      # extend the data as needed to maintain same number
      # of points in slope estimate
      ix   <- seq(min(10, length(x)))
      xfit <- fit$x[ix]
      yfit <- fit$y[ix]

      # regress the data with a second order polynomial
      # don't use nls() here because it chokes when data
      # is very close to linear
      ext.fit <- lm(y ~ x + x^2, data=data.frame(x=xfit, y=yfit))

      newx <- sort(min(xfit) - mean(diff(xfit)) * seq(dlag))
      newy <- predict(ext.fit, newdata=list(x=newx))


      slope <- diff(c(newy,fit$y),lag=dlag) / diff(c(newx,fit$x),lag=dlag)
      z[good] <- slope
      z
   	}, scale=scale, dlag=dlag, is.d1=is.d1)

  dstat[dstat == Inf] <- NA

  # adjust the slopes for the information dimension
  if (is.d1)
    slopes <- 1 / slopes

  names(slopes) <- dimension

  z <- list(stat=stat, scale=scale, dstat=dstat)
  oldClass(z) <- "chaoticInvariant"

  attr(z, "dimension")   <- dimension
  attr(z, "epsilon")     <- epsilon
  attr(z, "n.embed")     <- n.embed
  attr(z, "n.reference") <- n.reference
  attr(z, "n.neighbor")  <- n.neighbor
  attr(z, "tlag")        <- tlag
  attr(z, "olag")        <- olag
  attr(z, "resolution")  <- resolution
  attr(z, "series.name") <- series.name
  attr(z, "series")      <- series[seq(min(c(2048,length(series))))]
  attr(z, "xlab")        <- xlab
  attr(z, "ylab")        <- ylab
  attr(z, "metric")      <- metric
  attr(z, "invariant")   <- invariant
  attr(z, "fit")         <- fits
  attr(z, "slope")       <- slopes

  z
}

###
# eda.plot.chaoticInvariant
###

"eda.plot.chaoticInvariant" <- function(x, type="o", lty=1, cex=1, main.cex=0.7 * cex, gap=0.14, adj=1, main.line=0.5, ...)
{
  old.plt <- splitplot(2,2,1,gap=gap)
  on.exit(par(old.plt))

  xatt <- attributes(x)

  # Time History
  xlab  <- "Time"
  xdata <- xatt$series
  ylab  <- xatt$series.name

  plot(c(1,length(xdata)), range(xdata,na.rm=TRUE), type="n", xlab=xlab,
    ylab=ylab, cex=cex)
  lines(xdata, type="l", col=1, cex=cex)
  mtext("Time History", cex=main.cex, adj=adj, line=main.line)

  # 2-D Embedding
  splitplot(2,2,2,gap=gap)
  plot(embedSeries(xatt$series, dimension=2, tlag=xatt$tlag, series.name=ylab), add=TRUE, col=1)
  mtext("2D Embedding", cex=main.cex, adj=adj, line=main.line)

  # plot stat and slope curves
  splitplot(2,2,3,gap=gap)
  plot(x, type="stat", cex=cex, legend=FALSE, add=TRUE, adj=adj, main.line=main.line)
  splitplot(2,2,4,gap=gap)
  plot(x, type="slope", cex=cex, add=TRUE, adj=adj, main.line=main.line)

  invisible(NULL)
}

###
# plot.chaoticInvariant
###

"plot.chaoticInvariant" <- function(x, type="stat", lty=1,
  fit=TRUE, grid=FALSE, plot.type="b", legend=TRUE, add=FALSE,
  cex=1, main.cex=0.7*cex, main=NULL, adj=1, main.line=0.5, ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  xatt <- attributes(x)

  if (is.null(main))
    main <- properCase(xatt$invariant)

  type <- match.arg(type, c("stat","dstat","slope","entropy"))

  if (type == "stat"){
    xdata    <- x$scale
    ydata    <- x$stat
    ylab     <- xatt$ylab
    fit.data <- xatt$fit
    main     <- paste(main, "curves")
  }
  else if (type == "dstat"){
    xdata <- x$scale
    ydata <- x$dstat

    ylab  <- paste("derivative of log2 ", xatt$ylab, sep="")
    plot.type <- "o"
    fit   <- FALSE
  }
  else if (type == "slope"){

    xdata <- xatt$dimension
    ydata <- xatt$slope

    ylab <- ifelse1(is.element(xatt$invariant, "information dimension"),
      "d1(E)","d2(E)")

    plot(xdata, ydata, type=plot.type, xlab="Embedding Dimension, E",
      ylab=ylab, cex=cex, ylim=c(0,max(ydata)))

    fit    <- FALSE
    legend <- FALSE
  }
  else if (type == "entropy"){

    if (!is.element(xatt$invariant, c("correlation dimension","information dimension")))
      stop("Entropy plots are only available for correlation and information dimension analyses")
    is.d1 <- is.element(xatt$invariant, "information dimension")

    xdata <- x$scale
    ydata <- x$stat

    # create boxplots
    ydata2 <- abs(diff(t(ydata)))

    if (is.d1)
      ydata2 <- ydata2 * t(ydata[,seq(numCols(ydata)-1)])
    ydata2 <- data.frame(ydata2)

    names(ydata2) <- round(xdata,2)
    p <- boxplot(ydata2, plot=FALSE)

    # try to isolate only the most populous median entropy values
    # by forming a NCLASS bin histogram and limiting the
    # data to only those values that correspond to the
    # dominate mode
    p.median <- p$stats[3,]
    phist    <- hist(p.median, nclass=5, plot=FALSE)
    bestbar  <- rev(order(phist$counts))[1]
    best.entropy.range <- phist$breaks[bestbar + (0:1)]
    best.scales <- which(p.median >= min(best.entropy.range) &
      p.median <= max(best.entropy.range) & p.median > 0)
    outliers <- which(is.element(p$group, best.scales))

    # now that we have obtained the scale over which the data exhibits
    # the most populous median entropy values, subsample the boxplot
    # statistics accordingly
    p$stats <- p$stats[,best.scales]
    p$n     <- p$n[best.scales]
    p$conf  <- p$conf[,best.scales]
    p$names <- p$names[best.scales]
    p$out   <- p$out[outliers]
    p$group <- p$group[outliers]
    bxp(p,outline=FALSE, whisklty=1, ylim=range(p$stats), srt=90)

    # add median value of median of entropies vector (each element
    # of this vector corresponds to the median value at a particular
    # scale and over all dimensions)
    median.median.entropy <- median(p$stats[3,])
    abline(h=median.median.entropy, lty=4, lwd=3, xpd=FALSE)
    mtext(paste(round(median.median.entropy,3)), side=4, at=median.median.entropy,
      line=0.5, adj=0.5, col=2, srt=-90)

    # add labels
    mtext(xatt$xlab, side=1, line=3)
    mtext(paste(ifelse1(is.d1,"Kolmogorov-Sinai", "Correlation"), "Entropy"), side=2, line=3)

    return(invisible(median.median.entropy))
  }

  iy <- seq(numCols(ydata))

  if (is.element(type, c("stat","dstat")))
    matplot(x=xdata,y=ydata,pch=iy, col=iy, type=plot.type, lty=lty,xlab=xatt$xlab, ylab=ylab, cex=cex)

  if (fit)
   for (i in iy)
     abline(fit.data[[i]], lty=1, col=i, xpd=FALSE)

  if (grid)
    gridOverlay(lty=4)

  if (nchar(main))
    mtext(main, cex=main.cex, adj=adj, line=main.line)

  if (legend)
  {
  	legend("topleft",
          title="SLOPES",
	      legend = paste(paste("DIM", xatt$dimension),round(xatt$slope, 3), sep=": "),
	      lty=lty,
	      pch=iy,
	      col=iy,
	      text.col=iy)
   }
  invisible(NULL)
}

###
# print.chaoticInvariant
###

"print.chaoticInvariant" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)

  if (charmatch("maximal Lyapunov exponent", xatt$invariant, nomatch=FALSE))
    slope <- apply(matrix(xatt$slope, ncol=length(unique(xatt$dimension))), MARGIN=2, mean)
  else
    slope <- xatt$slope

  main <- paste(properCase(xatt$invariant), " for ", xatt$series.name, sep ="")

  z <- list(
    "Embedding points"=xatt$n.embed,
    "Embedding dimension(s)"=unique(xatt$dimension),
    "Epsilon(s)"=ifelse1(is.missing(xatt$epsilon), NULL, round(xatt$epsilon, 3)),
    "Time lag"=xatt$tlag,
    "Oribital lag"=xatt$olag,
    "Scale points/octave"=ifelse1(is.missing(xatt$n.resolution), NULL, xatt$resolution),
    "Reference points"=ifelse1(is.missing(xatt$n.reference), NULL, xatt$reference),
    "Distance metric"=paste("L-", xatt$metric, sep=""),
    "Invariant estimate(s)"=round(slope, 3))

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# print.summary.chaoticInvariant
###

"print.summary.chaoticInvariant" <- function(x, ...)
{
  NextMethod("print")
  invisible(x)
}

################################################
##
## Class fractalBlock
## Constructor function: fractalBlock
## Methods:
##
##   eda.plot.fractalBlock
##   plot.fractalBlock
##   print.fractalBlock
##
################################################

###
# Constructor: fractalBlock
###

"fractalBlock" <- function(domain, estimator, exponent, exponent.name, scale, stat, stat.name, detrend, overlap,
  data.name, sum.order, series, logfit, sdf=NULL)
{
  # check input argument class
  checkScalarType(domain,"character")
  checkScalarType(estimator,"character")
  checkScalarType(exponent.name,"character")
  checkScalarType(exponent,"numeric")
  if (!is.null(scale))
    checkVectorType(scale,"numeric")
  if (!is.null(stat))
    checkVectorType(stat,"numeric")
  checkScalarType(stat.name,"character")
  if (!is.null(detrend))
    checkScalarType(detrend,"character")
  if (!is.null(overlap))
    checkScalarType(overlap,"numeric")
  checkScalarType(data.name,"character")
  checkScalarType(sum.order,"integer")
  checkVectorType(series,"numeric")
  if (!is.null(logfit) && !is.element(class(logfit),c("lm","lms","lts")))
    stop("logfit must be a member of class \"lm\", \"lms\", or \"ltsreg\"")
  if (!is.null(overlap) && ((overlap < 0) | (overlap >= 1)))
    stop("Overlap factor must be in the range [0,1)")
  if (!is.null(sdf) && !is(sdf,"SDF"))
    stop("sdf must be an object of class \"SDF\"")
  if (any(scale < 0.0))
    stop("Negative scale(s) not allowed")
  domain <- match.arg(lowerCase(domain),c("time","frequency"))

  # pack primary objects into a list
  z <- exponent
  names(z) <- exponent.name
  oldClass(z) <- "fractalBlock"

  # assign attributes
  attr(z, "domain")        <- domain
  attr(z, "exponent.name") <- exponent.name
  attr(z, "scale")         <- scale
  attr(z, "scale.ratio")   <- ifelse1(length(scale) > 1, scale[2]/scale[1], NA)
  attr(z, "stat")          <- stat
  attr(z, "stat.name")     <- stat.name
  attr(z, "estimator")     <- estimator
  attr(z, "detrend")       <- detrend
  attr(z, "overlap")       <- overlap
  attr(z, "n.sample")      <- length(series)
  attr(z, "data.name")     <- data.name
  attr(z, "sum.order")     <- sum.order
  attr(z, "series")        <- series
  attr(z, "logfit")        <- logfit
  attr(z, "sdf")           <- sdf

  return(z)
}

###
# eda.plot.fractalBlock
###

"eda.plot.fractalBlock" <- function(x, cex=1, cex.main=1, col=2, adj.main=1, ...)
{
  old.plt <- splitplot(2,1,1)
  on.exit(par(old.plt))

  xatt <- attributes(x)

  direction <- ifelse1(xatt$sum.order > 0, "cumulative summation",
    xatt$sum.order < 0, "difference", "")

  tit.add <- ifelse1(xatt$sum.order != 0,
    paste(":", ordinal(abs(xatt$sum.order)), "order", direction), "")

  plot(x=time(xatt$series), y=xatt$series, col=col,
    xlab="TIME", ylab=xatt$data.name, cex=cex)
  title(paste("Time History", tit.add,sep=""), cex=cex.main, adj=adj.main)

  splitplot(2,1,2)
  plot(x, add=TRUE)

  tit.add <- ifelse1(!is.null(xatt$overlap) && xatt$overlap > 0,
    paste(": ", round(xatt$overlap * 100, 2), "% overlap", sep="" ), "")

  detrend.str <- ifelse1(is.null(xatt$detrend),"", paste(", Detrending: ", xatt$detrend, sep=""))
  mtext(paste(xatt$estimator, tit.add, detrend.str, sep=""), cex=cex.main, adj=adj.main, line=0.5)

  invisible(NULL)
}

###
# plot.fractalBlock
###

"plot.fractalBlock" <- function(x, pch=18, col=c("black","red"), lty=c(1,1),
   grid=list(lty=2, col=16, density=3), key=TRUE, add=FALSE, cex=1, ...)
{

  # define local functions
  rescale <- function(x, xmin, xmax){

    dx <- xmax-xmin
    x  <- x - min(x)
    x / max(x) * dx + xmin
  }

  lty  <- rep(lty,2)[1:2]
  col  <- rep(col,2)[1:2]
  xatt <- attributes(x)

  if (is.null(xatt$scale) || is.null(xatt$stat)){
  	cat("No data to plot\n")
  	return(invisible(NULL))
  }

  # for spectral regression data, the log has already
  # been taken, so don't do so below
  if (is.element(xatt$domain,"time")){
    xdata <- log(xatt$scale)
    ydata <- log(xatt$stat)
  }
  else if (is.element(xatt$domain,"frequency")){
    xdata <- xatt$scale
    ydata <- xatt$stat
  }

  fit   <- attr(x,"logfit")
  exponent <- as.numeric(x)

  plot(x=xdata,y=ydata, type="b", lty=lty[1], col=col[1], pch=pch,
    xlab="log(scale)", ylab=paste("log(",xatt$stat.name,")",sep=""), cex=cex)
  lines(x=xdata,y=fit$fitted.values, type="l", lty=lty[2], col=col[2])

  # add key if requested
  if (key){

  		legend("topright",
  		  legend=c(xatt$stat.name, paste(xatt$exponent.name,"=",round(exponent,3))),
  		  lty=lty, col=col, pch=c(18,-1))

  }

  invisible(NULL)
}

###
# print.fractalBlock
###

"print.fractalBlock" <- function(x,  justify="left", sep=":", n.digits=5, ...)
{
  xatt <- attributes(x)
  main <- paste(properCase(xatt$estimator), " for ", xatt$data.name, sep ="")
  direction <- ifelse1(xatt$sum.order > 0, "cumulative summation",
    xatt$sum.order < 0, "difference", "")

  z <- list(
    as.numeric(x),
    "Domain"=properCase(xatt$domain),
    "Statistic"=xatt$stat.name,
    "Length of series"=xatt$n.sample,
    "Block detrending model"=xatt$detrend,
    "Block overlap fraction"=xatt$overlap,
    "Scale ratio"=ifelse1(is.element(xatt$domain,"time"),xatt$scale.ratio, NULL),
    "Preprocessing"=ifelse1(xatt$sum.order != 0, paste(ordinal(abs(xatt$sum.order)), "order", direction), NULL))
  names(z)[1] <- paste(xatt$exponent.name, "estimate")

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  if (is.element(xatt$domain,"time")){
  	mat <- rbind(xatt$scale,xatt$stat)
    dimnames(mat) <- list(c("Scale",xatt$stat.name), rep("",length(xatt$scale)))
    print(mat,digits=n.digits)
  }
  else if (is.element(xatt$domain,"frequency") && !is.null(xatt$sdf)){
  	cat("\n")
    print(xatt$sdf)
  }

  invisible(x)
}

################################################
##
## Class FNN
## Constructor function: FNN
## Methods:
##
##   plot.FNN
##   print.FNN
##
################################################

###
# plot.FNN
###

"plot.FNN" <- function(x, add=FALSE, xlab="Embedding Dimension", ylab="FNN %", acol="blue", rcol="red", ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  xatt <- attributes(x)
  E    <- xatt$dimension

  plot(c(0,E),c(0,100),type="n", ylab="",xlab="", ...)
  lines(x[1,], type="b", pch="r", col=rcol)
  lines(x[2,], type="b", pch="a", col=acol)
  gridOverlay()
  title(xlab=xlab, ylab=ylab)
  text(rep(E-0.3,2),c(97,92),labels=c("r: rtol", "a: atol"),col=c(rcol,acol), adj=1)
}

###
# print.FNN
###

"print.FNN" <- function(x, justify="left", sep=":", digits=3, ...)
{
  xatt <- attributes(x)
  main <- paste("False Nearest Neighbors for", xatt$data.name)

  z <- list(
    "Series points"=xatt$n.sample,
    "Embedding dimension(s)"=seq(xatt$dimension),
    "Time lag"=xatt$tlag,
    "Oribital lag"=xatt$olag,
    "Neighbor tolerance (rtol)"=xatt$rtol,
    "Attractor tolerance (atol)"=xatt$atol)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)
  cat("Test results (%):\n")
  print(oldUnclass(x)[,], digits=digits)

  invisible(x)
}

################################################
##
## Class embedSeries
## Constructor function: embedSeries
## Methods:
##
##   [.embedSeries
##   as.matrix.embedSeries
##   plot.embedSeries
##   print.embedSeries
##
################################################

###
# [.embedSeries
###

#setMethod("[", "embedSeries",
#function(x, i, j, ..., drop = TRUE)
#{
#  if(missing(drop))
#    "[.embedSeries"(x, i, j, ...)
#  else
#    "[.embedSeries"(x, i, j, ..., drop = drop)
#})


"[.embedSeries" <- function(x, i, j, ..., drop=FALSE)
{
	narg <- nargs() - !missing(drop)
  if (narg < 3)
    stop("Incorrect indexing format. Should be of the form x[i,j]")
	if (missing(i))
	  i <- seq(nrow(x))
	if (missing(j))
	  j <- seq(ncol(x))

  oldUnclass(x)[i, j, drop=drop]
}

###
# as.matrix.embedSeries
###

"as.matrix.embedSeries" <- function(x, ...) x[,]

###
# plot.embedSeries
###

"plot.embedSeries" <- function(x, tlag=NULL, dimension=NULL,
			series.name=NULL, pch=".", add=FALSE, cex=1, col="black", ...)
{
  "embedPlot" <- function(x, tlag=NULL, dimension=NULL,
    series.name=NULL, pch=".",  cex=1, col="black", ...)
  {
    if (is.null(series.name))
      series.name <- deparseText(substitute(x))

    if (is.matrix(x)){
      z <- x
      if (is.null(dimension))
        dimension <- ncol(z)
    }
    else{
      z <- embedSeries(x, dimension=dimension, tlag=tlag, series.name=series.name)
      dimension <- ncol(z)
    }

    if (dimension == 1){
      plot(z[,1], rep(0,length(z[,1])), xlab=dimnames(z)[[2]][1], ylab="",
        pch=pch, cex=cex, col=col, ...)
    }
    else if (dimension == 2){
      xlab  <- dimnames(z)[[2]][1]
      ylab  <- dimnames(z)[[2]][2]
      xdata <- z[,1]
      ydata <- z[,2]

      plot(range(xdata), range(ydata), xlab=xlab, ylab=ylab, type="n")
      points(xdata, ydata, pch=pch, cex=cex, col=col, ...)
    }
    else
      invisible(scatterplot3d::scatterplot3d(as.matrix(z), ...))

    invisible(NULL)
  }

  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  embedPlot(x, tlag=tlag, dimension=dimension,
    series.name=series.name, pch=pch, cex=cex, col=col, ...)
  invisible(NULL)
}

###
# print.embedSeries
###

"print.embedSeries" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)
  main <- paste("Embedding matrix for ", xatt$series.name, sep ="")

  z <- list(
    "Number of points"=xatt$n.embed,
    "Embedding dimension(s)"=xatt$emb.dim,
    "Epsilon(s)"=ifelse1(is.null(xatt$epsilon), NULL, round(xatt$epsilon, 3)),
    "Time lag"=xatt$tlag)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

################################################
##
## Class lyapunov
## Constructor function: lyapunov
## Methods:
##
##   plot.lyapunov
##   print.lyapunov
##   print.summary.lyapunov
##   summary.lyapunov
##
################################################

###
# plot.lyapunov
###

"plot.lyapunov" <- function(x, add=FALSE, grid=TRUE, cex=1,
  xlab="SCALE",ylab="Local Lyapunov Exponents", whisklty=1, col=2:8, range=0,
  boxwex=0.2, ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  n.scale  <- numCols(x[[1]])
  scales   <- names(x[[1]])
  ix       <- seq(length(scales))
  boxwidth <- 1.5 / n.scale
  ylim     <- range(unlist(x))

  plot(c(0, n.scale + 1), ylim, xlab=xlab, ylab=ylab,
    cex=cex, type="n", axes=FALSE)

  box()
  axis(2, cex=cex)

  if (grid)
    gridOverlay(density=5)

  x.mean <- summary(x)$mean

  if (length(col) < length(x))
    col <- c(col, rep(col[length(col)], length(x)-length(col)))

  for (i in seq(along=x)){
    par(new=TRUE)
    p <- boxplot(x[[i]], ylim=ylim, boxcol=col[i],
      range=range, whisklty=whisklty, boxwex=boxwex, cex=cex, ..., plot=TRUE)
    lines(seq(length(x[[i]])), x.mean[[i]], col=col[i], lwd=2)
  }

  invisible(NULL)
}

###
# print.lyapunov
###

"print.lyapunov" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)
  main <- paste("Local Lyapunov Spectrum for", xatt$data.name)

  z <- list(
    "Series points"=xatt$n.sample,
    "Sampling interval"=xatt$sampling.interval,
    "Embedding dimension"=xatt$dimension,
    "Local dimension"=xatt$local.dimension,
    "Time lag"=xatt$tlag,
    "Orbital lag"=xatt$olag,
    "Reference point indices"=xatt$reference,
    "Jacobian, neighborhood size"=xatt$n.reference,
    "Jacobian, distance metric"=paste("L-", xatt$metric, sep=""),
    "Jacobian, polynomial order"=xatt$polynomial.order,
    "Scales"=xatt$scales)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)
  invisible(z)
}

###
# print.summary.lyapunov
###

"print.summary.lyapunov" <- function(x, ...)
{
  NextMethod("print")
  invisible(x)
}

###
# summary.lyapunov
###

"summary.lyapunov" <- function(object, ...){

  x.mean <- data.frame(lapply(object, function(x) colMeans(as.matrix(x))))
  x.var  <- data.frame(lapply(object, colVars))
  x.median <- data.frame(lapply(object, colMedians))

  nms <- 1:length(object)

  names(x.mean) <- nms
  names(x.var)  <- nms
  names(x.median) <- nms

  z <- list(mean=x.mean, var=x.var, median=x.median)
  # oldClass(z) <- c("summary.lyapunov","list")
  oldClass(z) <- c("summary.lyapunov")

  z
}

################################################
##
## Class spaceTime
## Constructor function: spaceTime
## Methods:
##
##   as.matrix
##   eda.plot.spaceTime
##   plot.spaceTime
##   print.spaceTime
##
################################################

###
# as.matrix.spaceTime
###

"as.matrix.spaceTime" <- function(x, mode="any", names=TRUE, ...)
{
  checkScalarType(mode,"character")
  checkScalarType(names,"logical")
  y <- x
  attributes(y) <- NULL
  dims <- dim(x)
  nrow <- dims[1]
  ncol <- dims[2]
  dimnms <- dimnames(x)
  y <- matrix(as.vector(y, mode=mode), nrow=nrow, ncol=ncol)
  if(names)
    dimnames(y) <- dimnms
  y
}

###
# eda.plot.spaceTime
###

"eda.plot.spaceTime" <- function(x, add=FALSE, density=10, type="l", cex=1.2, ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  olags <- attr(x,"olags")

  tsp.median <- apply(x, MARGIN = 1, median)
  tsp.dens <- density(tsp.median)

  old.plt <- splitplot(2,1,1)
  on.exit(par(old.plt))
  plot(tsp.dens, xlab="Median Spatial Separation", ylab="Density", type=type, cex=cex, zero=FALSE, ...)

  splitplot(2,1,2)
  plot(olags, tsp.median, xlab="Orbital Lag", ylab="Median Spatial Separation", type=type, cex=cex, ...)

  most.popular <- tsp.dens$x[rev(order(tsp.dens$y))[1]]

  yrange <- most.popular - c(0, sqrt(var(tsp.median)))
  abline(h=yrange, xpd=FALSE)

  xlim <- par("usr")[1:2]
  xp <- c(xlim, rev(xlim))
  yp <- rep(yrange,each=2)
  polygon(xp, yp, density=density)

  ilow <- min(which(tsp.median >= yrange[2]))
  rest <- tsp.median[-seq(ilow)]
  ihigh <- min(which(
    rest >= yrange[1] | rest <= yrange[2])) + ilow - 1

  low  <- olags[ilow]
  high <- olags[ihigh]

  if (!length(low))
    low <- olags[1]
  if (!length(high))
    high <- olags[1]

  mtext(ifelse1(low != high,
   paste("olag =", low,"to", high),
   paste("olag=", low, sep="")),
     side=3, adj=1, line=1, cex=cex)

  invisible(c(low,high))
}

###
# plot.spaceTime
###

"plot.spaceTime" <- function(x, lty=1, color=seq(8),
  ylab="Spatial Separation", xlab="Orbital Lag",
  add=FALSE, cex=1, main.cex=0.7*cex, main=NULL, pch=".", ...)
{
  xatt <- attributes(x)
  xx   <- xatt$olags
  yy   <- as.matrix(x)
  rows <- numRows(x)
  cols <- numCols(x)
  iy   <- ifelse1(cols <= 10, seq(cols), floor(seq(1,cols,length=10)))

  old.plt <- splitplot(1,1,1)
  on.exit(par(old.plt))
  old.mar <- par(mar=c(5,5,4,5))
  on.exit(par(old.mar))
  old.xpd <- par(xpd=TRUE)
  on.exit(par(old.xpd))

  matplot(x=xx, y=yy, col=color, lty=lty, type="l", pch=pch,
    xlab=xlab, ylab=ylab, add=add, cex=cex)

  em <- c(strwidth("m"), strheight("m"))

  text(x=par("usr")[2] + em[1], y=yy[rows, iy],
    labels=paste("p =", round(iy * xatt$probability,4)),
    col=rep(color, ceiling(cols/length(color)))[iy],adj=0, cex=0.8)

  invisible(NULL)
}

###
# print.spaceTime
###

"print.spaceTime" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)
  main <- paste("Space-Time Separation Analysis for", xatt$series.name)

  z <- list(
    "Embedding points"=xatt$n.embed,
    "Embedding dimension"=xatt$dimension,
    "Time lag"=xatt$tlag,
    "Probability"=xatt$probability,
    "Range of time separations"=paste(range(xatt$olags), collapse=" to "))

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}


################################################
##
## Class FDSimulate
## Constructor function: FDSimulate
## Methods:
##
##   plot.FDSimulate
##   print.FDSimulate
##
################################################

###
# print.FDSimulate
###

"print.FDSimulate" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)
  main <- paste("Time Varying FD Process Simulation")

  z <- list(
    "Range delta"=paste(range(xatt$delta), collapse=" to "),
    "Number of unique deltas"=length(unique(xatt$delta)),
    "Range innovations variance"=paste(range(xatt$innov), collapse=" to "),
    "Number of unique innov. var."=length(unique(xatt$innov)),
    "Method"= xatt$method)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# plot.FDSimulate
###

"plot.FDSimulate" <- function(x, simulation=TRUE, delta=TRUE, innovations.var=TRUE,
  lty=1, color=seq(8), xlab="Time Index", add=FALSE, ...)
{

  xatt <- attributes(x)

  z <- list(delta=xatt$delta, "innovations\nvariance"=xatt$innov, "tvfd\nsimulation"=asVector(x))
  plots <- c(delta, innovations.var, simulation)
  if (!any(plots))
    stop("Must specify at least one variable to plot: simulation, delta, or innovations.var")
  z <- z[plots]

  stackPlot(x=seq(along=z[[1]]), y=z, lty=lty, col=color, xlab=xlab)

  invisible(NULL)
}

"asVector" <- function(x) if (inherits(x, "signalSeries")) x@data else as.vector(x)

"eda.plot" <- function (x, ...) UseMethod("eda.plot")