################################################
## IFULTOOLS plot functions
##
##   autoKey
##   autoText
##   gridOverlay
##   scaleZoom
##   sparsestQuadrant
##   splitplot
##   stackPlot
##   stringSize
##
###################################################

###
# autoKey
###

"autoKey" <- function(x, y=NULL, args.=NULL, nquadrant=5)
{
  # automatically places a key in the sparsest
  # location of the space spanned by the
  # the x and y input coordinates. The argument 'args.' is a list
  # of arguments that define the key and 'nquadrant'
  # is the number of equi-sized quadrants to use in dividing
  # the space in both the x- and y-direction. e.g., if
  # nquadrant=3, the x-y plane is partitioned into a
  # 3x3 uniform grid and the position of the sparsest
  # quadrant (as defined by sparsetQuadrant())
  # is used to place the key. this function
  # does not take into consideration lines connecting
  # the data points and only considers the spatial distribution
  # of the points in the x-y plane. futhermore, this function
  # assumes that par("usr") is approximately c(range(x),range(y)),
  # i.e., that the user will not explicitly make space outside
  # span of the data to place the key. in this case, the user
  # should explicitly place the key as usual.

  checkScalarType(nquadrant,"integer")
  if (nquadrant < 1)
    stop("nquadrant must be positive")

  # obtain position list of sparset quadrant
  # arg checks are performed in sparsestQuadrant()
  pos <- sparsestQuadrant(x, y, nquadrant=nquadrant)

  # place the key
  args <- mergeList(args., pos)
  newargs <- c(args[c("x","y")], list(legend=""))
  do.call("legend",args=newargs)

  invisible(NULL)
}

###
# autoText
###

"autoText" <- function(x, y=NULL, text.="", cex=1, col=1, nquadrant=5)
{
  checkScalarType(text.,"character")
  checkScalarType(nquadrant,"integer")
  checkScalarType(cex,"numeric")
  checkScalarType(col,"numeric")
  if (col < 0)
    stop("col must be non-negative")
  if (nquadrant < 1)
    stop("nquadrant must be positive")

  pos <- sparsestQuadrant(x, y, nquadrant=nquadrant)
  textlist <- list(x=pos$x, y=pos$y, labels=text.[1], cex=cex, col=col,
    adj=pos$corner[1])
  do.call("text",args=textlist)
  invisible(NULL)
}

###
# gridOverlay
###

"gridOverlay" <- function(lty="dotted", col="gray", density=3, ...)
{
  checkScalarType(density,"numeric")
  if (lty < 0)
    stop("lty must be non-negative")
  if (col < 0)
    stop("col must be non-negative")
  if (density < 1)
    stop("grid density must be positive")

  xaxp   <- par("xaxp")
  yaxp   <- par("yaxp")
  usr    <- par("usr")
  dxtick <- (xaxp[2] - xaxp[1]) / xaxp[3] / density;
  dytick <- (yaxp[2] - yaxp[1]) / yaxp[3] / density;

  xtick <- seq(xaxp[1], usr[2], by=dxtick)
  ytick <- seq(yaxp[1], usr[4], by=dytick)

  # make sure we are not replicating an axis line
  nxtick <- length(xtick)
  nytick <- length(ytick)

  if (nxtick > 1){

    xtol <- abs(diff(xtick[1:2])) / 10
    if (abs(xtick[1] - usr[1]) < xtol)
      xtick <- xtick[-1]
    if (abs(xtick[length(xtick)] - usr[2]) < xtol)
      xtick <- xtick[-length(xtick)]
  }

  if (nytick > 1){
    ytol <- abs(diff(ytick[1:2])) / 10
    if (abs(ytick[1] - usr[3]) < ytol)
      ytick <- ytick[-1]
    if (abs(ytick[length(ytick)] - usr[4]) < ytol)
      ytick <- ytick[-length(ytick)]
  }

  abline(v=xtick, h=ytick, lty=lty, col=col, ...)

  invisible(NULL)
}

###
# scaleZoom
###

"scaleZoom" <- function(x, y, z=NULL, zoom=NULL, logxy="y",
  xy.linked=FALSE, xlab="x", ylab="y", log.base=2)
{
  # check input argument types and lengths
  checkScalarType(logxy,"character")
  checkScalarType(xlab,"character")
  checkScalarType(ylab,"character")
  checkScalarType(log.base,"integer")
  if (log.base <= 1)
    stop("log.base must be > 1")

  if (!is.null(zoom)){

    # check zoom vector
    x.range <- range(x)
    y.range <- range(y)

    if (zoom[1] < x.range[1] || zoom[2] > x.range[2])
      stop("x coordinates of zoom vector are out of range")
    if (zoom[3] < y.range[1] || zoom[4] > y.range[2])
      stop("y coordinates of zoom vector are out of range")

    ix.zoom <- which(x >= zoom[1] & x <= zoom[2])
    if (xy.linked){
      x <- x[ix.zoom]
      y <- y[ix.zoom]
    }

    iy.zoom <- which(y >= zoom[3] & y <= zoom[4])

    if (xy.linked){
      x <- x[iy.zoom]
      y <- y[iy.zoom]
    }
    else{
	  x <- x[ix.zoom]
	  y <- y[iy.zoom]
    }

    if (!is.null(z))
      z <- z[ix.zoom, iy.zoom]
  }
  else
    ix.zoom <- iy.zoom <- seq(x)

  if(is.element(logxy, "y")){
    y <- logb(y, base=log.base)
    ylab <- paste("log", log.base, "(", ylab, ")", sep="")
  }
  if(is.element(logxy, "x")){
    x <- logb(x, base=log.base)
    xlab <- paste("log", log.base, "(", xlab, ")", sep="")
  }

  list(x=x, y=y, z=z, ix=ix.zoom, iy=iy.zoom, xlab=xlab, ylab=ylab)
}

###
# sparsestQuadrant
###

"sparsestQuadrant" <- function(x, y=NULL, nquadrant=5)
{
  # partitions the space spanned by the input x-y data coordinates
  # into a uniform (nquadrant, nquadrant) grid and finds
  # the "sparsest quadrant": defined as the quadrant
  # that contains a minimum population (pmin) of x-y data points
  # and is and farthest (in an L-inf sense) from quadrants
  # containing a population > pmin.

  # define local functions
  "flat2rect" <- function(x, ncols)
     cbind(((x - 1) %% ncols), ((x - 1) %/% ncols)) + 1

  # ensure that x is a numeric vector
  if (!isVectorAtomic(x) || !is.numeric(x))
    stop("x must be a numeric vector")

  # define default for y coordinate if missing
  if (is.null(y)){

    y <- x

    if (inherits(x,"ts"))
      x <- time(x)
    else if (is(x, "signalSeries"))
      x <- positions(x)[]
    else
      x <- seq(along=x)
  }

  # coerce x and y into vectors
  x <- as.vector(x)
  y <- as.vector(y)

  # check input argumensts
  if (length(x) != length(y))
    stop("x and y vectors must be the same length")
  if (nquadrant < 2)
    stop("nquadrant must be exceed unity")

  # cut data into uniform divisions in each direction
  xcut <- cut(x, nquadrant)
  ycut <- cut(y, nquadrant)

  # form contingency table
  population <- table(xcut, ycut)

  # obtain smallest population value (minp)
  # and flattened index locations of those
  # elements in the contingency table that
  # that exceed minp
  minp <- min(population)
  nz   <- which(population > minp)

  if (length(nz)){

    # calculate rectangular coordinates (row,col)
    # of population values which exceed minimum
    nz <- flat2rect(nz, nquadrant)

    # calculate rectangular coordinates of
    # population values at minimum
    z <- flat2rect(which(population == minp), nquadrant)

    # for each minimum value, calculate the distance
    # (in an L-Inf sense) to the nearest minp value
    Linf <- unlist(apply(z, MARGIN=1,
      function(x, nz){
        min(rowMaxs(abs(sweep(nz,MARGIN=2,x))))
      }, nz=nz))

    # choose the quadrant with the maximum L-Inf distance,
    # yielding the "most remote" quadrant
    best  <- z[order(Linf)[length(Linf)],]
    xquad <- best[1]
    yquad <- best[2]
  }
  else{

    # all quadrants have equal density.
    # might as well choose the (1,1) block
    # of the contingency table
    xquad <- yquad <- 1
  }

  # convert most remote quadrant to x- and y-ranges
  zx <- levels(xcut)[xquad]
  zy <- levels(ycut)[yquad]

  token  <- ","
  zx <- substring(zx,2,nchar(zx)-1)
  zy <- substring(zy,2,nchar(zy)-1)
  xrange <- as.numeric(unlist(strsplit(zx, token)))
  yrange <- as.numeric(unlist(strsplit(zy, token)))
 
  # map ranges to key() corner and x,y values
  # use linear interpolation to locate legend
  # reference point and placement on the graph
  x <- approx(x=c(1,nquadrant), y=xrange, xout=xquad)$y
  y <- approx(x=c(1,nquadrant), y=yrange, xout=yquad)$y
  corner <- approx(x=c(1,nquadrant), y=0:1, xout=c(xquad, yquad))$y

  return(list(x=x, y=y, corner=corner))
}

###
# splitplot
###

"splitplot" <- function(nrows, ncols, i=1, new=as.logical(i > 1 && i <= nrows*ncols), gap=0.15)
{
  if (i == 1)
   frame()

  # check input arguments
  nrows <- as.integer(nrows)
  ncols <- as.integer(ncols)
  i <- as.integer(i)
  checkScalarType(new,"logical")
  if (nrows < 1)
    stop("subplot number of rows (nrow) must be a positive integer")
  if (ncols < 1)
    stop("subplot number of columns (ncol) must be a positive integer")
  checkRange(i, c(1, nrows*ncols))

  # calculate corresponding subplot row and column
  irow <- ((i - 1) %/% ncols) + 1
  icol <- ((i - 1) %% ncols) + 1

  # define normalized width of plot
  deltax <- (1 - (ncols + 1)*gap) / ncols #+ ifelse(ncols == 1, gap/nrows,0)
  deltay <- (1 - (nrows + 1)*gap) / nrows

  # calculate x and y plot coordinates
  x1 <- icol * gap  + (icol - 1) * (deltax + gap/ncols)
  x2 <- x1 + deltax
  y1 <- 1 - (irow * (gap  + deltay)) + (nrows-irow)*gap/nrows/2
  y2 <- y1 + deltay

  old.par <- par(plt=c(x1,x2,y1,y2))
  old.par <- c(old.par, new=FALSE)
  par(new=new)

  invisible(old.par)
}

###
# stackPlot
###

"stackPlot" <- function(x, y=NULL, xlty=NULL,
  bty="n", lty=1, col=1:8, lwd=1, rescale=TRUE, add=FALSE, cex=1, xaxs="r", xpd=TRUE,
  yaxis=list(add=TRUE, ndigit=3, col=1:8, lty=1, lwd=3, side="left", cex=1),
  xlab=list(text="", cex=1, srt=0, col=1),
  ylab=list(text=NULL, cex=1, srt=0, col=1:8, side="right"),
  main=list(text="", cex=1, srt=0, col=1, adj=0.5), ylim=NULL)
{
  # save plot parameters and restore upon exit
  if (!add){
   frame()
   old.plt <- splitplot(1,1,1)
   on.exit(par(old.plt))
   old.xpd <- par(xpd=xpd)
   on.exit(par(old.xpd))
   on.exit(par(c(old.plt, list(new=FALSE))))
  }

  if (is.null(y)){
  	if (is(x,"signalSeries")){
      y <- x@data
      ylab <- list(text=x@units, side="left")
      xlab <- list(text=x@units.position)
      main <- list(text=x@title)
      x <- as(positions(x),"numeric")
      xaxs <- "r"
  	}
  	else
  	  stop("y is missing with no default")
  }

  if (!is.list(ylab) && is.character(ylab))
    ylab <- list(text=ylab)
  if (!is.list(xlab) && is.character(xlab))
    xlab <- list(text=xlab)
  if (!is.list(main) && is.character(main))
    main <- list(text=main)

  # merge input argument lists
  xlab  <- mergeList(list(text="", cex=1, srt=0, col=1), xlab)
  ylab  <- mergeList(list(text=NULL, cex=1, srt=0, col=1:8, side="right"), ylab)
  main  <- mergeList(list(text="", cex=1, srt=0, col=1, adj=0.5), main)
  yaxis <- mergeList(list(add=TRUE, ndigit=3, col=1:8, lty=1, lwd=3, side="right", cex=1), yaxis)

  # coerce possible single row matrix
  # into a vector so that subsequent
  # data.frame() operation creates
  # a single-column data.frame.
  if (isVectorAtomic(y))
    y <- as.vector(y)

  y <- data.frame(y, check.names=FALSE)
  ynms <- names(y)

  if (!isVectorAtomic(x))
    stop("x must be a vector as defined by isVectorAtomic()")
  checkVectorType(xlab$text,"character")
#  checkVectorType(xlab$col,"integer")
  checkVectorType(xlab$srt,"integer")
#  checkVectorType(ylab$col,"integer")
  checkVectorType(ylab$srt,"integer")
  if (is.character(ylab$side))
    ylab$side <- match.arg(ylab$side, c("left","right"))
  else if (is.integer(ylab$side))
    ylab$side <- match.arg(ylab$side, c(2,4))
  else
    stop("ylab$side must be an object of class \"integer\" (2,4) ",
      "or class \"character\" (\"left\", \"right\")")
  checkVectorType(yaxis$add,"logical")
  checkVectorType(yaxis$cex,"integer")
  checkVectorType(yaxis$ndigit,"integer")
  checkVectorType(yaxis$lty,"integer")
  checkVectorType(yaxis$lwd,"integer")
  if (is.character(yaxis$side))
    yaxis$side <- switch(match.arg(yaxis$side, c("left","right")), left=2,right=4)
  else if (is.integer(yaxis$side))
    yaxis$side <- match.arg(yaxis$side, c(2,4))
  else
    stop("yaxis$side must be an object of class \"integer\" (2,4) ",
      "or class \"character\" (\"left\", \"right\")")
  if (!is.null(ylim)){
    checkVectorType(ylim,"numeric")
    if (length(ylim) != 2)
      stop("ylim must be a vector of two elements")
    rescale <- FALSE
  }

  # now that we know we have rectangular data, ensure
  # that x and y components have a consistent length
  if (numRows(x) != numRows(y))
    stop("x and y-data components have unequal lengths")

  # develop y-axis labels.
  # ensure labels are unique
  # for those named by as.data.frame()
  # for lack of an existing name, replace
  # with a suitable alternative
  if (!is.null(ylab$text)){
    if (!is.character(ylab$text))
      stop("ylab$text must be an object of class \"character\"")
   iy <- intersect(seq(along=ynms), seq(along=ylab$text))
   ynms[iy] <- ylab$text[iy]
  }
  ynms[grep("^X\\.", ynms)] <- NA

  # calculate y-axis parameters
  nplot <- numCols(y)
  yranges.original <- colRanges(y)

  # rescale the data if requested
  if (rescale){
    yy <- scale(y)[,]

    # if any y column is a constant,
    # all NAs will be returned by scale().
    # in this case return the original data
    ibad <- as.vector(apply(yy, MARGIN=2, function(x) all(is.na(x))))
    if (any(ibad)){
      for (i in which(ibad))
        yy[,i] <- y[,i]
    }

    y <- matrix(yy, ncol=numCols(y))

    yranges <- colRanges(y)
  }
  else if (!is.null(ylim))
    yranges <- matrix(rep(sort(ylim),nplot), ncol=nplot)
  else
    yranges <- yranges.original

  ydeltas <- diff(yranges)
  ygap    <- stdev(as.vector(unlist(y)), na.rm=TRUE) * logb(nplot,base=2)
  ybias   <- c(0, cumsum(ydeltas[,-nplot] + ygap)) + yranges[2,]


  # bias the y-data to separate the plots
  ybiased <- sweep(y, MARGIN=2, ybias, FUN="-")
  ycenter <- colMeans(as.matrix(ybiased))

  # initialize variables
  xlim <- range(x)

  if (!is.null(ylim)){

    ylim <- sort(ylim)

    # make sure that specified y-range spans data
    if (!all(ylim[1] <= yranges.original[1,] & ylim[2] >= yranges.original[2,]))
      stop("Specified ylim range does not span the range of all input series")

    ylim2 <- c(-(abs(diff(range(ylim))) * nplot + (nplot-1)*ygap), 0)
  }
  else
    ylim2 <- range(ybiased)

  # plot data
  plot(xlim, ylim2, type="n", yaxt="n", xlab=xlab$text, ylab="",
    lwd=min(lwd+2,3), bty=bty, cex=cex, xaxs=xaxs, xpd=xpd)
  matlines(x, ybiased, col=col, lty=lty, lwd=lwd, xpd=xpd)
  if (nchar(main$text))
    title(main$text, adj=main$adj, col=main$col, cex=main$cex, srt=main$srt)

  # add y-labels
  xgap <- 0.05 * diff(range(x))
  if (ylab$side == "left" || ylab$side == 2){
	  text.xpos <- min(min(x) - xgap, par("usr")[2])
      text.adj  <- 1
  }
  else{
	  text.xpos <- min(max(x) + xgap, par("usr")[2])
      text.adj  <- 0
  }

  text(rep(text.xpos, nplot), ycenter, ynms,
    cex=ylab$cex, col=ylab$col, srt=ylab$srt, adj=text.adj )

  # add x-lines if requested
  if (!is.null(xlty))
    abline(h=colMaxs(ybiased), lty=xlty, lwd=lwd)

  # add yaxis if requested
  if (yaxis$add){
	  yaxis.adj <- ifelse1(yaxis$side == 2, 1, 0)

	  for (i in seq(nplot)){

	    if (is.null(ylim)){
	      at <- range(ybiased[,i])
	      labels <- round(yranges.original[,i], yaxis$ndigit)
	    }
	    else{
	      at <- sort(ylim) - max(ylim) - (i - 1) * (ygap + abs(diff(range(ylim))))
	      labels <- round(ylim, yaxis$ndigit)
	    }

	    axis(side=yaxis$side, at=at, labels=labels,
        adj=yaxis.adj, cex=yaxis$cex, lwd=yaxis$lwd,
        lty=yaxis$lty, srt=0, las=2)
	  }
  }

  invisible(NULL)
}

###
# stringSize
###

"stringSize" <- function(x, srt=0, adj=0.5, cex=1)
{
  checkScalarType(x,"character")
  checkScalarType(cex,"numeric")
  checkScalarType(srt,"numeric")
  checkScalarType(adj,"numeric")
  if (cex <= 0)
    stop("cex must be positive")
  if (adj < 0 || adj > 1)
    stop("adj must be on the interval [0,1]")

  # calculates the relative x and
  # y spans of the input character string x
  # in the current par("usr") coordinates.

  # define local functions
  equivalentEms <- function(x)
  {

    # number of ems needed to equal the width of 10 of the following letters.
    # this is only a gross approximation using the default font for Windows
    ems <- list(a=6,b=6,c=6,d=6,e=6,f=3,g=6,h=6,i=3,j=3,k=6,l=3,m=10,
      n=6,o=6,p=6,q=6,r=4,s=6,t=4,u=6,v=4.5,w=8,x=6,y=6,z=6,
      A=8,B=8,C=8,D=8,E=8,F=7,G=9,H=8,I=2.5,J=5.5,K=8,L=7.5,M=10,N=8,
      O=9,P=8,Q=9,R=8,S=8,T=7.5,U=8,V=8,W=12,X=7.5,Y=8,Z=7.5)

    standard <- unlist(ems[strsplit(x,"")[[1]]]) / 10
    nspecial <- nchar(x) - length(standard)

    return(sum(standard) + nspecial / 2)
  }

  if (!is.character(x) || length(x) > 1)
    stop("x must be a single character string")
  if (adj < 0 || adj > 1)
    stop("adj must be on [0,1]")

  srt  <- srt * pi / 180
  em.size <- c(strwidth("m"), strheight("m"))
  span <- equivalentEms(x) * em.size * c(cos(srt), sin(srt))
  sgn  <- sign(span)

  # adjust for orthogonal srt and cex
  # e.g., if srt=90, then the centered x-width will be 1em
  span <- as.list(pmax(abs(span), em.size) * cex)
  names(span) <- c("x","y")

  # adjust for adj
  xright <- approx(x=c(0,1), y=c(span$x,0), xout=adj)$y
  xleft  <- xright - span$x
  yright <- approx(x=c(0,1), y=c(span$y,0), xout=adj)$y
  yleft  <- yright - span$y

  # cobine estimates
  xout <- c(xleft, xright)
  yout <- c(yleft, yright)

  # adjust for sign
  if (sgn[1] < 0)
    xout <- -rev(xout)
  if (sgn[2] < 0)
    yout <- -rev(yout)

  return(list(x=xout, y=yout))
}
