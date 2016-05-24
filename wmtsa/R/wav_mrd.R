#############################################################
## WMTSA package multiresolution decomposition functionality
##
##  Functions:
##
##    wavMRDSum
##
##  Constructor Functions and methods:
##
##    wavMRD
##
##      [.wavMRD
##      [<-.wavMRD
##      [[.wavMRD
##      as.matrix.wavMRD
##      boxplot.wavMRD
##      crystal.names.wavMRD
##      dotchart.wavMRD
##      eda.plot.wavMRD
##      plot.wavMRD
##      print.wavMRD
##      print.summary.wavMRD
##      reconstruct.wavMRD
##      stack.plot.wavMRD
##      summary.wavMRD
##
#############################################################

###
# wavMRD (constructor)
###

"wavMRD" <- function(x, level=NULL, osc=NULL)
{
  dict     <- x$dictionary
  n.sample <- dict$n.sample
  n.level  <- dict$n.level

  # put modwt crystals into a matrix the crystals MUST be in rows
  xform     <- lapply(x$data, function(x) matrix(x, nrow=1))
  is.dwt    <- is(x, "wavTransform") && is.element(x$xform, c("dwt","modwt"))
  data.name <- wavTitle(x)

  # obtain the wavelet and scaling filters
  wavelet.filter <- dict$analysis.filter$high
  scaling.filter <- dict$analysis.filter$low
  filters        <- list(wavelet.filter,scaling.filter)

  # map the transform to the transform type
	supported.types <- c("modwt","modwpt","dwt","dwpt")
	type <- match(x$xform, supported.types) - 1

	if (is.missing(type))
	  stop("Unsupported transform type. Must be one of the following: ",
	    paste(supported.types, collapse=", "))

  # in the case where both the level and osc
  # arguments are null, develop a full
  # decomposition. if the transform is a (mo)dwpt,
  # then select the first level by default
  if (is.null(level) && is.null(osc)){
    if (is.dwt){
      level <- c(1:n.level,n.level)
      osc   <- c(rep(1,n.level), 0)
    }
    else{
      level <- rep(n.level, 2^n.level)
      osc   <- seq(0, 2^n.level-1)
    }
  }
  else if (is.null(level) && !is.null(osc)){
    stop(paste("The level argument must be specified",
	  "when local node(s) are specified"))
  }
  else if (!is.null(level) && is.null(osc)){
    if (is.dwt){
      level <- unique(sort(level))
      J     <- length(level)
      osc   <- rep(1,J)

      if (level[J] == n.level){
	    level <- c(level,n.level)
	    osc   <- c(osc,0)
      }
    }
    else{
      level    <- unique(sort(level))
      newlevel <- NULL
      for (j in level){
	    osc      <- c(osc, seq(0, 2^j-1))
	    newlevel <- c(newlevel, rep(j, 2^j))
      }
      level <- newlevel
    }
  }

  # check the level(s)
  if (is.null(level))
    level <- 1
  if (any(level > n.level))
    stop("Level(s) exceeds maximum level in transform")

  # check the osc index
  if (is.null(osc)){
    # if a node is specified then a single
    # level should be specified
    level <- level[1]
    if (is.dwt && (x$xform == "dwt")){
      if (level == dict$n.levels)
	osc <- 0:1
      else
	osc <- 1
    }
    else
      osc <- seq(0,2^level - 1)
  }

  # sort unique nodes
  if (length(level) == 1)
    osc <- sort(unique(osc))

  # form detail matrix
  detail <- matrix(ncol=n.sample, nrow=length(osc))

  # calculate the crystal detail sequence
  crystals <- rep("", length(osc))

  for (n in seq(length(osc))){

    j <- min(n, length(level))

    detail[n,] <-
      as.vector(
        itCall("RS_wavelets_transform_packet_detail",
          xform, filters, as.integer(level[j]), as.integer(osc[n]), as.integer(type)))
    #
          #COPY=rep(FALSE,5),
          #CLASSES=c("list","list","integer","integer","integer"),
          #PACKAGE="ifultools"))

    if (is.dwt)
      crystals[n] <- paste(switch(osc[n]+1,"S","D"), level[j], sep="")
    else
      crystals[n] <- paste("W", level[j], ".", osc[n], sep="")
  }

  # name the details
  dimnames(detail) <- list(crystals, paste("t=",seq(0,n.sample-1),sep=""))
  detail <- t(detail)

  # store name
  series.name <- wavTitle(x)
  dict <- x$dictionary
  dict$series.name <- series.name

  # store attributes
  attr(detail, "crystal.names") <- crystals
  attr(detail, "dictionary")    <- dict
  attr(detail, "type")          <- "crystal"
  attr(detail, "crystal.type")  <- class(x)
  attr(detail, "xform")         <- x$xform
  attr(detail, "data.name")     <- data.name

  # cast in mrd class
  oldClass(detail) <- "wavMRD"

  detail
}

###
# [.wavMRD
##

"[.wavMRD" <- function(x,  ...,drop=TRUE)
{
  i <- ..1
  xattr <- attributes(x)
  y <- as.matrix(x)[, i, drop=FALSE]

  oldClass(y) <- oldClass(x)

  attr(y, "dictionary") <- xattr$dictionary
  attr(y, "tspar") <- xattr$tspar
  attr(y, "type") <- xattr$type
  attr(y, "crystal.type") <- xattr$crystal.type
  attr(y, "crystal.names") <- dimnames(y)[[2]]
  y
}

###
# [<-.wavMRD
##

"[<-.wavMRD" <- function(x, i, value)
{
  xattr <- attributes(x)
  if(is.numeric(value)){
    y <- as.matrix(x)
    y[, i] <- value
    attributes(y) <- xattr
    return(y)
  }
  if(!is.null(value))
    stop("Invalid value for assignment")

  # value is NULL; so return all but x[i]
  if(is.numeric(i))
    x[-i]
  else if(is.logical(i))
    x[!i]
  else if(is.character(i))
    x[!match(dimnames(as.matrix(x))[[2]],i,nomatch=0)]
  else
    stop("Invalid subscript for assignment")
}

###
# [[.wavMRD
##

"[[.wavMRD" <- function(x, ...)
{
  xattr <- attributes(x)
  y     <- as.vector(x[..1])
  if(!is.null(xattr$tspar)) attr(y, "tspar") <- xattr$tspar
  else if(!is.null(xattr$tsp)) attr(y, "tsp") <- xattr$tsp
  y
}


###
# as.matrix.wavMRD
##

"as.matrix.wavMRD" <- function(x, mode="any", names=TRUE, ...)
{
  y <- x
  attributes(y) <- NULL
  dims   <- dim(x)
  nrow   <- dims[1]
  ncol   <- dims[2]
  dimnms <- dimnames(x)

  y      <- matrix(as.vector(y, mode=mode), nrow=nrow, ncol=ncol)
  if (names)
    dimnames(y) <- dimnms
  y
}

###
# boxplot.wavMRD
##

"boxplot.wavMRD" <- function(x, data=TRUE, ...)
{
  n   <- ncol(x)
  nms <- rev(crystal.names(x))
  if(is.logical(data) && data){
    xlist <- vector("list", n+1)
    xlist[[1]] <- reconstruct.wavMRD(x)
    for(i in 1:n) xlist[[n+2-i]] <- x[i]
    nms <- c("Data", nms)
  }
  else{
    xlist <- vector("list", n)
    for(i in 1:n) xlist[[n+1-i]] <- x[i]
  }
  boxplot.default(xlist, names=nms, ...)
  invisible(NULL)
}

###
# crystal.names.wavMRD
###

"crystal.names.wavMRD" <- function(x, ...)
  attr(x, "crystal.names")

###
# dotchart.wavMRD
###

"dotchart.wavMRD" <- function(x, data=TRUE, xlab="Energy(100%)", new=FALSE, ...)
{
  pns <- rev(crystal.names(x))
  da  <- reconstruct.wavMRD(x)
  de  <- sum(da^2)
  xe  <- rev(apply(oldUnclass(x)^2, 2, "sum"))/de
  if(is.logical(data) && data){
    xe <- c(1, xe)
    names(xe) <- c("Data", pns)
  }
  else
    names(xe) <- pns

  dotchart(xe, xlab=xlab, xlim=c(0, 1), new=new, ...)
  invisible(NULL)
}

###
# eda.plot.wavMRD
###

"eda.plot.wavMRD" <- function(x, data=TRUE, cex=.6, mex=.6, ...)
{
  old.plt <- splitplot(1,2,1)
  on.exit(par(old.plt))
  boxplot(x, data=data, cex=cex, ..., srt=90, adj=1)
  splitplot(1,2,2)
  dotchart.wavMRD(x, data=data, cex=cex, new=FALSE, ...)
  invisible(NULL)
}

###
# plot.wavMRD
###

"plot.wavMRD" <- function(x, n.top=15, vgap=.05, col=1, show.sum=NULL, add=FALSE, sort.energy=FALSE, ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  # attenuate clipping so that the text labels will fit
  old.xpd <- par(xpd=NA)
  on.exit(par(old.xpd))

  if (is.null(show.sum))
    show.sum <- all(dim(x) > 1)

  if (is.complex(x))
    stop("cannot plot complex data yet")

  nms <- crystal.names(x)
  xnc <- ncol(x)

  x  <- as.matrix(x)
  tt <- seq(nrow(x))

  # extract top crystals
  n.top <- min(n.top, length(nms))
  itop  <- seq(n.top)

  ixtop <- rev(order(colSums(x)))[itop]
  if (!sort.energy)
   ixtop <- sort(ixtop)

  x <- x[, ixtop, drop=FALSE]
  top.nms <- dimnames(x)[[2]]
  if (show.sum){
    x <- cbind(rowSums(x), x)
    dimnames(x) <- list(dimnames(x)[[1]],c("sum",top.nms))
  }

  nms <- dimnames(x)[[2]]

  nc <- ncol(x)

  if (length(col) < nc)
    col <- rep(col, nc)

  # scale all to same height
  for(j in 1:nc)
    x[,j] <- x[,j] - min(x[,j])

  x    <- x / ((1 + vgap) * max(x))
  yloc <- (nc:1) - 1

  plot(tt, x[,1], type="n", axes=FALSE, ylim=c(0, nc), xlab="", ylab="", ...)
  for(j in 1:nc)
    lines(tt, x[,j]+yloc[j], col=col[j], ...)

  text(rep(par("usr")[1], length(nms)), colMeans(as.matrix(x))+yloc, nms, adj=1)
  axis(side=1, at=pretty(tt), outer=FALSE, ...)

  invisible(yloc)
}

###
# print.wavMRD
###

"print.wavMRD" <- function(x, justify="left", sep=":", ...)
{
  dict        <- attr(x, "dictionary")
  series.name <- attr(x, "data.name")
  type        <- attr(x, "crystal.type")
  pns         <- crystal.names(x)

  if(type == "wp" || type == "dwt" || type == "dwtconv")
    main <- paste("Wavelet Decomposition of", series.name)
  else if(type == "wavTransform")
    main  <- paste(upperCase(attr(x, "xform")), "Multiresolution Decomposition of", series.name)
  else
  	main  <- paste("Local Cosine Decomposition of", series.name)

  crystals <- ifelse1(length(pns) > 12,
    paste(paste(pns[1:3], collapse=" "), " ... ",
      paste(pns[seq(length(pns)-2,length=3)], collapse=" "), " (", length(pns), " crystals)", sep=""),
    pns)

  z <- c(as.list(dict),
    list("Signal Components"=crystals))

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# print.summary.wavMRD
###

"print.summary.wavMRD" <- function(x, digits=max(2, .Options$digits-4), ...)
{
  cat("\nCorrelation Matrix:\n")
  print(round(x$cor, digits=digits), ...)
  cat("\nVariances:\n")
  print(round(x$var, digits=digits), ...)
  cat("\nStatistics for Components:\n")
  print(round(x$crystal.stat, digits=digits), ...)
  invisible(x)
}

###
# reconstruct.wavMRD
###

"reconstruct.wavMRD" <- function(x, ...) apply(x, 1, "sum")

###
# wavStackPlot.wavMRD
###

"wavStackPlot.wavMRD"  <- function(x, n.top=15, vgap=.05, col=1, ...)
  plot.wavMRD(x, n.top=n.top, vgap=vgap, col=col, ...)

###
# summary.wavMRD
###

"summary.wavMRD" <- function(object, ...)
{
  nm   <- crystal.names(object)
  m    <- length(nm)
  data <- reconstruct.wavMRD(object)
  allE <- sum(data^2)
  smat <- matrix(0, m, 9)
  dimnames(smat) <- list(nm, c("Min", "1Q", "Median", "3Q", "Max",
			"Mean", "SD", "MAD", "Energy.%"))
  x <- as.matrix(object)
  for(i in 1:m){
    xi           <- x[,i]
    smat[i, 1:5] <- quantile(xi, c(0, .25, .5, .75, 1))
    smat[i, 6]   <- mean(xi)
    smat[i, 7]   <- sqrt(var(xi))
    smat[i, 8]   <- mad(xi)
    smat[i, 9]   <- sum(xi^2)/allE
  }
  x             <- cbind(data, x)
  cor           <- cor(x)
  var           <- diag(var(x))
  dimnames(cor) <- list(c("", nm), c("Data", nm))
  names(var)    <- c("Data", nm)
  out           <- list(cor=cor, var=var, crystal.stat=smat)

  oldClass(out) <- "summary.wavMRD"
  out
}

###
# wavMRDSum
###

"wavMRDSum" <- function(x, wavelet="s8",
  levels=1, xform="modwt", reflect=TRUE,
  keep.smooth=TRUE, keep.details=TRUE)
{
  # check input arguments
  if (is.element(class(x),c("named","signalSeries")))
    x <- as.vector(x)
  if (is.complex(x))
    stop("Time series must be real-valued")
  if (!isVectorAtomic(x))
    stop("Time series must be a vector")
  if (length(x) < 2)
    stop("Time series must contain more than one point")
  if (anyMissing(x))
    stop("Time series contains NA values")
  checkScalarType(xform, "character")
  checkScalarType(wavelet, "character")
  checkVectorType(levels,"integer")
  checkScalarType(reflect,"logical")
  checkScalarType(keep.smooth,"logical")
  checkScalarType(keep.details,"logical")

  # initialize variables
  xform <- match.arg(lowerCase(xform), c("dwt","modwt"))
  decimated <- xform == "dwt"

  # apply logic
  if (!(keep.smooth | keep.details))
    stop("One of keep.smooth or keep.details must be TRUE")

  # obtain the wavelet and scaling filters. the wavelet argument
  # is checked here
  filters <- wavDaubechies(wavelet=wavelet,
    normalized=!decimated)[c("wavelet","scaling")]

  # initialize length parameters
  L <- length(filters$wavelet)
  N <- length(x)
  n.level <- max(levels)

  if (n.level < 1)
    stop("Number of wavelet transform decomposition levels must be positive")
  if ( (reflect && n.level > ilogb((N-1)/(L-1) + 1, base=2) ) ||
        (!reflect && n.level > ilogb(N, base=2) ) )
     stop("Number of wavelet transform decomposition levels exceeds maximum")

  if (reflect){

		# reflect the last Lj points of the series, where
		# Lj = (2^n.level - 1)(L - 1) + 1, and L=length(filters$wavelet)
		Lj <- (2^n.level - 1) * (L - 1) + 1
		ix <- seq(N-1, max(N-Lj,1))

		x <- c(x[rev(seq(along=ix))],x, x[ix])
  }

  # calculate a discrete wavelet transform
  xform <- ifelse1(decimated, wavDWT, wavMODWT)
  W <- xform(x, n.level=n.level, wavelet=wavelet, keep.series=FALSE)

  # specify the oscillation indices ala a DWPT:
  # set to NULL to have wavMRD() automatically
  # fill in the appropriate values in the case that
  # keep.details=TRUE. If keep.details=FALSE, set
  # the oscillation index to 0 (corresponding to the
  # first crystal of the last DWPT level) and (to be sure)
  # re-specify the level.
  if (keep.details)
  	osc <- NULL
  else{
    osc <- 0
    levels <- max(levels)
  }

  # calculate the details and smooths.
  # ix an iy are the indices of the rows and columns, respectively,
  # to extract from the MRD matrix (the columns of the MRD matrix
  # contain the details and smooth). The column indices are obtained
  # based on the following truth table:
  #
  # keep.details  keep.smooth  iy
  # ------------  -----------  -------------
  # TRUE          TRUE         1:(J+1)
  # TRUE          FALSE        1:J
  # FALSE         TRUE         1
  #
  # where J is the number of levels.
  ix  <- ifelse1(reflect, seq(from=length(ix) + 1, length=N), seq(N))
  J   <- length(levels)
  iy  <- seq(ifelse1(keep.details & keep.smooth, J + 1,
    keep.details & !keep.smooth, J,
    !keep.details & keep.smooth, 1))
  mrd <- as.matrix(wavMRD(W, level=levels, osc=osc))[ix,iy]

  # form smooth approximation
  ifelse1(numCols(mrd) > 1, rowSums(mrd), mrd)
}
