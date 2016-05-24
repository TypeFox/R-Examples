#############################################################
## WMTSA package wavelet transform coefficient functionality
##
##  Functions:
##
##    wavIndex
##
##  Constructor Functions and methods:
##
##    wavBoundary
##
##      [.wavBoundary
##      plot.wavBoundary
##      print.wavBoundary
##      print.summary.wavBoundary
##      summary.wavBoundary
##
#############################################################

###
# wavBoundary (constructor)
##

"wavBoundary" <- function(x)
{
  # define local functions

  "crystal.energy" <- function(x) sum(x^2)

  "crystals.sJ2d1" <- function(str){
    sindex <- which(as.logical(regexpr("s",str) > 0))
    ifelse1(sindex == 1, str, rev(str))
  }

  "extract.coefficients" <- function(W, low, high, interior=TRUE)
  {
    crystals <- names(W)

    obj <- lapply(crystals,
      function(crystal, W, low, high, interior){

       # grab the appropriate crystal
       WW <- W[[crystal]]

       # make an an index vector locating the interior coefficients
       interior.indices <- seq(from=low[crystal], to=high[crystal])
       all.indices      <- seq(from=1, to=length(WW))

       if (interior)
         indices.to.grab <- interior.indices
       else
         indices.to.grab <- all.indices[setdiff(all.indices,interior.indices)]

       y <- WW[indices.to.grab]

       if (length(indices.to.grab) && length(y))
         names(y) <- paste(indices.to.grab)

       invisible(y)
     },
     W=W, low=low, high=high, interior=interior)

    names(obj) <- crystals
    obj
  }

  # obtain the interior wavelet and scaling coefficient ranges.
  # wavIndex() also performs a check on the input class.
  separation <- wavIndex(x)

  # obtain the wavelet transform data (with
  # crystal names attached to the data),
  # and series information
  wave.class <- class(x)

  if (is(x,"wavBoundary"))
    return(x)
  else if (is(x,"wavTransform")){
    series    <- x$series
    dict      <- x$dictionary
    crystals  <- wavSortCrystals(names(x$data), reverse= TRUE)
    WW        <- x$data[crystals]
    xform     <- x$xform
  }

  # extract non-boundary (interior) coefficients
  # and boundary coefficients
  low  <- separation$interior$low
  high <- separation$interior$high

  interior <- extract.coefficients(WW, low, high, interior=TRUE)
  boundary <- extract.coefficients(WW, low, high, interior=FALSE)

  y <- list(boundary=boundary,
	    interior=interior,
	    all=WW,
	    dictionary=dict,
	    series=series,
	    crystals=crystals,
	    xform=xform,
	    index=separation)

  # add energy
  y$boundary$energy <- unlist(lapply(boundary, crystal.energy))
  y$interior$energy <- unlist(lapply(interior, crystal.energy))

  # add length
  y$boundary$length <- separation$boundary$length[crystals]
  y$interior$length <- separation$interior$length[crystals]

  # add total enery and length objects
  y$all$energy <- y$boundary$energy + y$interior$energy
  y$all$length <- separation$all$length[crystals]

  # add names to each coefficient in the all object
  # which identifies it as a boundary or interior coefficient
  for (j in 1:length(crystals)){

    crystal                 <- crystals[j]
    coeff.labels            <- rep("b", length(y$all[[crystal]]))
    nb.index                <- low[crystal]:high[crystal]
    coeff.labels[nb.index]  <- "nb"
    names(y$all[[crystal]]) <- coeff.labels
  }

  # record shift
  y$shifted <- x$shifted
  oldClass(y) <- "wavBoundary"

  y
}

###
# print.wavBoundary
##

"print.wavBoundary" <- function(x, justify="left", sep=":", ...)
{
  series      <- x$series
  series.name <- wavTitle(x)
  dict        <- x$dictionary
  nms         <- wavSortCrystals(x$crystals)

  main <- paste("Boundary coefficient separation of",
	  upperCase(x$xform), "series:", series.name)

  z <- list(
    "Wavelet"=dict$wavelet,
    "Length of series"=dict$n.sample,
    "Number of levels"=dict$n.levels,
    "Boundary correction rule"=dict$boundary,
    "Filtering technique"=ifelse1(dict$conv,"Convolution","Correlation"),
    "Zero phase shifted"=x$shifted )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)
  cat("Crystals: ")
  if (length(nms) > 15)
    cat(paste(nms, collapse=" "),
	" ... (",length(nms)," bases)\n",sep="", fill=TRUE)
  else cat(nms, fill=TRUE)

  invisible(x)
}

###
# summary.wavBoundary
##

"summary.wavBoundary" <- function(object, ...)
{
  crystals <- wavSortCrystals(object$crystals)

  row.names <- c("Interior","Boundary","[Ratio]","Total")

  l.nonbound <- object$interior$length[crystals]
  l.bound    <- object$boundary$length[crystals]
  l.ratio    <- round(l.nonbound / l.bound)
  l.total    <- object$all$length[crystals]

  cmat <- rbind(l.nonbound, l.bound, l.ratio, l.total)
  dimnames(cmat) <- list(row.names,crystals)

  e.nonbound <- object$interior$energy[crystals]
  e.bound    <- object$boundary$energy[crystals]
  e.ratio    <- round(e.nonbound / e.bound)
  e.total    <- object$all$energy[crystals]

  emat <- rbind(e.nonbound, e.bound, e.ratio, e.total)
  dimnames(emat) <- list(row.names,crystals)

  obj <- list(emat=emat,
	cmat=cmat)

  oldClass(obj) <- "summary.wavBoundary"

  return(obj)
}

###
# print.summary.wavBoundary
##

"print.summary.wavBoundary" <- function(x, digits=max(2, .Options$digits - 4), ...)
{
  cat("\nNumber of coefficients:\n\n")
  print(x$cmat, digits=digits, ...)

  cat("\nEnergy of Wavelet Boundary \nand Interior Crystals:\n\n")
  print(x$emat, digits=digits, ...)
  cat("\n")

  invisible(x)
}

###
# plot.summary.wavBoundary
##

"plot.wavBoundary" <- function(x, x.axis=TRUE, y.axis=TRUE, type="l", plot=TRUE, xlab=NULL,
  title=NULL, bars=FALSE, vgap=.05, grid=FALSE, times=NULL, grid.lty=par("lty")+1,
  same.scale=NULL, zerocenter=FALSE, zeroline=FALSE, col=c("red","black"),
  bar=list(lty=1, col=c("red","black"), lwd=2), ...)
{
  series         <- x$series
  series.omitted <- as.logical(length(series@data)==0)
  series.name    <- wavTitle(x)
  pos            <- series@positions
  position.units <- series@units.position
  data.units     <- series@units
  is.shifted     <- x$shifted

  if (is.null(xlab))
    xlab <- ifelse1(is.null(series.name), "Position",
      ifelse1(length(position.units), position.units, "Position"))

  if (nchar(series.name) > 5)
    series.name <- paste(substring(series.name,1,5), "...")

  wavelet     <- x$dictionary$wavelet
  xform   <- upperCase(x$xform)

  # create title
  if (is.null(title)){
    title <- paste(xform, "(", wavelet, ")")
    if (!is.null(series.name))
      title <- paste(title, ":", series.name, sep="")
  }

  # obtain time vector for plotting

  combined.length <- length(x$boundary$d1) + length(x$interior$d1)
  if (is.null(times)){

    times <- ifelse1(series.omitted, seq(length=x$dictionary$n.sample),
      seq(from=pos@from, by=pos@by, length=pos@length))

    if (is.null(times))
      times <- seq(length=combined.length)
  }

  # obtain crystals names, boundary
  # crystals, and interior crystals
  crystals  <- x$crystals
  if (is.null(crystals))
    crystals <- as.character(1:n)

  alldata  <- x$all[crystals]
  nonbound <- x$interior[crystals]
  n <- length(crystals)

  if (length(type)!=n)
    type <- rep(type[1], n)
  ycenters <- mi <- rep(0., n)

  # compute yrange
  for(i in 1:n){

    yi <- alldata[[i]]
    mi[[i]] <- ifelse1(zerocenter, 2*max(abs(yi), na.rm=TRUE)*(1+vgap),
      diff(range(yi, na.rm=TRUE))*(1+vgap))

    if (mi[[i]] < .Machine$double.eps*100)   # constant case
      mi[[i]] <- diff(range(c(0, yi), na.rm=TRUE))*(1+vgap)
  }

  if (length(same.scale)==1 && is.logical(same.scale)){
    if (same.scale) same.scale <- 1:n
    else same.scale <- numeric(0)
  }
  if (length(same.scale)){
    maxmi <- max(unlist(mi[same.scale]), na.rm=TRUE)
    for(i in same.scale)
      mi[[i]] <- maxmi
  }

  # compute xrange
  tt.delta <- times[2] - times[1]
  tt.start <- times[1] - tt.delta/2
  tt.range <- times[length(times)] - times[1] + tt.delta

  # plot axes
  if (plot){
    xlim <- tt.start + c(0, tt.range) + c(-.03, .03+bars*.03)*tt.range
    ylim <- c(-n, 0)

    plot(xlim, ylim, type="n", axes=FALSE,
			     xlab="", ylab="",xlim=xlim, ylim=ylim, col=1)
    if (x.axis)
      axis(side=1, at=pretty(times), line=1, srt=0, col=1, ...)
    if (y.axis)
      axis(side=2, at=-(1:n)+.5, labels=crystals, tick=FALSE, line=1, srt=0, col=1, ...)
  }

  # TODO: replace this with wavStackPlot functionality
  # plot the data

  for(i in 1:n){

    yi <- x$all[[i]]
    ni <- length(yi)
    tt <- ((0:(ni-1))/ni + 1/(2*ni))*tt.range + tt.start

    if (zerocenter){
      ycenters[[i]] <- -i+.5
      yi <- yi/mi[[i]]+ycenters[[i]]
      yzero <- TRUE
    }
    else if (all(yi>0, na.rm=TRUE) || all(yi<0, na.rm=TRUE)){
      ycenters[[i]] <- -i
      my <- min(yi, na.rm=TRUE)
      if (all(yi==my)) my <- 0
      yi <- (yi-my)/mi[[i]] + ycenters[[i]] + vgap/2
      yzero <- FALSE
    }
    else if (all(yi==0, na.rm=TRUE)){
      ycenters[[i]] <- -i+.5
      yi <- yi+ycenters[[i]]
      yzero <- TRUE
    }
    else{
      ycenters[[i]] <- -i + .5 - mean(range(yi, na.rm=TRUE)/mi[[i]])
      yi <- yi/mi[[i]] + ycenters[[i]]
      yzero <- TRUE
    }
    if (plot){
      if (type[i]=="h")
	      segments(tt, rep(ycenters[[i]], ni), tt, as.vector(yi),col=col[i])
      else if (type[i]=="s"){
	      tt <- tt - (tt[2]-tt[1])/2
	      tt <- c(tt, tt[length(tt)]+(tt[2]-tt[1]))
	      tt <- as.vector(rbind(tt, tt))
	      yi <- as.vector(rbind(c(yi[2], yi), c(yi, yi[length(yi)-1])))

	      lines(tt, yi, col=col[i], ...)
      }
      else {

        # plot the data
        lines(tt, yi, type=type[i], col=bar$col[2], lwd=1, ...)

        # draw lines separating the boundary from the
        # interior coefficients ...
        left.interior  <- as.integer(names(nonbound[[i]][1]))
        right.interior <- as.integer(names(nonbound[[i]][length(nonbound[[i]])]))

        # ... add left interior coefficient line
        x12 <- tt[left.interior]
        y1  <- max(yi)
        y2  <- min(yi)

        if (!is.missing(x12) && !is.missing(y1) && !is.missing(y2))
          segments(x12, y1, x12, y2, lty=bar$lty, col=bar$col[1], lwd=bar$lwd)

        # add right interior coefficient line if not already
        # on far right edge of the signal
        if (is.shifted){
          x12 <- tt[right.interior]
          if (!is.missing(x12) && !is.missing(y1) && !is.missing(y2))
            segments(x12, y1, x12, y2, lty=bar$lty, col=bar$col[1], lwd=bar$lwd)
        }

        # add text to help the user locate the interior
        # and boundary coefficients
        if (i == 1 && length(x12) > 0) {
          separation <- 0.1
        if (!is.shifted){
          mtext('<-- b', side=3, at=x12, col=bar$col[2], adj=1 + separation)
          mtext('nb -->', side=3, at=x12, col=bar$col[1], adj=- separation)
        }
        else{
          mtext('<-- nb -->', side=3, at=mean(tt[c(left.interior,right.interior)]), col=bar$col[1], adj=0.5)
        }
     }

      }
      if (yzero && zeroline)
	     segments(tt[1], ycenters[[i]], tt[length(tt)], ycenters[[i]])
    }
  }
  if (plot && bars){
    mi <- 1/(mi)
    mi <- mi/max(mi)
    xloc <- times[length(times)] + .03*tt.range
    segments(xloc, -(1:n)+.5+mi/2, xloc, -(1:n)+.5-mi/2, lty=bar$lty, col=bar$col[1], lwd=bar$lwd)
  }
  if (plot && grid)
    abline(h=-(1:(n-1)), lty=grid.lty)

  title(main=title, xlab=xlab)

  invisible(ycenters)
}

###
# wavIndex
###

"wavIndex" <- function(x)
{
  # develop local functions
  "make.crystal.names" <- function(levels, wavelet=TRUE)
    paste(ifelse1(wavelet, "d", "s"),levels,sep="")

  wave.class <- class(x)
  supported.classes <- c("wavBoundary","wavTransform")

  if (!any(wave.class == supported.classes))
    stop("Class of input currently not supported")

  dict <- x$dictionary

  N      <- dict$n.sample
  L      <- length(dict$analysis.filter$low)
  J      <- dict$n.levels
  levels <- seq(length=J)

  # create crystal names
  wavelet.crystals <- make.crystal.names(levels, wavelet=TRUE)
  scaling.crystals <- make.crystal.names(levels, wavelet=FALSE)
  all.crystals     <- c(wavelet.crystals, scaling.crystals)

  # obtain index of boundary coefficients (Lj) and
  # maximum number of coefficients per scale (Nj)
  # from d1, ..., sJ order
  if (wave.class == "wavTransform"){

    if (!is.element(x$xform, c("dwt","modwt")))
      stop("Only the DWT and MODWT are supported for boundary/interior coefficient identification.")

    type <- mutilsTransformType(x$xform)


    separation <- itCall("RS_wavelets_transform_coefficient_boundaries",
			as.integer(J), as.integer(L), as.integer(N), type)
        #,
			#COPY=c(FALSE,FALSE,FALSE,FALSE),
			#CLASSES=c("integer","integer","integer","integer"),
                        #PACKAGE="ifultools")
    interior.low    <- rep(separation[[1]], 2)
    interior.high   <- rep(separation[[2]], 2)
    interior.length <- rep(separation[[3]], 2)
    boundary.length <- rep(separation[[4]], 2)
    all.length      <- interior.length + boundary.length
  }
  else if (wave.class == "wavBoundary"){
    Nj              <- x$all$length[wavelet.crystals]
    boundary.length <- x$boundary$length[wavelet.crystals]

    # calculate the length of interior coefficients
    # and total length of each scale
    interior.low    <- rep(pmin(Nj, boundary.length + 1), 2)
    interior.high   <- rep(Nj, 2)
    interior.length <- rep(Nj - boundary.length, 2)
    all.length      <- interior.length + boundary.length
  }

  # give names to the interior and all crystals
  names(interior.low)    <- all.crystals
  names(interior.high)   <- all.crystals
  names(interior.length) <- all.crystals
  names(boundary.length) <- all.crystals
  names(all.length)      <- all.crystals

  # form interior, boundary, and all lists. save into one list
  interior <- list(low=interior.low, high=interior.high, length=interior.length)
  boundary <- list(length=boundary.length)
  all      <- list(length=all.length)
  obj      <- list(interior=interior, boundary=boundary, all=all)

  # obtain the approximate zero phase shifts
  obj$shift <- wavZeroPhase(wavelet=dict$wavelet, levels=levels)

  obj
}
