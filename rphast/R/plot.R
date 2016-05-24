##' Is this a track?
##' @param x An object to test
##' @param ... ignored
##' @return A logical indicating whether x is an object of type track
##' @keywords plot
##' @export
##' @author Melissa J. Hubisz
is.track <- function(x, ...) {
  if (is.null(attr(x, "class"))) return(FALSE)
  attr(x, "class") == "track"
}


##' Get the coordinate range of a list of RPHAST results
##' @param ... a list of tracks
##' @param na.rm logical, indicating if \code{NA}'s should be omitted
##' @return a numeric vector of length two giving the minimum and maximum
##' coordinates in any wig or feature track in the list.  MSA tracks are
##' *only* used if there are no wig or feature tracks.
##' @keywords plot
##' @method range track
##' @export range.track
##' @export
##' @author Melissa J. Hubisz
range.track <- function(..., na.rm=FALSE) {
  l <- list(...)
  if (length(l)==1 && !is.track(l[[1]])) 
    l <- l[[1]]
  currrange <- c()
  msarange <- c()
  for (i in 1:length(l)) {
    if (!is.track(l[[i]])) stop("element ", i, " is not a track")
    if (l[[i]]$type == "wig") {
      currrange <- c(currrange, range(l[[i]]$data$coord, na.rm=na.rm))
    } else if (l[[i]]$type != "msa") {
      currrange <- c(currrange, range.feat(l[[i]]$data, na.rm=na.rm))
    } else { #msa track
      msarange <- c(msarange, coord.range.msa(l[[i]]$data, l[[i]]$refseq))
    }
  }
  if (length(currrange) == 0L) return(range(msarange))
  range(currrange)
}



##' Smooth a wig plot in rphast
##' @param coord The x coordinates of un-smoothed plot
##' @param score The scores cooresponding to the x coordinates (should be same length as coord)
##' @param numpoints The number of points to use in the new plot
##' @return A data frame with numpoints rows and columns "coord" and "score"
##' with smoothed values.  If \code{length(coord) <= numpoints}, it will contain the
##' original data
##' @export
##' @author Melissa J. Hubisz
smooth.wig <- function(coord, score, numpoints=300) {
  if (length(coord) != length(score)) stop("smooth.wig expects length(coord) == length(score)")
  if (length(coord) <= numpoints) return(data.frame(coord=coord, score=score))

  xlim <- range(coord)
  windowSize <- (xlim[2]-xlim[1])/(numpoints-1)
  newcoord <- seq(from=xlim[1], to=xlim[2], by=windowSize)

  newscore <- sapply(newcoord, function(x) {
    f <- (coord >= (x-windowSize) &
          coord <= (x+windowSize))
    if (sum(f)==0) return(NA)
    sum(score[f])/sum(f)
  })
  f <- !is.na(newscore)
  data.frame(coord=newcoord[f], score=newscore[f])
}


##' Make browser-like plot in rphast
##' @param x a list of tracks, created by the as.track.wig or as.track.feat
##' @param doLabels Logical.  Whether to plot the label above each plot.  Will be
##' recycled to the length of x.  Does not affect printing of shortLabels.
##' @param labels Labels to appear directly above each plot.
##' @param cex.axis The character expansion factor for axis annotations.
##' @param cex.labels The character expansion factor for the labels
##' @param cex.shortLabels The character expansion factor for the shortLabels
##' @param relWigSize The relative size of wig plots compared to feature plots
##' @param relMsaSize The relative size of msa plots compared to feature plots
##' @param xlim The range of the x coordinate to be plotted.  If \code{NULL} (the default), will
##' use the entire range represented in the resultList.
##' @param xlab The label for the x axis
##' @param ylab The label for the y axis
##' @param blankSpace The amount of vertical blank space between each plot.  This should be a single numeric
##' value between 0 and 1, representing the total fraction of the plot occupied by blank space.
##' @param axisDigits The number of digits to use on the y-axis for wig plots.
##' @param labelSpace The total fraction of vertical space given to plot labels.
##' @param belowLabelSpace The amount of space between a label and the plot it corresponds to,
##' in fractions of a character width.
##' @param lmar The size of the left margin (in number of lines)
##' @param ... Other options to be passed to \code{plot}.  See \link{par}.
##' @seealso \code{plotPhast}, which may be easier to use but less flexible
##' @keywords plot
##' @method plot track
##' @export plot.track
##' @export
##' @author Melissa J. Hubisz
plot.track <- function(x,
                       doLabels=TRUE,
                       cex.axis=1.0,
                       cex.labels=1.0,
                       cex.shortLabels=0.75,
                       relWigSize=5,
                       relMsaSize=5,
                       xlim=NULL,
                       xlab="coord", ylab="", blankSpace=0.25, axisDigits=3,
                       labelSpace=min(length(x)*0.05, 0.25),
                       belowLabelSpace=0.2, lmar=4, ...) {
  if (is.track(x)) tracks <- list(x)
  else tracks <- x
  numresult <- length(tracks)

  doLabels <- rep(doLabels, length.out=numresult)
  check.arg(doLabels, "doLabels", "logical", null.OK=FALSE,
            min.length=numresult, max.length=numresult)
  numlabels <- sum(doLabels)
  
  check.arg(relWigSize, "relWigSize", "numeric", null.OK=FALSE)
  check.arg(relMsaSize, "relMsaSize", "numeric", null.OK=FALSE)
  check.arg(xlim, "xlim", "numeric", null.OK=TRUE, min.length=2L, max.length=2L)
  check.arg(xlab, "xlab", "character", null.OK=FALSE)
  check.arg(ylab, "ylab", "character", null.OK=FALSE)
  check.arg(blankSpace, "blankSpace", "numeric", null.OK=FALSE)
  if (blankSpace < 0 || blankSpace > 1) stop("blankSpace should be between 0 and 1")
  check.arg(axisDigits, "axisDigits", "integer", null.OK=FALSE)

  resultType <- character(length(tracks))
  for (i in 1:length(tracks)) resultType[i] <- tracks[[i]]$type
  plotScale <- as.numeric(ifelse(resultType=="wig", relWigSize,
                                 ifelse(resultType=="msa", relMsaSize, 1)))
  plotScale <- plotScale/sum(plotScale)

  if (is.null(xlim))
    xlim <- range.track(tracks)

  plot(c(0), c(0), type="n", xlim=xlim, ylim=c(0, 1), xlab=xlab,
       ylab=ylab, yaxt="n", bty="n", cex.axis=cex.axis, ...)
  
  maxy <- 1
  if (numlabels > 0L) {
    check.arg(belowLabelSpace, "belowLabelSpace", "numeric", null.OK=FALSE)
    check.arg(labelSpace, "labelSpace", "numeric", null.OK=FALSE)
    cex.labels <- rep(cex.labels, length.out=numresult)
    check.arg(cex.labels, "cex.labels", "numeric", null.OK=FALSE, min.length=NULL, max.length=NULL)
    if (labelSpace < 0 || labelSpace > 1) stop("labelSpace should be between 0 and 1")
    labelSize <- sum(plotScale)*labelSpace/numlabels
    plotScale <- plotScale*(1.0-labelSpace)
  }
  cex.shortLabels <- rep(cex.shortLabels, length.out=numresult)
  check.arg(cex.shortLabels, "cex.shortLabels", "numeric", null.OK=FALSE,
            min.length=NULL, max.length=NULL)
  plotScale <- plotScale*(1-blankSpace)
  blankSpace <- blankSpace/numresult
  
  for (i in 1:numresult) {
    el <- tracks[[i]]
    if (doLabels[i]) {
      miny <- maxy - labelSize
      text(x=mean(xlim), y=miny, el$name, pos=3, offset=belowLabelSpace, cex=cex.labels[i])
      maxy <- miny
    }
    miny <- maxy - plotScale[i]
    yrange <- c(miny, maxy)

    if (resultType[[i]] == "wig") {
      coord <- el$data$coord
      f <- (coord >= xlim[1] & coord <= xlim[2])
      coord <- coord[f]
      score <- el$data$score[f]
      if (el$smooth) {
        if (is.null(el$numpoints)) smoothData <- smooth.wig(coord, score)
        else smoothData <- smooth.wig(coord, score, el$numpoints)
        coord <- smoothData$coord
        score <- smoothData$score
      }
      if (is.null(el$ylim)) {
        oldrange <- range(score)
        span <- oldrange[2] - oldrange[1]
        if (span >= 10) numDigits <- 0
        numDigits <- 1 - floor(log10(span))
        oldrange <- round(oldrange, digits=numDigits)
      }
      else oldrange <- el$ylim
      newscore <- (score - oldrange[1])*(maxy - yrange[1])/(oldrange[2]-oldrange[1]) + yrange[1]
      lines(coord, newscore, col=el$col)
      axis(side=2, at=yrange, labels=FALSE)
      mtext(format(oldrange[1], digits=axisDigits), side=2, line=0.5, at=yrange[1], las=1, cex=cex.axis)
      mtext(format(oldrange[2], digits=axisDigits), side=2, line=0.5, at=yrange[2], las=1, cex=cex.axis)
      if (!is.null(el$horiz.line)) {
        horiz.line <- (el$horiz.line - oldrange[1])*(maxy - yrange[1])/(oldrange[2]-oldrange[1]) + yrange[1]
        for (i in 1:length(el$horiz.line)) 
          abline(h=horiz.line, lty=el$horiz.lty[i], col=el$horiz.col[i])
      }
    } else if (resultType[[i]] == "msa") {
      plot.msa(el$data, xlim=xlim, ylim=yrange, add=TRUE, pretty=el$pretty, refseq=el$refseq, nuc.text=el$nuc.text, nuc.text.pos=el$nuc.text.pos, nuc.text.col=el$nuc.text.col)
    } else if (resultType[[i]] == "feat") {
      plot.feat(el$data,
                y=mean(yrange), height=(yrange[2]-yrange[1]),
                add=TRUE, col=el$col, fill.col=el$col)
    } else if (resultType[[i]] == "gene") {
      plot.gene(el$data, y=mean(yrange), height=(yrange[2]-yrange[1]),
                add=TRUE, col=el$col, arrow.density=el$arrow.density)
    } else stop("don't know track type ", resultType[[i]])
    if (!is.null(el$short.label))
      mtext(el$short.label, side=2, line=0.5, at=mean(yrange),
            las=1, cex=cex.shortLabels[i])
    maxy <- miny-blankSpace
  }

}


##' Create a wig track
##' @param wig A "wig" object (Must have elements wig$coord and wig$score which should both
##' be numeric vectors).  coord/score may be passed directly instead.
##' @param name The name of the track (a character string)
##' @param coord (Alternative to wig) A numeric vector of coordinates (to be used for x-axis)
##' @param score (Alternative to wig) A numeric vector of scores (y-axis coords), should be
##' same length as coord.
##' @param short.label An optional character string to be displayed in left
##' hand margin of track
##' @param col The color to be used to plot this track.
##' @param ylim The limits to be used on the y-axis.  If NULL use entire range
##' of score.
##' @param smooth A logical value indicating whether to perform smoothing when plotting
##' this track
##' @param numpoints (Only used if \code{smooth==TRUE}).  An integer value indicating how many
##' points to display in the smoothed wig.
##' @param horiz.line If non-NULL, draw horizontal lines on the display at the given y coordinates
##' @param horiz.lty If horiz.line is defined, use this line type.
##' @param horiz.col If horiz.line is defined, use this color
##' @return An object of type \code{track} which can be plotted with the plot.track
##' function
##' @keywords plot
##' @export
##' @author Melissa J. Hubisz
as.track.wig <- function(wig=NULL, name, coord=NULL, score=NULL, short.label=NULL,
                         col="black", ylim=NULL, smooth=FALSE, numpoints=250,
                         horiz.line=NULL, horiz.lty=2, horiz.col="black") {
  if (is.null(wig)) {
    if (is.null(coord) || is.null(score))
      stop("If wig not provided, coord and score must both be provided")
  } else {
    if (!(is.null(coord) && is.null(score)))
      stop("If wig is provided, coord and score should not be provided")
    if (ncol(wig) != 2 || names(wig)[1] != "coord") 
      stop("wig should have two columns and the first should be named \"coord\"")
    coord <- wig$coord
    score <- wig[,2]
  } 
  rv <- list()
  attr(rv, "class") <- "track"
  rv$data <- data.frame(coord=coord, score=score)
  rv$name <- name
  rv$short.label <- short.label
  rv$col <- col
  if (!is.null(ylim))
    rv$ylim <- ylim
  rv$smooth <- smooth
  if (rv$smooth)
    rv$numpoints <- numpoints
  if (!is.null(horiz.line)) {
    rv$horiz.line <- horiz.line
    rv$horiz.lty <- rep(horiz.lty, length.out=length(horiz.line))
    rv$horiz.col <- rep(horiz.col, length.out=length(horiz.line))
  }
  rv$type="wig"
  rv
}


##' Create an alignment track
##' @param x An object of type \code{msa}
##' @param name The name of the track (a character string)
##' @param refseq A character string identifying the sequence whose
##' coordinate range to use in the plot.  A value of \code{NULL} implies
##' the frame of reference of the entire alignment.
##' @param short.label An optional character string to be displayed in the
##' left-hand margin of the track
##' @param pretty If \code{TRUE}, display bases in the non-reference species
##' which are identical to the reference species as a dot.
##' @param nuc.text If not NULL, can be a vector of character strings.  Each
##' character string should be the same length as the MSA with respect to refseq.
##' @param nuc.text.pos If nuc.text is not NULL, can be either "top" or "bottom"
##' to indicate where to place nuc.text relative to the alignment.  Will be recycled
##' to the length of nuc.text.
##' @param nuc.text.col If nuc.text is not NULL, color to be used for printing nuc.text.
##' Will be recycled to the length of nuc.text.
##' @return An object of type \code{track} which can be plotted with the
##' plot.track function
##' @keywords plot
##' @seealso plot.track
##' @note alignment plots will only be displayed if the plot is zoomed in
##' enough to show the alignment data.
##' @export
##' @author Melissa J. Hubisz
as.track.msa <- function(x, name, refseq=names.msa(x)[1],
                         short.label=NULL, pretty=TRUE,
                         nuc.text=NULL, nuc.text.pos="bottom",
                         nuc.text.col="black") {
  if (!is.msa(x)) stop("x is not msa object")
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(name, "name", "character", null.OK=FALSE)
  check.arg(short.label, "short.label", "character", null.OK=TRUE)
  check.arg(pretty, "pretty", "logical", null.OK=FALSE)
  rv <- list()
  attr(rv, "class") <- "track"
  rv$data <- x
  rv$name <- name
  rv$short.label <- short.label
  rv$pretty <- pretty
  rv$refseq <- refseq
  rv$nuc.text <- nuc.text
  rv$nuc.text.pos <- nuc.text.pos
  rv$nuc.text.col <- nuc.text.col
  rv$type <- "msa"
  rv
}


##' Create a features track
##' @param x An object of type \code{feat}
##' @param name The name of the track (a character string)
##' @param short.label An optional character string to be displayed in
##' left hand margin of track
##' @param col The color to use plotting this track (can be a single
##' color or a color for each element)
##' @param is.gene A logical value; if \code{TRUE}, extract and plot gene
##' information from features.  The features which will be plotted are the
##' ones with types "CDS", "exon", or "intron".  All others will be ignored.
##' @param arrow.density (Only used if \code{is.gene==TRUE}.  The number of
##' lines per inch used to denote strand in gene plots.
##' @return An object of type \code{track} which can be plotted with plot.track
##' function
##' @keywords plot
##' @export
##' @example inst/examples/as-track-feat.R
##' @author Melissa J. Hubisz
as.track.feat <- function(x, name, short.label=NULL, col="black",
                          is.gene=FALSE, arrow.density=10) {
  rv <- list()
  attr(rv, "class") <- "track"
  rv$data <- x
  rv$name <- name
  rv$short.label <- short.label
  rv$col <- col
  if (is.gene) {
    rv$type = "gene"
    rv$arrow.density=arrow.density
  } else {
    rv$type="feat"
  }
  rv
}


