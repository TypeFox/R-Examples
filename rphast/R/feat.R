.makeObj.feat <- function(isPointer) {
  feat <- list()
  if (isPointer) {
    attr(feat, "class") <- c("feat", "list")
  } else attr(feat, "class") <- c("feat", "data.frame")
  feat
}


##' Creates a copy of a features object
##'
##' If x is stored in R (as it is by default), then this is no different
##' than x2 <- x.  But if it is stored as a pointer to a structure in C,
##' then this is the only way to make an explicity copy of the features.
##' @title Features copy
##' @param x an object of type \code{feat}
##' @return a features object which can be modified independently from the
##' original object
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
copy.feat <- function(x) {
  if (is.null(x$externalPtr)) return(x)
  result <- .makeObj.feat(TRUE)
  result$externalPtr <- .Call.rphast("rph_gff_copy", x$externalPtr)
  result
}



##' Read a features object from a file
##'
##' The function will guess the format of the input file automatically.
##'
##' @title Read a Feature File (GFF, BED, or GenePred)
##' @param filename the name of the file (can be GFF, BED, GenePred, or wig:
##' rphast will auto-detect)
##' @param pointer.only Whether to store object by reference instead of a
##' data.frame
##' @return If \code{pointer.only==FALSE}, a data.frame with columns corresponding
##' to the GFF specification.  Otherwise, an object which is a pointer to
##' an object stored in C.
##' @seealso \code{\link{feat}} for more description of features objects.
##'
##' \code{\link{msa}} for more explanation of the pointer.only option.
##'
##' \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
##' for a detailed description of GFF file format.  The columns in features
##' objects mirror the GFF column definitions.
##' 
##' \url{http://genome.ucsc.edu/FAQ/FAQformat} for descriptions
##' of BED and GenePred formats.
##' @keywords GFF
##' @keywords Genepred
##' @keywords BED
##' @keywords features
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/read-feat.R
##' @export
read.feat <- function(filename, pointer.only=FALSE) {
  filename <- check.arg(filename, "filename", "character", null.OK=FALSE,
                        min.length=1L, max.length=1L)
  feat <- .makeObj.feat(TRUE)
  feat$externalPtr <- .Call.rphast("rph_gff_read", filename)
  if (!pointer.only) {
    feat <- from.pointer.feat(feat)
  }
  feat
}


##' Create a new features object
##'
##' See \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
##' for more detailed description of each parameter.
##' 
##' All arguments which are provided should be vectors of equal length.
##' 
##' If pointer.only==FALSE, the new object is a data frame, with
##' columns mirroring the GFF Specification
##' Otherwise, it is a list containing a single element, which is
##' a pointer to an object stored in C.
##' @title Features Objects
##' @param seqname a character vector containing the name of the sequence.  If
##' the features correspond to regions of an alignment, then seqname
##' should be the name of the sequence in the alignment that is used
##' as the frame of reference in the features.  To use the entire alignment
##' as a frame of reference, set seqname to "MSA".
##' @param src The source of the feature
##' @param feature The feature type name
##' @param start The start of the feature.  Sequence numbering begins at 1.
##' @param end The end of the feature.  This is the last coordinate included
##' in the feature.
##' @param score The feature score, or NA if there is no score.
##' @param strand A character string which is either "+", "-", or "." (if
##' strand is not available or relevant).
##' @param frame A 0, 1, or 2, which specifies whether the feature is in frame.
##' @param attribute A feature attribute (character string).
##' @param pointer.only Whether to store object as a pointer to an object
##' in C, rather than as a data.frame in R.
##' @return If pointer.only==FALSE, returns a data.frame whose format
##' mirrors the GFF specification.  Otherwise, returns a list with a single
##' object, which is a external pointer to a C structure representing a
##' features object.
##' @seealso \code{\link{read.feat}}
##'
##' \code{\link{msa}} for more details on the pointer.only option.
##' @keywords features
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/feat.R
##' @export
feat <- function(seqname="default", src=".", feature=".",
                 start, end, score=NULL,
                 strand=NULL, frame=NULL, attribute=NULL,
                 pointer.only=FALSE) {
  
  check.arg(start, "start", "integer", null.OK=FALSE, min.length=NULL,
             max.length=NULL)
  len <- length(start)
  seqname <- rep(seqname, length.out = len)
  check.arg(seqname, "seqname", "character", null.OK=FALSE, min.length=len,
            max.length=len)
  src <- rep(src, length.out = len)
  check.arg(src, "src", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  feature <- rep(feature, length.out = len)
  check.arg(feature, "feature", "character", null.OK=FALSE, min.length=len,
             max.length=len)
  check.arg(end, "end", "integer", null.OK=FALSE, min.length=len,
             max.length=len)
  if (!is.null(score)) score <- rep(score, length.out=len)
  check.arg(score, "score", "numeric", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(strand)) strand <- rep(strand, length.out=len)
  check.arg(strand, "strand", "character", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(frame)) frame <- rep(frame, length.out = len)
  check.arg(frame, "frame", "integer", null.OK=TRUE, min.length=len,
             max.length=len)
  if (!is.null(attribute)) attribute <- rep(attribute, length.out = len)
  check.arg(attribute, "attribute", "character", null.OK=TRUE,
             min.length=len, max.length=len)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

# NOTE: I'd like to impose this, but in some genepred files elements have length -1
#  if (sum(end - start < 0) > 0)
#    stop("Error: end should be >= start for all features")

  if (pointer.only) {
    if (!is.null(score)) score <- as.numeric(score)
    if (!is.null(strand)) strand <- as.character(strand)
    if (!is.null(frame)) frame <- as.integer(frame)
    if (!is.null(attribute)) attribute <- as.character(attribute)
    ptr <- .Call.rphast("rph_gff_new", as.character(seqname),
                        as.character(src), as.character(feature),
                        as.integer(start), as.integer(end),
                        score, strand, frame, attribute)
    x <- .makeObj.feat(TRUE)
    x$externalPtr <- ptr
  } else {
    x <- data.frame(seqname=seqname, src=src, feature=feature,
                    start=start, end=end)
    if (!is.null(score)) x <- cbind(x, score)
    if (!is.null(strand)) x <- cbind(x, strand)
    if (!is.null(frame)) x <- cbind(x, frame)
    if (!is.null(attribute)) x <- cbind(x, attribute)
    attr(x, "class") <- c("feat", "data.frame")
  }
  x
}


##' Take a set of features stored in R and return one stored by reference
##' @title Features To Pointer
##' @param x an object of type \code{feat} stored as a data frame in R
##' @return an object of type \code{feat} stored by reference as a pointer to an
##' object created in C.
##' @seealso \code{\link{feat}} for more details on features storage
##' options.
##' @author Melissa J. Hubisz
##' @example inst/examples/as-pointer-feat.R
##' @export
as.pointer.feat <- function(x) {
  if (!is.null(x$externalPtr))
    return(x)
  feat(x$seqname, x$src, x$feature,
       x$start,   x$end, x$score,
       x$strand,  x$frame,  x$attribute,
       pointer.only=TRUE)
}


##' Prints a features object.
##' @title Printing a features Object
##' @param x an object of type \code{feat}
##' @param ... further arguments to be passed to or from other methods
##' @keywords features
##' @seealso \code{\link{write.feat}}
##' @author Melissa J. Hubisz and Adam Siepel
##' @method print feat
##' @export
##' @export print.feat
print.feat <- function(x, ...) {
  cat(paste("Features object\n"))
  write.feat(x, NULL)
}


##' Write a features object to a file in GFF format.
##' @title Writing a features Object
##' @param x an object of type \code{feat}
##' @param file The name of the file to write to (will be overwritten)
##' @keywords features GFF
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/write-feat.R
##' @export
write.feat <- function(x, file) {
  check.arg(file, "file", "character", null.OK=TRUE)
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  invisible(.Call.rphast("rph_gff_print", file, x$externalPtr))
}


##' Get the number of rows in a features object
##' @title Number of Features
##' @param x An object of type \code{feat}
##' @return An integer containing the number of rows in each features object
##' @export
##' @export nrow.feat
##' @method nrow feat
##' @keywords features
##' @example inst/examples/nrow-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
nrow.feat <- function(x) {
  if (is.data.frame(x)) {
    attr(x, "class") <- "data.frame"
    nrow(x)
#    attr(x, "class") <- "data.frame"
#    NextMethod("nrow", x)
#    NextMethod()
#    UseMethod("nrow", x)
  } else {
    if (is.null(x$externalPtr)) stop("invalid object sent to nrow.feat")
    .Call.rphast("rph_gff_numrow", x$externalPtr)
  }
}


##' Get the number of columns in a features object
##' @title Number of Columns in Features
##' @param x An object of type \code{feat}
##' @return An integer containing the number of columns in the features object
##' @note If the features object is stored as a pointer in C, the number
##' of columns is always 9.
##' @export
##' @export ncol.feat
##' @method ncol feat
##' @keywords features
##' @example inst/examples/ncol-feat.R
##' @author Melissa J. Hubisz
ncol.feat <- function(x) {
  if (is.data.frame(x)) {
    attr(x, "class") <- "data.frame"
#    NextMethod("nrow")
    ncol(x)
  } else return(9)  # feat objects stored in C always have 9 columns
}


##' Prints a brief summary of a features object.
##' @title Features Summary
##' @param object an object of type \code{feat}
##' @param ... further arguments to be passed to or from other methods
##' @export
##' @export summary.feat
##' @method summary feat
##' @keywords features
##' @example inst/examples/summary-feat.R
##' @author Melissa J. Hubisz
summary.feat <- function(object, ...) {
  if (is.null(object$externalPtr)) {
    as <- "stored as data frame"
    object <- as.pointer.feat(object)
  } else as <- "stored as a pointer to a C structure"
  nrow <- nrow.feat(object)
  cat(paste("features object with", nrow, "rows", as))
  cat("\n")
}


##' Convert a features object to a data frame
##' @title Features to Data Frame
##' @param x an object of type \code{feat}
##' @param row.names optional names for each feature
##' @param optional logical, if \code{TRUE}, setting row names and 
##' converting column names (to syntactic names: see
##' \code{\link{make.names}} is optional.
##' @param ... additional arguments to be passed to other methods
##' @return a data frame containing features data
##' @seealso \code{\link{feat}} for a description of features data frames,
##' and \code{\link{as.pointer.feat}} for conversion in the other
##' direction.
##' @export
##' @export as.data.frame.feat
##' @method as.data.frame feat
##' @example inst/examples/as-data-frame-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
as.data.frame.feat <- function(x, row.names=NULL, optional=FALSE, ...) {
  if (is.data.frame(x)) return(x)
  if (!is.null(x$externalPtr)) {
    x <- .Call.rphast("rph_gff_dataframe", x$externalPtr)
  }
  attr(x, "class") <- "list"
  as.data.frame(x, row.names, optional, ...)
}


##' Convert a features object from C memory (external pointer) to R memory
##' @param x A features object stored as a pointer to C memory
##' @return A features object (inheriting from the data.frame class) stored
##' in R memory
##' @export
##' @author Melissa J. Hubisz
from.pointer.feat <- function(x) {
  x <- as.data.frame.feat(x)
  attr(x, "class") <- c("feat", "data.frame")
  x
}



##' Get the dimensions of a features object
##' @title Feature dimensions
##' @param x an object of type \code{feat}
##' @return An integer vector of length two containing the number of
##' rows and number of columns in the features object.
##' @export
##' @export dim.feat
##' @method dim feat
##' @keywords features
##' @author Melissa J. Hubisz
##' @example inst/examples/dim-feat.R
dim.feat <- function(x) {
  c(nrow.feat(x), ncol.feat(x))
}


##' Get the range of a features object
##' @title Features range
##' @param ... Objects of type \code{feat}
##' @param na.rm Whether to remove values of NA before calculating range.
##' @return A vector of size 2 indicating minimum and maximum coord in
##' the features object
##' @method range feat
##' @export range.feat
##' @export
##' @keywords features
##' @author Melissa J. Hubisz
range.feat <- function(..., na.rm=FALSE) {
  feats <- list(...)
  mins <- numeric(length(feats))
  maxs <- numeric(length(feats))
  for (i in 1:length(feats)) {
    x <- feats[[i]]
    if (is.data.frame(x)) {
      r <- range(c(x$start, x$end), na.rm=na.rm)
    } else {
      r <- c(.Call.rphast("rph_gff_minCoord", x$externalPtr),
             .Call.rphast("rph_gff_maxCoord", x$externalPtr))
    }
    mins[i] <- r[1]
    maxs[i] <- r[2]
  }
  c(min(mins, na.rm=na.rm),
    max(maxs, na.rm=na.rm))
}


##' plot features
##' @title Features plot
##' @param x an object of type \code{feat}
##' @param y the location of the plot on the y axis
##' @param height the height of the boxes
##' @param plottype either "r" for rectangles or "a" for arrows, "b"
##' for arrows within rectangles, or "l" for line segments only.
##' @param arrow.density If plottype=="a" or "b", then this gives the density
##' of arrows in arrows per inch.  Otherwise it gives the density of
##' shading lines in the rectangles, and a value of \code{NULL} implies
##' no shading lines.
##' @param angle angle (in degrees) of the shading lines or arrows.
##' @param col color to draw the boxes/lines/arrows with.
##' @param fill.col Color to fill the rectangles with.  If \code{NULL}
##' then do not fill.
##' @param lty line type for lines, arrows, borders, and shading
##' @param lwd line width for lines, arrows, borders and shading
##' @param add if \code{TRUE}, add to existing plot
## @param labels whether to label the boxes.  If a vector of strings gives
## the lables for each element.  If \code{TRUE}, use x$feature for labels.
##' @param xlim A numerical vector of length 2 giving the range for the x-axis.
##' @param ylim A numerical vector of length 2 giving the range for the y-axis.
##' @param ... graphical parameters to be passed to \code{plot}.
##' @method plot feat
##' @export plot.feat
##' @export
##' @keywords features plot
##' @example inst/examples/plot-feat.R
##' @author Melissa J. Hubisz
plot.feat <- function(x, y=0, height=1, plottype="r",
                      arrow.density=5,
                      angle=30,
                      col="black",
                      fill.col=if (plottype == "r") col else NULL,
                      lty=par("lty"),
                      lwd=par("lwd"),
                      add=FALSE,
                      #                     labels=FALSE,
                      xlim=range.feat(x),
                      ylim=c(y-height*3/4, y+height*3/4), ...) {

  if (!is.null(x$externalPtr))
    x <- from.pointer.feat(x)
  if (is.null(x$start) || is.null(x$end))
    stop("invalid features object")
  if (plottype=="a" && is.null(x$strand))
    stop("cannot plot feature arrows without strand data")

  if (!add) 
    plot(c(0), c(0), type="n", xlim=xlim, ylim=ylim, xlab="Coordinate", ylab="", ...)

  f <- x$end >= xlim[1] & x$start <= xlim[2]
  if (sum(f)==0) return(invisible(NULL))
  x <- x[f,]
  
  if (plottype=="r" || plottype=="b") {
    rect(x$start-0.5, y-height/2, x$end+0.5, y+height/2,
         col=fill.col, border=col, lty=lty,
         lwd=lwd)
  }
  if (plottype == "l")
    segments(x$start-0.5, y, x$end+0.5, col=col, lty=lty, lwd=lwd)
  if (plottype == "a" || plottype == "b") {
    usr <- par("usr")
    pin <- par("pin")
    upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L])/pin
    arrow.density <- upi[1L]/arrow.density  #convert from lines per inch
                                            # to coordinates per line
    f <- x$strand == "+"
    angle <- angle*pi/180 #convert to radians
    xlen <- height/2/tan(angle)/upi[2]*upi[1]

    if (sum(f) > 0) {
      xplus <- x[f,]
      for (i in 1:nrow.feat(xplus)) {
        start <- xplus[i,]$start-0.5
        end <- xplus[i,]$end+0.5
        xlen <- height/2/tan(angle)/upi[2]*upi[1]
        x0 <- seq(from = end + arrow.density*floor(xlen/arrow.density),
                  to = start-xlen, by=-arrow.density)
#        x0 <- seq(from=end+xlen, to=start-xlen, by=-arrow.density)
        x1 <- x0 - xlen
        y0top <- rep(y, length(x0))
        y0bottom <- rep(y, length(x0))
        y1top <- rep(y + height/2, length(x0))
        y1bottom <- rep(y - height/2, length(x0))

        f <- x0 > end
        if (sum(f) > 0L) {
          y0top[f] <- y  + height/2*(x0[f] - end)/xlen
          y0bottom[f] <- y - height/2*(x0[f]-end)/xlen
          x0[f] <- end
        }
        f <- x1 < start
        if (sum(f) > 0L) {
          y1top[f] <- y + height/2*(1.0-(start - x1[f])/xlen)
          y1bottom[f] <- y - height/2*(1.0-(start - x1[f])/xlen)
          x1[f] <- start
        }
        f <- x1 < x0
        if (sum(f) > 0L) {
          segments(x0[f], y0top[f], x1[f], y1top[f], col=col, lty=lty, lwd=lwd)
          segments(x0[f], y0bottom[f], x1[f], y1bottom[f], col=col, lty=lty, lwd=lwd)
        }
      }
    }
    
    f <- x$strand == "-"
    if (sum(f) > 0) {
      xneg <- x[f,]
      for (i in 1:nrow.feat(xneg)) {
        start <- xneg[i,]$start-0.5
        end <- xneg[i,]$end+0.5
        x0 <- seq(from=start - arrow.density*floor(xlen/arrow.density),
                  to = end+xlen, by=arrow.density)
#        x0 <- seq(from=start-xlen, to=end+xlen, by=arrow.density)
        x1 <- x0 + xlen
        y0top <- rep(y, length(x0))
        y0bottom <- rep(y, length(x0))
        y1top <- rep(y + height/2, length(x0))
        y1bottom <- rep(y - height/2, length(x0))

        f <- x1 > end
        if (sum(f) > 0L) {
          y1top[f] <- y  + height/2*(1.0-(x1[f] - end)/xlen)
          y1bottom[f] <- y - height/2*(1.0-(x1[f]-end)/xlen)
          x1[f] <- end
        }
        f <- x0 < start
        if (sum(f) > 0L) {
          y0top[f] <- y + height/2*(start - x0[f])/xlen
          y0bottom[f] <- y - height/2*(start - x0[f])/xlen
          x0[f] <- start
        }
        f <- x0 < x1
        if (sum(f) > 0L) {
          segments(x0[f], y0top[f], x1[f], y1top[f], col=col, lty=lty, lwd=lwd)
          segments(x0[f], y0bottom[f], x1[f], y1bottom[f], col=col, lty=lty, lwd=lwd)
        }
      }
    }
  }
  invisible(NULL)
}


##' make gene plot
##' @title Gene plot
##' @param x An object of type \code{feat}
##' @param y the location of the plot on the y axis
##' @param height the height of the boxes
##' @param arrow.density The density of the arrows in arrows per inch
##' @param angle angle (in degrees) of the arrow heads
##' @param col color to use for plotting
##' @param lty line type for arrows, borders, and shading
##' @param lwd line width for arrows, borders and shading
##' @param add if \code{TRUE}, add to existing plot
## @param labels whether to label the boxes.  If a vector of strings gives
## the lables for each element.  If \code{TRUE}, use x$feature for labels.
##' @param xlim A numerical vector of length 2 giving the range for the x-axis.
##' @param ylim A numerical vector of length 2 giving the range for the y-axis.
##' @param ... graphical parameters to be passed to \code{plot}.
##' @method plot gene
##' @export plot.gene
##' @export
##' @keywords features plot
##' @example inst/examples/plot-gene.R
##' @author Melissa J. Hubisz
plot.gene <- function(x, y=0, height=1,
                      arrow.density=5,
                      angle = 30,
                      col="black",
                      lty=par("lty"),
                      lwd=par("lwd"),
                      add=FALSE,
                      #                     labels=FALSE,
                      xlim=range.feat(x),
                      ylim=c(y-height*3/4, y+height*3/4), ...) {
  if (!is.null(x$externalPtr))
    x <- from.pointer.feat(x)
  if (is.null(x$start) || is.null(x$end))
    stop("invalid features object")

  seqname <- unique(x$seqname)
  if (length(seqname) > 1L)
    stop("feature has multiple sequences, can only plot 1")

  if (!add) 
    plot(c(0), c(0), type="n", xlim=xlim, ylim=ylim, ...)

  f <- x$end >= xlim[1] & x$start <= xlim[2]
  if (sum(f) == 0) return(invisible(NULL))
  x <- x[f,]

  f <- x$feature == "intron"
  if (sum(f) > 0L) 
    plot.feat(x[f,], plottype="l", y=y, lty=1, col=col, add=TRUE)
  fexon <- x$feature == "exon"
  fcds <- x$feature == "CDS"
  if (sum(fexon) > 0L) {
    if (sum(fcds) > 0L) {
      region.bounds <- feat(seqname=seqname, start=min(x$start), end=max(x$end))
      exon <- coverage.feat(x[fexon,], x[fcds,], region.bounds,
                            get.feats=TRUE,
                            not=c(FALSE, TRUE, FALSE))
    } else exon <- x[fexon,]
    if (nrow.feat(exon) > 0L)
      plot.feat(exon, plottype="a", y=y, height=height/5, add=TRUE, col=col,
                angle=angle, arrow.density=arrow.density)
  }
  if (sum(fcds) > 0L)
    plot.feat(x[fcds,], plottype="b", y=y, height=height, add=TRUE, col=col,
              angle=angle, arrow.density=arrow.density)
}



##' Features kernel density
##' @param x An object of type \code{feat}
##' @param type a character string, denoting the value to compute
##' the density for.  Currently the only valid types are "length"
##' and "score"
##' @param ... additional arguments to be passed to \code{density}
##' @return A kernel density object as defined by \code{\link{density}}
##' @export
##' @export density.feat
##' @method density feat
##' @keywords features
##' @importFrom stats density
##' @author Melissa J. Hubisz
density.feat <- function(x, type="length", ...) {
  if (type == "length") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call.rphast("rphast_gff_lengths", x$externalPtr)
    } else {
      vals <- x$end - x$start
    }
  } else if (type == "score") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call.rphast("rphast_gff_getScores", x$externalPtr)
    } else {
      vals <- x$score
    }
  } else stop("unknown type (should be \"length\" or \"score\")")
  density(vals, ...)
}

##' plot histogram of feature lengths
##' @param x an object of type \code{feat}
##' @param type a character string, denoting the value to make the histogram with.
##' Currently the only valid types are "length" or "score"
##' @param ... additional arguments to be passed to \code{hist}
##' @method hist feat
##' @export hist.feat
##' @export
##' @keywords features
##' @author Melissa J. Hubisz
hist.feat <- function(x, type="length", ...) {
  if (type == "length") {
    if (!is.null(x$externalPtr)) {
      starts <- .Call.rphast("rph_gff_starts", x$externalPtr)
      ends <- .Call.rphast("rph_gff_ends", x$externalPtr)
      vals <- ends - starts + 1
    } else {
      vals <- x$end - x$start + 1
    }
  } else if (type == "score") {
    if (!is.null(x$externalPtr)) {
      vals <- .Call.rphast("rph_gff_scores", x$externalPtr)
    } else {
      vals <- x$score
    }
  } else stop("unknown type (should be \"length\" or \"score\")")
  hist(vals, ...)
}

##' Feature overlap
##'
##' Creates a features object containing all the features from one set which overlap
##' features from another.
##' @param x An object of type \code{feat} containing features to select
##' @param filter An object of type \code{feat} which determines which elements of x to select
##' @param numbase The number of bases of overlap between x and filter required to choose
##' a record.  Use NULL to ignore (but then min.percent must be defined)
##' @param min.percent The minimum percent that a record must overlap with the combined records in filter
##' in order to be chosen
##' @param overlapping If \code{FALSE}, choose records with less than numbase overlapping bases,
##' and less than min.percent fraction overlap if min.percent is not \code{NULL}
##' @param get.fragments If \code{FALSE}, entire records are selected from x based on whether
##' they meet selection criteria.   
##' If \code{TRUE}, return only the fragments of x that overlap
##' with filter.  In this case, the same fragments may be output multiple times, if they are
##' selected by multiple entries in filter.  numbase and min.percent apply in either case.
##' When this option is used, the return value is a list with two gffs.  The first (named frags)
##' contains the overlapping fragments, and the second (filter.frags) contain the fragments
##' from filter which selected the overlapping fragments. 
##' @param pointer.only If \code{TRUE}, the return object will only be a
##' pointer to an object stored in C (useful for very large features; advanced use only).
##' @return an object of type \code{feat} containing the selected entries from x (unless
##' get.fragments==TRUE, then it returns a list with two \code{feat} objects; see
##' get.fragments).
##' @note If either x or filter are feature objects stored as a pointer to C memory,
##' then this function may reorder the elements in these objects, but leave them
##' otherwise unchanged.
##' @export
##' @keywords features
##' @example inst/examples/overlap-feat.R
##' @author Melissa J. Hubisz
overlap.feat <- function(x, filter, numbase=1, min.percent=NULL,
                         overlapping=TRUE, get.fragments=FALSE, pointer.only=FALSE) {
  check.arg(numbase, "numbase", "integer", null.OK=TRUE)
  check.arg(min.percent, "min.percent", "numeric", null.OK=TRUE)
  if (is.null(numbase) && is.null(min.percent)) stop("one of numbase or min.percent should not be NULL")
  if (!is.null(numbase) && numbase < 0) stop("numbase should be at least 1 (if it isn't NULL)")
  if (!is.null(min.percent) && (min.percent < 0 || min.percent > 1))
    stop("min.percent should be NULL or in the range (0,1)")
  check.arg(overlapping, "overlapping", "logical", null.OK=FALSE)
  check.arg(get.fragments, "get.fragments", "logical", null.OK=FALSE)

  if (is.null(x$externalPtr)) 
    x <- as.pointer.feat(x)
  if (is.null(filter$externalPtr))
    filter <- as.pointer.feat(filter)

  rv <- list()
  rv[[1]] <- .Call.rphast("rph_gff_overlapSelect",
                          x$externalPtr, filter$externalPtr,
                          numbase, min.percent, !overlapping,
                          get.fragments)
  if (get.fragments) {
    rv <- rphast.simplify.list(rv)
    if (!pointer.only) {
      rv$frags <- from.pointer.feat(rv$frags)
      rv$filter.frags <- from.pointer.feat(rv$filter.frags)
    }
    return(rv)
  }
  feat.rv <- .makeObj.feat(TRUE)
  feat.rv$externalPtr <- rv[[1]]
  if (!is.null(rv) && !pointer.only) {
    feat.rv <- from.pointer.feat(feat.rv)
  }
  if (nrow.feat(feat.rv) == 0) return(NULL)
  feat.rv
}


##' Get inverse features
##' @param x An object of type \code{feat}
##' @param region.bounds An object of type \code{feat} which defines
##' the boundaries of all relevant chromosomes in the first argument
##' @param pointer.only If \code{TRUE}, return a pointer to a structure stored
##' in C (advanced use only).
##' @return An object of type \code{feat} which contains all regions in
##' region.bounds that are not in the first argument
##' @note If x is stored as a pointer to C memory, then its elements will
##' be sorted by this function.  region.bounds will not be changed.
##' @export
##' @keywords features
##' @author Melissa J. Hubisz
inverse.feat <- function(x, region.bounds, pointer.only=FALSE) {
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  if (is.null(region.bounds$externalPtr))
    region.bounds <- as.pointer.feat(region.bounds)
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_inverse",
                                 x$externalPtr,
                                 region.bounds$externalPtr)
  if (!pointer.only)
    rv <- from.pointer.feat(rv)
  rv
}


##' Features coverage
##' @param ... objects of type \code{feat}
##' @param or if \code{TRUE}, get the coverage of union of feat arguments.
##' or is \code{FALSE} by default, which takes the intersection.
##' @param not If not \code{NULL}, a vector of logicals the same length
##' as the number of features
##' provided (or will be recycled to this length).
##' For each value which is \code{TRUE}, then any base *not* included in this
##' feature will be counted.  (The negation is done before any other operation).
##' If \code{NULL}, do not negate any features.  There must be at least one
##' feature which is not negated (so that boundaries can be established).
##' @param get.feats if \code{TRUE}, return an object of type \code{feat}
##' representing the intersection (or union of \code{or==TRUE}) of the
##' features
##' @param pointer.only (Only used if \code{get.feats==TRUE}).  If \code{TRUE},
##' the features object returned will be stored as a pointer to an object
##' in C.
##' @return The number of bases covered by the feat arguments, or the
##' combined feat object if \code{get.feats==TRUE}.
##' @note Any features object passed into this function which is stored as a
##' pointer to an object stored in C may be reordered (sorted) by this function.
##' @export
##' @keywords features
##' @example inst/examples/coverage-feat.R
##' @author Melissa J. Hubisz
coverage.feat <- function(..., or=FALSE, not=FALSE, get.feats=FALSE,
                          pointer.only=FALSE) {
  check.arg(or, "or", "logical", null.OK=FALSE)
  check.arg(get.feats, "get.feats", "logical", null.OK=FALSE)
  if (get.feats) 
    check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  featlist <- list(...)
  check.arg(not, "not", null.OK=TRUE, min.length=1L,
            max.length=length(featlist))
  if (is.null(not)) {
    not <- rep(FALSE, length(featlist))
  } else {
    not <- rep(not, length.out=length(featlist))
  }
  if (sum(not == FALSE) == 0L)
    stop("at least one feature must have not==FALSE")
  region.bounds <- featlist[[which(not==FALSE)[1]]]
  for (i in 1:length(featlist)) {
    x <- featlist[[i]]
    if (not[i]) x <- inverse.feat(x, region.bounds)
    if (is.null(x$externalPtr)) 
      x <- as.pointer.feat(x)
    featlist[[i]] <- x$externalPtr
  }
  if (get.feats) {
    rv <- .makeObj.feat(TRUE)
    rv$externalPtr <- .Call.rphast("rph_gff_featureBits", featlist,
                                   or, get.feats)
    if (pointer.only) return(rv)
    return(from.pointer.feat(rv))
  }
  .Call.rphast("rph_gff_featureBits", featlist, or, get.feats)
}


##' Add UTRs to features
##' @param x An object of type \code{feat}.  CDS regions must be present with type "CDS", and
##' the transcript_id must be indicated in the attribute field.
##' @return An object of type \code{feat}, with all the entries of the original object, but
##' also with UTR annotations.
##' @note If x is stored as a pointer to an object stored in C, then UTRs will be added
##' to x.
##' @export
##' @keywords features
##' @example inst/examples/add-UTRs-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
add.UTRs.feat <- function(x) {
  if (is.null(x$externalPtr)) {
    x <- as.pointer.feat(x)
    getDataFrame <- TRUE
  } else getDataFrame <- FALSE
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_add_UTRs", x$externalPtr)
  if (getDataFrame) return(from.pointer.feat(rv))
  rv
}

##' Add introns to features
##' @param x An object of type \code{feat}.  CDS regions must be present with type "CDS", and
##' the transcript_id must be indicated in the attribute field.
##' @return An object of type \code{feat}, with all the entries of the original object, but
##' also with intron annotations.
##' @note If x is stored as a pointer to an object stored in C, introns will be
##' added to x.
##' @export
##' @keywords features
##' @example inst/examples/add-introns-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
add.introns.feat <- function(x) {
  if (is.null(x$externalPtr)) {
    x <- as.pointer.feat(x)
    getDataFrame <- TRUE
  } else getDataFrame <- FALSE
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_add_introns", x$externalPtr)
  if (getDataFrame) return(from.pointer.feat(rv))
  rv
}

##' Add start/stop codon, 3'/5' splice signals to features
##' @param x An object of type \code{feat}.  CDS regions must be present with type "CDS", and
##' the transcript_id must be indicated in the attribute field.
##' @return An object of type \code{feat}, with all the entries of the original object, but
##' also with stop codons, start, codons, 3' splice, and 5' splice sites annotated.
##' @note \itemize{
##' \item{If x is stored as a pointer to an object stored in C, signals will be
##' added to x.}
##' \item{Does not correctly handle case of splice site in middle of start
##' or stop codon.}}
##' @export
##' @keywords features
##' @example inst/examples/add-signals-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
add.signals.feat <- function(x) {
  if (is.null(x$externalPtr)) {
    x <- as.pointer.feat(x)
    getDataFrame <- TRUE
  } else getDataFrame <- FALSE
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_add_signals", x$externalPtr)
  if (getDataFrame) return(from.pointer.feat(rv))
  rv
}


##' Fix start and stop signals
##' @param x An object of type \code{feat}.  CDS regions must be present with type
##' "CDS", and the transcript_id must be indicated in the attribute field.
##' Start and stop codons should have feature type "start_codon" and "stop_codon"
##' (as produced by addSignals.feat).
##' @return An object of type \code{feat}, in which CDS regions are ensured to
##' include start codons and exclude stop codons, as required by the GTF2 standard.
##' @note \itemize{
##' \item{If x is stored as a pointer to an object stored in C, signals will be
##' added to x.}
##' \item{Assumes at most one start_codon and at most one stop_codon per transcript.}
##' }
##' @export
##' @keywords features
##' @example inst/examples/fix-start-stop-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
fix.start.stop.feat <- function(x) {
  if (is.null(x$externalPtr)) {
    x <- as.pointer.feat(x)
    getDataFrame <- TRUE
  } else getDataFrame <- FALSE
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_fix_start_stop", x$externalPtr)
  if (getDataFrame) return(from.pointer.feat(rv))
  rv
}


##' concatenate feature objects
##' @param ... objects of type \code{feat} to be combined into a single object
##' @return An object of type \code{feat} containing entries from all
##' given features
##' @export
##' @keywords features
##' @author Melissa J. Hubisz and Adam Siepel
rbind.feat <- function(...) {
  feat <- .makeObj.feat(TRUE)
  featlist <- list(...)
  idx <- 1
  p.featlist <- list()
  for (i in 1:length(featlist)) {
    currfeat <- featlist[[i]]
    if (!is.null(currfeat) && nrow.feat(currfeat) > 0L) {
      if (is.null(currfeat$externalPtr))
        currfeat <- as.pointer.feat(currfeat)
      p.featlist[[idx]] <- currfeat$externalPtr
      idx <- idx+1
    }
  }
  if (idx == 1) return(NULL)
  feat$externalPtr <- .Call.rphast("rph_gff_append", p.featlist)
  from.pointer.feat(feat)
}

##' Split features by length
##' @param x An object of type \code{feat}
##' @param f The maximum length of features in new object.  Can be a
##' vector giving a different length for each row, or a single numeric
##' value.  Values will be recycled to the same length
##' as \code{nrow.feat(x)}.
##' @param drop A logical value saying whether to drop "left-over" elements
##' which do not have exactly length f.
##' @param start.from A character string, current valid values are "left"
##' (start split at smallest coordinate for each feature), or "right"
##' (start splitting at the last coordinate and work down).  Values will
##' be recycled to the length of \code{nrow.feat(x)}
##' @param pointer.only If \code{TRUE}, return an object which is a pointer
##' to a features object stored in C (advanced use only).
##' @param ... Currently not used (for S3 compatibility).
##' @return An object of type \code{feat} with the same features as x but
##' with all features of length > max.length broken into segments
##' (starting from the first position in feature).  The last piece
##' of each split segment may be smaller than max.length
##' @method split feat
##' @export split.feat
##' @export
##' @keywords features
##' @author Melissa J. Hubisz
split.feat <- function(x, f, drop=FALSE, start.from="left",
                       pointer.only=FALSE, ...) {
  featSize <- nrow.feat(x)
  check.arg(f, "f", "integer", null.OK=FALSE,
            min.length=1, max.length=NULL)
  max.length <- f
  check.arg(drop, "drop", "logical", null.OK=FALSE)
  check.arg(start.from, "start.from", "character", null.OK=FALSE,
            min.length=1L, max.length=nrow.feat(x))
  if (sum(start.from!="left" & start.from!="right") > 0L)
    stop("start.from invalid value (should be \"left\" or \"right\")")
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  splitFeat <- .makeObj.feat(TRUE)
  splitFeat$externalPtr <- .Call.rphast("rph_gff_split", x$externalPtr,
                                        max.length, drop, 
                                        ifelse(start.from=="left", 0, 1))
  if (!pointer.only)
    splitFeat <- from.pointer.feat(splitFeat)
  splitFeat
}


##' Sort a GFF
##' @param x An object of type \code{feat}
##' @param decreasing Set to TRUE to sort from highest to lowest coordinates
##' @param ... Currently not used
##' @return An object of type \code{feat} sorted primarily by
##' seqname, then by start position, then by end position.
##' @method sort feat
##' @export sort.feat
##' @export
##' @note If x is stored as a pointer to an object in C, the object will
##' be modified to the return value.
##' @keywords features
##' @author Melissa J. Hubisz and Adam Siepel
sort.feat <- function(x, decreasing = FALSE, ...) {
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_sort", x$externalPtr)
  rv <- from.pointer.feat(rv)
  if (decreasing) 
    rv <- rv[dim(rv)[1]:1,]
  rv
}
  

##' Composition of features with respect to annotations
##' @param x An object of type \code{feat}.
##' @param annotations An object of type \code{feat} containing
##' some annotations.
##' @return A data frame with two columns and a row for
##' each type of element in the annotations.  The second
##' column gives the fraction of 
##' x which fall in the corresponding annotation type.
##' Given non-overlapping annotations which cover the
##' entire range of interest, the second column should sum
##' to 1 (otherwise not).
##' @export
##' @note If x or annotations are passed to this function as pointers to objects
##' stored in C, they will be sorted after the function call.
##' @keywords features
##' @author Melissa J. Hubisz
composition.feat <- function(x, annotations) {
  if (!is.null(annotations$externalPtr))
    annotations <- from.pointer.feat(annotations)
  if (!is.null(x$externalPtr))
    x <- from.pointer.feat(x)
  annTypes <- unique(annotations$feature)
  rv <- list()
  for (anntype in annTypes) {
    annfeat <- annotations[annotations$feature == anntype,]
    rv[[anntype]] <- coverage.feat(x, annfeat)/coverage.feat(x)
  }
  data.frame(type=names(rv), composition=as.numeric(rv), stringsAsFactors=TRUE)
}


##' Enrichment of features with respect to annotation types
##' @param x An object of type \code{feat}
##' @param annotations An object of type \code{feat} containing
##' some annotations.
##' @param region.bounds An object of type \code{feat} representing
##' the boundary coordinates of the regions of interest (such as chromosome
##' boundaries).  All elements from the first two arguments should fall
##' entirely within region.bounds.
##' @return A data frame with two columns and a row for each type of
##' element in the annotations.  The second column gives the 
##' fold-enrichment
##' of x across the corresponding annotation type,
##' which is equal to the
##' fraction of x which fall within the annotation type,
##' divided by the fraction of the entire region covered by the
##' annotation type.
##' @export
##' @note If any of the arguments to this function are passed as pointers to
##' objects stored in C, they will be sorted after this function call.
##' @keywords features
##' @author Melissa J. Hubisz
enrichment.feat <- function(x, annotations, region.bounds) {
  if (!is.null(annotations$externalPtr))
    annotations <- from.pointer.feat(annotations)
  if (!is.null(x$externalPtr))
    x <- from.pointer.feat(x)
  annTypes <- unique(annotations$feature)
  rv <- list()
  totalNumBase <- coverage.feat(region.bounds)
  for (anntype in annTypes) {
    annfeat <- annotations[annotations$feature == anntype,]
    rv[[anntype]] <- as.numeric(coverage.feat(x, annfeat, region.bounds))*as.numeric(totalNumBase)/(as.numeric(coverage.feat(annfeat, region.bounds))*as.numeric(coverage.feat(x, region.bounds)))
  }
  data.frame(type=names(rv), enrichment=as.numeric(rv), stringsAsFactors=TRUE)
}


##' Remove overlapping genes
##' @param x An object of type \code{feat}, usually read from a genepred
##' file.  Should have attributes labelled "transcript_id" which identify
##' features belonging to the same gene.
##' @param incomparables Not currently used (present for S3 compatibility).
##' @param ... Not currently used.
##' @return An object of type \code{feat} with overlapping genes removed.
##' If the features are scored, then the feature with the highest score
##' is kept; otherwise the feature with the longest length.  If
##' x is a pointer to an object stored in C, the return value will also
##' be a pointer (and x will be altered to the return value).
##' @note \itemize{
##' \item{Long UTRs can have undesirable effects; may want to filter these
##' out first.}
##' \item{If x is a pointer to an object in C, it will be modified (to
##' the return value).}}
##' @method unique feat
##' @export unique.feat
##' @export
##' @keywords features
##' @note VERY IMPORTANT: this function is not currently implemented to look
##' at chromosomes (ie the seqname field of the feature).  Therefore any
##' genes which have overlapping coordinates REGARDLESS OF THE CHROMOSOME
##' will be pruned to a single "non-overlapping" gene.  To get around this,
##' first subset the features by chromosome and call uniq.feat on each
##' subset.
##' @note Also, this algorithm considers genes to be overlapping even if
##' they are on different strands.  If this is undesirable, then subset
##' the features by strand as well as chromosome.
##' @author Melissa J. Hubisz and Adam Siepel
unique.feat <- function(x, incomparables=FALSE, ...) {
  if (is.null(x$externalPtr)) {
    x <- as.pointer.feat(x)
    getDataFrame <- TRUE
  } else getDataFrame <- FALSE
  rv <- .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_nonOverlapping_genes", x$externalPtr)
  if (getDataFrame)
    rv <- from.pointer.feat(rv)
  rv
}


##' Extract value from tag-value formatted attribute in features object
##' @param x A features object of type \code{feat}.  The attribute field
##' should be in tag-value format (either GFF 2 standard; ie,
##' "tag1 val1a val1b; tag2 val2 ; ...", or, GFF 3 standard; ie,
##' "tag1=val1a,val1b;tag2=val2; ...".
##' where vals are in quotes if they are strings.  
##' @param tag The tag whose values are to be extracted.
##' @return If there is at most one relevant value for each feature,
##' a character vector of the same length as x will be returned, containing
##' the value for each feature, or NA where the tag does not exist for that
##' feature.  If some elements have multiple values, then the return value
##' will be a list with the same length as x, each element being a character
##' vector containing the values for the corresponding element of x (or
##' NA for no value).
##' @export
##' @keywords features GFF
##' @example inst/examples/tagval-feat.R
##' @author Melissa J. Hubisz
tagval.feat <- function(x, tag) {
  check.arg(tag, "tag", "character", null.OK=FALSE)
  if (is.null(x$externalPtr)) {
    if (is.null(x$attribute)) return (rep(NA, nrow(x)))
    x <- as.pointer.feat(x)
  }
  rv <- rphast.simplify.list(.Call.rphast("rph_gff_one_attribute",
                                          x$externalPtr, tag))
  maxlen <- max(sapply(rv, length))
  if (maxlen == 1L) rv <- as.character(rv)
    f <- rv==""
  if (sum(f) > 0L)
    rv[f] <- NA
  rv
}


##' Extract value from tag-value formatted attributes
##' @param x A vector of character strings in tag-val format (as
##' described in the GFF 2 standard; ie, "tag1 val1a val1b; tag2 val2 ; ...",
##' or in the GFF 3 format "tag1=val1a,val1b; tag2=val2; ..."),
##' where vals are in quotes if they are strings.
##' @param tag The tag whose values are to be extracted.
##' @return If there is at most one value per tag for each element of x,
##' a character vector of the same length as x will be returned, containing
##' the value for each element, or NA if the tag does not exist for that
##' element.  If some elements have multiple values, then the return value
##' will be a list with the same length as x, each element being a character
##' vector containing the values for the corresponding element of x (or
##' NA for no value).
##' @export
##' @keywords GFF
##' @example inst/examples/tagval.R
##' @author Melissa J. Hubisz
tagval <- function(x, tag) {
  check.arg(x, "x", "character", null.OK=FALSE, min.length=1L, max.length=NULL)
  tagval.feat(x=feat(start=rep(1, length(x)),
                end=rep(2, length(x)), attribute=x),
              tag=tag)
}


##' Combine adjacent features with the same "feature" field
##' @param x An object of type feat
##' @return A features object in which adjacent features are
##' combined into one longer feature.
##' @note If x is stored as a pointer to a C structure, then
##' x will be modified to the return value.
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
flatten.feat <- function(x) {
  isPointer <- !is.null(x$externalPtr)
  if (!isPointer) x <- as.pointer.feat(x)
  sort.feat(x)
  .Call.rphast("rph_gff_flatten", x$externalPtr)
  if (!isPointer) return(from.pointer.feat(x))
  x
}


##' Read a wig file
##'
##' Reads fixed or variable step wig files.  Stores them as a
##' features object.
##' @param file The file to read
##' @param pointer.only If \code{TRUE}, store as a pointer to a
##' C structure
##' @return A GFF object representing data in wig file
##' @export
read.wig <- function(file, pointer.only=FALSE) {
  read.feat(file, pointer.only=pointer.only)
}


##' Write a fixedStep wig file
##' @title Writing a wig file
##' @param chrom A character vector giving chromosome name for each point.
##' Will be recycled to length(start)
##' @param start An integer vector giving start coordinate for each point.
##' @param score A numeric vector giving score at each point
##' Will be recycled to length(start)
##' @param span An integer giving span (ie, length) of each element (all
##' elements must have the same length, so only a single value is allowed).
##' @param file The name of the file to write to (will be overwritten).
##' A value of NULL implies write to console.
##' @param append Whether to append to the file.  If FALSE, file will be
##' overwritten.
##' @keywords wig fixedStep
##' @author Melissa J. Hubisz
##' @export
##' @example inst/examples/write-wig.R
write.wig <- function(chrom, start, score, span=1, file=NULL, append=FALSE) {
  chrom <- check.arg(chrom, "chrom", "character", null.OK=FALSE, max.length=NULL)
  start <- check.arg(start, "start", "integer", null.OK=FALSE, max.length=NULL)
  score <- check.arg(score, "score", "numeric", null.OK=FALSE, max.length=NULL)
  span <- check.arg(span, "span", "integer", null.OK=FALSE)
  file <- check.arg(file, "file", "character", null.OK=TRUE)
  append <- check.arg(append, "append", "logical", null.OK=FALSE)
  if (span <= 0) stop("span should be integer >= 1")

  if (length(chrom) > length(start))
    stop("length(chrom) > length(start)")
  if (length(score) > length(start))
    stop("length(score) > length(start)")
  
  invisible(.Call.rphast("rph_wig_print",
                         feat(seqname=chrom, start=start, end=start+span-1,
                              score=score, pointer.only=TRUE)$externalPtr,
                         file,
                         append))
}


##' Write a features object in fixedStep wig format
##' @param x An object of type feat
##' @param file The name of the file to write to.  A value of NULL
##' implies write to console.
##' @param append If TRUE, append to the file.  Otherwise overwrite.
##' @keywords wig fixedStep GFF features
##' @author Melissa J. Hubisz
##' @note Wig format only contains chromosome, coordinates, and score.  Any
##' other data will be lost.
##' @note This function will quit with an error if the elements of x are
##' not all the same length (as required by fixedStep wig format).
##' @note If x is stored as a pointer to a C structure, the elements
##' will be sorted by this function.
##' @example inst/examples/write-wig-feat.R
##' @export
write.wig.feat <- function(x, file=NULL, append=FALSE) {
  append <- check.arg(append, "append", "logical", null.OK=FALSE)
  file <- check.arg(file, "file", "character", null.OK=TRUE)
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  invisible(.Call.rphast("rph_wig_print", x$externalPtr, file, append))
}
           

##' Convert coordinates from one frame of reference to another
##'
##' Converts coordinates of features in a GFF according to a multiple
##' alignment.  Will map from the coordinate system of any sequence
##' to any other sequence; can also map to or from the coordinate system
##' of the entire alignment.
##' @param x A features object; if the from parameter is -1, the first
##' column should indicate which frame of reference is used for that row
##' (ie, the species name).
##' For coordinates in the frame of reference of the entire alignment, set
##' the first column to "MSA".
##' @param align An msa object containing the alignment
##' @param from A single character string or integer, used to indicate the
##' current frame of reference of the features.  If the rows of the features are
##' not all in the same frame of reference, set this to -1 and indicate the
##' frame of reference in the 1st column of the features.  Otherwise, it can
##' be specified here as a single integer (from 0 to nrow.msa(align)), with
##' 0 indicating the frame of reference of the entire alignment, and 1
##' indicating the 1st species, 2 the second, etc.  Or it can be a single
##' character string giving the name of the species, with "MSA" indicating
##' the entire alignment.
##' @param to A single character string or integer, used to indicate the frame
##' of reference to convert to.  This is specified in the same way as the
##' "from" argument, above, except that -1 is not an option.
##' @return A features object with elements in the frame of refernece
##' indicated by the "to" argument.
##' @note Ignores any offset in MSA.  All coordinates should start with the
##' first position in the alignment as 1.
##' @note If the endpoints of an element have gaps in the "to" species,
##' the elements will be truncated
##' @example inst/examples/convert-coords-feat.R
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
convert.coords.feat <- function(x, align, from=-1, to=1) {
  x <- copy.feat(x)
  if (is.null(x$externalPtr))
    x <- as.pointer.feat(x)
  if (!is.msa(align)) stop("align should be of type msa")
  if (length(from) != 1L) {
    cat("from=", from, "\n")
    cat("length(from)=", length(from), "\n")
    stop("from should be character or integer of length 1")
  }
  is.wholenumber <-
         function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (is.character(from)) {
    if (from=="MSA") from <- 0
    else {
      if (!is.element(from, names.msa(align)))  {
        stop("argument from does not match any names in the alignment")
      }
      from <- which(names.msa(align)==from)
    }
  } else if (is.numeric(from)) {
    if (from < -1 || from > nrow.msa(align) || !is.wholenumber(from))
      stop("argument from should be in the range -1, .., ", nrow.msa(align))
  } else 
    stop("from should be a character or integer of length 1")
  # now from is an integer of length 1

  # do same for to
  if (length(to) != 1L) 
    stop("to should be character or integer of length 1")
  if (is.character(to)) {
    if (to=="MSA") to <- 0
    else {
      if (!is.element(to, names.msa(align))) 
        stop("argument to does not match any names in the alignment")
      to <- which(names.msa(align)==to)
    }
  } else if (is.numeric(to)) {
    if (to < 0 || to > nrow.msa(align) || !is.wholenumber(to))
      stop("argument to should be in the range 0, .., ", nrow.msa(align))
  } else 
    stop("to should be a character or integer of length 1")
  # now to is an integer of length 1


  
  rv <-  .makeObj.feat(TRUE)
  rv$externalPtr <- .Call.rphast("rph_gff_convert_coords",
                                 x$externalPtr,
                                 as.pointer.msa(align)$externalPtr,
                                 to)
  from.pointer.feat(rv)
}
