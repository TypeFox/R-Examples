##' Constructor for a \code{dataset} object.
##'
##' @title Constructor for a \code{dataset} object.
##' @param o An \code{outline} object.
##' @param dataset The name of the dataset
##' @param Ds A list of data point sets, with each set being a 2
##' column matrix of X and Y coordinates of data point locations. Each
##' item in the list should be named. Elements with these names should
##' also be in the \code{cols} argument (see below).
##' @param Ss A list of landmarks. These do not need to be named. If
##' any elements are  named, the names should map onto an element in
##' the \code{cols} argument. Any elements that are named \code{""}
##' will be plotted using the default colour.
##' @param cols A list of colours in which to plot datapoints and landmarks.
##' @param raw A place to put raw data in whatever format is desired.
##' @param Gs  A list of grouped point sets, with each set being a 3
##' column matrix of X and Y coordinates and the value Z of the
##' variable at that point.  Each item in the list should be
##' named. Elements with these names should also be in the \code{cols}
##' argument.
##' @return A \code{dataset} object.
##' @author David Sterratt
Dataset <- function(o, dataset, Ds, Ss, cols, raw, Gs=NULL) {
  d <- o
  class(d) <- addClass("dataset", o)
  d$dataset <- dataset
  d$Ds <- Ds
  d$Ss <- Ss
  if (is.null(names(Ss))) {
    names(d$Ss) <- rep("", length(d$Ss))
  }
  d$cols <- cols
  d$raw <- raw
  d$Gs <- Gs
  return(d)
}

##' @export
nameLandmark.dataset <- function(d, i, name) {
  if (!is.na(i)) {
    new.names <- names(d$Ss)
    ## If this name already exists, replace it with ""
    j <- getLandmarkID(d, name)
    if (!is.na(j)) {
      new.names[j] <- ""
    }
    new.names[i] <- name
    names(d$Ss) <- new.names
  }
  return(d)
}

getLandmarkID <- function(d, name) {
  id <- which(names(d$Ss) == name)
  if (length(id) == 1) {
    return(id)
  } else {
    return(NA)
  }
}

##' @title Get IDs of groups of data within a dataset
##' @param r \code{\link{Dataset}} object
##' @return Array of IDs
##' @author David Sterratt
##' @method getIDs dataset
##' @export
getIDs.dataset <- function(r) {
  return(names(r$Ds))
}

##' @title Flat plot of Dataset
##' @param x \code{\link{Dataset}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param datapoints If \code{TRUE}, display data points.
##' @param grouped If \code{TRUE}, dipslay grouped data.
##' @param landmarks If \code{TRUE}, dipslay landmarks.
##' @param ids IDs of groups of data within a dataset, returned using
##' \code{\link{getIDs}}.
##' @param ... Graphical parameters to pass to plotting functions
##' @method flatplot dataset
##' @author David Sterratt
##' @export
flatplot.dataset <- function(x, axt="n", ylim=NULL,
                             datapoints=TRUE,
                             grouped=FALSE,
                             landmarks=TRUE,
                             ids=getIDs(x),
                             ...) {
  ## This will call projection.reconstructedOutline(), and hence
  ## Outline(), but without drawing a grid.  The grid will be drawn
  ## later, after all the data has appeared.
  NextMethod(grid=FALSE)

  if (datapoints) {
    with(x, {
      for (id in ids) {
        if (!is.null(Ds[[id]])) {
          points(Ds[[id]][,1], Ds[[id]][,2],
                 col=cols[[id]], pch=20)
        }
      }
    })
  }
  
  if (grouped) {
    Gs <- x$Gs
    for (id in ids) {
      if (!is.null(Gs[[id]])) {
        if (nrow(Gs[[id]]) > 0) {
          text(Gs[[id]][,1], Gs[[id]][,2], Gs[[id]][,3],
               col=x$cols[[id]])
        }
      }
    }
  }

  if (landmarks) {
    with(x, {
      if (length(Ss) > 0) {
        for (i in 1:length(Ss)) {
          name <- names(Ss)[i]
          col <- ifelse(is.null(name) || (name==""), "default", name)
          lines(Ss[[i]][,1], Ss[[i]][,2],
                col=cols[[col]])
        }
      }
    })
  }
  ## This will call flatplot.reconstructedOutline() and will draw a
  ## grid but not add an image. Thus the grid appears over the data.
  NextMethod(add=TRUE,
             image=FALSE)
}

all.equal.dataset <- function(target, current, ...) {
  return((all.equal(target$Ds, current$Ds) == TRUE) &
         (all.equal(target$Gs, current$Gs) == TRUE) &
         (all.equal(target$Ss, current$Ss) == TRUE))
}
