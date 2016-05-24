ijroi.checkDatadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.roi")))
}

##' Read a retinal dataset in IJROI format. Each dataset is a folder
##' containing a file called \code{outline.roi} that specifies the
##' outline in X-Y coordinates. It may also contain a file
##' \code{datapoints.csv}, containing the locations of data points;
##' see \code{\link{read.datapoints}} for the format of this file. The
##' folder may also contain a file \code{od.roi} specifying the
##' coordinates of the optic disc.
##' 
##' @title Read a retinal dataset in IJROI format
##' @param dataset Path to directory containing \code{outline.roi}
##' @return A \code{\link{RetinalDataset}} object
##' @author David Sterratt
ijroi.read.dataset <- function(dataset) {
  ## Read the raw data
  roi <- RImageJROI::read.ijroi(file.path(dataset, "outline.roi"))
  out <- roi$coords

  ## Read scale
  scale <- read.scale(dataset)
  
  ## If there is an image, read it
  im <- read.image(dataset)

  ## ImageJ ROI format plots has the coordinate (0, 0) in the top
  ## left.  We have the coordinate (0, 0) in the bottom left. We need
  ## to transform P so that the outline appears in the correct
  ## orientation.
  P <- out
  offset <- ifelse(is.null(im), max(P[,2]), nrow(im))
  P[,2] <- offset - P[,2]
  
  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour. There are no datapoints in this
  ## format, but we may have landmarks.
  Ds <- list()
  cols <- list(OD="blue",
               default="orange")
  dat <- read.datapoints(dataset)
  Ds <- c(Ds, dat$Ds)
  cols <- c(cols, dat$cols)
  Ds <- lapply(Ds, function(P) {cbind(P[,1], offset - P[,2])})

  Ss <- list()

  ## Read in an Optic Disc. FIXME: this should actually be marked as
  ## the Optic Disc
  od.file <- file.path(dataset, "od.roi")
  if (file.exists(od.file)) {
    roi <- RImageJROI::read.ijroi(od.file)
    out <-  roi$coords
    out[,2] <- offset - out[,2]
    Ss[["OD"]] <- out
  }
  
  ## Create forward and backward pointers
  o <- Outline(P, scale, im)
  o <- simplify.outline(o)
  
  ## Check that P is more-or-less closed
  ## if (vecnorm(P[1,] - P[nrow(P),]) > (d.close * diff(range(P[,1])))) {
  ##    stop("Unable to find a closed outline.")
  ## }

  d <- Dataset(o, dataset, Ds, Ss, cols=cols, raw=list(outline=out))
  a <- AnnotatedOutline(d)
  a <- RetinalDataset(a)
  return(a)
}

