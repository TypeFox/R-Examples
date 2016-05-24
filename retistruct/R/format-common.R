read.scale <- function(dataset) {
  ## If there is a scale file, read it
  scale <- NA
  scfile <- file.path(dataset, "scale.csv")
  if (file.exists(scfile)) {
    print("Reading scale file")
    sc <- read.csv(file.path(dataset, "scale.csv"))
    scale <- sc[1,1]
    if (!is.numeric(scale)) {
      stop("Scale file has not been read correctly. Check it is in the correct format.")
    }
  } else {
    warning("Scale file does not exist. Scale bar will not be set.")
  }
  return(scale)
}

read.image <- function(dataset) {
  im <- NULL
  imfile <- file.path(dataset, "image.png")
  if (file.exists(imfile)) {
    message("Reading image")
    im <- as.raster(png::readPNG(imfile))
  }
  return(im)
}

##' Read data points from a file \code{dataponts.csv} in the directory
##' \code{dataset}. The CSV should contain two columns for every
##' dataset. Each pair of columns must contain a unique name in the
##' first cell of the first row  and a valid colour in the second
##' cell of the first row. In the remaining rows, the X coordinates of
##' data points should be in the first column and the Y coordinates
##' should be in the second column.
##'
##' @title Read data points in CSV format
##' @param dataset Path to directory containing \code{dataponts.csv}
##' @return List containing
##' \item{\code{Ds}}{List of sets of datapoints. Each set comprises a 2-column matrix and each set is named.}
##' \item{\code{cols}}{List of colours for each dataset. There is one element that corresponds to each element of \code{Ds} and which bears the same name.}
##' @author David Sterratt
read.datapoints <- function(dataset) {
  datfile <- file.path(dataset, "datapoints.csv")
  Ds <- list()
  cols <- c()
  if (file.exists(datfile)) {
    message("Reading datapoints")
    ## Read file. stringsAsFactors=FALSE prevents conversion to factors
    dat <- read.csv(file.path(datfile), stringsAsFactors=FALSE)

    ## Go through pairs of columns
    while(ncol(dat) >= 2) {
      ## Extract first two columns
      d <- dat[,1:2]
      dat <- dat[,-(1:2)]
      names <- colnames(d)

      ## Convert strings to numeric. Suppress warnings as sapply
      ## complains about coercion to NA
      suppressWarnings({d <- sapply(d, as.numeric, USE.NAMES=FALSE)})
      ## Force conversion to matrix, necessary when the data has only
      ## one row
      d <- matrix(d, ncol=2)
      
      ## Any strings (e.g. empty ones) that don't convert will be
      ## converted to NA. Get rid of these.
      d <- na.omit(d)
      attr(d, "na.action") <- NULL

      ## Add to lists with appropriate names
      
      D <- list(d)
      names(D) <- names[1]
      Ds <- c(Ds, D)

      col <- list(names[2])
      names(col) <- names[1]
      cols <- c(cols, col)
    }
  }
  return(list(Ds=Ds, cols=cols))
}
