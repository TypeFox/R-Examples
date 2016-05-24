`addMisclass` <-
function(x, misclass=0, bands.missed=0, parents=FALSE, parent.cols=c(1,2), seed)
{
  
  ## Description: Adds missing data to objects of class autoMarker or
  ## autoCross as specified.

  ## Arguments:
  ## x : object of class simAutoMarker or simAutoCross, or a matrix with
  ##     dominant markers scored as 0 or 1 
  ## misclass:      proportion misclassified specified as for na.proportion
  ##                (Default: 0)
  ## bands.missed: proportion of bands that are not scored when they are actually
  ##               present. Note this is applied to correctly specified markers after
  ##               markers are misclassified (Default: 0)
  ## parents:      if TRUE then misclassify parental alleles, otherwise
  ##               misclassify offspring marker alleles
  ##               (Default: FALSE)
  ## parent.cols: for object of simAutoClass the columns containg
  ##              parental markers
  ## seed:     random number generator (RNG) state for random number
  ##           which will be set at start to reproduce results

  ## Values:
  ## x: returns object of class autoMarker or autoCross, or a matrix with
  ##    dominant markers scored as 0 or 1 with extra components
  ##    misclass.info: list with five elements
  ##                proportion: numeric proportion misclassified
  ##                index: indicates which markers were set as misclassified
  ##                bands.proportion: numeric proportion of marker bands missed
  ##                bands.index: indicates which markers bands were missed
  ##                call: matches arguments when function called
  ##                time.generated: time/date when misclassifieds added
  ##                seed: seed for random number generation


  if (!missing(seed)) {
    set.seed(seed)
  }

  ## allow for simAutoMarkers, simAutoCross and matrix x

  if (class(x) == "simAutoMarkers" | class(x) == "simAutoCross" ) {
    if (class(x) == "simAutoMarkers") {
      markers <- x$markers
    } else {
      if (parents) {
        markers <- x$markers[,parent.cols]
      } else {
        markers <- x$markers[,-parent.cols]
      }
    }
  } else {
    if (is.matrix(x)) {
      markers <- x
    } else {
      stop("x should be a matrix or of class simAutoMarkers or simAutoCross")
    }
  }

  ## add in measurement error/misclassification if set

  if (misclass != 0) {
    if (misclass < 0 | misclass > 1)
      stop("Error: 'misclass' should be between 0 and 1")
    ## choose misclassified 
    s <- sample( 1:length(markers), misclass*length(markers))
    ## swap them
    misclass.index <- cbind(row(markers)[s],col(markers)[s])
    markers[misclass.index] <- (markers[misclass.index] == 0) + 0
    misclass.info <- list(proportion=misclass, index=misclass.index)
  } else {
    misclass.info <- list(proportion=misclass)
  }

  ## add in bands.missed to markers correctly classified as 1's

  if (bands.missed != 0) {
    if (bands.missed < 0 | bands.missed > 1)
      stop("Error: 'bands.missed' should be between 0 and 1")
    m <- markers
    m[misclass.index] <- 3
    ok <- setdiff(1:length(markers),s)  # those elements not misclassified
    index.ok <- cbind(row(markers)[ok],col(markers)[ok])
    ones <- ok[markers[index.ok]==1]    # ok markers that re 1
    ## choose bands.missed
    sb <- sample( ones, bands.missed*length(ones))
    bands.index <- cbind(row(markers)[sb],col(markers)[sb])
    ## set them to zero
    markers[bands.index] <- 0
    misclass.info$bands.index <- bands.index
    misclass.info$bands.proportion=bands.missed
  } else {
    misclass.info <- list(bands.proportion=bands.missed)
  }
  

  ## set up markers for simAutoMarkers or simAutoCross if necessary
  
  if (class(x) == "simAutoMarkers" | class(x) == "simAutoCross") {
    res <- x
    if (class(x) == "simAutoMarkers") {
      res$markers <- markers
      res$seg.ratios <- segregationRatios(markers)
      res$misclass.info <- misclass.info
    } else {
      parType <- rep( 1:3, times=x$no.parType)
      if (parents) {
        res$markers[,parent.cols] <- markers
        res$p01$parent.markers <- markers[parType==1,]
        res$p10$parent.markers <- markers[parType==2,]
        res$p11$parent.markers <- markers[parType==3,]
        res$parent.misclass.info <- misclass.info
      } else {
        res$markers[,-parent.cols] <- markers
        res$p01$markers <- markers[parType==1,]
        res$p10$markers <- markers[parType==2,]
        res$p11$markers <- markers[parType==3,]
        res$p10$seg.ratios <- segregationRatios(res$p10$markers)
        res$p01$seg.ratios <- segregationRatios(res$p01$markers)
        res$p11$seg.ratios <- segregationRatios(res$p11$markers)
        res$misclass.info <- misclass.info
      }
    }
    if (!missing(seed)) res$misclass.info$seed <- seed
    res$misclass.info$call <- match.call()
    res$misclass.info$time.generated <- date()
    
  } else {
    res <- markers
  }
  return(res)
}

