`addMissing` <-
function(x, na.proportion=0, parent.cols=c(1,2), seed)
{
  
  ## Description: Adds missing data to objects of class autoMarker or
  ## autoCross as specified.

  ## Arguments:
  ## x : object of class simAutoMarker or simAutoCross, or a matrix with
  ##     dominant markers scored as 0 or 1 
  ## na.proportion: proportion missing at random or a list
  ##                with two components indiv and marker each containing
  ##                c(prop. markers missing, prop. missing)
  ##                (Default: 0)
  ## parent.cols:  columns containing parental markers (etc) not altered
  ##               only used if object of class simAutoCross
  ## seed:     random number generator (RNG) state for random number
  ##           which will be set at start to reproduce results

  ## Values:
  ## x: returns object of class autoMarker or autoCross, or a matrix with
  ##    dominant markers scored as 0 or 1 with extra component na.proportion
  ##    which has the following elements
  ##    na.proportion: proportion missing at random or a list
  ##                   with two components indiv and marker each containing
  ##                   c(prop. markers missing, prop. missing)
  ##    time.generated: time/date when data set generated + when missing added
  ##    seed:  random number generator seed which could be used to
  ##           reproduce results (I hope)
  ##    call:  matches arguments when function called

  if (!missing(seed)) {
    set.seed(seed)
  }

  ## allow for simAutoMarkers, simAutoCross and matrix x

  if (class(x) == "simAutoMarkers" | class(x) == "simAutoCross") {
    if (class(x) == "simAutoMarkers") {
      markers <- x$markers
    } else {
      markers <- x$markers[,-parent.cols]
    }
  } else {
    if (is.matrix(x)) {
      markers <- x
    } else {
      stop("x should be a matrix or of class simAutoMarkers or simAutoCross")
    }
  }
  n.individuals <- ncol(markers)
  n.markers <- nrow(markers)
      
  ## drop markers if na.proportion set
  
  if (mode(na.proportion) == "numeric") {
    if (na.proportion != 0) {
      if (na.proportion < 0 | na.proportion > 1)
        stop("Error: na.proportion should be between 0 and 1")
      ## choose missing and set to NA
      s <- sample( 1:length(markers), na.proportion*length(markers))
      markers[cbind(row(markers)[s],col(markers)[s])] <- NA
    }
  } else {
    if (mode(na.proportion) == "list") {
      s.rows <- sample( 1:n.markers , na.proportion$marker[1]*n.markers)
      s.cols <- sample( 1:n.individuals , na.proportion$indiv[1]*n.individuals)
      for (i in s.rows){
        markers[i, sample(1:n.individuals,
                          na.proportion$indiv[2]*n.individuals)] <- NA
      }
      for (i in s.cols){
        markers[sample(1:n.markers,
                       na.proportion$marker[2]*n.markers), i] <- NA
      }
    }
  }
  
  if (class(x) == "simAutoMarkers" | class(x) == "simAutoCross") {
    res <- x
    if (class(x) == "simAutoMarkers") {
      res$markers <- markers
      res$seg.ratios <- segregationRatios(markers)
    } else {
      res$markers[,-parent.cols] <- markers
      res$seg.ratios <- segregationRatios(res$markers)
      parType <- rep( 1:3, times=x$no.parType)
      res$p10$markers <- markers[parType==1,]
      res$p01$markers <- markers[parType==2,]
      res$p11$markers <- markers[parType==3,]
      res$p10$seg.ratios <- segregationRatios(res$p10$markers)
      res$p01$seg.ratios <- segregationRatios(res$p01$markers)
      res$p11$seg.ratios <- segregationRatios(res$p11$markers)
    }
    res$na.proportion <- list(na.proportion=na.proportion,
                              time.generated <- date(),call=match.call())
    if (!missing(seed)) res$na.proportion$seed <- seed
  } else {
    res <- markers
  }

  return(res)
}

