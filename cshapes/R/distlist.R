distlist <- function(date, type="mindist", tolerance=0.1, useGW=T) {

  # check input
  if (!inherits(date, "Date")) {
    stop("date is not of type Date")
  }
  
  if (date < as.Date("1946-1-1") | date > as.Date("2015-6-30")) {
    stop("Specified date is out of range")
  }
  
  if (!(type %in% c("mindist", "capdist", "centdist"))) {
  	stop("Wrong type argument. Possible values: mindist, capdist, centdist")
  }
  
  if (tolerance<0) {
  	stop("Tolerance must be >=0")
  }
  
  # minimum distance
  if (type=="mindist") {
  	dmat <- distmatrix(date, type="mindist", tolerance, useGW)
  }
      
  # capital distance
  if (type=="capdist") {
  	dmat <- distmatrix(date, type="capdist", useGW=useGW)
  }
  
  # centroid distance
  if (type=="centdist") {
  	dmat <- distmatrix(date, type="centdist", useGW=useGW) 
  }
 	
  dimension <- ncol(dmat)
  dimlabels <- colnames(dmat)
  
  ccode1 <- as.vector(sapply(dimlabels, function (x) rep(x, dimension)))
  ccode2 <- rep(dimlabels, dimension)
  resultframe <- data.frame(cbind(as.numeric(ccode1), as.numeric(ccode2), as.numeric(dmat)))
  colnames(resultframe) <- c("ccode1", "ccode2", type)
  resultframe
}

