##
## scale.R - scale estimators 
##
## Scale estimation functions...
##
## Authors:
##   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>
##   Dipl. Stat. Karen Schettlinger <schettlinger@statistik.uni-dortmund.de>
##
## Bugs:
## 20061017 (OME):
## * The correction coefficients need to be recalculated since the
##   implementation of 'SN', 'QN' and 'MAD' changed.
##

match.scale.estimator <- function(name) {
  if(!is.character(name))
    name <- as.character(name)
  scale <- switch(name,
                  QN=rf.QN,
                  SN=rf.SN,
                  MAD=rf.MAD,
                  LSH=rf.LSH)
  if (is.null(scale)) {
    ## No match found. Throw error message.
    stop("value of 'scale' must be one of 'MAD', 'LSH', 'QN' or 'SN'")
  }
  return(scale)
}

rf.MAD <- function(d, n=length(d)) {
  if (n <= nrow(sizecorrection))
    return(mad(d, constant=sizecorrection$MAD[n]))
  else
    return(mad(d))
}

rf.SN <- function(d, n=length(d)) {
  if (n <= nrow(sizecorrection))
    return(Sn(d, constant=sizecorrection$SN[n]))
  else
    return(Sn(d))
}

rf.QN <- function(d, n=length(d)) {
  if (n <= nrow(sizecorrection))
    return(Qn(d, constant=sizecorrection$QN[n]))
  else
    return(Qn(d))
}

rf.LSH <- function(d, n=length(d)) {
  m <- ceiling((n - 1)/2)
  x <- sort(d)
  res <- min(abs(x[(m+1):n] - x[1:(n-m)]))

  if (n <= nrow(sizecorrection))
    return(res * sizecorrection$LSH[n])
  else
    return(res)
}
