# 
# vcovPC.R
# function to produce panel-corrected variance covariance matrices for use in coeftest
# 

vcovPC <- function(x, ...){
  UseMethod("vcovPC")
}

vcovPC.default <- function(x, groupN, groupT, pairwise=FALSE, ...){
  mf <- match.call()
  if (is.null(groupN) | is.null(groupT)){
    stop("You must provide groupN or groupT. Call: ", mf, ".")
  }
  pc <- pcse(object=x, groupN=groupN, groupT=groupT, pairwise=pairwise)
  return(pc$vcov)
}
