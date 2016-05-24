#' Returns the upper bound value for a given boxnumber
#' @param boxnumber An integer boxnumber 
#' @param boxsize single numeric value, used as the boxsize
#' @param log TRUE, if logarithmic scales should be used
box2upper <- function(boxnumber, boxsize=1, log=FALSE) {
  # FIXME improve sanity checks
  if(sum(floor(boxnumber)!=boxnumber,na.rm=T)>0)
    stop("Error: Only integer values allowed as boxnumber!")
  if (log) {
    myexp <- function(x){
      exp(x)
    }
  } else {
    myexp <- function(x){
      x
    }
  }
  
  myexp((boxnumber+1)*boxsize)
}
