
setYears <- function(object,nm=NULL) {
  getYears(object) <- nm
  return(object)
}