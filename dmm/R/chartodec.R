chartodec <-
function(cvec,ndec){
#  chartodec()  -  char vector cvec to decimal number (or NA) vector dvec
#      ie - use function on one trait (= one col of df)
  if(!is.character(cvec)){
    stop("chartodec: argument cvec must be of type character\n")
  }
  dvec <- as.numeric(cvec)/(10^ndec)
  return(dvec)
}
