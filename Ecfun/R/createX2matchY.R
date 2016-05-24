createX2matchY <- function(x, y){
  cix <- classIndex(x)
  ciy <- classIndex(y)
  ciz <- max(cix, ciy)
  if(ciz<2) return(NULL)
  chz <- index2class(ciz)
  do.call(chz, list(length(y)))
}

