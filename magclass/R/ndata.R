ndata <- function(x,fulldim=FALSE) {
  if(fulldim==FALSE){
    return(dim(x)[3])
  } else {
    return(fulldim(x)[[1]][3])
  }
}