unwrap <- function(x,sep=".") {
  if(!is.magpie(x)) stop("Input is not a MAgPIE object. unwrap works only for MAgPIE objects")
  dim <- fulldim(x,sep)
  if(any(duplicated(getNames(x)))) stop("Malformed MAgPIE object. Duplictaed names detected!")
  if(prod(dim[[1]])!=prod(dim(x))) stop("Malformed MAgPIE object. Different number of entries in original and unwrapped object! (prod(dim(in))!=prod(dim(out)))")
  reorder <- dimnames(wrap(array(NA,dim[[1]],dim[[2]]),list(1,2,NA),sep=sep))[[3]]
  if(!is.null(reorder)) x <- x[,,reorder]
  return(array(as.vector(x),dim[[1]],dim[[2]]))
}




