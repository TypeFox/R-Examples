
magpiesort <- function(x) {
  if(!is.magpie(x)) stop("Input is not a MAgPIE object!")
  if(any(dim(x)==0)) return(x)
  if(dim(x)[1]==1) {
    spatial_order <- 1
  } else if(length(grep("\\.[0-9]*$",dimnames(x)[[1]]))==dim(x)[1]) {
    spatial_order <- order(as.numeric(gsub("^[A-Z]+\\.","",dimnames(x)[[1]])))
  } else {
    spatial_order <- order(dimnames(x)[[1]])
  }
  if(!is.null(dimnames(x)[[2]])) temporal_order <- order(dimnames(x)[[2]])
  else temporal_order <- 1:dim(x)[2]
  return(x[spatial_order,temporal_order,])
}