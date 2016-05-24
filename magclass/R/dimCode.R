
dimCode <- function(dim,x,missing=0){
  #function to translate dim to dim code
  if(is.character(dim)) {
    dnames <- dim
    dim <- match(dim,getSets(x),nomatch=0)
    if(length(dim)>length(dnames)) stop('One or more elements were found more than once in x!')
    names(dim) <- dnames
    
    #translate sub-datadimensions to 3.1, 3.2,...
    dim[dim>=3] <- 3 + (dim[dim>=3]-2)/10
  }
  if(any(dim>=4) | any(dim<1)) {
    if(missing=="stop") stop("illegal dimension. Use either dimension 1, 2, or 3, or if you want to address subdimensions use 3.1, 3.2, ...")
    dim[dim>=4 | dim<1] <- missing
  }
  return(dim)
}