getYears <- function(x,as.integer=FALSE) {
  if(as.integer) {
    return(as.integer(substring(dimnames(x)[[2]],2)))  
  } else {
    return(dimnames(x)[[2]])
  }
}

"getYears<-" <- function(x,value) {
  if(!is.null(value)) if(length(value)!=nyears(x)) stop("Wrong number of years supplied!")
  if(nyears(x)==0) return(x)
  if(is.null(value) & nyears(x)!=1) stop("Setting years to NULL is not possible as the number of years is not 1!")
  if(is.null(value)) {
    tmp <- list(NULL,NULL,NULL)
    if(!is.null(dimnames(x)[[1]])) tmp[[1]] <- dimnames(x)[[1]]
    if(!is.null(dimnames(x)[[3]])) tmp[[3]] <- dimnames(x)[[3]]
    names(tmp) <- names(dimnames(x))
    dimnames(x) <- tmp
  } else {
    if(all(is.numeric(value))) value <- gsub(" ","0",format(value,width=4))
    if(all(nchar(value)==4)) value <- paste("y",value,sep="")
    if(any(nchar(value)!=5) | any(substr(value,1,1)!="y")) stop("Wrong year format. Please supply either integer values or years in the format y0000!")
    dimnames(x)[[2]] <- value
  }
  return(x)
}
