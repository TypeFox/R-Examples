#' Coerces data.frames to xpssFrame objects
#'
#' Function to check if an object is an xpssFrame, or coerce it to one if possible.
#' @name as.xpssFrame
#' @usage as.xpssFrame(x)
#' is.xpssFrame(x)
#' @aliases as.xpssFrame
#' @param x any (non-empty) R object. 
#' @author Bastian Wiessner
#' @return \code{as.xpssFrame} returns an xpssFrame, with the essential attributes \code{DO_IF}, \code{FILTER}, \code{SPLIT_FILE}, \code{TEMPORORY} and \code{WEIGHTS} for the actual dataset and additional variable.label and varname attributes for every variable in the dataset. 
#' \cr \cr \code{is.xpssFrame} returns TRUE if its argument is a xpssFrame (that is, has "xpssFrame" amongst its classes) and FALSE otherwise. 
#' @examples
#' 
#' data(fromXPSS)
#' is.xpssFrame(fromXPSS) 
#' temp <- data.frame(x=1:5, y = 2:6, z=c("a","b","c","d","e"))
#' is.xpssFrame(temp)
#' temp <- as.xpssFrame(temp)
#' is.xpssFrame(temp)
#' @export


as.xpssFrame  <- function(x) {
  
  if("xpssFrame" %in% class(x))
  {
    stop("Data is already a xpssFrame object")
  }
  
  if(!(is.element("data.frame",class(x)))){
    x <- as.data.frame(x)
  }
  
  class(x) <- c("xpssFrame","data.frame")
  attr(x,"FILTER") <- FALSE
  attr(x,"TEMPORARY") <- FALSE  
  attr(x,"SPLIT_FILE") <- FALSE
  attr(x,"DO_IF") <- FALSE
  attr(x,"SELECT_IF") <- FALSE
  attr(x,"WEIGHTS") <- FALSE
  
  
  
  
  for(i in 1:length(x)) {
    attr(x[[i]],"varname") <- paste(names(x[i]),sep="")  
    attr(x[[i]],"variable.label") <- paste(names(x[i]),sep="") 
    if(is.element("factor",lapply(x,class))){
      var <- names(which(lapply(x,class) == "factor"))
      for(i in 1:length(var)) {
        eval(parse(text =paste("x$",var[[i]]," <- as.character(x$",var[[i]],")",sep="")))  
      }
    }
  }
  
  
  return(x)
}

#' @rdname as.xpssFrame
#' @aliases as.xpssFrame
#' @keywords internal
#' @export

is.xpssFrame  <- function(x) {
  if(is.element("xpssFrame",class(x)) && (!is.null(attributes(x)$FILTER)) && (!is.null(attributes(x)$TEMPORARY)) && (!is.null(attributes(x)$SPLIT_FILE)) && (!is.null(attributes(x)$DO_IF)) && (!is.null(attributes(x)$SELECT_IF)) && (!is.null(attributes(x)$WEIGHTS))) {
    return(T)
  } else{
    return(F)
  }
}
