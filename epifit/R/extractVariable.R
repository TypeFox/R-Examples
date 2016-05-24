##' Extract variables according to mode from data.frame.
##'
##' This function extract variables which match specified mode from data.frame, and make a new data frame.
##' @param data a data.frame from which numeric variables are extracted.
##' @param mode a character specifying object type. Object modes of \sQuote{numeric}, \sQuote{character}, \sQuote{factor}, and \sQuote{logical} are supported.
##' @return a data.frame which includes only specified mode of variables.
##' @examples
##' df <- data.frame(id=seq(1,10), str=letters[1:10], fac=factor(seq(1,10)), stringsAsFactors=FALSE)
##' extractVariable(df)
##' extractVariable(df, mode="character")
##' extractVariable(df, mode="factor")
##' @export
extractVariable <- function(data=NULL, mode="numeric"){

  if(is.null(data)||!is.data.frame(data))
    stop("data is not specified or not data.frame")

  idx <- rep(TRUE, ncol(data))
  
  funcname <- paste("is.", mode, sep="")
  if(exists(funcname, mode="function", envir=.BaseNamespaceEnv)){
    func <- get(funcname, mode="function", envir=.BaseNamespaceEnv)
  } else {
    stop("Invalid mode")
  }
  
  for(i in 1:ncol(data)){
    if(!func(data[,i])){
      idx[i] = FALSE
    }
  }
  
  return(data[,idx,drop=FALSE])
}
