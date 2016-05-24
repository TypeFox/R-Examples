##' List incompatible values when converted into numeric values.
##'
##' List incompatible values when converted into numeric values. Character or factor variables in data.frame are scanned to check whether the values in it is compatible, and incompatible values are returned as list. This function is intended for check before use of \code{convertFromFactor}.
##' @title List incompatible values when converted into numeric values.
##' @param data a data.frame to search for incomatible values when converted into numeric values
##' @return a list which contains incompatible values for each variable in data.frame when converted into numeric variable.
##' @seealso \code{\link{convertFromFactor}}, \code{\link{showContents}}.
##' @examples
##' a <- factor(rnorm(5))
##' b <- c("a", "b", "c", "d", "e")
##' c <- c("1", "2", "3", "4", NA)
##' d <- c("1", "2", "3", "4", ".")
##' dat <- data.frame(a,b,c,d)
##' listNumericIncompatibility(dat)
##' @export
listNumericIncompatibility <- function(data=NULL){
  
  if(!is.data.frame(data))
    stop("data argument must be data.frame")

  varlist <- colnames(data)

  result <- list()

  for(var in varlist){

    content <- data[[var]]

    if(is.factor(content)){
      level <- levels(content)
      for(i in 1:length(level)){
        tryCatch(
          {as.numeric(level[i])},
          warning=function(e){
            result[[var]] <<- c(result[[var]], level[i])
          }
          )
      }
    } else if(is.character(content)){
      level <- levels(as.factor(content))
      for(i in 1:length(level)){
        tryCatch(
          {as.numeric(level[i])},
          warning=function(e){
            result[[var]] <<- c(result[[var]], level[i])
          }
          )
      }
    }
  }
  return(result)
}
