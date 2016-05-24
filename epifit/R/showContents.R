##' Show variable information in data.frame.
##'
##' This function shows contents in data.frame object, and typical use is to list modes of variables in data.frame. In addition to variable mode, a variety of functions can be used such as \code{min}, \code{max} and \code{mean} for summarizing variable information in data.frame.
##' @title Show variable information in data.frame
##' @param data a data.frame whose contents to be shown
##' @param functions a vector of functions to be applied to each variable in data.frame. Default added function is \code{mode}, which is almost same functionality as storage.mode except that factor is separated from integer.
##' @param na.action a character which specifies missing-data filter function. This is applied to data.frame supplied by data argument. Default is \sQuote{na.pass}.
##' @param na.omit.for.variable a logical value specifying whether to apply na.omit for each variable. This option is useful to supply NA free data to each function. When set to TRUE (default), the number of observation may differ for each variable, and it is recommended to check the number of obsevations by show.nobs option.
##' @param show.nobs a logical valye specifying whether to add the number of observations for each function as nobs column in result object.
##' @return a data.frame which contains the results of each function.
##' @seealso \code{\link{countNA}}
##' @examples
##' dat <- data.frame(a = rnorm(100))
##' dat <- transform(dat, b = ifelse(a > -2, a, NA))
##' dat <- transform(dat, c = ifelse(a > 0, "positive", "negative"))
##' showContents(dat)
##' showContents(dat, functions=c(modes, storage.mode, min, max, median))
##' @export
showContents <- function(data=NULL, functions=c(modes), na.action="na.pass", na.omit.for.variable=TRUE, show.nobs=TRUE){

  if(!is.data.frame(data)){
    stop("data argument must be data.frame")
  }

  nfunc <- length(functions)

  if(nfunc < 1){
    stop("No function is specified")
  }

  funcname <- match.call()
  if(length(funcname) < 3){
    funcname <- "modes"
  } else {
    funcname <- funcname[[3]]
    funcname <- as.character(funcname[2:length(funcname)])
  }

  varlist <- colnames(data)

  nobs <- numeric(length(varlist))

  if(length(varlist) < 1){
    stop("data is empty")
  }

  if (na.action == "na.pass") {
    data <- na.pass(data)
  }
  else if (na.action == "na.fail") {
    data <- na.fail(data)
  }
  else if (na.action == "na.omit") {
    data <- na.omit(data)
  }
  else if (na.action == "na.exclude") {
    data <- na.exclude(data)
  }
  else {
    stop("invalid na.action specification")
  }

  # each list element stores for result from each function
  result <- list()

  # i = 1 as special case
  if(na.omit.for.variable){
    content <- na.omit(data[[varlist[1]]])
  } else {
    content <- data[[varlist[1]]]
  }

  nobs[1] <- length(content)
  
  for(j in 1:nfunc){
    res <- NA
    tryCatch(
      res <- functions[[j]](content),
      error = function(e) {warning("Function (", funcname[j], ") is not applicable in some variables")}
      )
    
    if(length(res) > 1){
      warning("Element-wise function (", funcname[j], ") cannot be used")
      result[[j]] <- NA
    } else {
      result[[j]] <- res
    }
  }

  if(length(varlist) > 1){
    for(i in 2:length(varlist)){
      if(na.omit.for.variable){
        content <- na.omit(data[[varlist[i]]])
      } else {
        content <- data[[varlist[i]]]
      }

      nobs[i] <- length(content)

      for(j in 1:nfunc){
        res <- NA
        tryCatch(
          res <- functions[[j]](content),
          error = function(e) {}
          )
        
        if(length(res) > 1){
          result[[j]] <- c(result[[j]], NA)
        } else {
          result[[j]] <- c(result[[j]], res)
        }
      }
    }
  }
  
  if(show.nobs){
    result <- data.frame(nobs, as.data.frame(result))
    colnames(result) <- c("nobs", funcname)
  } else {
    result <- as.data.frame(result)
    colnames(result) <- funcname
  }
  rownames(result) <- varlist
  return(result)
}
