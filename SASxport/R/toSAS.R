
#' @export

toSAS <- function(x, format, format.info=NULL)
  UseMethod("toSAS")

#' @export

toSAS.numeric <- function(x, format=SASformat(x), format.info=NULL)
  {
    retval <- as.numeric(x)
    attr(retval, "SASformat")=format
    retval
  }

#' @export

toSAS.logical <- function(x, format=SASformat(x), format.info=NULL)
  {
    retval <- as.character(x)
    attr(retval, "SASformat")=format
    retval
  }
  
#' @export

toSAS.character <- function(x, format=SASformat(x), format.info=NULL)
  {
    retval <- as.character(x)
    attr(retval, "SASformat")=format
    retval
  }

#' @export

toSAS.factor <- function(x, format=SASformat(x), format.info=NULL)
  {
    finfo <- process.formats(format.info)
    if( (length(format>0)) && (format %in% names(finfo)) )
      {
        labels <- finfo[[SASformat(x)]]$label
        values <- finfo[[SASformat(x)]]$value
        retval <- values[match( x, labels)]
      }
    else
      {
        retval <- as.character(x)
      }
    attr(retval, "SASformat")=format
    retval
  }

#' @export

toSAS.POSIXt <- function( x, format="DATETIME16.", format.info=NULL)
  {
    sasBaseSeconds <- as.numeric(ISOdatetime(1960,1,1,0,0,0))
    retval <- unclass(as.POSIXct(x))  - sasBaseSeconds  # sasBaseSeconds is negative
    attr(retval,"SASformat") <- format
    retval
  }

#' @export

toSAS.Date <- function(x, format="DATE9.", format.info=NULL )
  {
    sasBase <- as.Date(strptime("01/01/1960", "%m/%d/%Y", tz="GMT")) # days
    retval <- as.numeric( as.Date(x) - sasBase)
    attr(retval, "SASformat") <- format
    retval
  }

#' @export

toSAS.default <- function(x, format=SASformat(x), format.info=NULL)
  {
    retval <- as.character(x)
    attr(retval, "SASformat") <- format
    retval
  }

#' @export

toSAS.chron <- toSAS.POSIXt
