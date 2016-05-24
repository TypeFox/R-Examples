getN <- function(model, what="CAc", series=NULL, digits=NULL)
{
  ## 1  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  what <- match.arg(what, c("CAc","CAs","CLc","CLs"))
  x <- model[[what]]
  if(is.null(x))
    stop("Element '", what, "' not found")

  ## 2  Extract series
  if(is.null(series))
    series <- unique(x$Series)
  if(length(series) > 1)
  {
    output <- lapply(series, function(s) getN(model=model, what=what, series=s, digits=digits))
    names(output) <- series
  }
  else
  {
    ok.series <- x$Series %in% series; if(!any(ok.series)) stop("Please check if the 'series' argument is correct")
    x <- x[!is.na(x$Obs) & ok.series,]

    ## 3  Create output
    output <- tapply(x$SS, x$Year, mean)
    output <- drop(as.matrix(output))
    if(!is.null(digits))
      output <- round(output, digits=digits)
  }

  return(output)
}
