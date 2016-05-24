estN.int <- function(P, Phat)
{
  numerator <- sapply(rownames(P), function(t) sum(Phat[t,]*(1-Phat[t,])))
  denominator <- sapply(rownames(P), function(t) sum((P[t,]-Phat[t,])^2))
  nhat <- as.numeric(numerator / denominator)

  return(nhat)
}



estN <- function(model, what="CAc", series=NULL, init=NULL, FUN=mean, ceiling=Inf, digits=0)
{
  ## 1  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  what <- match.arg(what, c("CAc","CAs","CLc","CLs"))
  x <- model[[what]]
  if(is.null(x))
    stop("Element '", what, "' not found")
  x$Column <- if(substring(what,1,2)=="CA") x$Age else x$Length

  ## 2  Extract series
  if(is.null(series))
    series <- unique(x$Series)
  if(length(series) > 1)
  {
    output <- lapply(series, function(s)
                     estN(model=model, what=what, series=s, init=init, FUN=FUN, ceiling=ceiling, digits=digits))
    names(output) <- series
  }
  else
  {
    ok.series <- x$Series %in% series; if(!any(ok.series)) stop("Please check if the 'series' argument is correct")
    x <- x[!is.na(x$Obs) & ok.series,]

    ## 3  Calculate nhat
    P <- tapply(x$Obs, list(x$Year,x$Column), mean)     # pool sexes
    Phat <- tapply(x$Fit, list(x$Year,x$Column), mean)  # pool sexes
    nhat <- estN.int(P, Phat)

    ## 4  Prepare init
    years <- unique(x$Year)
    if(is.null(init))
      init <- tapply(x$SS, x$Year, unique)
    else if(class(init) == "scape")
      init <- tapply(init[[what]]$SS, init[[what]]$Year, unique)
    else if(identical(init,FALSE))
      init <- nhat
    if(length(init)!=1 && length(init)!=length(years))
      stop("'init' should be NULL, model, vector, FALSE, or 1")

    ## 5  Prepare output
    output <- as.numeric(FUN(nhat) * init / FUN(init))
    output <- pmin(output, ceiling)
    if(!is.null(digits))
      output <- round(output, digits=digits)
    if(length(output) == length(years))
      names(output) <- years
  }

  return(output)
}
