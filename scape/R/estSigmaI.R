estSigmaI <- function(model, what="s", series=NULL, init=NULL, FUN=mean, p=1, digits=2)
{
  ## 1  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  what <- match.arg(what, c("c","s"))
  x <- if(what=="c") model$CPUE else model$Survey
  if(is.null(x))
    stop("Element '", what, "' not found")

  ## 2  Extract series
  if(is.null(series))
    series <- unique(x$Series)
  if(length(series) > 1)
  {
    output <- lapply(series, function(s)
                     estSigmaI(model=model, what=what, series=s, init=init, FUN=FUN, p=p, digits=digits))
    names(output) <- series
  }
  else
  {
    ok.series <- x$Series %in% series; if(!any(ok.series)) stop("Please check if the 'series' argument is correct")
    x <- x[!is.na(x$Obs) & ok.series,]

    ## 3  Calculate sigmahat
    eps <- log(x$Obs) - log(x$Fit)
    rss <- sum(eps^2)
    n <- nrow(x)
    sigmahat <- sqrt(rss/(n-p))

    ## 4  Prepare init
    years <- unique(x$Year)
    if(is.null(init))
      init <- tapply(x$CV, x$Year, identity)  # use array to behave like estN
    else if(class(init) == "scape")
    {
      component <- if(what=="c") init$CPUE else init$Survey
      init <- component[!is.na(component$Obs),]
      init <- tapply(init$CV, init$Year, identity)
    }
    else if(identical(init,FALSE))
      init <- 1
    if(length(init)!=1 && length(init)!=length(years))
      stop("'init' should be NULL, model, vector, FALSE, or 1")

    ## 5  Prepare output
    output <- as.numeric(sigmahat * init / FUN(init))
    if(!is.null(digits))
      output <- round(output, digits=digits)
    if(length(output) == length(years))
      names(output) <- years
  }

  return(output)
}
