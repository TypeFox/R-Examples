"s" <-
  function(x, df = 4, spar = 1)
{
  scall <- deparse(sys.call())
  if(missing(df)){
    if(!missing(spar))df<-0
  }
    
  if(ncol(as.matrix(x)) > 1)
    stop(paste(
               "The default smoother is bivariate; you gave a matrix as an argument in ",
               scall, "\n"))
  if(!is.null(levels(x))) {
    if(inherits(x, "ordered"))
      x <- as.numeric(x)
    else stop("unordered factors cannot be used as smoothing variables"
              )
  }
  attr(x, "spar") <- spar
  attr(x, "df") <- df
  real.call <- substitute(gam.s(data[[scall]], z, w, spar = spar, df = df
                                ))
  attr(x, "call") <- real.call
  attr(x, "class") <- "smooth"
  a <- is.na(x)
  if(any(a))
    attr(x, "NAs") <- seq(along = x)[a]
  x
}
