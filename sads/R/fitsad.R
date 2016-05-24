fitsad <- function(x, sad=c("bs","gamma","geom","lnorm","ls","mzsm","nbinom","pareto","poilog","power","volkov", "weibull"), ...){ 
  dots <- list(...)
  sad <- match.arg(sad)
  fit <- get(paste("fit", sad, sep=""), mode = "function")
  if(!"trunc" %in% names(dots) && (sad %in% c("poilog", "geom", "nbinom"))) {
    #BUGFIX: some objects (ie, those created with zeroes=T) have some species with 0 individuals
    # we remove those species before fitting to avoid zero-truncated fits throwing an error
    if (sum(x==0) > 0) {
      warning("Removing zeroes from abundance vector")
      x <- x[x > 0]
  }
}
  do.call(fit, c(list(x = x), dots))
}
