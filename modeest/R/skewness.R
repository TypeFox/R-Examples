# Author : authors of package 'fBasics'
# Modifications : P. Poncet
skewness <-
function (x, ...) 
{
  UseMethod("skewness")
}

skewness.default <-
function(x,                                       # sample (the data)
         na.rm = FALSE,                           # should missing values be removed ?
         method = c("moment", "fisher", "bickel"),# method
         M = shorth(x),                           # mode, if Bickel's method is choosen
         ...)
{
########################################################
# Modified 'skewness' function, from package 'fBasics'
# This version now includes Bickel's measure of skewness
########################################################

  method <- match.arg(tolower(method), c("moment", "fisher", "bickel"))
  if (!is.numeric(x)) {
    stop("argument 'x' is must be numeric")
  }
  if (na.rm) x <- x[!is.na(x)]
  nx <- length(x)
  if (is.integer(x)) x <- as.numeric(x)
  if (method == "moment") {
    skewness <- sum((x - mean(x))^3/sqrt(var(x))^3)/nx
  }
  if (method == "fisher") {
    if (nx < 3) skewness <- NA
    else skewness <- ((sqrt(nx * (nx - 1))/(nx - 2)) * (sum(x^3)/nx))/((sum(x^2)/nx)^(3/2))
  }
  if (method == "bickel") {
    cdf.value <- (length(x[x < M]) + length(x[x == M])/2)/nx
    skewness <- 1-2*cdf.value
  }
  attr(skewness, "method") <- method
  return(skewness)
}
