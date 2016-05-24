
# currently "outliers" is not used by the functions in the package
# but it may be useful to define the input of some functions, 
# e.g. argument "mo" in functions "outliers.effects" and "outliers.regressors"

outliers <- function(type, ind, weight = NULL)
{
  if (is.null(weight))
    weight <- rep(1, length(type))

  m <- data.frame(type = factor(type, levels = c("IO", "AO", "LS", "TC", "SLS")), 
    ind = ind, coefhat = weight)

  #structure(m, class = "outliers")
  m
}
