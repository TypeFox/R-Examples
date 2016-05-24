.packageName <- "rv"


postsim <- function(fit)
{
  warning("'postsim' is defunct. Use 'posterior' instead.")
  UseMethod('postsim')
} 

