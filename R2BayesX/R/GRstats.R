GRstats <- function(object, term = NULL, ...)
{
  if((n <- length(object)) < 2L)
    stop("at least two models ar needed for calculation!")
  return(coda::gelman.diag(samples(object, model = NULL, term, coda = TRUE), ...))
}

