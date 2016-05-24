#####################################################
## indirect_sampling.r
##
## indirect sampling estimators, including
## multiplicity
##
## TODO -- handle NAs in code below (this is partially done)
## TODO -- make handling of column vars, etc, uniformly handled
##         when passed in as args. (ie, use get.var; now some fns
##         still don't use get.var)
## TODO -- perhaps an easy to use interface for nsum estimates
##         with bootstrap (and then use this in nsum.internal.validation)
## TODO -- when using defaults (for example, taking
##         popn size info from dataframe attributes,
##         should we print out a message to the screen?
##         or perhaps have a (default) verbose mode?)
## TODO -- think about code to get 45q15 from these data...

#####################################################
##' multiplicity.estimator
##'
##' compute multiplicity estimate of the population
##' size
##'
##' @param survey.data the dataframe with the survey results
##' @param mult.col the name or index of the column that contains, for each
##'             respondent, the individual value of the number known divided
##'             by the sum of the multiplicities (TODO MORE DETAIL)
##' @return the multiplicity estimate of the hidden
##'         population's size (as a prevalence)
##' @export
multiplicity.estimator <- function(survey.data,
                                   mult.col="mult")
{

  mult.vals <- try(subset(survey.data,
                          select=mult.col),
                   silent=TRUE)

  if( inherits(mult.vals, "try-error") ||
      ncol(mult.vals) != 1 ||
      ! is.numeric(mult.vals[,1]) ) {
    stop(paste(mult.col,
               "does not identify a valid column for computing multiplicities.\n"))
  }

  mult.vals <- as.numeric(mult.vals[,1])

  res <- mean(mult.vals)

  return(res)

}

#####################################################
##' indirect estimator
##' (generalized weight share method / gwsm)
##'
##' compute gwsm estimate of the population
##' size
##'
##' @param survey.data the dataframe with the survey results
##' @param gwsm.col the name or index of the column that contains, for each
##'             respondent, the individual value of the number known divided
##'             by the sum of the multiplicities (TODO MORE DETAIL)
##' @return the multiplicity estimate of the hidden
##'         population's size (as a prevalence)
##' @export
gwsm.estimator <- function(survey.data,
                           gwsm.col="mult")
{

  stop("gwsm.estimator is not yet implemented.")
  
  ## TODO -- this is just a skeleton; the stuff below has not
  ## been edited yet. think about the best way to set this up
  mult.vals <- try(subset(survey.data,
                          select=gwsm.col),
                   silent=TRUE)

  if( inherits(mult.vals, "try-error") ||
      ncol(mult.vals) != 1 ||
      ! is.numeric(mult.vals[,1]) ) {
    stop(paste(gwsm.col,
               "does not identify a valid column for computing multiplicities.\n"))
  }

  mult.vals <- as.numeric(mult.vals[,1])

  res <- mean(mult.vals)

  return(res)

}
