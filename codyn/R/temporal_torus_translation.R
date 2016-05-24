#' @title Cyclic Shift Permutations
#' @aliases temporal_torus_translation
#' @description Performs a user-specified function on a null ecological community created via cyclic shift permutations (Harms et al. 2001, Hallett et al. 2014).
#' The null community is formed by randomly selected different starting years for the time series of each species.
#' This generates a null community matrix in which species abundances vary independently but within-species autocorrelation is maintained.
#' The user-specified function must require a species x time input.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the replicate column. Defaults to \code{NA}, indicating no replicates (i.e., data are from a single plot).
#' @param FUN A function to calculate on the null community
#' @param bootnumber Number of null simulations to run
#' @return The cyclic_shift function returns an S3 object of class "cyclic_shift" and param "out".
#' The length of the "out" param  is the number of null iterations as specified by bootnumber.
#' If multiple replicates are specified, null values are averaged among replicates for each interation, but a different cyclic shift permutation is applied to each replicate within an interation.
#' @details The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' @examples
#' # Calculate a covariance matrix on a null community
#' data(knz_001d)
#' a1_cyclic <- cyclic_shift(subset(knz_001d, subplot == "A_1"),
#'                    time.var = "year",
#'                    species.var = "species",
#'                    abundance.var = "abundance",
#'                    FUN = cov,
#'                    bootnumber = 10)
#' @references
#' Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient." Ecology 95, no. 6 (2014): 1693-1700.
#'
#' Harms, Kyle E., Richard Condit, Stephen P. Hubbell, and Robin B. Foster. "Habitat Associations of Trees and Shrubs in a 50-Ha Neotropical Forest Plot." Journal of Ecology 89, no. 6 (2001): 947-59.
#' @importFrom stats confint
#' @export
cyclic_shift <- function(df, time.var,
                         species.var,
                         abundance.var,
                         replicate.var=NA,
                         FUN,
                         bootnumber){

  ## match arg to check method name
  ## switch to go from character to function

  # check time and abundance are numeric
  check_numeric(df, time.var, abundance.var)

  if (is.na(replicate.var)) {

    # check there unique species x time combinations
    check_single_onerep(df, time.var = time.var, species.var = species.var)

    # cyclic shift on one rep
    out <- cyclic_shift_onerep(df = df, time.var = time.var,
                               species.var = species.var,
                               abundance.var = abundance.var, FUN = FUN,
                               bootnumber = bootnumber)
  } else {

    # remove unused levels if replicate.var is a factor
    df[replicate.var] <- if (is.factor(df[[replicate.var]]) == TRUE) {
      factor(df[[replicate.var]])
    } else {
      df[replicate.var]
    }

    # check there unique species x time x replicate combinations
    check_single(df, time.var, species.var, replicate.var)

    # sort and apply cyclic shift to all reps
    df <- df[order(df[[replicate.var]]),]
    X <- split(df, df[replicate.var])
    lout <- lapply(X, cyclic_shift_onerep, time.var,
                   species.var, abundance.var, FUN, bootnumber)
    out <- do.call("rbind", lout)

    # take the mean value of each bootnumber interation across reps
    out <- colMeans(out)
  }

  # assign the output shift as an S3 object
  shift <- structure(list(out = out), class = "cyclic_shift")

  return(shift)
}


#' @title Confidence Intervals from a Cyclic Shift Permutation
#' @description Calculates confidence intervals for the S3 object produced by \code{cyclic_shift}
#' @param object An object of class \code{cyclic_shift}
#' @param parm which parameter is to be given a confidence interval. At present there is only one option: the mean of the null distribution. Defaults to "out", referring to the null distribution in objects of class \code{cyclic_shift}.
#' @param level the confidence level required.
#' @param ... further arguments to \code{quantile}
#' @return A dataframe with the following columns:
#' \itemize{
#'  \item{lowerCI: }{A numeric column with the lowest confidence interval value.}
#'  \item{upperCI: }{A numeric column with the highest confidence interval value.}
#'  \item{nullMean: }{A numeric column with the average value of the specified test statistic when calculated on a null community.}
#' }
#' @examples
#' # Calculate a covariance matrix on a null community
#' data(knz_001d)
#' a1_cyclic <- cyclic_shift(subset(knz_001d, subplot == "A_1"),
#'                    time.var = "year",
#'                    species.var = "species",
#'                    abundance.var = "abundance",
#'                    FUN = cov,
#'                    bootnumber = 10)
#'
#' # Return CI on a1_cyclic
#' confint(a1_cyclic)
#' @references
#' Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient." Ecology 95, no. 6 (2014): 1693-1700.
#'
#' Harms, Kyle E., Richard Condit, Stephen P. Hubbell, and Robin B. Foster. "Habitat Associations of Trees and Shrubs in a 50-Ha Neotropical Forest Plot." Journal of Ecology 89, no. 6 (2001): 947-59.
#' @export
confint.cyclic_shift <- function(object,
                                 parm = "out",
                                 level = 0.95, ...){

  # set lower and upper confidence intervals
  li <- (1 - level)/2
  ui <- 1 - li

  out <- object[[parm]]

  # calculate lower, upper CI and the null mean
  lowerCI <- stats::quantile(out, li, ...)
  upperCI <- stats::quantile(out, ui, ...)
  nullmean <- mean(out)
  output <- data.frame(lowerCI, upperCI, nullmean)
  row.names(output) <- NULL

  # results
  return(output)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################



#' A function to generate a community dataframe with a random start time for each species
#'
#' @param comdat A community dataframe
#' @return rand.comdat A randomized community dataframe
#' @importFrom permute shuffleSeries
shuffle_community <- function(comdat){

  # create empty matrix with same number of rows, columns as comdat
  nr <- nrow(comdat)
  nc <- ncol(comdat)
  rand.comdat <- matrix(NA, nrow = nr, ncol = nc)
  rand.start <- sample.int(nr, nc, replace = TRUE)

  # fill in the matrix with cyclic-shifted time series
  for (i in seq_len(nc)) {
    rand.comdat[, i] <- permute::shuffleSeries(comdat[,i])
  }
  rand.comdat <- as.data.frame(rand.comdat)
  names(rand.comdat) <- names(comdat)
  row.names(rand.comdat) <- row.names(comdat)

  # return the null community generated via cyclic shifts
  return(rand.comdat)
}


#' A function to calculate a non-S3 cyclic shift on one replicate
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param FUN A function to calculate on the null community
#' @param bootnumber The number of null model iterations returned
#' @return out A vector of  test statistics calculated on the null community
cyclic_shift_onerep <- function(df,
                                time.var,
                                species.var,
                                abundance.var,
                                FUN,
                                bootnumber){

  # check time and abundance are numeric
  check_numeric(df, time.var, abundance.var)

  # calculate the test statistic specified by FUN on null communities as many times specified by bootnumber
  out <- replicate(bootnumber, FUN(shuffle_community(transpose_community(df, time.var,  species.var, abundance.var))))

  # results
  return(out)
}

temporal_torus_translation <- function(df, time.var="year",
                                       species.var="species",
                                       abundance.var="abundance", FUN){
  .Deprecated("cyclic_shift")

  cyclic_shift(df = df, time.var = time.var,
                species.var = species.var,
                abundance.var = abundance.var,
                replicate.var = NA,
                FUN = FUN,
                bootnumber = 1)
}
