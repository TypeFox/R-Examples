#' @title Variance Ratio
#' @description Computes the ratio of the variance of aggregate species abundances
#' in a community to the sum of the variances of individual, component species. A
#' variance ratio = 1 indicates that species do not covary,  a variance ratio > 1
#' indicates predominately positive covariance among species and a variance
#' ratio < 1 indicates predominately negative covariance (Schluter 1984).
#'
#' Includes a cyclic shift null modeling option to test if the variance ratio significantly
#' differs from 1. The null community is created by randomly selecting different
#' starting points for each species' time series, which generates a community in
#' which species abundances vary independently but within-species autocorrelation
#' is maintained (Hallett et al. 2014). This randomization is repeated a user-specific
#' number of times and confidence intervals are reported for the resultant null
#' distribution of variance ratios. If the dataframe includes multiple replicates,
#' the variance ratios for the actual and null communities are averaged within each
#' iteration unless specified otherwise.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param bootnumber The number of null model iterations used to calculated confidence intervals
#' @param level The confidence level for the null mean
#' @param replicate.var The name of the (optional) replicate column
#' @param average.replicates If true returns the variance ratio and CIs averaged
#' @param li (deprecated) lower confidence interval
#' @param ui (deprecated) upper confidence interval
#' across replicates; if false returns the variance ratio and CI for each replicate. Defaults to true.
#' @return The variance_ratio function returns a data frame with the following attributes:
#' \itemize{
#'  \item{VR: }{A numeric column with the actual variance ratio value.}
#'  \item{lowerCI: }{A numeric column with the lowest confidence interval value.}
#'  \item{upperCI: }{A numeric column with the highest confidence interval value.}
#'  \item{nullmean: }{A numeric column with the average null variance ratio value.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance;
#' time.var, species.var and abundance.var are used to indicate which columns
#' contain those variables. If multiple replicates are included in the data frame,
#' that column should be specified with replicate.var. Each replicate should
#' reflect a single experimental unit - there must be a single abundance value
#' per species within each time point and replicate.
#'
#' Null model confidence intervals default to the standard lowest 2.5\% and
#' upper 97.5\% of the null distribution, typically these do not need to be change,
#' but they can be user-modified to set more stringent CIs.
#'  @references
#'  Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins,
#'  Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. (2014)
#'  "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient."
#'  Ecology 95, no. 6: 1693-1700. doi: 10.1890/13-0895.1
#'
#'  Schluter, Dolph. (1984) "A Variance Test for Detecting Species Associations,
#'  with Some Example Applications." Ecology 65, no. 3: 998-1005. doi:10.2307/1938071.
#' @importFrom stats aggregate confint
#' @export
#' @examples
#'  data(knz_001d)
#'
#'  # Calculate the variance ratio and CIs averaged within replicates
#'  # Here the null model is repeated once, for final use it is recommended to set a
#'  # large bootnumber (eg, 10000)
#'
#'  res_averagedreplicates <- variance_ratio(knz_001d,
#'                time.var = "year",
#'                species.var = "species",
#'                abundance.var = "abundance",
#'                bootnumber = 1,
#'                replicate = "subplot")
#'
#'  #Calculate the variance ratio and CIs for each replicate
#'
#'  res_withinreplicates <- variance_ratio(knz_001d,
#'                time.var = "year",
#'                species.var = "species",
#'                abundance.var = "abundance",
#'                bootnumber = 1,
#'                replicate = "subplot",
#'                average.replicates = FALSE)
variance_ratio <- function(df, time.var,
                           species.var,
                           abundance.var,
                           bootnumber,
                           replicate.var = NA,
                           average.replicates = TRUE,
                           level = 0.95, li, ui) {

  ## check for use of li, ui

  if ((!missing(li) | !missing(ui))) {
    warning("argument li and ui are deprecated; please use level instead.",
            call. = FALSE)
  }

  # check to make sure abundance is numeric data
  check_numeric(df, time.var, abundance.var)

  # if no replicates, calculate single variance ratio
  if  (is.na(replicate.var)) {

    ## check the structure of the data
    check_single_onerep(df, time.var, species.var)

    ## calculate observed variance ratio
    VR <- variance_ratio_longformdata(df, time.var, species.var, abundance.var)

    ## null models
    nullval <- cyclic_shift(df, time.var = time.var,
                            species.var = species.var,
                            abundance.var = abundance.var,
                            FUN = variance_ratio_matrixdata,
                            bootnumber = bootnumber)

    nullout <- confint(nullval)

    # bind actual value and CI
    output <- cbind(nullout, VR)

  } else {

    ## drop any unused levels
    df <- droplevels(df)

    # if multiple replicates, check all replicates have values
    check_single(df, time.var, species.var, replicate.var)

    # calculate average variance ratio across replicates
    if (average.replicates == TRUE) {
      ## check the data
      check_multispp(df, species.var, replicate.var)

      # split the dataset into replicates, according to replicate.var
      df <- df[order(df[[replicate.var]]),]
      X <- split(df, df[replicate.var])

      ## observed variance ratio
      VR <- mean(unlist(lapply(X, FUN = variance_ratio_longformdata, time.var, species.var, abundance.var)))

      ## Use cyclic_shift to calculate the null distribution
      nullval <- cyclic_shift(df = df, time.var = time.var,
                              species.var = species.var,
                              abundance.var = abundance.var,
                              replicate.var = replicate.var,
                              FUN = variance_ratio_matrixdata,
                              bootnumber = bootnumber)

      nullout <- confint(nullval)

      output <- cbind(nullout, VR)

    } else {

      # calculate the variance ratio for replicate
      check_multispp(df, species.var, replicate.var)
      df <- df[order(df[[replicate.var]]),]
      X <- split(df, df[replicate.var])

      ## workaround, necessary because you are not allowed to pass an argument called FUN to lapply
      cyclic_shift_nofun <- function(f = variance_ratio_matrixdata){
        function(...) {
          cyclic_shift(FUN = f, ...)
        }
      }

      ## also do separate null models (ie on split list)
      null_list <- lapply(X = X, FUN = cyclic_shift_nofun(),
                          time.var = time.var,
                          species.var = species.var,
                          abundance.var = abundance.var,
                          replicate.var = NA,
                          bootnumber = bootnumber)

      ## calculate confidence intervals for each
      null_intervals <- lapply(null_list, confint)

      ## combine as data.frames and preserve replicate var
      repnames <- lapply(names(null_intervals), as.data.frame)

      ## combine replicate names with intervals, then combine them all
      nullout <- do.call("rbind", Map(cbind, repnames, null_intervals))

      ## use the replicate.var name
      names(nullout)[1] <- replicate.var

      VR_list <- lapply(X, FUN = variance_ratio_longformdata,
                        time.var, species.var, abundance.var)

      ## combine as data.frames and preserve replicate var
      VRnames <- lapply(names(VR_list), as.data.frame)
      VR_df <- lapply(VR_list, as.data.frame)

      ## combine replicate names with intervals, then combine them all
      VRout <- do.call("rbind", Map(cbind, VRnames, VR_df))

      names(VRout) <- c(replicate.var, "VR")

      output <- merge(nullout, VRout, by = replicate.var)

    }
  }

  ## result
  row.names(output) <- NULL
  return(output)
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' A function to calculate the variance ratio
#'
#' @param comdat A community dataframe
#' @return var.ratio The variance ratio of the community
variance_ratio_matrixdata <- function(comdat){
    check_sppvar(comdat)
    all.cov <- stats::cov(comdat, use = "pairwise.complete.obs")
    col.var <- apply(comdat, 2, stats::var)
    com.var <- sum(all.cov)
    pop.var <- sum(col.var)
    var.ratio <- com.var/pop.var
    return(var.ratio)
}

#' A function to calculate the variance ratio from a longform dataframe
#'
#' @param df A dataframe containing time.var, replicate.var, species.var and abundance.var columns
#' @param time.var The name of the time.var column from df
#' @param species.var The name of the species.var column from df
#' @param abundance.var The name of the abundance.var column from df
#' @return var.ratio The variance ratio of the community
variance_ratio_longformdata <- function(df, time.var, species.var, abundance.var){
    com.use <- transpose_community(df, time.var, species.var, abundance.var)
    var.ratio <- variance_ratio_matrixdata(com.use)
    return(var.ratio)
}
