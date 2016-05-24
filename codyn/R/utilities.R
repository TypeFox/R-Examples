#' Convert from a longform abundance dataframe to a time by species dataframe.
#'
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return A dataframe of species abundances x time
transpose_community <- function(df, time.var, 
                                species.var, 
                                abundance.var) {
  df <- as.data.frame(df)
  
  # remove unused levels if species is a factor
  df[species.var] <- if(is.factor(df[[species.var]]) == TRUE){factor(df[[species.var]])} else {df[species.var]}
  
  # sort by time and species
  df <- df[order(df[[time.var]], df[[species.var]]),]
 
  # cast as a species x time dataframe; NAs to 0s
  comdat <- tapply(df[[abundance.var]], list(df[[time.var]], as.vector(df[[species.var]])), sum)
  comdat[is.na(comdat)] <- 0
  comdat <- as.data.frame(comdat)
  
  # results
  return(comdat)
}

#' check names of data frames
#'
#' @param given Vector of variable names as supplied by user
#' @param data Data frame containing variables
check_names <- function(given, data) {
    for (i in given) {
        assertthat::assert_that(assertthat::has_name(data, i))
    }
}

#' Utility function to warn users that either multiple records exist within replicates, or that data may be spanning mutiple replicates but no replicate.var has been specified
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
check_single_onerep <- function(df, time.var, species.var){
  if(max(table(df[[time.var]], df[[species.var]]))>1) warning("Either data span multiple replicates with no replicate.var specified or multiple records within years for some species") }

#' Utility function to ensure only a single record exists for a given species within one replicate, for one time point.
#' @param df A dataframe containing time.var, species.var, and replicate.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param replicate.var The name of the replicate column from df

check_single <- function(df, time.var, species.var, replicate.var){
  X <- split(df, df[[replicate.var]])
  checksingle <- lapply(X, FUN = function(xx) apply(table(xx[[species.var]], xx[[time.var]]), 2, function(x) any(x>1)))
  reptest <- unlist(lapply(checksingle, any))
  yrtest <- lapply(checksingle, which)

  if(any(unlist(checksingle))){
    if(length(names(reptest)[which(reptest)]) == 1){

    stop(paste("In replicate", names(reptest)[which(reptest)], "there is more than one record for species at the time point", unlist(lapply(yrtest, names))))
    }
      else  {
        toprint <- unlist(lapply(yrtest, names))
    stop("For the following replicates in the following time points, there are more than one records for species: \n", paste(names(toprint), collapse = "\t"), "\n", paste(toprint, collapse = "\t"))
      }
  }
}

#' Utility to check for numeric abundance and time variables
#'
#' @param df A dataframe containing time.var, species.var, and replicate.var columns
#' @param time.var The name of the time column from df
#' @param abundance.var The name of the replicate column from df

check_numeric <- function(df, time.var, abundance.var) {
  if(!is.numeric(df[[time.var]])) { stop("Time variable is not numeric") }
  if(!is.numeric(df[[abundance.var]])) { stop("Abundance variable is not numeric") }
  }

#' Utility function to stop calculations if only one species occurs in at least one replicate
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param species.var The name of the species column from df
#' @param replicate.var The name of the replicate column from df
check_multispp <- function(df, species.var, replicate.var){

  df <- droplevels(df)

  spptable <- table(df[[species.var]], df[[replicate.var]])
  if (min(colSums(spptable > 1)) < 2) {
    stop("One or more replicates consists of only a single species;
       please remove these replicates prior to calculations ")
  }
}

#' Utility function to stop calculations if the species never change in a replicate
#' @param comdat A community dataframe
#' @importFrom stats var
check_sppvar <- function(comdat){
  sppvar <- sum(apply(comdat, 2, var))
  if(sppvar == 0)
    stop("One or more replicates consist of species that never vary;
         please remove these replicates before calculation")}

