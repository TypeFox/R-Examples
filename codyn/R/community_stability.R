#' @title Community Stability
#' @description Calculates the stability of the overall community over time as the temporal mean / temporal standard deviation of aggregate species abundances (Tilman 1999).
#'
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the optional replicate column
#' @return The community_stability function returns a numeric stability value unless a replication column is specified in the input data frame.
#' If replication is specified, the function returns a data frame with the following attributes:
#' \itemize{
#'  \item{stability: }{A numeric column with the stability values.}
#'  \item{replicate.var: }{A column that shares the same name and type as the replicate.var column in the input data frame.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there should be a single community represented within each time point and replicate.
#' @references
#' Tilman, D. "The Ecological Consequences of Changes in Biodiversity: A Search for General Principles." Ecology 80, no. 5 (July 1999): 1455-74. doi:10.1890/0012-9658(1999)080[1455:TECOCI]2.0.CO;2.
#' @examples
#' data(knz_001d)
#' community_stability(knz_001d[knz_001d$subplot=="A_1",], 
#'                      time.var = "year", 
#'                      abundance.var = "abundance") # for one subplot
#' community_stability(knz_001d,
#'                      time.var = "year", 
#'                      abundance.var = "abundance",
#'                      replicate.var = "subplot") # across all subplots
#' @importFrom stats aggregate as.formula sd
#' @export

community_stability <- function(df, time.var,
                                abundance.var,
                                replicate.var = NA){
  
  # check time and abundance are numeric
  check_numeric(df, time.var, abundance.var)

  # remove zeros  
  df <- df[which(df[[abundance.var]] > 0),]

  if (is.na(replicate.var)) {

    #sum abundance within a year
    aggform <- as.formula(paste(abundance.var, "~", time.var, sep = ""))
    data2 <- stats::aggregate(aggform, data = df, sum)
    output <- stability_onerep(data2, abundance.var)

  } else {

    df <- df[order(df[[replicate.var]]),]
    df <- droplevels(df)

    # sum abundance within a replicate and year, sort
    aggform <- as.formula(paste(abundance.var, "~", replicate.var, "+", time.var, sep = ""))
    data2 <- aggregate(aggform, data = df, sum)
    data2 <- data2[order(data2[[replicate.var]]),]
    
    # caculate stability on all reps
    X <- split(data2, data2[replicate.var])
    out <- lapply(X, stability_onerep, abundance.var)
    reps <- unique(data2[replicate.var])
    output <- cbind(reps, do.call("rbind", out))
    names(output) = c(replicate.var, "stability")
    output <- subset(output, !is.na(output$stability))

  }

  # result
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



#' A function to calculate species synchrony over time within one replicate
#'
#' @param df A dataframe containing x column
#' @param x The column to calculate stability on
#' @return Stability of x, calculated as the mean/sd
stability_onerep <- function(df, x){

  assertthat::assert_that(assertthat::has_name(df, x))
  assertthat::assert_that(is.numeric(df[[x]]))

  sync_var <- df[[x]]

  stability <- mean(sync_var, na.rm = TRUE)/stats::sd(sync_var, na.rm = TRUE)

  return(stability)
}

