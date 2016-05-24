#' Provides a summary of missingness in a dataset.
#'
#' Generates a summary of the percentage of missing data in a dataset.  Provides
#' insight on the appropriateness of imputation methods.  For example, if 30\%
#' of data is missing, then perhaps this is too much to impute.
#'
#' @param aDataFrame A data.frame containing columns that will be assessed for
#'    missingness.  
#'
#' @return A data frame that summarizes percent missingness for each column of 
#'    a dataset.
#'
#' @examples
#' data(example_references_metagear)
#' impute_missingness(example_references_metagear)
#'
#' @export impute_missingness

impute_missingness <- function(aDataFrame) {

  columnMissingness <- colSums(is.na(aDataFrame)) / nrow(aDataFrame) * 100
  cat("\nSummary of missingness:\n\n")
  summary.columns <- data.frame(
    "COLUMN" = names(columnMissingness), 
    "PERCENT_MISSINGNESS" = columnMissingness,
    "IMPUTATIONS" = colSums(is.na(aDataFrame)))
  print(summary.columns, 
        row.names = FALSE, 
        digits = 2, 
        na.print = "", 
        quote = FALSE )
  cat("\nTotal missingness: ", 
      format(sum(is.na(aDataFrame)) / prod(dim(aDataFrame)) * 100, digits = 2),
      "% (" ,sum(is.na(aDataFrame)), 
      " imputations needed)\n\n", 
      sep = "")
  return (summary.columns)
}
