#' Delete cases with set amount of missing values
#'
#' @param df A \code{data.frame},
#' @param n Number of \code{NAs} allowed, defaults to \code{ncol(df) - 1}.
#'
#' @return A filtered version of the input \code{data.frame}.
#' @export
#' @note Adapted from \url{http://stackoverflow.com/a/30461945/409362}.
#' @examples
#' \dontrun{
#' df <- delete_na(df)
#' }
delete_na <- function(df, n = ncol(df) - 1) {
  log      <- apply(df, 2, is.na)
  logindex <- apply(log, 1, function(x) sum(x) <= n)

  return(df[logindex, ])
}

#' Convert all labels to factor variables
#'
#' @param df A \code{data.frame}
#'
#' @return An identical \code{data.frame} with labelled data converted to factors
#' @export
#' @importFrom sjmisc is_labelled
#' @importFrom haven as_factor
#' @examples
#' \dontrun{
#' data %<>% labels_to_factor
#' }
labels_to_factor <- function(df) {
  for (column in names(df)) {
    if (sjmisc::is_labelled(df[[column]])) {
      df[[column]] <- haven::as_factor(df[[column]])
    }
  }
  return(df)
}

#' Re-label a vector after subsetting
#'
#' @param x A vector with now unused labels
#' @return Identical vector with appropriate labels
#' @export
#' @importFrom sjmisc set_labels
#' @importFrom sjmisc get_labels
#' @examples
#' \dontrun{
#' x <- drop_labels(x)
#' }
drop_labels <- function(x) {
  sjmisc::set_labels(x, labels = sjmisc::get_labels(x)[unique(x)])
}

