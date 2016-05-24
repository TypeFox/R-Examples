#' Profile atomic vectors or data.frames.
#'
#' \code{profile} takes vectors and data.frames and returns a data.frame containing important descriptive statistics.
#'
#' @param .x a vector or data.frame to be profiled.
#' @return  a data.frame containing important descriptive statistics.
#' @examples
#' # Example
#' profile(mtcars)
#' profile(iris)
#' profile(state.name)
#' @export
profile <- function(.x) UseMethod("profile")

#' Profile numeric atomic vectors.
#'
#' \code{profile_numeric} takes numeric vectors and returns a data.frame containing important descriptive statistics.
#'
#' @param .x a numeric atomic vector or data.frame to be profiled.
#' @return  a data.frame containing important descriptive statistics.
#' @examples
#' # Example
#' profile_numeric(1:100)
#' @export
profile_numeric <- function(.x) {
  if(!is.numeric(.x)) stop(".x must be a numeric vector!")
  if(!is.atomic(.x)) stop(".x must be an atomic vector!")
  profile_functions <- list(
    .count_elements = function(.x, ...) length(.x),
    .count_uniques = function(.x, ...) length(unique(.x)),
    .percent_uniques = function(.x, ...) length(unique(.x)) / length(.x),
    .count_NULLs = function(.x, ...) sum(is.null(.x)),
    .percent_NULLs = function(.x, ...) sum(is.null(.x)) / length(.x),
    .count_NAs = function(.x, ...) sum(is.na(.x)),
    .percent_NAs = function(.x, ...) sum(is.na(.x)) / length(.x),
    .count_zeroes = function(.x, ...) sum(.x == 0),
    .percent_zeros = function(.x, ...) sum(.x == 0) / length(.x),
    .mean_value = function(.x, ...) mean(.x, ...),
    .sd_value = function(.x, ...) stats::sd(.x, ...),
    .q0_value = function(.x, ...) min(.x, ...),
    .q25_value = function(.x, ...) stats::quantile(.x, probs = 0.25, ...),
    .q50_value = function(.x, ...) stats::median(.x, ...),
    .q75_value = function(.x, ...) stats::quantile(.x, probs = 0.75, ...),
    .q100_value = function(.x, ...) max(.x, ...),
    .top_5_values = function(.x, ...) frequent_values(.x, n = 5, type = "most"),
    .bottom_5_values = function(.x, ...) frequent_values(.x, n = 5, type = "least")
  )
  return(as.data.frame(lapply(profile_functions,
                              function(.f) .f(.x, na.rm = TRUE)),
                       row.names = "1",
                       stringsAsFactors = FALSE))
}

#' Profile non-numeric atomic vectors.
#'
#' \code{profile_nonnumeric} takes non-numeric atomic vectors and returns a data.frame containing important descriptive statistics.
#'
#' @param .x a non-numeric atomic vector to be profiled.
#' @return  a data.frame containing important descriptive statistics.
#' @examples
#' # Example
#' profile(letters)
#' @export
profile_nonnumeric <- function(.x) {
  if(is.numeric(.x)) stop(".x must be a non-numeric vector!")
  if(!is.atomic(.x)) stop(".x must be an atomic vector!")
  profile_functions <- list(
    .count_elements = function(.x, ...) length(.x),
    .count_uniques = function(.x, ...) length(unique(.x)),
    .percent_uniques = function(.x, ...) length(unique(.x)) / length(.x),
    .count_NULLs = function(.x, ...) sum(is.null(.x)),
    .percent_NULLs = function(.x, ...) sum(is.null(.x)) / length(.x),
    .count_NAs = function(.x, ...) sum(is.na(.x)),
    .percent_NAs = function(.x, ...) sum(is.na(.x)) / length(.x),
    .count_zeroes = function(.x, ...) sum(.x == 0),
    .percent_zeros = function(.x, ...) sum(.x == 0) / length(.x),
    .mean_value = function(.x, ...) NA,
    .sd_value = function(.x, ...) NA,
    .q0_value = function(.x, ...) min(.x, ...),
    .q25_value = function(.x, ...) NA,
    .q50_value = function(.x, ...) NA,
    .q75_value = function(.x, ...) NA,
    .q100_value = function(.x, ...) max(.x, ...),
    .top_5_values = function(.x, ...) frequent_values(.x, n = 5, type = "most"),
    .bottom_5_values = function(.x, ...) frequent_values(.x, n = 5, type = "least")
  )
  return(as.data.frame(lapply(profile_functions,
                              function(.f) .f(.x, na.rm = TRUE)),
                       row.names = "1",
                       stringsAsFactors = FALSE))
}


#' @describeIn profile Method for numeric.
#' @export
profile.numeric <- function(.x) return(profile_numeric(.x))

#' @describeIn profile Method for character.
#' @export
profile.character <- function(.x) return(profile_nonnumeric(.x))

#' @describeIn profile Method for data.frames.
#' @export
profile.data.frame <- function(.x) {
  call <- match.call()
  .df <- do.call(rbind, lapply(.x, profile))
  .metadata <- data.frame(.column_name=rownames(.df),
                          .column_class=sapply(.x, class),
                          .column_type=sapply(.x, typeof),
                          stringsAsFactors = FALSE,
                          row.names = as.character(1:nrow(.df)))
  .profiled_dataset <- cbind(.metadata, .df)
  rownames(.profiled_dataset) <- NULL
  return(.profiled_dataset)
}

#' @describeIn profile Method for default.
#' @export
profile.default <- function(.x) return(profile_nonnumeric(as.character(.x)))
