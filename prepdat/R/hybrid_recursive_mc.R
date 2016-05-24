#' Hybrid-recursive Outlier Removal Procedure with Moving Criterion
#'
#' @description Hybrid-recursive outlier removal procedure with moving
#'   criterion according to Van Selst & Jolicoeur (1994).
#' @param exp_cell Numeric vector on which the outlier removal method takes
#'   place. If experimental cell has 4 trials or less it will result in
#'   \code{NA}.
#' @references Grange, J.A. (2015). trimr: An implementation of common response
#'   time trimming methods. R Package Version 1.0.0.
#'   \url{https://cran.r-project.org/package=trimr}
#'
#' Van Selst, M., & Jolicoeur, P. (1994). A solution to the effect of sample
#' size on outlier elimination. \emph{The quarterly journal of experimental
#' psychology, 47}(3), 631-650.
#' @return A vector with the mean of \code{exp_cell} after removing outliers,
#'  percent of trials removed, and total number of trials in \code{exp_cell} before
#'  outlier removal.
hybrid_recursive_mc <- function(exp_cell) {

  # Call non_recursive_mc function
  nrmc <- non_recursive_mc(exp_cell)

  # Call modified_recursive_mc function
  mrmc <- modified_recursive_mc(exp_cell)

  # Hybrid procedure: get the mean of non-recursive and the modified-recursive
  # procedures
  hrmc <- mean(c(nrmc[1], mrmc[1]))

  # Get the mean percent of trials removed
  percent_removed <- mean(c(nrmc[2], mrmc[1]))

  # Get total number of trials before trimming
  total_trials <- length(exp_cell)

  # Return
  return_data <- c(hrmc, percent_removed, total_trials)
  return(return_data)
}
