#' Non-recursive Outlier Removal Procedure with Moving Criterion
#'
#' @description Non-recursive outlier removal procedure with moving criterion
#'    according to Van Selst & Jolicoeur (1994).
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
#'  percent of trials removed, number of trials removed in the procedure,and
#'  total number of trials in \code{exp_cell} before outlier removal.
non_recursive_mc <- function(exp_cell) {

  # Load creiterion cutoffs according to Table 4 in vanSelst & Jolicoeur (1994)
  linear_interpolation <- linear_interpolation

  # Get sample size of data
  sample_size <- length(exp_cell)

  # The procedure takes place only if there are at least 4 trials
  if (sample_size >= 4) {

    # If sample is greater than 100, use SDs for 100.
    if (sample_size > 100) {sample_size <- 100}

    # Look up sample size in the table and find the required standard deviation
    # this looks in the "non_recursive" column, obviously
    cutoff <- linear_interpolation$non_recursive[sample_size]

    ## Now use these values to complete the trimming
    # Calculate SD of exp_cell according to denominator n
    sd_exp_cell <- sd(exp_cell) * sqrt((length(exp_cell) - 1) / (length(exp_cell)))
    # Find max SDs
    sd_max <- cutoff * sd_exp_cell
    # Find maximum cutoff value
    maxcutoff <- sd_max + mean(exp_cell)
    # Find minimum cutoff value
    mincutoff <- mean(exp_cell) - sd_max
    # Find trials to include based on min and max cutoffs
    included_trials <- exp_cell < maxcutoff & exp_cell > mincutoff

    ## Compute final data (mean, number of trials removes and percent removed)
    final_data <- exp_cell[included_trials]
    # Number of trails removed
    num_removed <- length(exp_cell[!included_trials])
    # Percent removed
    percent_removed <- 100 - ((length(final_data) / length(exp_cell)) * 100)
    # Mean
    final_data <- mean(final_data)
    # Total number of trials
    total_trials <- length(exp_cell)

  } else {
    #If there 4 RTs or less, then return NA in final_data
    final_data <- NA
    percent_removed <- 0
    num_removed <- 0
    total_trials <- length(exp_cell)
  }

  # Return
  return_data <- c(final_data, percent_removed, num_removed, total_trials)
  return(return_data)

}




