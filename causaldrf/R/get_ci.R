#' This calculates an upper and lower bound from bootstrap matrix
#'
#' This function takes a matrix containing the bootstrapped coefficients
#' from a parametric ADRF estimator and returns upper and lower 95 percent
#' confidence lines.
#'
#' @param grid_val is the vector of grid values on \code{treat} axis
#' @param coef_mat contains the bootstrapped parameter estimates.
#' @param degree is 1 for linear and 2 for quadratic outcome model
#'
#' @return
#' \code{get_ci} returns upper and lower 95 percent confidence lines.
#'
#' @usage
#'
#' get_ci(grid_val,
#'        coef_mat,
#'        degree)
#'
#' @export
#'



get_ci <- function(grid_val,
                   coef_mat,
                   degree){


  # gets the ci from bootstrapped parameter estimates

  if (degree == 1){


  grid_ci <- matrix( numeric(length(grid_val) * nrow(coef_mat) ), nrow = nrow(coef_mat))


  for ( i in 1:nrow(coef_mat)){
    grid_ci[i, ] <- coef_mat[i, 1] + coef_mat[i, 2] * grid_val
  }
  param_sorted <- apply(coef_mat, 2, sort, decreasing=F)

  # sort each column of grid_ci
  ci_sorted <- apply(grid_ci, 2, sort, decreasing=F)

  # get the 97.5% and 2.5% bands

  upper_ci <- ci_sorted[ceiling(0.975 * nrow(coef_mat)), ]
  lower_ci <- ci_sorted[floor(0.025 * nrow(coef_mat)), ]

  ci_est <- rbind(upper_ci, lower_ci)


  } else if (degree == 2){


    grid_ci <- matrix( numeric(length(grid_val) * nrow(coef_mat) ), nrow = nrow(coef_mat))


    for ( i in 1:nrow(coef_mat)){
      grid_ci[i, ] <- coef_mat[i, 1] + coef_mat[i, 2] * grid_val + coef_mat[i, 3] * grid_val^2
    }
    param_sorted <- apply(coef_mat, 2, sort, decreasing=F)

    # sort each column of grid_ci
    ci_sorted <- apply(grid_ci, 2, sort, decreasing=F)

    # get the 97.5% and 2.5% bands

    upper_ci <- ci_sorted[ceiling(0.975 * nrow(coef_mat)), ]
    lower_ci <- ci_sorted[floor(0.025 * nrow(coef_mat)), ]

    ci_est <- rbind(upper_ci, lower_ci)


  }
  return(list(ci_est))

}

