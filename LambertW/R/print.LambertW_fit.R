#' @rdname LambertW_fit-methods
#' @description
#' \code{print.LambertW_fit} prints only very basic information about
#' \eqn{\widehat{\theta}} (to prevent an overload of data/information in the
#' console when executing an estimator).
#' @export
print.LambertW_fit <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("Estimation method: ", x$method, "\n")
  cat("Input distribution: ", x$distname, "\n")
  cat(ifelse(x$use.mean.variance,
             "mean-variance",
             "location-scale"), 
      "Lambert W x F type ('h' same tails; 'hh' different tails; 's' skewed): ", 
      x$type, "\n")
  cat("\n Parameter estimates:\n")
  if (x$method == "IGMM") {
    x$params.hat <- x$tau
  }
  print(x$params.hat)
  if (x$method == "IGMM") {
    param.txt <- "gamma"
    if (x$type == "h") {
      param.txt <- "delta"
    } else if (x$type == "hh") {
      param.txt <- "pair (delta_l, delta_r)"
    }
    if (is.na(x$sub.iterations)) {
      cat(paste("\n Obtained after", x$iterations, 
                "iterations for mu_x and sigma_x, and \n on average", 
                round(x$sub.iterations/x$iterations, 2), 
                "iterations to find the optimal", param.txt, "in each run."))
    } else {
      cat(paste("\n Obtained after", x$iterations, "iterations."))
    }
  }
  cat("\n")
} 