print.moult <- function(x, digits = max(3, getOption("digits") - 3), ...)      # see print.lm
 {  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
        0.85)), "", sep = "\n")
    
    if (!x$converged) {
      cat("model did not converge\n")
    }
    else {
    
      cat("\n Coefficients for moult duration:\n")
      print(x$coefficients$duration, digits = digits)

      cat("\n Coefficients for mean start date of moult:\n")
      print(x$coefficients$mean, digits = digits)

      cat("\n Coefficients for standard deviation in start date of moult:\n")
      print(x$coefficients$sd, digits = digits)

    }

    invisible(x)
 }