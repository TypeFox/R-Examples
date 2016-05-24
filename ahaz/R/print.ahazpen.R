"print.ahazpen" <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    ## Purpose: Print 'ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahazpen' object
    ##   digits: number of digits in output
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

    cat("* No. predictors:            ", format(x$nvars, width = 7,digits = digits),  "\n")
    cat("* No. observations:          ", format(x$nobs, width = 7,digits = digits), "\n")
    cat("* Max no. predictors in path:", format(max(x$df), width = 7,digits = digits),"\n")
    cat("* Penalty parameter lambda:\n")
    cat("    -No. grid points:", format(length(x$lambda), width = 7,digits = digits), "\n")
    cat("    -Min value:      ",  format(min(x$lambda), width = 7,digits = digits), "\n")
    cat("    -Max value:      ", format(max(x$lambda), width = 7,digits = digits))
    cat("\n\n")
  
    invisible(x)
  }
