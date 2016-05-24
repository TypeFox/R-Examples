"print.ahaz" <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    ## Purpose: Print 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    cat("\nCall: ", deparse(x$call), "\n\n")
  
    cat("Coefficients:\n")
    c<-coef(x)
    print.default(format(c, digits = digits), 
                  print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}
