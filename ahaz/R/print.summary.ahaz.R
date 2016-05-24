"print.summary.ahaz" <- function(x, digits=max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"), ...)
  {
    ## Purpose: Print summary of 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'summary.ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    cat("\nCall:\n")
    dput(x$call)
    cat("\n")    
    cat("  n =", x$nvars,"\n\n")
    if(any(is.na(x$coefficients))){
       cat("Coefficients: (some undefined due to singularities) \n\n")
     }
       else {    
         if(x$univar)
           cat("Univariate coefficients:\n")
         else
           cat("Coefficients:\n")
     }
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars)
    cat("\n")
    if(!is.null(x$waldtest))
      cat("Wald test = ", format(round(x$waldtest["test"], digits=digits)), "  on ",
          x$waldtest["df"], " df,", "   p=", format(x$waldtest["pvalue"],digits=digits),
          "\n\n", sep = "")
  }
