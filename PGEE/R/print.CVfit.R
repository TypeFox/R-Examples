print.CVfit <-
function(x, digits = NULL, quote = FALSE, prefix = "", ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
    
    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")
  
    cat(sprintf("\n%d-fold CV results:\n", x$fold))
    print(cbind("lambda"=x$lam.vect, "Cv"=x$cv.vect))
    cat("\nOptimal tuning parameter:\n")
    optimalTuning <- c("Best lambda"=x$lam.opt)
    print(optimalTuning) 
   
    # return object invisibly
    invisible(x)
}
