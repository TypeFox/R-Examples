
# version du 13.06.13


print.flexrsurv=function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts)) 
      cat("  [contrasts: ", apply(cbind(names(co), co), 
                                  1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  
  if(!is.null(x$optim.control$reltol)){
    digitsll <- round(abs(log10(x$optim.control$reltol)))
    if ( digitsll < 0) {
      digitsll <- digits
    }
  } else {
      digitsll <- digits
  }

  cat("\nLog-likelihood:\t ", format(x$loglik, digits = digitsll), "\n")
  
  invisible(x)
  
  
}
