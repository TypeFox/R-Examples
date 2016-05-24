`print.summary.mlcm` <-
function(x, 
	digits = max(3, getOption("digits") - 4), ...){
  cat("\n Maximum Likelihood Conjoint Measurement\n")
  cat("\nLink:\t")
  cat(x$link)
  cat("\t\tModel:\t")
  cat(x$model)	
  cat("\n\nPerceptual Scale:\n")
    print.default(format(x$pscale, digits = digits), 
    	quote = FALSE, ...)
  if(x$method == "formula"){
    cat("\nformula:\t", as.character(x$formula))
    cat("\n\np:\t", x$par)
    }
  cat("\n")  
  cat("\nStandard Errors:\n") 
  print.default(format(x$se, digits = digits),
   		quote = FALSE, ...)
  cat("\nsigma:\t")
    cat(format(x$sigma, digits = digits))
    cat("\t\tlogLik:\t")
    cat(format(x$logLik, digits = digits))
    cat("\t\tAIC:\t")
    cat(format(x$aic, digits = digits))
    cat("\n\n")	
    invisible(x)
}

