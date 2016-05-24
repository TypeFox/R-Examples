################################################################################
##  print.manyglm prints the  manyglm object in a nice way        	          ##
## so far this is identical with the print.glm method 	                      ##
## maybe we want to add some more specific things later on    			          ##
################################################################################

print.manyglm <- function( x, digits = max(3, getOption("digits") - 3), dig.tst=max(2, min(5, digits-1)), ... )  {

  cat("\nCall: ", deparse(x$call), "\n")
	
#	family <- x$family$family
#		if(	substr(family[1],1, 12)=="quasipoisson") {
#		print(quasipoisson(link=x$family$link))
#	} else if (substr(	family[1], 1, 17)=="Negative Binomial"){
#		print(negative.binomial(link=x$family$link))
#	} else 	
  print(x$family)  

  if (x$show.coef==TRUE) {
     if (length(coef(x))) {
        cat("Coefficients")
#       if (is.character(co <- x$contrasts)) 
#          cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(round(x$coefficients, digits=dig.tst), digits = digits), print.gap = 2, quote = FALSE)
      }
      else cat("No coefficients\n\n")   
  }      
  if (x$show.fitted==TRUE) {
        cat("Fitted Values\n")
        print.default(format(round(x$fitted.values, digits=dig.tst), digits = digits), print.gap = 2, quote = FALSE)
  }
  if (x$show.residuals=="pearson") {
        cat("Standardized Pearson Residuals\n")
        print.default(format(round(x$Pearson.residuals, digits=dig.tst), digits = digits), print.gap = 2, quote = FALSE)
  }
  else if (x$show.residuals=="PIT") {
        cat("Uniform PIT Residuals\n")
        print.default(format(round(x$PIT.residuals, digits=dig.tst), digits = digits), print.gap = 2, quote = FALSE)
  }

  if (x$family=="quasipoisson" | x$family=="negative.binomial"){
      p <- length(x$phi)
      cat("\nNuisance Parameter(s) phi estimated by the", x$theta.method, "method.\n")
      print.default(format(round(x$phi,digits=dig.tst)),print.gap=2,quote=FALSE)
  }		

  cat("\nDegrees of Freedom:", NROW(x$y)-1, "Total (i.e. Null);", x$df.residual, "Residual\n")

  if (nchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")

  table <- rbind(round(x$two.loglike, digits=dig.tst), round(x$deviance, digits=dig.tst), round(x$aic, digits=dig.tst))
  rownames(table) <- c("2*log-likelihood:", "Residual Deviance:", "AIC:")
  colnames(table) <- colnames(x$coefficients)
  cat("\n")
  print.default(format(table, digits = digits), print.gap = 2, quote = FALSE)

#  if(!is.null(x$ave.goodness.fit)){
#    cat("\n")
#    cat("Deviance statistics averaged across all variables:")
#    ave.goodness     <- matrix(nrow=2, ncol=2)
#    ave.goodness[1,] <- c(format( x$ave.goodness.fit[1],
#      digits = max(5, digits + 1)) ,
#      paste(" on", format(x$df.null), "degrees of freedom")  )
#    ave.goodness[2,] <- c(format( x$ave.goodness.fit[2],
#      digits = max(5, digits + 1)),
#      paste(" on", format(x$df.residual), "degrees of freedom")  )
#    dimnames(ave.goodness) <-
#      list( c("Null deviance:", "Residual deviance:"), c("","") )
#    print.default(ave.goodness, quote =FALSE, right = TRUE, na.print = "NA",...)
#    cat( c("AIC:" , format( x$ave.goodness.fit[3],
#      digits = max(4, digits + 1)))) # ,quote =FALSE, right = FALSE, na.print = "NA",...)
#    cat( c("       -2 x log-likelihood:", format( x$ave.goodness.fit[4],
#      digits = max(4, digits + 1))) )
#    cat("\n")
#  }
  invisible(x)
}


