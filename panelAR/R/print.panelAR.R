print.panelAR <- function(x, digits=max(3,getOption("digits")-3),...){
	if(x$call$autoCorr=="none"){
		autoCorr.Method <- "no autocorrelation"
		} else{
			autoCorr.Method <- "AR(1) Prais-Winsten correction"
		}
	panelCorr.Method  <- switch(x$call$panelCorrMethod,none="homoskedastic variance",phet="panel heteroskedasticity-robust standard errors",pwls="panel weighted least squares",pcse="panel-corrected standard errors",parks="Parks-Kmenta FGLS")

	cat(paste("\nPanel Regression with ",autoCorr.Method, " and ",panelCorr.Method,"\n", sep = "")) 
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
        
    if(any(x$aliased)){
    	coef <- rep(NA,length(x$aliased))
    	names(coef) <- names(x$aliased)
    	coef[!x$aliased] <- coef(x)
    } else{
    	coef <- coef(x)
    }
    cat("Coefficients:\n")
    print.default(format(coef, digits = digits), print.gap = 2, 
    quote = FALSE)
    cat("\n")
    invisible(x)
}
