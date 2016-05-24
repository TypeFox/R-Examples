print.summary.panelAR <- function(x,digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"),...){
	if(x$call$autoCorr=="none"){
		autoCorr.Method <- "no autocorrelation"
		} else{
			autoCorr.Method <- "AR(1) Prais-Winsten correction"
		}
	panelCorr.Method  <- switch(x$call$panelCorrMethod,none="homoskedastic variance",phet="panel heteroskedasticity-robust standard errors",pwls="panel weighted least squares",pcse="panel-corrected standard errors",parks="Parks-Kmenta FGLS")
	
	table.structure <- cbind(c("Total obs.:","Number of panels:","Number of times:"),
	c(x$panelStructure$N,x$panelStructure$N.panel,x$panelStructure$N.time),
	c("Avg obs. per panel","Max obs. per panel","Min obs. per panel"),
	c(round(x$panelStructure$N.avg,digits),x$panelStructure$N.max,x$panelStructure$N.min))
	dimnames(table.structure) <- list(rep("",nrow(table.structure)),rep("",ncol(table.structure)))

	cat(paste("\nPanel Regression with ",autoCorr.Method, " and ",panelCorr.Method,"\n", sep = ""))
	cat(paste("\n",ifelse(x$panelStructure$balanced,"Balanced","Unbalanced")," Panel Design:",sep=""))
	print.default(table.structure, quote=F,print.gap=1)
	
	if (any(x$aliased)) {
            cnames <- names(x$aliased)
            coefs <- matrix(NA, length(x$aliased), 4, dimnames = list(cnames, 
                colnames(x$coefficients)))
            coefs[!x$aliased, ] <- x$coefficients
        } else{
        	coefs <- x$coefficients
    }
	cat("\nCoefficients:\n")
	printCoefmat(coefs,digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    if(!is.null(x$r2)){cat(paste("\nR-squared: ",round(x$r2,4),sep=""))}
    cat(paste("\nWald statistic: ",round(x$wald["value"],4),", Pr(>Chisq(",x$wald["df"],")): ",round(x$wald["Pr(>Chisq)"],4),"\n",sep=""))
	}