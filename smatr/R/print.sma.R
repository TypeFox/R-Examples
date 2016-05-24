print.sma <- function(x, ..., coefbygroup=FALSE){

	obj <- x

	#cat("\nCall:\n", deparse(obj$call), "\n\n", sep = "")
	
	sectionline <- function()cat(rep("-",60),"\n",sep="")
	
	Method <- switch(obj$method,
		OLS = "Ordinary Least-Squares Regression",
		MA = "Major Axis",
		SMA = "Standardized Major Axis")
	
	cat("Call:",paste(deparse(obj$call),collapse="\n"),"\n\n")
	cat("Fit using", Method, "\n\n")

	if(obj$log != ''){
	cat("These variables were log-transformed before fitting:",obj$log,"\n\n")
	cat("Confidence intervals (CI) are at ",100*(1-obj$alpha),"%\n\n",sep="")
	}
	
	sectionline()
	
	ngroups <- length(obj$coef)
	
	if(ngroups == 1){
		if(obj$intercept)cat("Coefficients:\n")
		if(!obj$intercept)cat("Coefficients (no intercept included):\n")
		tcc <- t(obj$coef[[1]])
		rownames(tcc)[1] <- "estimate"
		print(tcc)
		cat("\n")
		
		cat("H0 : variables uncorrelated\n")
		cat("R-squared :", obj$r2[[1]], "\n")
		cat("P-value :", format.pval(obj$pval[[1]]), "\n\n")
		
		if(obj$slopetestdone){
		sectionline()
		cat("H0 : slope not different from",obj$slopetest[[1]]$test.value,"\n")
		cat("Test statistic : r=",signif(obj$slopetest[[1]]$r,4),"with",obj$n[[1]]-2,"degrees of freedom under H0\n")
		cat("P-value :",format.pval(obj$slopetest[[1]]$p),"\n")
		}
	} else {
		
		# if(obj$intercept)cat("Coefficients for common slope:\n")
		# if(!obj$intercept)cat("Coefficients for common slope (no intercept included):\n")
		# tcc <- t(obj$commoncoef)
		# rownames(tcc)[1] <- "estimate"
		# print(tcc)
		# cat("\n")
		
		if(obj$method != "OLS"){
		cat("Results of comparing lines among groups.\n\n")
		
		#if(obj$gt == "slopecom"){
			cat("H0 : slopes are equal.\n")
			cat("Likelihood ratio statistic :", signif(obj$commoncoef$LR,4),"with",obj$commoncoef$df,"degrees of freedom\n")
			cat("P-value :",format.pval(obj$commoncoef$p),"\n")
		#}
		
		sectionline()
		if(obj$gt == "shiftcom"){
			cat("\n")
			cat("H0 : no shift along common axis.\n")
			cat("Wald statistic:",signif(obj$gtr$stat,4),"with",obj$gtr$df,"degrees of freedom\n")
			cat("P-value :",format.pval(obj$gtr$p),"\n")
			sectionline()
		}
		if(obj$gt == "elevcom"){
			cat("\n")
			cat("H0 : no difference in elevation.\n")
			cat("Wald statistic:",signif(obj$gtr$stat,4),"with",obj$gtr$df,"degrees of freedom\n")
			cat("P-value :",format.pval(obj$gtr$p),"\n")
			sectionline()
		}
		cat("\n")
		} else cat("Cannot perform common slope test with method == \"OLS\" , use lm() instead.\n\n")
		
		
		if(!all(is.na(obj$commonslopetestval))){
		
			x <- obj$commonslopetestval
			cat("H0 : common slope not different from",x$b,"\n")
			cat("Likelihood ratio statistic =",signif(x$LR,4),"with",x$df,"degrees of freedom under H0\n")
			cat("P-value :",format.pval(x$p),"\n\n")
		
		}
		
		# Multiple comparisons.
		
		if(obj$multcompdone != "none"){
			cat("Results of multiple comparisons among groups.\n\n")
			cat("Test for pair-wise difference in",obj$multcompdone,":\n")
			#browser()
			dfr <- obj$multcompresult
			dfr[,3] <- format.pval(obj$multcompresult[,3])
			print(dfr[,1:4])
			cat("\n")
			
			if(obj$multcompmethod == "default")
				cat("No adjustment of P values.\n")
			if(obj$multcompmethod == "Bonferroni")
				cat("Sidak adjustment of P values.\n")
			sectionline()
		}
		
		
		
		if(coefbygroup){
		cat("Coefficients by group in variable \"",obj$groupvarname,"\"\n\n",sep="")
		for(i in 1:ngroups){
			cat("Group:",names(obj$coef)[i],"\n")
			tcc <- t(obj$coef[[i]])
			rownames(tcc)[1] <- "estimate"
			print(tcc)
			cat("\n")
			cat("H0 : variables uncorrelated.\n")
			cat("R-squared :", obj$r2[[i]], "\n")
			cat("P-value :", format.pval(obj$pval[[i]]), "\n")
			cat("\n")
			
			if(obj$slopetestdone){
				cat("H0 : slope not different from",obj$slopetest[[1]]$test.value,"\n")
				cat("Test statistic: r=",signif(obj$slopetest[[i]]$r,4),"with",obj$n[[i]]-2,"degrees of freedom under H0\n")
				cat("P-value :",format.pval(obj$slopetest[[i]]$p),"\n\n")
			}
		}
		
		cat("\n")
		} else cat("Use the summary() function to print coefficients by group.\n")

		
		
			
	}
	
	if(obj$elevtestdone){
	sectionline()
	if(ngroups == 1){
		cat("H0 : elevation not different from",obj$elevtest[[1]]$test.value,"\n")
		cat("Test statistic: t=",signif(obj$elevtest[[1]]$t,4),"with",obj$n[[1]]-2,"degrees of freedom under H0\n")
		cat("P-value :",format.pval(obj$elevtest[[1]]$p),"\n")
	} else {
		cat("Warning : elev.test ignored in fit to multiple groups.\n")
	
	}
	} 
	
	
}







