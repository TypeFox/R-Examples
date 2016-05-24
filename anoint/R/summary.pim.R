pim.summary <- function (object,...) 
{
	
	#TAKEN FROM SUMMARY.GLM
	
	if(!object@exact&all(dim(object@boot.pim)==c(1,1))){
		print(object@boot.pim)
		warning("Coef. table requires bootstrap resamples")
	}
	else if(object@exact&all(dim(object@boot.pim)==c(1,1))){
		# FIT WITH APPROXIMATE VCOV, NO VARIANCE FOR RESPONSIVENESS	

		cov.mat <- coef(object)
		theta <- cov.mat[names(cov.mat)=="theta"]
		cov.mat <- cov.mat[names(cov.mat)!="theta"]
		cov.mat <- c(cov.mat[names(cov.mat)!="Treatment"],cov.mat[names(cov.mat)=="Treatment"])
		cov.p <- length(cov.mat)
		df <- nrow(object@formula@data)-cov.p
		se <- sqrt(diag(vcov(object)))
		tvalue <- cov.mat/se
		pvalue <- 2*pt(-abs(tvalue),df)
        coef.table <- cbind(cov.mat,se,tvalue,pvalue)
		
        dimnames(coef.table) <- list(names(cov.mat),c("Estimate","SE","t value","p-value"))

        print(coef.table)    
        cat("\ntheta: ",theta,"\n")    
        warning("Std. error does not account for responsiveness parameter so could be underestimated.")
           
 		invisible(coef.table)		
	}
	else{		
		cov.mat <- coef(object)
		cov.p <- length(cov.mat)
		df <- nrow(object@formula@data)-cov.p
		se <- sqrt(diag(vcov(object)))
		tvalue <- cov.mat/se
		pvalue <- 2*pt(-abs(tvalue),df)
        coef.table <- cbind(cov.mat,se,tvalue,pvalue)

        dimnames(coef.table) <- list(names(cov.mat),c("Estimate","SE","t value","p-value"))
        
        print(coef.table)
             
 invisible(coef.table)
 }
}

setMethod("summary","pim",pim.summary)
