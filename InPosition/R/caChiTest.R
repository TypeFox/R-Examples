caChiTest <- function(DATA,res,critical.value=2){
	
	##cols
	j.vals.to.test <- res$ExPosition.Data$cj * matrix(res$ExPosition.Data$eigs,nrow(res$ExPosition.Data$cj),ncol(res$ExPosition.Data$cj),byrow=TRUE) * sum(DATA)
	j.df <- nrow(DATA) - 1
	#p.vals
	j.p.vals <- 1-pchisq(j.vals.to.test,j.df)
	j.signed.vals <- sign(res$ExPosition.Data$fj) * sqrt(j.vals.to.test)
	rownames(j.signed.vals) <- colnames(DATA)
	j.significant.vals <- j.p.vals < (2*(1-pnorm(critical.value)))
	rownames(j.significant.vals) <- rownames(j.signed.vals)
	
	##rows
	i.vals.to.test <- res$ExPosition.Data$ci * matrix(res$ExPosition.Data$eigs,nrow(res$ExPosition.Data$ci),ncol(res$ExPosition.Data$ci),byrow=TRUE) * sum(DATA)
	i.df <- ncol(DATA) - 1
	#p.vals
	i.p.vals <- 1-pchisq(i.vals.to.test,i.df)
	i.signed.vals <- sign(res$ExPosition.Data$fi) * sqrt(i.vals.to.test)
	rownames(i.signed.vals) <- colnames(DATA)
	i.significant.vals <- i.p.vals < (2*(1-pnorm(critical.value)))
	rownames(i.significant.vals) <- rownames(i.signed.vals)	
	
	##omni
	omni.val <- sum(res$ExPosition.Data$eigs * sum(DATA))
	omni.df <- (nrow(DATA)-1) * (ncol(DATA)-1)
	omni.p <- 1-pchisq(omni.val,omni.df)
	
	return(list(j.sig.vals=j.significant.vals, j.signed.vals=j.signed.vals, j.p.vals=j.p.vals, i.sig.vals=i.significant.vals, i.signed.vals=i.signed.vals, i.p.vals=i.p.vals, omni.val=omni.val,omni.p=omni.p))
	
}