runPvalues <-
function(data, addit = FALSE, covariable=NULL, family=binomial) {
	padd <- NULL  
	dataA <- ff(data, model="add")
	if(addit==FALSE) {
		dataD <- ff(data, model="dom")
		dataR <- ff(data, model="rec")
		dataC <- ff(data, model="codom")
		pdom <- prec <- pcodom <- NULL
		for(i in 2:ncol(data)){
			pdom <- c(pdom, pvalFmla(dataD, i, covariable, family))
      prec <- c(prec, pvalFmla(dataR, i, covariable, family))
      pcodom <- c(pcodom, pvalFmla(dataC, i, covariable, family))
    	padd <- c(padd, pvalFmla(dataA, i, covariable, family))
		}
		pvalors <- data.frame(pdom, pcodom, prec, padd)
    		pvalors$min <- apply(pvalors,1,min,na.rm=TRUE)
    		return(pvalors$min)
	}
	if (addit==TRUE) { 
		for(i in 2:ncol(data)){
			padd <- c(padd, pvalFmla(dataA, i, covariable, family)) 
		}
	} 	
  	return(padd)
}
