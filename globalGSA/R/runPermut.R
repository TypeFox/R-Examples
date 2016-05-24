runPermut <-
function(data, addit=FALSE, covariable=NULL, family=binomial) {
	TraitR <- sample(data[,1])
   	dataR <- data.frame(TraitR, data[,2:ncol(data)])
    	return(runPvalues(dataR, addit, covariable, family))
  }
