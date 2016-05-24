# Internal function to convert approximate covariance matrix of 'variance' parameter estimates to natural scale

refactor_apVar<-function(apVar, cov_class){
	
	npars<-dim(apVar)[1]
	
	if(cov_class=="covBM"){
		for(i in 1:(npars-2)){
						apVar[npars-1, i]<- 2*apVar[npars,i] + apVar[npars-1,i]
						}
						apVar[npars-1, npars-1]<- 4*apVar[npars,npars] + 4*apVar[npars,npars-1] + apVar[npars-1,npars-1]
						apVar[npars-1, npars]<- 2*apVar[npars,npars] + apVar[npars-1,npars]
		for(i in 1:npars){
						if(i!=npars-1) {apVar[i, npars-1]<-apVar[npars-1, i]}
						}
		attr(apVar,"Pars")["corStruct"]<- as.numeric(attr(apVar,"Pars")["corStruct"] + 2*attr(apVar,"Pars")["lSigma"])
		}
		
	if(cov_class=="covFracBM" || cov_class=="covIOU"){
		for(i in 1:(npars-3)){
						apVar[npars-2, i]<- 2*apVar[npars,i] + apVar[npars-2,i]
						}
						apVar[npars-2, npars-2]<- 4*apVar[npars,npars] + 4*apVar[npars,npars-2] + apVar[npars-2,npars-2]
						apVar[npars-2, npars-1]<- 2*apVar[npars,npars-1] + apVar[npars-2,npars-1]
						apVar[npars-2, npars]<- 2*apVar[npars,npars] + apVar[npars-2,npars]
		for(i in 1:npars){
						if(i!=npars-2) {apVar[i, npars-2]<-apVar[npars-2, i]}
						}
		attr(apVar,"Pars")["corStruct1"]<- as.numeric(attr(apVar,"Pars")["corStruct1"] + 2*attr(apVar,"Pars")["lSigma"])
		}		
		
	return(apVar)
	}
