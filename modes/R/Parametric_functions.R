##	Parametric Functions



#########################
#	Ashman's D		#
#########################


Ashmans_D<- function(mu1, mu2, sd1, sd2,...){
	D= sqrt(2)*(abs(mu1-mu2))/sqrt((sd1^2)+(sd1^2))
	return(D)
}


#####################################
#	Bimodality Separation		#
#####################################


bimodality_separation<-function(mu1,mu2,sd1,sd2,...){
	S=0.5*(mu1-mu2)/(sd1+sd2)
	return(S)
}


