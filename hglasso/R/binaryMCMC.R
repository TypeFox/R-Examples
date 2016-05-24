

######################################################
# Simulate Binary Network Data
######################################################
binaryMCMC<-function(n,Theta,burnin,skip,trace=FALSE){
	
	
	  
p <- nrow(Theta)	
X <- matrix(numeric(n*p),nrow=n)
temp <- sample(c(1,0),p,replace=TRUE)

######################################################
# samples from burn in period
######################################################
for(t in 1:burnin){
	for(j in 1:p){
		temp[j]<-rbinom(1,1,exp(Theta[j,j]+sum((temp*Theta[j,])[-j]))/(1+exp(Theta[j,j]+sum((temp*Theta[j,])[-j]))))
	}
	if(t%%10000==0){
	}		 	
}

######################################################
# Samples obtained by skipping 500 samples 
######################################################
k<-1
samples<-1
while(samples<=n){
	for(j in 1:p){
		temp[j]<-rbinom(1,1,exp(Theta[j,j]+sum((temp*Theta[j,])[-j]))/(1+exp(Theta[j,j]+sum((temp*Theta[j,])[-j]))))		
	}
	if(k%%skip==1){
		X[samples,] <- temp
		samples<-samples+1
		if(trace==TRUE){print(samples)}
	}
	k<-k+1
}

return(X)
}



