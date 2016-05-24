`garsim` <-
function(n,phi,X=matrix(0,nrow=n),beta=as.matrix(0), sd=1,
 family= "gaussian", transform.Xbeta= "identity", link= "identity",
 minimum=0, zero.correction= "zq1", c=1, theta=0){
	# Initialize the parameters
	y<-matrix(0,nrow=n)
	m<-matrix(0,nrow=n)
	lambda<-matrix(0,nrow=n)
	y.transformed<-matrix(0,nrow=n)
	phi<-as.matrix(phi)
	p<-nrow(phi)
	if(transform.Xbeta== "identity"){
		Xbeta<- as.matrix(as.matrix(X)%*%as.vector(beta))
	}
	if(transform.Xbeta== "exponential"){
		Xbeta<- exp(as.matrix(as.matrix(X)%*%as.vector(beta)))
	}
	
	if(link== "identity"){
		if(family== "gaussian"){
			# Loop
			for(t in 1:n){ 
				if(t<(p+1)){
					m[t]<- Xbeta[t]
					y[t]<-rnorm(1, m[t], sd=sd)
				}
				else{
					yphi<-t(y[(t-p):(t-1)])%*%rev(phi)
					Xbetaphi<-t(-Xbeta[(t-p):(t-1)])%*%rev(phi)
					m[t]<- yphi + Xbetaphi + Xbeta[t]
					y[t]<-rnorm(1, m[t], sd=sd)
				}
			}
		}
		if(family== "poisson"){
			for(t in 1:n){
				if(t<(p+1)){
					m[t]<- Xbeta[t]
					y[t]<-rpois(1,m[t])
				}
				else{
					yphi<-t(y[(t-p):(t-1)])%*%rev(phi)
					Xbetaphi<-t(-Xbeta[(t-p):(t-1)])%*%rev(phi)
					m[t]<- max(yphi + Xbetaphi + Xbeta[t] ,minimum) #The max() is #to ensure m[t]>0, m[t] is not allowed to be <0, and for m[t]=0 all following simulated data will be #0.
					y[t]<-rpois(1,m[t])
				}
			}
		}
		if(family== "negative.binomial"){
			#require(MASS)
			for(t in 1:n){
				if(t<(p+1)){
					m[t]<- Xbeta[t]
					y[t]<-rnegbin(1,m[t], theta)
				}
				else{
					yphi<-t(y[(t-p):(t-1)])%*%rev(phi)
					Xbetaphi<-t(-Xbeta[(t-p):(t-1)])%*%rev(phi)
					m[t]<- max(yphi + Xbetaphi + Xbeta[t] ,minimum) 
					#The max() is #to ensure m[t]>0, m[t] is #not allowed to be <0, and for m[t]=0 all following simulated #data will be #0.
					y[t]<-rnegbin(1,m[t],theta)
				}
			}
		}
	}
	
	if(link== "log"){
		if(family== "gaussian"){
			stop(cat("Gaussian model with log link is not implemented, "))
		} 
		if(family== "poisson"){
			if(zero.correction== "zq1"){
				for(t in 1:n){
					if(t<(p+1)){
						m[t]<- Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rpois(1,lambda[t])
						y.transformed[t]<-max(c, y[t])
					}
					else{
						yphi<-t(log(y.transformed[(t-p):(t-1)]))%*%rev(phi)
						Xbetaphi<-t(-Xbeta[(t-p):(t-1)])%*%rev(phi)
						m[t]<- yphi + Xbetaphi + Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rpois(1,lambda[t])
						y.transformed[t]<-max(c, y[t])
					}
				}
			}
			if(zero.correction== "zq2"){
				for(t in 1:n){
					if(t<(p+1)){
						m[t]<- Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rpois(1,lambda[t])
						y.transformed[t]<- y[t]+c
					}
					else{
						yphi<-t(log(y.transformed[(t-p):(t-1)]))%*%rev(phi)
						Xbetaphi<-t(-log(exp(Xbeta[(t-p):(t-1)])+c))%*%rev(phi)
						m[t]<- yphi + Xbetaphi + Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rpois(1,lambda[t])
						y.transformed[t]<-y[t]+ c
					}
				}
			}
		}
		
		if(family== "negative.binomial"){
			#require(MASS)
			if(zero.correction== "zq1"){
				for(t in 1:n){
					if(t<(p+1)){
						m[t]<- Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rnegbin(1,lambda[t], theta)
						y.transformed[t]<-max(c, y[t])
					}
					else{
						yphi<-t(log(y.transformed[(t-p):(t-1)]))%*%rev(phi)
						Xbetaphi<-t(-Xbeta[(t-p):(t-1)])%*%rev(phi)
						m[t]<- yphi + Xbetaphi + Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rnegbin(1,lambda[t], theta)
						y.transformed[t]<-max(c, y[t])
					}
				}
			}
			if(zero.correction== "zq2"){
				for(t in 1:n){
					if(t<(p+1)){
						m[t]<- Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rnegbin(1,lambda[t], theta)
						y.transformed[t]<- y[t]+c
					}
					else{
						yphi<-t(log(y.transformed[(t-p):(t-1)]))%*%rev(phi)
						Xbetaphi<-t(-log(exp(Xbeta[(t-p):(t-1)])+c))%*%rev(phi)
						m[t]<- yphi + Xbetaphi + Xbeta[t]
						lambda[t]<-exp(m[t])
						y[t]<-rnegbin(1,lambda[t], theta)
						y.transformed[t]<-y[t]+ c
					}
				}
			}
		}
	}
	return(y)
}

