BDSSanal <- function(par,times,ttype,rho=0,sprob,root=0,survival=1,maxpar=1000){ 
	if (length(ttype)==1 && ttype==0) {ttype<-times*0+1}  
	l<- par[1]
	m<-(1-sprob)*par[2]
	psi<-(sprob)*par[2]

	transmission<-times[which(ttype==1)]
	sampling<-times[which(ttype==0)]
	
	if (root==1){
		transmission<-c(transmission, max(transmission))
	}	
	extant<-length(transmission)-length(sampling)

	if (min(l,m,psi)<0  || max(l,m,psi)>maxpar ){lik<- - 10^100 } else {   
		lik<- -(root+1)*log(2*l)-survival*(root+1)*log(1-p0sersamp(max(transmission),l,m,psi,rho))
		if (extant>0){
			lik<-lik+extant*log(rho)}
		for (i in 1:length(transmission)){
			lik<- lik+log(2*l/qfunc(transmission[i],l,m,psi,rho))}
		if (length(sampling)>0){   
			for (i in 1:length(sampling)){
				lik<- lik+log(psi*qfunc(sampling[i],l,m,psi,rho))}}}
	lik <- lik - (length(transmission)-1-root)*log(2)
	-lik
}

c1 <- function(l,m,psi){sqrt((l-m-psi)^2 + 4*l*psi) }  
c2 <- function(l,m,psi,rho=0){  - (l-m-psi-2*l*rho)/c1(l,m,psi) }
qfunc <- function(t,l,m,psi,rho=0){  1/4*((1-c2(l,m,psi,rho))*exp(-t*c1(l,m,psi)/2) + (1+c2(l,m,psi,rho))*exp(t*c1(l,m,psi)/2)  )^2 }
p0sersamp <- function(t,l,m,psi,rho=0){  
	co2<-c2(l,m,psi,rho)
	co1<-c1(l,m,psi)
	res<-l+m+psi+co1*(exp(-co1*t)*(1-co2)-(1+co2))/(exp(-co1*t)*(1-co2)+(1+co2))	
	res<-res/(2*l)
	res
	}

