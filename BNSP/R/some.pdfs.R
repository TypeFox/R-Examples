cut.pois<-function(y,rate,E){
    qnorm(ppois(y,E*rate))
}

dnorm.pois<-function(x,y,mu,Sigma,rate,E){
    dnorm(x,mean=mu[1],sd=sqrt(Sigma[1,1]))*(
    pnorm(cut.pois(y,rate,E),mean=(mu[2]+Sigma[1,2]*(x-mu[1])/Sigma[1,1]),sd=sqrt(Sigma[2,2]-Sigma[1,2]^2/Sigma[1,1]))-
    pnorm(cut.pois(y-1,rate,E),mean=(mu[2]+Sigma[1,2]*(x-mu[1])/Sigma[1,1]),sd=sqrt(Sigma[2,2]-Sigma[1,2]^2/Sigma[1,1])))
}

