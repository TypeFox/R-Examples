SampleNMM <- function(Y){
    nmm <- function(theta,y){
        ## Normal-Uniform mixture model
        A<-theta[1]
        mu<-theta[2]
        sigma<-theta[3]
        psigma<-log(dchisq(sigma^2,1))
        L<-sum(log(A/(max(y,na.rm=TRUE)-min(y,na.rm=TRUE))+(1-A)*dnorm(y,mu,sigma)))+psigma
        #L=sum(log(A+(1-A)*dnorm(y,mu,0.02)))
        #print(c(theta,L))
        return(-L)
    }
    a<-optim(c(0.5,0.5,0.5),nmm,y=Y)
    A<-a$par[1]
    mu<-a$par[2]
    sigma<-a$par[3]
    if(A>1)A<-1 else if(A<0)A<-0
    if(mu>1)mu<-1 else if (mu<0)mu<-0
    psigma<-log(dchisq(sigma^2,1))
    L<-sum(log(A/(max(Y,na.rm=TRUE)-min(Y,na.rm=TRUE))+(1-A)*dnorm(Y,mu,sigma)))+psigma
    BIC<--2*L+3*log(length(Y))
    a$par<-c(A,mu,sigma)
    return(list(BIC=BIC,fit=a$par))
}