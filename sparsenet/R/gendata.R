gendata=function(N,p,nonzero,rho,snr=3,alternate=TRUE){
 if(nonzero>0){
     beta=rev(seq(from=1,to=nonzero))
     if(alternate)beta=beta*(rep(c(1,-1),nonzero/2+1)[seq(nonzero)])
     beta=c(beta,rep(0,p-nonzero))

 }
     else beta=rep(0,p)
     x=matrix(rnorm(N*p),N,p)
     if(rho>0){
         z=rnorm(N)
         x=x+sqrt(rho/(1-rho))*z
         x=x*sqrt(1-rho)
     }
     signal.sd=sqrt((1-rho)*sum(beta^2)+rho*(sum(beta)))
     noise.sd=signal.sd/snr
     y=x%*%beta+rnorm(N)*noise.sd
     list(x=x,y=y,beta=beta,signal.sd=signal.sd,noise.sd=noise.sd,rho=rho,snr=snr)
 }



