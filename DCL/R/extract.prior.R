## The 3CL parameter estimation ##################### 
extract.prior<-function(Xtriangle,Npaid,Ntriangle,Plots=TRUE,n.cal=NA,
                        Fj.X=NA,Fj.N=NA,Fj.Npaid=NA)
{
  Xtriangle<-as.matrix(Xtriangle)
  Ntriangle<-as.matrix(Ntriangle)
  Npaid<-as.matrix(Npaid)
  m<-nrow(Xtriangle)
  # m is the dimension of the triangle and d=m-1 is the maximum delay
  
  ## 1. Estimate the cl parameters
  clm.N<-clm(Ntriangle,n.cal=n.cal,Fj=Fj.N)
  alpha.N<-clm.N$alpha
  beta.N<-clm.N$beta
  Nhat<-as.matrix(clm.N$triangle.hat)
  
  clm.X<-clm(Xtriangle,n.cal=n.cal,Fj=Fj.X) 
  alpha.X<-clm.X$alpha
  beta.X<-clm.X$beta
  Xhat<-as.matrix(clm.X$triangle.hat)
  
  clm.R<-clm(Npaid,n.cal=n.cal,Fj=Fj.Npaid) 
  beta.R<-clm.R$beta
  alpha.R<-clm.R$alpha    
  
  #### 3. The inflation parameters.
  ##  Inflation in the development (j) i.e. mu_{j-l,l}= mu_j =beta.X/beta.R 
  inflat.j <- beta.X/beta.R
  inflat.j[beta.X==beta.R]<-NA
  
  ## The probability of zero-claims = alpha.R/alpha.N
  Qi<-1-(alpha.R/alpha.N)
  
  if (Plots==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
    # c(bottom, left, top, right)
    
    plot(1:m,inflat.j,type='b',lwd=2,lty=1,col=4,pch=19,cex=0.8,
         xlab='development period',ylab='',cex.main=1,
         main='Severity inflation',
         xaxt='n')
    axis(1:m,as.character(0:(m-1)))
    grid(ny=5,nx=NA)
    
    plot(1:m,Qi,type='b',lwd=2,col=4,pch=19,ylab='',cex=0.8,
         xlab='underwriting period',xaxt='n',cex.main=1,
         main='Probability of zero-claims')
    axis(1:m,as.character(1:m))
    grid(ny=5,nx=NA)
    par(mfrow=c(1,1))
  }  
  
  return(list(inflat.j=inflat.j,Qi=Qi))
}
