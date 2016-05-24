EMVS.probit=function(y,x,epsilon=.0005,v0s=.025,nu.1=100,nu.gam=1,a=1,b=ncol(x),beta.initial=NULL,
                     sigma.initial=1,theta.inital=.5,temp=1,p=ncol(x),n=nrow(x)){  #set reasonable defauly for v0
  
  if(length(beta.initial)==0){
    beta.initial=rep(NaN,times=ncol(x))
    while(sum(is.nan(beta.initial))>0)
    {
      beta.initial=CSDCD.logistic(p,n,x,ifelse(y==1,1,-1),rep(.0001,p),75)
    }
  }
  
  scrap=0
  cat("\n")
  
  L=length(v0s)     
  cat("\n","Running Probit across v0's","\n")
  cat(rep("",times=(L+1)),sep="|")
  cat("\n")
  
  intersects=numeric(L)    # intersection points between posterior weighted spike and slab
  log_post=numeric(L)      # logarithm of the g-function models associated with v0s
  sigma.Vec=numeric(L)
  theta.Vec=numeric(L)
  log_post=numeric(L) 
  index.Vec=numeric(L)
  beta.Vec=matrix(0,L,p)    # L x p matrix of MAP beta estimates for each spike  
  p.Star.Vec=matrix(0,L,p)    # L x p matrix of conditional posterior inclusion probabilities
  
  for (i in (1:L)){
    bail.count=0
    nu.0=v0s[i]
    beta.Current=beta.initial
    beta.new=beta.initial
    sigma.EM=sigma.initial
    theta.EM=theta.inital
    
    eps=epsilon+1
    iter.index=1
    while(eps>epsilon && iter.index<20){ #
      
      if(bail.count>=3)
      {
        cat("\n","Iteration scrapped!","\n",sep="")
        list=list(betas=NULL,intersects=NULL,sigmas=NULL,
                  niters=NULL,posts=NULL,thetas=NULL,v0s=NULL,scrap=1)
        return(list)
        
      }
      
      #print(iter.index)
      ############### E STEP #######################
      
      d.Star=rep(NA,p)
      p.Star=rep(NA,p)
      
      for(j in 1:p){
        
        gam.one=dnorm(beta.Current[j],0,sqrt(nu.1))**temp*theta.EM**temp
        gam.zero=dnorm(beta.Current[j],0,sqrt(nu.0))**temp*(1-theta.EM)**temp
        
        p.Star[j]=gam.one/(gam.one+gam.zero)
        d.Star[j]=((1-p.Star[j])/nu.0)+(p.Star[j]/nu.1)
      }
      
      
      ############# LATENT VARIABLE FOR PROBIT ##################
      
      M=rep(NA,n)
      y.Star=rep(NA,n)
      for(j in 1:n){
        M[j]=ifelse(y[j]==0,-dnorm(-x[j,]%*%beta.Current,0,1)/pnorm(-x[j,]%*%beta.Current,0,1),
                    dnorm(-x[j,]%*%beta.Current,0,1)/(1-pnorm(-x[j,]%*%beta.Current,0,1)))
        y.Star[j]=x[j,]%*%beta.Current+M[j]
        
      }
      
      
      ############### M STEP #######################
      
      d.Star.Mat=diag(d.Star,p)
      beta.old=beta.new
      
      ######## Sherman-Morrison-Woodbury Formula ####################
      #beta.EM[i,]=(solve(v.Star.Mat)-(solve(v.Star.Mat)%*%t(x))%*%solve(diag(1,n)+x%*%solve(v.Star.Mat)%*%t(x))%*%(x%*%solve(v.Star.Mat)))%*%(t(x)%*%y.Star)
      #beta.Current=solve(t(x)%*%x+d.Star.Mat)%*%t(x)%*%y.Star
      beta.Current=CSDCD.random(p,n,x,y.Star,d.Star)
      #print(beta.Current[1:10])
      
      if(sum(is.nan(beta.Current))==0)
      {
        
        sigma.EM=1
        theta.EM=(sum(p.Star)+a-1)/(a+b+p-2)
        
        eps=max(abs(beta.new-beta.Current))
        
        beta.new=beta.Current
        iter.index=iter.index+1
        
      }
      if(sum(is.nan(beta.Current))>0){
        beta.Current=beta.old
        bail.count=bail.count+1
      }
    }
    
    p.Star.Vec[i,]=p.Star
    beta.Vec[i,]=beta.new
    sigma.Vec[i]=sigma.EM
    theta.Vec[i]=theta.EM
    
    index.Vec[i]=iter.index
    
    index=p.Star>0.5
    
    c=sqrt(nu.1/v0s[i])
    
    w=(1-theta.Vec[i])/theta.Vec[i]  
    
    if (w>0){
      intersects[i]=sigma.Vec[i]*sqrt(v0s[i])*sqrt(2*log(w*c)*c^2/(c^2-1))}else{
        intersects[i]=0}
    
    cat("|",sep="")
    
  }
  
  list=list(betas=beta.Vec,intersects=intersects,sigmas=sigma.Vec,
            niters=index.Vec,posts=p.Star.Vec,thetas=theta.Vec,v0s=v0s,scrap=0)
  return(list)
  
}