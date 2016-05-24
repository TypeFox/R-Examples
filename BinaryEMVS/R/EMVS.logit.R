EMVS.logit=function(y,x,epsilon=.0005,v0s=5,nu.1=1000,nu.gam=1,lambda.var=.001,a=1,b=ncol(x),
                    beta.initial=rep(1,p),sigma.initial=1,theta.inital=.5,temp=1,p=ncol(x),n=nrow(x),SDCD.length=50){
  
  if(length(beta.initial)==0){
    beta.initial=rep(1,p)
  }
  L=length(v0s)   
  
  cat("\n")
  
  cat("\n","Running Logit across v0's","\n")
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
    nu.0=v0s[i]
    
    beta.Current=beta.initial
    beta.new=beta.initial
    sigma.EM=sigma.initial
    theta.EM=theta.inital
    
    eps=epsilon+1
    iter.index=1
    while(eps>epsilon && iter.index<20){
      d.Star=rep(NA,p)
      p.Star=rep(NA,p)
      for(j in 1:p){
        
        gam.one=dnorm(beta.Current[j],0,sigma.EM*sqrt(nu.1))**temp*theta.EM**temp
        gam.zero=dnorm(beta.Current[j],0,sigma.EM*sqrt(nu.0))**temp*(1-theta.EM)**temp
        
        p.Star[j]=gam.one/(gam.one+gam.zero)
        d.Star[j]=((1-p.Star[j])/nu.0)+(p.Star[j]/nu.1)
      }
      
      #cat("max p.Star", max(p.Star),"\n")
      #cat("d.Star.EM:  ", d.Star[1:5],"\n")
      ############### M STEP #######################
      
      d.Star.Mat=diag(d.Star,p)
      beta.Current=rep(NA,p)
      count.while=0
      while(is.na(min(beta.Current))){
        beta.Current=CSDCD.logistic(p,n,x,y,d.Star,SDCD.length)
        count.while=count.while+1
        #cat("This is count.while:",count.while,"\n")
      }
      
      ######## VARIANCE FORUMULA IS DIFFERENT FROM CONTINUOUS AND PROBIT CASE ###########
      #sigma.EM[i]=sqrt((sum(log(1+exp(-y*x%*%beta.EM[i,])))+sum((sqrt(d.Star.Mat)%*%beta.EM[i,])**2)+lambda.var*nu.gam)/(n+p+nu.gam))
      sigma.EM=sqrt((sum((sqrt(d.Star.Mat)%*%beta.Current)**2)+lambda.var*nu.gam)/(n+p+nu.gam+2))
      theta.EM=(sum(p.Star)+a-1)/(a+b+p-2)
      eps=max(abs(beta.new-beta.Current))
      #print(eps)
      beta.new=beta.Current
      iter.index=iter.index+1  
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
            niters=index.Vec,posts=p.Star.Vec,thetas=theta.Vec,v0s=v0s)
  return(list)
  
}