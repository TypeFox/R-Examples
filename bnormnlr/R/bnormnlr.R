library(mvtnorm)
require(numDeriv)
bnlr<-function(y,f1,f2,f1g=NULL,f2g=NULL,x,z,bta0,gma0,b=rep(0,length(bta0)),B=diag(10^6,length(bta0)),
               g=rep(0,length(gma0)),G=diag(10^6,length(gma0)),Nc){
  ###########
  ##
  ###########
  if(!(is.function(f1)) | !(is.function(f2))){stop("f1 and f2 must be functions")}
  if(!is.null(f1g)){if(!(is.function(f1))){stop("f1g must be a function")}}
  if(!is.null(f2g)){if(!(is.function(f2))){stop("f2g must be a function")}}
  
  if(is.null(dim(x))){if(length(y)!=length(x)){stop("length(y)!=length(x)")}}else
  {if(dim(x)[1]!=length(y)){stop("dim(x)[1]!=length(y)")}}
  if(is.null(dim(z))){if(length(y)!=length(z)){stop("length(y)!=length(z)")}}else
  {if(dim(z)[1]!=length(y)){stop("dim(z)[1]!=length(y)")}}
  
  if(!is.null(f1)){if(length(f1(bta0,x))!=length(y)){stop("f1(bta0,x) != length(y)")}}
  if(!is.null(f2)){if(length(f2(gma0,z))!=length(y)){stop("f2(gma0,z) != length(y)")}}
  if(!is.null(f1g)){if(!all(dim(f1g(bta0,x))==c(length(y),length(bta0)))){stop("f1g(bta0,x) is not of appropriate dimensions")}}
  if(!is.null(f2g)){if(!all(dim(f2g(gma0,z))==c(length(y),length(gma0)))){stop("f2g(gma0,z) is not of appropriate dimensions")}}
  
  f1 <- match.fun(f1);f2 <- match.fun(f2)
  if(is.null(f1g)){f1g<-function(bta,x){jacobian(f1,x=bta,cov = x)}}
  if(is.null(f2g)){f2g<-function(bta,x){jacobian(f2,x=bta,cov = x)}}
  B.inv<-solve(B)
  B.inv_b<-solve(B)%*%b
  G.inv<-solve(G)
  G.inv_g<-solve(G)%*%g
  pb<-length(bta0);pg<-length(gma0)
  Dev<-rep(NA,Nc)
  bta.c<-matrix(rep(NA,Nc*pb),ncol=pb)
  gma.c<-matrix(rep(NA,Nc*pg),ncol=pg)
  accept.gma<-0
  accept.bta<-0
  #############
  ###
  ##############
  bta.c[1,]<-bta0
  gma.c[1,]<-gma0
  Dev[1]<- -2*sum(dnorm(y,f1(bta0,x),sqrt(f2(gma0,z)),log=T))
  ################
  ##gibbs
  ################
  for(i in 2:Nc){
    #######bta|y,gma
    X_hat<-f1g(bta.c[i-1,],x)
    mu<-f1(bta.c[i-1,],x)
    sgma<-f2(gma.c[i-1,],z)
    y1_hat<-drop(X_hat%*%bta.c[i-1,]+y-mu)
    
    B.<-solve(B.inv+t(X_hat)%*%(X_hat/sgma))
    b.<-drop(B.%*%(B.inv_b+t(X_hat)%*%(y1_hat/sgma)))
    bta.p<-drop(rmvnorm(1,b.,B.))#####proposal
    
    X_hat.p<-f1g(bta.p,x)
    mu.p<-f1(bta.p,x)
    y1_hat.p<-drop(X_hat.p%*%bta.p+y-mu.p)
    B.p<-solve(B.inv+t(X_hat.p)%*%(X_hat.p/sgma))
    b.p<-drop(B.p%*%(B.inv_b+t(X_hat.p)%*%(y1_hat.p/sgma)))
    
    r1<-exp(dmvnorm(bta.p,b,B,log=T)-dmvnorm(bta.c[i-1,],b,B,log=T)+
              sum(dnorm(y,mu.p,sqrt(sgma),log =T)-dnorm(y,mu,sqrt(sgma),log =T))+
              dmvnorm(bta.c[i-1,],b.p,B.p,log=T)-dmvnorm(bta.p,b.,B.,log=T))
    ###evaluate criteria
    if(runif(1)<r1){
      accept.bta<-accept.bta+1
      bta.c[i,]<-bta.p
    } else{
      bta.c[i,]<-bta.c[i-1,]
    }
    #######gma|y,bta
    Z_hat<-f2g(gma.c[i-1,],z)
    mu<-f1(bta.c[i,],x)
    y2_hat<-drop(Z_hat%*%gma.c[i-1,]-sgma+(y-mu)^2)
    
    G.<-solve(G.inv+t(Z_hat)%*%(.5*Z_hat/sgma^2))
    g.<-drop(G.%*%(G.inv_g+t(Z_hat)%*%(.5*y2_hat/sgma^2)))
    gma.p<-drop(rmvnorm(1,g.,G.))#####proposal
    
    Z_hat.p<-f2g(gma.p,z)
    sgma.p<-f2(gma.p,z)
    y2_hat.p<-drop(Z_hat.p%*%gma.p-sgma.p+(y-mu)^2)
    G.p<-solve(G.inv+t(Z_hat.p)%*%(.5*Z_hat.p/sgma.p^2))
    g.p<-drop(G.p%*%(G.inv_g+t(Z_hat.p)%*%(.5*y2_hat.p/sgma.p^2)))
    
    r2<-exp(dmvnorm(gma.p,g,G,log=T)-dmvnorm(gma.c[i-1,],g,G,log=T)+
              sum(dnorm(y,mu,sqrt(sgma.p),log=T)-dnorm(y,mu,sqrt(sgma),log=T))+
              dmvnorm(gma.c[i-1,],g.p,G.p,log=T)-dmvnorm(gma.p,g.,G.,log=T))
    
    if(runif(1)<r2){
      accept.gma<-accept.gma+1
      gma.c[i,]<-gma.p
    } else{
      gma.c[i,]<-gma.c[i-1,]
    }
    sgma<-f2(gma.c[i,],z)
    Dev[i]<- -2*sum(dnorm(y,mu,sqrt(sgma),log=T)) 
  }
  chains<-data.frame(bta.c,gma.c,Dev)
  colnames(chains)<-c(paste("bta",1:length(bta0),sep=""),paste("gma",1:length(gma0),sep=""),"Dev")
  return(list(chains=chains,accept.bta=accept.bta,accept.gma=accept.gma,y=y,x=x,z=z,f1=f1,f2=f2))
}

chainsum<-function(chains,q=c(0.025,.5,.975),burn=NULL){
  if(length(dim(chains))!=2){stop("chains must have 2 dim")}
  if(is.null(burn)){
    chns.sum<-sapply(1:dim(chains)[2],function(x){c(mean=mean(chains[,x],na.rm = T),quantile(chains[,x], probs = q, na.rm = T))})
    colnames(chns.sum)<-colnames(chains) 
  }else
  {
    if(all(burn>0)){
      chns.sum<-sapply(1:dim(chains)[2],function(x){c(mean=mean(chains[-burn,x],na.rm = T),quantile(chains[-burn,x], probs = q, na.rm = T))})
      colnames(chns.sum)<-colnames(chains)
    }else
    {stop("burn must be positive")}
  }
  return(chns.sum)
}

infocrit<-function(model,burn){
  chns.sum<-chainsum(model$chains,burn=burn)
  btas<-colnames(chns.sum)[grepl("bta", colnames(chns.sum))]
  gmas<-colnames(chns.sum)[grepl("gma", colnames(chns.sum))]
  
  pd<-chns.sum["mean","Dev"]+2*sum(dnorm(model$y,model$f1(chns.sum["mean",btas],model$x),
                                         sqrt(model$f2(chns.sum["mean",gmas],model$z)),log=T))
  DIC<-chns.sum["mean","Dev"]+pd
  AIC<-chns.sum["mean","Dev"]+2*length(c(btas,gmas))
  BIC<-chns.sum["mean","Dev"]+length(c(btas,gmas))*log(length(model$y))
  return(c(pd=pd,DIC=DIC,AIC=AIC,BIC=BIC))
}
