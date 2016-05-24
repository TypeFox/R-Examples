######################################### Ordinal Procedure
kfweOrd<-function(p,k=1,alpha=.01,ord=NULL,alpha.prime=alpha,J=qnbinom(alpha,k,alpha.prime),disp=TRUE,GD=FALSE){

  if(!is.null(ord)){# sort by ord
    o <- order(ord,decreasing=T)
  }
  else{o <- 1:length(p)}
  ps <- p[o]

  if(GD) alpha1 <- k*alpha.prime/(J+k)
  else alpha1 <- alpha.prime
  
  u<-cumsum(ps>alpha1)

  if(sum(u<=J)>0){
  	h <- rep(0,length(p))
  	h[1:max(which(u<=J))] <- 1
  	h[ps>alpha1] <- 0

        if(sum(h)<k)  h[(h==0)&&(ps<=alpha1)][1:min(sum((h==0)&&(ps<=alpha1)),k-1-sum(h))]=1
  	h[o] <- h
  }
  else{h <- rep(0,length(p))}
  
  
  
  if(disp==T) 	cat(paste("Ordered k-FWER procedure\n ",length(p)," tests, k=", k, ", alpha=",alpha,", individual alpha threshold=",round(alpha1,digits=7),"\n ",J," jumps allowed","\n ",sum(h)," rejections\n\n",sep=""))
  return(h==1)}
  
  
################kFWE by Guo e Romano - (Step Down and Single Step) procedure
kfweGR<-function(p,k=1,alpha=.01,disp=TRUE,SD=TRUE,const=10,alpha.prime=getAlpha(k=k,s=length(p),alpha=alpha,const=const)) {
  
  if(is.null(alpha.prime)) alpha.prime=getAlpha(k=k,s=length(p),alpha=alpha,const=50)
   
  rej <- rep(0,length(p))
  rej[p<=alpha.prime] <- 1
  n.rej <- sum(p<=alpha.prime)
  
  sd=((n.rej>=k)&(SD))
  while (sd){
    alpha.prime=getAlpha(k=k,s=length(p)-n.rej+k-1,alpha=alpha,const=const)
    rej[p<=alpha.prime] <- 1
    sd=n.rej<sum(rej)
    n.rej=sum(p<=alpha.prime)
  }
  
  if(disp) 	cat(paste("Guo and Romano k-FWER ",switch(SD,"Step Down ",""),"procedure\n ",length(p)," tests, k=", k, ", alpha=",alpha, "\n ",round(alpha.prime,digits=7)," individual alpha threshold\n ",n.rej," rejections\n\n",sep=""))
  return(rej==1)
}
  
############kFWE by Lehamann and Romano
kfweLR <- function(p,k=1,alpha=0.01,disp=TRUE) {
  s <- length(p)
  sdconst <- rep(1,s)
  
  sdconst[1:min(k,s)] <- k*alpha/s
  if(s>k) sdconst[(k+1):s] <- k*alpha/(s+k-((k+1):s))
  
  ps <- sort(p)
  u <- ps<sdconst
  res <- 0
  if(any(u)) {
  	w <- min(which(!u))-1
  	res <- ps[w]}
  p[which(p>res)] <- 1
  p[p<=alpha] <- 0
  h=(!p)
  if(disp) 	cat(paste("Lehmann e Romano k-FWER Step Down procedure\n ",length(p)," tests, k=", k, ", alpha=",alpha, "\n ",sum(h)," rejections\n\n",sep=""))
  return(h==1)
}


############### other functions

getAlpha <- function(s,k=1,alpha=.01,const=10){
  start<-1E-8
  stop<-1
  
  delta<-1
  while(delta>1E-8){
  	alphas<-start+((0:const)/(const))*(stop-start)
  	temp<-round(pbinom(k-1,s,alphas,lower.tail=F),digits=7)-alpha
  	temp2<-max(which(temp<=0))
  	alpha.prime<-alphas[temp2]
  	start<-alpha.prime
  	stop<-alphas[temp2+1]
  	delta<-abs(temp[temp2])
  }
  
  return(alpha.prime)
}

