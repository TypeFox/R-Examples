MAMSE=function(x,surv=FALSE,ub=NULL,lb=0){

 if(!is.list(x)){ stop("You must provide a list of samples") }

 if(is.null(dim(x[[1]]))){ 
   if(min(sapply(x,is.numeric))==0 | sum(sapply(x,dim)==NULL)>0){ 
     stop("All univariate samples must be numeric vectors.")
   }
   return(MAMSEunipo(x)) 
 }
 if(surv==TRUE){
   if(sum(sapply(x,dim)[2,]!=2)>0){ stop("The sample from each population must be in two columns: (1) value (2) indicator (1 or TRUE) = observed.") }
   if(is.null(ub)){ stop("You need to specify an upper bound (ub).") }
   return(MAMSEsurvpo(x,lb=lb,ub=ub))
 }
 if(min(sapply(x,is.numeric))==0){
   stop("Samples must be numerical matrices (or data.frames).")
 } 
 if(var(sapply(x,dim)[2,])>0){
   stop("All samples must have the same number of dimensions.")
 }
 
  y=lapply(x,function(z){apply(z,2,ranked)})
  return(MAMSEmultipo(y))

}

MAMSEunipo=function(x){

  if(length(x)==1) return(1)
  a=MAMSEuni(x)
  if(min(a)<0){
      y=x
      i=(1:length(x))[a<0]
      for(j in i){y[[j]]=NULL}
      a[i]=0
      a[-i]=MAMSEunipo(y)
  }
  a
}




MAMSEuni=function(x){

# New proposition to determine weights. Weights are provided
# for inference on population 1.
# Uses a counting measure on data from population 1
#
# Input:  A list of datasets (vectors)
# Output: The vector of weights


  m=length(x)

  n=sapply(1:m,function(i){length(x[[i]])})
  y=sort(x[[1]])
  N=n[1]
  
  comp=function(x,y,nx,ny){
 
  # x and y must be sorted; 
  # nx and ny are the length of the datasets (vectors) x and y)
  # Used for the calculation of MAMSE weights.
  #  
  # Output: b[i] = Fx(y[i]), the empirical CDF of X evaluated at y[i]
 
     .C("comp",
        as.numeric(x),
        as.numeric(y),
        as.integer(nx),
        as.integer(ny),
        out=integer(ny),
	PACKAGE="MAMSE")$out

  }

b=t(cbind((1:n[1]),sapply(2:m,function(i,x,y,n){comp(sort(x[[i]]),y,n[i],n[1])},n=n,x=x,y=y)))/n
  A=matrix(apply(b,2,function(b,n,m){
(b[1]-b[-1])%*%t(b[1]-b[-1])+
b[1]*(1-b[1])/n[1]+
diag(b[-1]*(1-b[-1])/n[-1],nrow=m-1)},n=n,m=m),(m-1)^2,n[1])
 A=matrix(apply(A,1,mean),m-1,m-1)
 d=rep(mean(b[1,]*(1-b[1,]))/n[1],m-1)
 
if(m>2){
 for(i in 1:(m-2)){
   for(j in (i+1):(m-1)){
     if(sum(A[i,]==A[j,])==m-1){
       add=rep(0,m-1)
       add[i]=1
       add[j]=-1
       A[j,]=add
       d[j]=0
     }
   }
 }
}
 l=solve(A,d)
 c(1-sum(l),l)
}

# Version for survival functions

MAMSEsurv=function(x,ub,lb){

# Quickly revised post-comp

# New proposition to determine weights. Weights are provided
# for inference on population 1.
# Uses the empirical measure on data from population 1
#
# Input:  x=list, each element=matrix: column 1 is data, column 2 is censoring
# Output: The vector of weights

  m=length(x)

  x=lapply(x,function(y){y[sort.list(y[,1]),]})
  n=unlist(lapply(x,function(y){dim(y)[1]}))
  
  Dt=x[[1]][x[[1]][,2]==1,1]
  w1=dS(x[[1]])
  w1=w1[x[[1]][,2]==1]
  w1=w1[Dt>=lb & Dt<=ub]
  Dt=Dt[Dt>=lb & Dt<=ub]
  N=length(Dt)
  
  Fi=matrix(unlist(lapply(x,KM,t=Dt)),N,m)
  sig=matrix(0,N,m)
  for(i in 1:m){
    sig[,i]=wvar(x[[i]],Dt,Fi[,i])
  }


  b=Fi
  A=matrix(sapply(1:length(Dt),function(i,b,n,m,sig){
   (b[i,1]-b[i,-1])%*%t(b[i,1]-b[i,-1])+
   sig[i,1]+
   diag(sig[i,-1],nrow=m-1)
   },n=n,m=m,sig=sig,b=b),(m-1)^2,length(Dt))
 A=matrix(apply(t(A)*w1,2,sum),m-1,m-1)

 d=rep(sum(sig[,1]*w1),m-1)
 
if(m>2){
 for(i in 1:(m-2)){
   for(j in (i+1):(m-1)){
     if(sum(A[i,]==A[j,])==m-1){
       add=rep(0,m-1)
       add[i]=1
       add[j]=-1
       A[j,]=add
       d[j]=0
     }
   }
 }
}
 l=solve(A,d)
 c(1-sum(l),l)
}

MAMSEsurvpo=function(x,ub,lb=0){

  m=length(x)
  if(m==1){ return(1) }
  
  if(sum(((x[[1]][,1]>=lb & x[[1]][,1]<=ub))&(x[[1]][,2]==1))<2){ 
    warning("Too few data points from Population 1 fall in the interval of interest.")
    return(c(1,rep(0,length(x)-1)))
  }
  
  MX=unlist(lapply(x,function(y){max(y[y[,2]==1,1])}))
  MN=unlist(lapply(x,function(y){min(y[y[,2]==1,1])}))

  a=rep(0,m)
  if(min(MX[-1])<ub||max(MN[-1])>=ub){
      y=x
      rem=(MX[-1]<ub)|(MN[-1]>=ub)
      for(j in sort((2:m)[rem],decreasing=TRUE)){y[[j]]=NULL}
      a[c(FALSE,rem)]=0
      a[c(TRUE,!rem)]=MAMSEsurvpo(y,ub,lb)
      return(a)
  }
 
  a=MAMSEsurv(x,ub,lb)
  if(min(a)<0){
      y=x
      rem=a[-1]<0
      for(j in sort((2:m)[rem],decreasing=TRUE)){y[[j]]=NULL}
      a[c(FALSE,rem)]=0
      a[c(TRUE,!rem)]=MAMSEsurvpo(y,ub,lb)
  }
  a
}

dS=function(x){

#  x=matrix: column 1 is data (must be sorted), column 2 is censoring
#  Output, weight of each data point in the K-M estimate
#  Note: cumsum(w) gives the CDF.

  n=length(x[,1])
  w=rep(1/n,n)

  for(i in 1:(n-1)){
    if(x[i,2]==0) {
      w[(i+1):n]=w[(i+1):n]+w[i]/(n-i)
      w[i]=0
    }
  }  
  if(w[n]==0){w[n]=0}
  w 
}

KM=function(x,time){
 # x = matrix: column 1 is data (must be sorted), column 2 is censoring
 # t = vector of times at which to evaluate KM
 # Output: CDF based on K-M estimate of the survival function
 
 csum(cbind(x[,1],dS(x)),time)

}

WKME=function(x,ub,lb=0,time=NULL,boot=NULL,REP=1000){

# x = sample, 2 cloumns, values + indicator
# boot = if bootstrap intervals are required, level \in (0,1)
# time = points on the line where to evaluate the curve
# OUT:
#  x= sorted values 
#  weight= weigth given to each value
#  km = values of KM at points time
#  time= sorted vector of times
#  CI = the value of the pointwise CI at "time"

  m=length(x)
  x=lapply(x,function(x){x[sort.list(x[,1]),]})
  weight=lapply(x,dS)
  km=NULL
  MAMSEKM=NULL
  a=list(NULL,NULL)
  if(!is.null(time)){ 
      w=MAMSEsurvpo(x,ub=ub)
      km=matrix(unlist(lapply(x,KM,t=time)),length(time),m)
      MAMSEKM=km%*%w
  }
      
  if(!is.null(boot)){
    if(boot<=0 | boot>=1){ stop("Parameter boot must be in (0,1) and represents the level of the CI.") }
    a=bootx(x,time=time,lb=lb,ub=ub,REP=REP,boot=boot)
  }
  return(list(x=lapply(x,function(y){y[,1]}),weight=weight,kme=km,time=time,kmeCI=a[[1]],wkme=MAMSEKM,wkmeCI=a[[2]]))
}


bootx=function(x,time,lb,ub,REP=1000,boot=0.95){

  m=length(x)
  n=unlist(lapply(x,length))/2
  
  Xw=lapply(x,dS)
  
  y=x
  for(i in 1:m){y[[i]][,2]=!y[[i]][,2]}
  Yw=lapply(y,dS)
  
  z=as.list(1:m)
  for(i in 1:m){
    Xw[[i]]=c(Xw[[i]],max(0,1-sum(Xw[[i]])))
    Yw[[i]]=c(Yw[[i]],max(0,1-sum(Yw[[i]])))
    z[[i]]=c(y[[i]][,1],Inf)
  }
  
 
samX=lapply(as.list(1:m),function(i,n,REP,prob){matrix(z[[i]][sample(1:(n[i]+1),REP*n[i],replace=TRUE,prob=prob[[i]])],n[i],REP)},n=n,REP=REP,prob=Xw)
 
samY=lapply(as.list(1:m),function(i,n,REP,prob){matrix(z[[i]][sample(1:(n[i]+1),REP*n[i],replace=TRUE,prob=prob[[i]])],n[i],REP)},n=n,REP=REP,prob=Yw)

  KMN=NULL
  KMW=NULL

  Z=as.list(1:m)
  for(i in 1:REP){
    for(j in 1:m){
      Z[[j]]=cbind(pmin(samX[[j]][,i],samY[[j]][,i]),samX[[j]][,i]<=samY[[j]][,i])
    }
    Z=lapply(Z,function(x){x[sort.list(x[,1]),]})
 
      w=MAMSEsurvpo(Z,lb=lb,ub=ub)

      km=matrix(unlist(lapply(Z,KM,t=time)),length(time),m)
      MAMSEKM=km%*%w
      KMN=c(KMN,km[,1])
      KMW=c(KMW,MAMSEKM)
    
  }
  
  KMN=matrix(KMN,length(time),REP)
  KMW=matrix(KMW,length(time),REP)
  
  CIN=apply(KMN,1,quantile,probs=c((1-boot)/2,1-(1-boot)/2))
  CIW=apply(KMW,1,quantile,probs=c((1-boot)/2,1-(1-boot)/2))  
  
  list(CIN,CIW)  
}
 

csum=function(x,y){
 # x=matrix, col 1 = time (must be sorted), col 2 = values
 # y=vector of times at which to calculate partial cusums
 # output: sum of the x[,2] such that x[,1] <= y
 
  i=1;j=1
  lx=length(x[,1])
  ly=length(y)
  out=rep(0,ly)
  
  while(i<=ly){
    while(j<=lx && x[j,1]<=y[i]){
      out[i]=out[i]+x[j,2]
      j=j+1
    }
    if(i<ly){out[i+1]=out[i]}
    i=i+1
  }
  out
}

wvar=function(x,time,Fi=NULL){
  # x = matrix: column 1 is data, column 2 is censoring
  # Fi= KM calculated at Dt[1]. Passed as an argument because it
  #     is already calculated in MAMSEsurv.
  # Dt= times at which to calculate var
  # output: weights for each death time in the sum used 
  #         for calculating the variance (only the 
  #         Greenwood formula, not \sigma^2).
  
  Dt=x[x[,2]==1,1]
  allt=x[,1]
  w=rep(0,length(Dt))
  if(is.null(Fi)){Fi=KM(x,time)}
  
  for(i in 1:length(Dt)){
    atrisk=sum(allt>Dt[i])
    w[i]=ifelse(atrisk>0,1/(sum(allt>=Dt[i])*sum(allt>Dt[i])),0)
  }
  
  v=sapply(time,function(y,Dt,w){sum(w[Dt<y])},Dt=Dt,w=w)
  (1-Fi)^2*v
}

# Programs for positively constrained weights
# -------------------------------------------


MAMSEmultipo=function(x){

  if(length(x)==1) return(1)
  a=MAMSEmulti(x)
  if(min(a)<0){
      y=x
      i=(1:length(x))[a<0]
      for(j in i){y[[j]]=NULL}
      a[i]=0
      a[-i]=MAMSEmulti(y)
  }
  a
}


# General programs
# ----------------


MAMSEmulti=function(x){

# MAMSE for multivariate data

  m=length(x)

  n=sapply(1:m,function(i){dim(x[[i]])})
  if(max(n[2,])>min(n[2,])){stop("Populations have data of different dimensions.")}
  d=n[2,1]
  n=n[1,]
  N=n[1]^d


  compmulti=function(x,y,nx,ny,p){

  # x and y are data sets nx*p and ny*p;
  # Calculates the empirical functions using comp.c
  #
  # Output: b[i] = Fx(y[i]), the empirical CDF of X evaluated at y[i]


     .C("compmulti",
        as.numeric(c(x,0)),
        as.numeric(c(y,0)),
        as.integer(nx),
        as.integer(ny),
        as.integer(p),
        out=integer(ny),
	PACKAGE="MAMSE")$out/nx

  }


  gridof=function(x){

  # Creates the grid of points of the cross-product of the sample x.

    d=dim(x)
    if(length(d)>2){stop("This is not a sample.")}
    n=d[1]
    d=d[2]
    if(d==1){stop("This sample is univariate. Use the appropriate function.")}

    out=matrix(0,n^d,d)
  
    out[,1]=rep(sort(x[,1]),n^(d-1))
    for(i in 2:d){
      out[,i]=rep(sort(x[,2]),each=n^(i-1))
    }
    out
  }

  
  b=matrix(0,m,N)
  X=gridof(x[[1]])

  for(i in 1:m){
    b[i,]=compmulti(x[[i]],X,n[i],N,d)
  }

    A=matrix(apply(b,2,function(b,n,m){
(b[1]-b[-1])%*%t(b[1]-b[-1])+
b[1]*(1-b[1])/n[1]+
diag(b[-1]*(1-b[-1])/n[-1],nrow=m-1)},n=n,m=m),(m-1)^2,N)
 A=matrix(apply(A,1,mean),m-1,m-1)
 d=rep(mean(b[1,]*(1-b[1,]))/n[1],m-1)
 
if(m>2){
 for(i in 1:(m-2)){
   for(j in (i+1):(m-1)){
     if(sum(A[i,]==A[j,])==m-1){
       add=rep(0,m-1)
       add[i]=1
       add[j]=-1
       A[j,]=add
       d[j]=0
     }
   }
 }
}
 l=solve(A,d)
 c(1-sum(l),l)
} 
 
ranked=function(x){ rank(x)/(length(x)+1) }

