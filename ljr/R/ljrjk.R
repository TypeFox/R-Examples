ljrjk <- function(j,k,y,n,tm,X=NULL,ofst=0,R=1000,alpha=0.05)
{
if (sum(duplicated(tm))>0){
 cat("Error: duplicated observation times\n")  
}
else{
 if ((is.null(X)==FALSE)&&(is.matrix(X)==FALSE))
  X=as.matrix(X)
 if (is.unsorted(tm)==TRUE){
  o=order(tm)
  y=y[o]
  n=n[o]
  tm=tm[o]
  if (is.null(X)==FALSE)
   X=as.matrix(X[o,])
  if ((length(ofst)>1)||(ofst[1]!=0))
   ofst=ofst[o]
 }
 N=length(y)
 m=ncol(X)
 if (is.null(m)) 
  m=0
 or.ofst=ofst
 if (length(ofst)==1)
  ofst=as.double(rep(ofst,N))
 else
  ofst=as.double(ofst)
 if (j==k){
  k=j+1
  cat("Warning: j should not equal k!\n")
 }
 if (j>k){
  i=j
  j=k
  k=i
 }
 j=as.integer(j)
 k=as.integer(k)
 out1 <- .C("testjk",j,k,as.double(y),as.double(n),as.double(tm),as.double(X),ofst,N,m,as.integer(R),p=double(1),PACKAGE="ljr")
 cat("Testing H0: ",j," joinpoint(s) vs. H1: ",k," joinpoints\n")
 cat("p-value=",out1$p,"\n")
 if (out1$p>alpha){
  cat("Null hypothesis is not rejected\n\n")
  out2 <- .C("ljrk",j,as.double(y),as.double(n),as.double(tm),as.double(X),ofst,beta=double(m+1),gamma=double(j+1),tau=double(j),N,m,ans=double(1),PACKAGE="ljr")
  cat('Model:\n')
  cat('y~Binom(n,p) where p=invlogit(eta)\n')
  if (m>0){
   if (is.null(dimnames(X)[[2]])==FALSE)
    m.variables=dimnames(X)[[2]]
   else
    m.variables=paste("X",1:ncol(X),sep="")
  }
  else
   m.variables=NULL
  m.variables=c('Intercept',m.variables)
  t.variables='t'
  if (length(out2$tau)>0)
   for (i in 1:length(out2$tau))
    t.variables=c(t.variables,paste('max(t-tau',i,',0)',sep=""))
  if ((or.ofst[1]==0)&(length(or.ofst)==1))
   cat('eta=b0')
  else
   cat('eta=ofst+b0')
  if (length(m.variables)>1)
   for (i in 2:length(m.variables))
    cat(paste(paste('+b',i-1,sep=""),'*',m.variables[i],sep=""))
  for (i in 0:(length(t.variables)-1))
   cat(paste(paste('+g',i,sep=""),'*',t.variables[i+1],sep=""))
  cat("\n\n")
  m.coef=c(out2$beta,out2$gamma)
  f1=data.frame(Variables=c(m.variables,t.variables),Coef=m.coef,row.names=c(paste('b',0:(length(m.variables)-1),sep=""),paste('g',0:(length(t.variables)-1),sep="")))
  print(f1)
  ret1=m.coef
  names(ret1)=c(m.variables,t.variables)
  return(list(Coef=ret1,wlik=out2$ans,pval=out1$p))
 } 
 else{
  cat("Null hypothesis is rejected\n\n")
  out2 <- .C("ljrk",k,as.double(y),as.double(n),as.double(tm),as.double(X),ofst,beta=double(m+1),gamma=double(k+1),tau=double(k),N,m,ans=double(1),PACKAGE="ljr")
  cat('Model:\n')
  cat('y~Binom(n,p) where p=invlogit(eta)\n')  
  if (m>0){
   if (is.null(dimnames(X)[[2]])==FALSE)
    m.variables=dimnames(X)[[2]]
   else
    m.variables=paste("X",1:ncol(X),sep="")
  }
  else
   m.variables=NULL
  m.variables=c('Intercept',m.variables)
  t.variables='t'
  for (i in 1:length(out2$tau))
   t.variables=c(t.variables,paste('max(t-tau',i,',0)',sep=""))
  if ((or.ofst[1]==0)&(length(or.ofst)==1))
   cat('eta=b0')
  else
   cat('eta=ofst+b0')
  if (length(m.variables)>1)
   for (i in 2:length(m.variables))
    cat(paste(paste('+b',i-1,sep=""),'*',m.variables[i],sep=""))
  for (i in 0:(length(t.variables)-1))
   cat(paste(paste('+g',i,sep=""),'*',t.variables[i+1],sep="")) 
  cat("\n\n")
  m.coef=c(out2$beta,out2$gamma)
  f1=data.frame(Variables=c(m.variables,t.variables),Coef=m.coef,row.names=c(paste('b',0:(length(m.variables)-1),sep=""),paste('g',0:(length(t.variables)-1),sep="")))
  j.labels=NULL
  if (length(out2$tau)>0)
   for (i in 1:length(out2$tau))
    j.labels=c(j.labels,paste('tau',i,'=',sep=""))
  f2=data.frame(j.labels,Joinpoint=out2$tau)
  if (length(out2$tau)>0)
   names(f2)=c(' ','  ') 
  else 
   names(f2)=NULL
  print(f1)
  cat('\nJoinpoints:\n')
  print(f2)
  ret1=m.coef
  names(ret1)=c(m.variables,t.variables)
  ret2=out2$tau
  names(ret2)=j.labels
  return(list(Coef=ret1,Joinpoint=ret2,wlik=out2$ans,pval=out1$p))
 }
}
}
