ljr1 <- function(y,n,tm,X=NULL,ofst=0,summ=TRUE)
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
   if (is.null(X)) 
    m=0
   else
    m=ncol(X)
   or.ofst=ofst
   if (length(ofst)==1)
    ofst=as.double(rep(ofst,N))
   else
    ofst=as.double(ofst)
   out <- .C("ljr1",as.double(y),as.double(n),as.double(tm),as.double(X),ofst,beta=double(m+1),gamma=double(2),tau=double(1),N,m,ans=double(1),PACKAGE="ljr")
   if (summ==TRUE){
    cat('Model:\n')
    cat('y~Binom(n,p) where p=invlogit(eta)\n')
   }
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
   for (i in 1:length(out$tau))
    t.variables=c(t.variables,paste('max(t-tau',i,',0)',sep=""))
   if (summ==TRUE)
    if ((or.ofst[1]==0)&(length(or.ofst)==1))
     cat('eta=b0')
    else
     cat('eta=ofst+b0')
   if (summ==TRUE){
    if (length(m.variables)>1)
     for (i in 2:length(m.variables))
      cat(paste(paste('+b',i-1,sep=""),'*',m.variables[i],sep=""))
     for (i in 0:(length(t.variables)-1))
      cat(paste(paste('+g',i,sep=""),'*',t.variables[i+1],sep=""))
    cat("\n\n")
   }
   m.coef=c(out$beta,out$gamma)   
   f1=data.frame(Variables=c(m.variables,t.variables),Coef=m.coef,row.names=c(paste('b',0:(length(m.variables)-1),sep=""),paste('g',0:(length(t.variables)-1),sep="")))
   j.labels=NULL
   for (i in 1:length(out$tau))
    j.labels=c(j.labels,paste('tau',i,'=',sep=""))
   f2=data.frame(j.labels,Joinpoint=out$tau)
   names(f2)=c(' ','  ')
   if (summ==TRUE){
    print(f1)
    cat('\nJoinpoints:\n')
    print(f2)
   }
   ret1=m.coef
   names(ret1)=c(m.variables,t.variables)
   ret2=out$tau
   names(ret2)=j.labels
   return(list(Coef=ret1,Joinpoint=ret2,wlik=out$ans))
  }
}
