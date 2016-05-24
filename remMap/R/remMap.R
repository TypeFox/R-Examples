##################################
###### 12-05-07: R package "remMap"
###### R functions for fitting multivariate regression with MAP panelty: 


remMap<-function(X.m, Y.m,lamL1, lamL2, phi0=NULL, C.m=NULL)
{
##########################################################
## X.m:  n (sample size) by p (number of predictors) data matrix; 
## Y.m:  n (sample size) by q (number of responses) data matrix; 
## C.m:  p (number of predictors) by q (number of responses) data matrix; 
## lamL1: numeric scale; specifying l1 penalty parameter. 
## lamL2: numeric scale; specifying l2 penalty parameter. 
## phi0: NULL or numeric matrix (p by q)

  n=nrow(X.m)
  p=ncol(X.m)
  q=ncol(Y.m)

  X.v=as.vector(t(X.m))
  Y.v=as.vector(t(Y.m))

  if(is.null(C.m))
  {
    C.v=rep(1, p*q)
  } else
  {
    C.v=as.vector(t(C.m))
  }

  iter.count=0
  RSS=0
  Edebug=rep(0, n*q)

  Phi.m=matrix(0, nrow=p, ncol=q)
  Phi.v=as.vector(Phi.m)
  Phi.vdebug=Phi.v
  ###### penalty parameters

  lambda1=lamL2 #### group lasso penalty
  lambda2=lamL1 #### lasso penatly
 
  #dyn.load("remMap.so")

  ###### begin estimation
  if(is.null(phi0))
  {
   junk=.C("MultiRegGroupLasso",
          as.integer(n),
          as.integer(p),
          as.integer(q),
          as.double(X.v),
          as.double(Y.v),
          as.integer(C.v),
          as.double(lambda1[1]),
          as.double(lambda2[1]),
          phi.out=as.double(Phi.v),
          phi.debug=as.double(Phi.vdebug),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)
          )  
  } else {
       phi.ini.v=as.vector(t(phi0))  
       junk=.C("MultiRegGroupLassoIni",
          as.integer(n),
          as.integer(p),
          as.integer(q),
          as.double(X.v),
          as.double(Y.v),
          as.integer(C.v),
          as.double(lambda1),
          as.double(lambda2),
          as.double(phi.ini.v),
          phi.out=as.double(Phi.v),
          phi.debug=as.double(Phi.vdebug),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)
          )  
  }


  phi.result=matrix(junk$phi.out, nrow=p, byrow=T)
  E=matrix(junk$Edebug, nrow=n, ncol=q, byrow=T)
  rss.v=apply(E^2,2,sum)
  iter.count=junk$n.iter

  return(list(phi=phi.result, rss.v=rss.v))
}

############################################################################
############################################################################

remMap.df<-function(X.m, Y.m, lamL1.v, lamL2.v, C.m=NULL)
{
##########################################################
## X.m:  n (sample size) by p (number of predictors) data matrix; 
## Y.m:  n (sample size) by q (number of responses) data matrix; 
## C.m:  p (number of predictors) by q (number of responses) data matrix; 
## lamL1.v: numeric scale/vector; specifying l1 penalty parameters. 
## lamL2.v: numeric scale/vector; specifying l2 penalty parameters. 

  n=nrow(X.m)
  p=ncol(X.m)
  q=ncol(Y.m)

  lam1.v=lamL2.v  #### group lasso penalty
  lam2.v=lamL1.v  #### 

  n1=length(lam1.v)
  n2=length(lam2.v)
  degree=matrix(0, nrow=n1, ncol=n2)
 
  X.v=as.vector(t(X.m))
  Y.v=as.vector(t(Y.m))
  degree.v=as.vector(t(degree)) 

  if(is.null(C.m))
  {
    C.v=rep(1, p*q)
  } else
  {
    C.v=as.vector(t(C.m))
  }

  debug=0   
  #dyn.load("remMap.so")

  junk=.C("MultiRegGroupLassoDegree",
          as.integer(n),
          as.integer(p),
          as.integer(q),
          as.double(X.v),
          as.double(Y.v),
          as.integer(C.v),
          as.integer(n1),
          as.double(lam1.v),
          as.integer(n2),
          as.double(lam2.v),
          degree=as.double(degree),
          debug=as.double(debug)
          )  

  result=t(matrix(junk$degree, nrow=n1, ncol=n2, byrow=T))
  return(result)
}


##############################################################################
############################## performing model selection using BIC
##############################################################################

remMap.BIC<-function(X.m, Y.m,lamL1.v, lamL2.v, C.m=NULL)
{
##########################################################
## X.m:  n (sample size) by p (number of predictors) data matrix; 
## Y.m:  n (sample size) by q (number of responses) data matrix; 
## C.m:  p (number of predictors) by q (number of responses) data matrix; 
## lamL1.v: numeric scale/vector; specifying l1 penalty parameters. 
## lamL2.v: numeric scale/vector; specifying l2 penalty parameters. 


  n=nrow(X.m)
  p=ncol(X.m)
  q=ncol(Y.m)

 ###### sort penalty parameter from large to small
 k1=length(lamL1.v) 
 lambda1.v=sort(lamL1.v)[k1:1]

 k2=length(lamL2.v)
 lambda2.v=sort(lamL2.v)[k2:1]
  
 ####### calculate degree freedom  
 df.m=remMap.df(X.m, Y.m, lamL1.v=lambda1.v, lamL2.v=lambda2.v, C.m=C.m) 
 ####### variables to record the result
 bic<-matrix(-1,k1,k2)
 phi=NULL

 ####### fit the model
 for(i in 1:k1)
  {
     phi.old=NULL 
     if(i>1) phi.old=phi.last
       for(j in 1:k2)        
         {          
          #print(paste(i,j))               
          cur.lam1=lambda1.v[i]          
          cur.lam2=lambda2.v[j]
          temp=remMap(X.m, Y.m, lamL1=cur.lam1, lamL2=cur.lam2, phi0=phi.old, C.m=C.m)          
    
          rss<-temp$rss.v             
          phi[[(i-1)*k2+j]]<-list(lam1=cur.lam1, lam2=cur.lam2, phi=temp$phi) 
          phi.old=temp$phi

          bic[i,j]<-sum(log(rss))*n+df.m[i,j]*log(n)            
          if(j==1) phi.last=temp$phi  
         }
  }
  rownames(bic)=paste("lamL1=", round(lambda1.v,3))
  colnames(bic)=paste("lamL2=", round(lambda2.v,3))
  
  result=list(BIC=bic, phi=phi) 
  return(result)
}



##############################################################################
############################## performing model selection using CV
##############################################################################

remMap.CV<-function(X, Y,lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)
{
##########################################################
## X.m:  n (sample size) by p (number of predictors) data matrix; 
## Y.m:  n (sample size) by q (number of responses) data matrix; 
## C.m:  p (number of predictors) by q (number of responses) data matrix; 
## lamL1.v: numeric scale/vector; specifying l1 penalty parameters. 
## lamL2.v: numeric scale/vector; specifying l2 penalty parameters. 


  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)

 ###### sort penalty parameter from large to small
 k1=length(lamL1.v) 
 lambda1.v=sort(lamL1.v)[k1:1]

 k2=length(lamL2.v)
 lambda2.v=sort(lamL2.v)[k2:1]
  
 l1.index<-as.vector(matrix(lambda1.v, nrow=k1, ncol=k1, byrow=FALSE))
 l2.index<-as.vector(matrix(lambda2.v, nrow=k2, ncol=k2, byrow=TRUE))
 l.index<-rbind(l1.index,l2.index)
 
 ############################### set seed and generate cross validation seperations
 set.seed(seed)

 index.cv<-NULL
 ran.order=sample(1:n, n, replace=F)
 f.n=floor(n/fold)
 for (f in 1:(fold-1))
    {  
     index.cv[[f]]<-ran.order[(f-1)*f.n+1:f.n] 
    } 
 index.cv[[fold]]<-ran.order[((fold-1)*f.n+1):n]
 
 ################################ begin cross validation
 phi.cv<-NULL 
 rss.cv.cv<-NULL 
 ols.cv.cv<-NULL   
 for(f in 1:fold)
  {  #print(paste("fold=",f))   
     index.cur.cv<-index.cv[[f]]   
     X.m<-X[-(index.cur.cv),]   
     Y.m<-Y[-(index.cur.cv),]   
     X.t<-X[index.cur.cv,]   
     Y.t<-Y[index.cur.cv,]   

     ols.cv.c<-array(0,dim=c(k1,k2,q))
     rss.cv.c<-array(0,dim=c(k1,k2,q))   
     phi.c<-NULL   
     phi.old<-NULL   
     for(i in 1:k1)
     {      
       phi.old<-NULL
       if(i>1) phi.old=phi.last      
       for(j in 1:k2)       
       {          
          #print(paste(i,j))               
          cur.lam1=lambda1.v[i]          
          cur.lam2=lambda2.v[j]
          temp=remMap(X.m, Y.m, lamL1=cur.lam1, lamL2=cur.lam2, phi0=phi.old, C.m=C.m)          

          rss<-temp$rss.v             
          phi.c[[(i-1)*k1+j]]<-list(lam1=cur.lam1, lam2=cur.lam2, phi=temp$phi) 
          phi.old=temp$phi
          try.adj=abs(phi.old)>1e-6                
          ols.cv.c[i,j,]<-apply(as.matrix(1:q),1,OLS.CV, Y.m=Y.m, X.m=X.m, Y.t=Y.t, X.t=X.t, Coef.m=try.adj)                           
          rss.cv.c[i,j,]<-RSS.CV(X.t,Y.t,phi.old)
          
          if(j==1) phi.last=temp$phi  
        } ### end of j iter     
      }### end of i iter
     ##save     
     phi.cv[[f]]<-phi.c        
     rss.cv.cv[[f]]<-rss.cv.c    
     ols.cv.cv[[f]]<-ols.cv.c    
   }###end fold loop
   ################################ analyze cross validation result
   ols.cv.sum<-array(0,dim=c(k1,k2,q))    
   for (f in (1:fold)){         
        ols.cv.sum<-ols.cv.sum+ols.cv.cv[[f]]        }
   ols.cv<-apply(ols.cv.sum,c(1,2),sum)   
   rownames(ols.cv)=paste("lamL1=", round(lambda1.v,3))
   colnames(ols.cv)=paste("lamL2=", round(lambda2.v,3))
 
             
   rss.cv.sum<-array(0,dim=c(k1,k2,q))    
   for (f in (1:fold)){         
        rss.cv.sum<-rss.cv.sum+rss.cv.cv[[f]]        }
   rss.cv<-apply(rss.cv.sum,c(1,2),sum)
   rownames(rss.cv)=paste("lamL1=", round(lambda1.v,3))
   colnames(rss.cv)=paste("lamL2=", round(lambda2.v,3)) 
    
    ################################## return CV results
   result=list(ols.cv=ols.cv, rss.cv=rss.cv,phi.cv=phi.cv, l.index=l.index)
   return(result)
}
  


##########################################################
##################  private functions
##########################################################

OLS.CV<-function(Y.m,X.m,Y.t,X.t,Coef.m,i)
{
  ##print(i)
  index.cur<-Coef.m[,i]
  Y.c<-matrix(Y.m[,i],ncol=1)
  Y.t.c<-matrix(Y.t[,i],ncol=1)
  n<-nrow(Y.m)
  
  if(sum(index.cur)==0){
      rss.t<-sum(Y.t.c^2)
      return(rss.t)
     } 
      
   X.c<-matrix(X.m[,index.cur],ncol=sum(index.cur))
   X.t.c<-matrix(X.t[,index.cur],ncol=sum(index.cur))
    
   if(sum(index.cur)>n){
      X.c<-X.c[,1:n]
      X.t.c<-X.t.c[,1:n] 
     }
          
    beta<-solve(t(X.c)%*%X.c)%*%t(X.c)%*%Y.c
    rss.t<-sum((Y.t.c-X.t.c%*%beta)^2)
   
    return(rss.t)
}


RSS.CV<-function(X.t,Y.t, phi.est)
{    
    Y.est<-X.t%*%phi.est       
    resi<-Y.t-Y.est
    rss<-apply(resi^2,2,sum)
    return(rss)
}

