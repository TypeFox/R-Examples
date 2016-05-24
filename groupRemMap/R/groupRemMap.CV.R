group.remmap.cv <-
function(X,Y,G,lam1.v, lam2.v,gamma=0.5, C.m=NULL,fold=10, seed=1)
{

  ## X.m:  n (sample size) by p (number of predictors) data matrix; 
  ## Y.m:  n (sample size) by q (number of responses) data matrix; 
  ## C.m:  p (number of predictors) by q (number of responses) data matrix; 
  ## lam1.v: numeric scale/vector; specifying l1 penalty parameters. 
  ## lam2.v: numeric scale/vector; specifying bridge penalty parameters. 

  
  X=as.matrix(X);
  Y=as.matrix(Y);
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)

 ###### sort penalty parameter from large to small
 k1=length(lam1.v) 
 lambda1.v=sort(lam1.v,decreasing=T)

 k2=length(lam2.v)
 lambda2.v=sort(lam2.v,decreasing=T)
  
 l1.index<-rep(lambda1.v,k2)
 l2.index<-as.vector(sapply(lambda2.v,function(v){rep(v,k1)}))
 l.index<-rbind(l1.index,l2.index)
 l.col=ncol(l.index)
 
 ############################### set seed and generate cross validation seperations
 set.seed(seed)

 cv.id=sample(rep(1:fold,length=n));
 
 ################################ begin cross validation
 #phi.cv<-NULL 
 rss.cv.cv<-NULL 
 ols.cv.cv<-NULL   

 phi.old=NULL;

 for(f in 1:fold)
 { 
     #print(paste("fold=",f))   
     index.cur.cv<-which(cv.id==f)   

     X.m<-as.matrix(X[-(index.cur.cv),])  
     Y.m<-as.matrix(Y[-(index.cur.cv),])  
     X.t<-as.matrix(X[index.cur.cv,])   
     Y.t<-as.matrix(Y[index.cur.cv,])   

     ols.cv.c<-matrix(0,nrow=l.col,ncol=q);
     rss.cv.c<-matrix(0,nrow=l.col,ncol=q);
     phi.c<-NULL;


     for(i in 1: l.col)
     {

          cur.lam1=l.index[1,i];          
          cur.lam2=l.index[2,i];
          #print(sprintf("penalty settting=%d, lam1=%.2f, lam2=%.2f",i, cur.lam1, cur.lam2))                     
          temp=group.remmap(X=X.m, Y=Y.m,G=G,lam1=cur.lam1,lam2=cur.lam2,gamma=gamma, phi0=phi.old, C.m=C.m);  
          #print(temp$iter.count) ;     

          #phi.c[[(i-1)*k1+j]]<-list(lam1=cur.lam1, tau=cur.lam2, phi=temp$phi);
          phi.cur=temp$phi
          try.adj=abs(phi.cur)>0;
          ols.cv.c[i,]<-sapply(1:q,function(v){OLS.CV(Y.m=Y.m, X.m=X.m, Y.t=Y.t, X.t=X.t, Coef.m=try.adj,v)});                           
          rss.cv.c[i,]<-RSS.CV(X.t,Y.t,phi.cur);
          #phi.old=phi.cur;

    }   ## end of i iter;

     #phi.cv[[f]]<-phi.c        
     rss.cv.cv[[f]]<-rss.cv.c    
     ols.cv.cv[[f]]<-ols.cv.c       


  } ## end of fold iter



   ################################ analyze cross validation result
   ols.cv.sum<-matrix(0,nrow=l.col,ncol=q);   
   for (f in (1:fold))
   {         
        ols.cv.sum<-ols.cv.sum+ols.cv.cv[[f]];  
   }
   ols.cv<-apply(ols.cv.sum,1,sum);   
   names(ols.cv)=sapply(1:l.col,function(v){paste("lam1=",l.index[1,v],",","lam2=",l.index[2,v],sep="")});

 
             
   rss.cv.sum<-matrix(0,nrow=l.col,ncol=q);     
   for (f in (1:fold))
   {         
        rss.cv.sum<-rss.cv.sum+rss.cv.cv[[f]]       
   }

   rss.cv<-apply(rss.cv.sum,1,sum)
   names(rss.cv)=sapply(1:l.col,function(v){paste("lam1=",l.index[1,v],",","lam2=",l.index[2,v],sep="")});

    
    ################################## return CV results
   #result=list(ols.cv=ols.cv, rss.cv=rss.cv,phi.cv=phi.cv, l.index=l.index)
    result=list(ols.cv=ols.cv, rss.cv=rss.cv, l.index=l.index);
   return(result)
}


##########################################################
##################  private functions
##########################################################

OLS.CV<-function(Y.m,X.m,Y.t,X.t,Coef.m,v)
{
  ##print(i)
  index.cur<-Coef.m[,v]
  Y.c<-matrix(Y.m[,v],ncol=1)
  Y.t.c<-matrix(Y.t[,v],ncol=1)
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

