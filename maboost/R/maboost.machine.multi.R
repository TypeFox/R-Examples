"maboost.machine.multi" <-
  function(x,y,test.x,test.y,iter=50,smoothfactor=1,nu=0.1,bag.frac=0.5,random.feature=TRUE,random.cost=TRUE
           ,lossObj,oldObj=NULL,na.action=na.action,...){
    kapstat<-function(tab=diag(2) ){
      if(dim(tab)[1]==1){
        return(0)
      }
      
        rs<-apply(tab,2,sum)
        cs<-apply(tab,1,sum)
        N<-sum(rs)
        c_cs= names(cs)
        E<-sum(rs[c_cs]*cs)/N^2
        O<-sum(tab[cbind(c_cs,c_cs)])/N
        return( (O-E)/(1-E) )
      
    }
    tmp<-function(i){
      ai=rep(0,K)
      for (l in 1 : K){
        
        ai[l]<-sample(which(y==l),1)
      }
      
      return(c(sample(setdiff(1:n,ai),n-val-K,FALSE),ai))
    }
    
    n=dim(x)[1]
    num_feat=dim(x)[2]
    num_zero=rep(0,iter)
    S_W=1;
    M=n;
    fit=list()
    y<-as.factor(y)
    dat<-data.frame(y=y,x)
    K=lossObj$K;
    Ctree=lossObj$Ctree;
    verbose=lossObj$verbose;
    W=matrix(1/((K-1)*M),M,K-1);
    w=rep(1/M,M)
    oobm.mat<-matrix(0,nrow=n,ncol=K)
    fits=matrix(0,M,K)
    
    atmp=alpha=vector(length=iter)
    oobm.err<-rep(0,iter)
    train.err<-rep(0,iter)
    train.kap<-rep(0,iter)
    start=0
    if(!is.null(test.y)){
      fit=oldObj$model$trees
      test.err<-rep(0,iter)
      test.kap<-rep(0,iter)
      test.n<-dim(test.x)[1]
      fits.test<-matrix(0,test.n,K)
    }
    if(!is.null(oldObj)){
      fit=oldObj$model$trees
      
        W=oldObj$model$lw;
        w=apply(W,1,sum);
      
      
      oobm.mat=oldObj$model$oob.str$oobm.mat
      fits=oldObj$model$F[[1]]
      start=oldObj$iter
      alpha[1:start]<-oldObj$model$alpha
      train.err[1:(start)]<-oldObj$model$err[,1]
      train.kap[1:(start)]<-oldObj$model$err[,2]
      oobm.err[1:(start)]<-oldObj$model$oob.str$oobm.err
      nu=oldObj$nu
      bag.frac=oldObj$bag.frac
      if(!is.null(test.y)){
        test.err[1:start]<-oldObj$model$err[,3]
        test.kap[1:start]<-oldObj$model$err[,4]
        fits.test<-oldObj$model$F[[2]]
      }
    }  
    val<-floor(bag.frac*n)
    a<-NULL
    if ((n-val-K)>0){
      if(val<n){
        a<-sapply(1:iter,tmp)
      } 
      if(is.vector(a)){
        a<-t(as.matrix(a))
      }
    }
    start<-start +1
    wfun=lossObj$wfun
    coefs=lossObj$coefs
    method=lossObj$method
    dcal=lossObj$dcal;
    predict.type=lossObj$predict.type
    
    f1<-f2<-0
    ind_feature=c(1,(1:num_feat)+1 )
    feat_cost=rep(1,num_feat);
    for (m in start:iter){
      xval=1:n
      if(!is.null(a)){
        xval=a[,m]
      }
      
      if(random.feature)  {ind_feature=c(1,sample(1:num_feat,floor(num_feat^.5),replace=F)+1 );feat_cost = rep(1,floor(num_feat^.5))}
      if(random.cost)     feat_cost=3*runif(length(ind_feature)-1)+1;
            
      
      if(Ctree){
        temp = w[xval];
       
        temp[temp==0]=temp[temp==0]+10^-12;
        ind = sample( 1:length(xval) , length(xval), replace=T, prob=temp);
        
       
        #fit[[m]] =C5.0( y~., data=dat[ind,ind_feature],trials=1,control=C5.0Control(minCases=round(M/(K+20)),CF=.1) );
        fit[[m]] =C5.0( y~., data=dat[ind,ind_feature],trials=1
                    ,na.action=na.pass,control=C5.0Control(...) );
        
      }
      else{
      
        fit[[m]] =rpart(y~.,data=dat[xval,ind_feature],weights=w[xval]
        ,method=method,x=FALSE,y=TRUE
        ,na.action=na.action,control=rpart.control(cost=feat_cost,...));
      }
      f<-predict.type(fit[[m]],dat)
     
      D=dcal(f,y)
     
      alpha[m]=nu*coefs(D,W,m)
      
      fits[cbind( (1:M),f )] = fits[cbind( (1:M),f )] + alpha[m];
      
      
      W=wfun(W,D,alpha[m])
      S_W=sum(W)
      w=apply(W,1,sum)/S_W; # this wasn't necessary if we could pass a cost matrix to weaklearner
      num_zero[m]=length(which(W==0));
      
      tab<-table(apply(fits,1,which.max),y)
     
      XX=row.names(tab)
     
      train.err[m]<- 1-sum(tab[cbind(XX,XX)])/n
      train.kap[m]<-1-kapstat(tab)
      if(verbose)  print(c('err',round(train.err[m],8),'eta',round(alpha[m],8),'sum_w',round(sum(w),8),'num_zero',num_zero[m],'max_w*M',round(max(w)*n*(K-1),4) ))
      
      indx<- setdiff(1:n,xval)
      
      
      if(length(indx)>0){
        btmp<-apply(fits[indx,],1,which.max)
        oobm.mat[cbind(indx,btmp)]<-  oobm.mat[cbind(indx,btmp)]+1
        
        oobm.err[m]<-length(which(apply(oobm.mat,1,which.max)!=y))/length(y)
      }
      if(is.null(test.y)){
        next
      }
      fit1<-predict.type(fit[[m]],test.x)
      fits.test[cbind( (1:test.n),fit1 )] = fits.test[cbind( (1:test.n),fit1 )] + alpha[m];
      tab<-table(apply(fits.test,1,which.max),test.y)
      XX=row.names(tab)
      test.err[m]<- 1- sum(tab[cbind(XX,XX)])/test.n
      test.kap[m]<-1-kapstat(tab)
      
    }
    
    
    errs<-cbind(train.err,train.kap)
    ans<-list()
    ans[[1]]=fits
    if(!is.null(test.y)){
      errs<-cbind(errs,test.err,test.kap)
      ans[[2]]=fits.test
    }
    obj=list(trees=fit,alpha=alpha,F=ans,errs=errs,oob.str=list(oobm.err=oobm.err,oobm.mat=oobm.mat),lw=W,lossObj=lossObj,num_zero=num_zero)
    return(obj)
  }