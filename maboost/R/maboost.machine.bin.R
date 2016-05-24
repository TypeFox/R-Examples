"maboost.machine.bin" <-
  function(x,y,test.x,test.y,iter=50,smoothfactor=1,nu=0.1
           ,bag.frac=0.5,random.feature=TRUE,random.cost=TRUE
           ,lossObj,oldObj=NULL,na.action=na.action,...){
    kapstat<-function(tab=diag(2) ){
      if(dim(tab)[1]==1){
        return(0)
      }
      if(dim(tab)[1]==dim(tab)[2]){
        rs<-apply(tab,2,sum)
        cs<-apply(tab,1,sum)
        N<-sum(rs)
        E<-sum(rs*cs)/N^2
        O<-sum(diag(tab))/N
        return( (O-E)/(1-E) )
      }else{
        return(0.5)
      }
    }
    tmp<-function(i){
      a1<-sample(which(y==1),1)
      a2<-sample(which(y==-1),1)
      ind<-c(a1,a2)
      return(c(sample(setdiff(1:n,ind),n-val-2,FALSE),ind))
    }
    
    n=dim(x)[1]
    num_feat=dim(x)[2]
    num_zero=rep(0,iter);
    fit=list()
    y<-as.numeric(y)
    dat<-data.frame(y=y,x)
    w_n=rep(1/n,n) 
    w=rep(1/n,n)
    oobm.mat<-matrix(0,nrow=n,ncol=2)
    fits=rep(0,n)
    ind_feature=c(1,(1:num_feat)+1 )
    feat_cost=rep(1,num_feat);
    
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
      fits.test<-rep(0,test.n)
    }
    if(!is.null(oldObj)){
      fit=oldObj$model$trees
      w=oldObj$model$lw
      w_n=w/sum(w)
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
    if(val<n){
      a<-sapply(1:iter,tmp)
    } 
    if(is.vector(a)){
      a<-t(as.matrix(a))
    }
    start<-start +1
    wfun=lossObj$wfun
    coefs=lossObj$coefs
    method=lossObj$method
    predict.type=lossObj$predict.type
    dcal=lossObj$dcal;
    Ctree=lossObj$Ctree;
    verbose=lossObj$verbose
    if(Ctree) dat$y=as.factor(dat$y);

    f1<-f2<-0
    for (m in start:iter){
      xval=1:n
      if(!is.null(a)){
        xval=a[,m]
      }
      if(random.feature)  {ind_feature=c(1,sample(1:num_feat,floor(num_feat^.5),replace=F)+1 );feat_cost = rep(1,floor(num_feat^.5))}
      if(random.cost)     feat_cost=3*runif(length(ind_feature)-1)+1;
      
      
      if(Ctree){
        temp = w_n[xval];
        
        temp[temp==0]=temp[temp==0]+10^-12;
       
        ind = sample( 1:length(xval) , length(xval), replace=T, prob=temp);
        
        
        #fit[[m]] =C5.0( y~., data=dat[ind,ind_feature],trials=1,control=C5.0Control(minCases=round(M/(K+20)),CF=.1) );
        fit[[m]] =C5.0( y~., data=dat[ind,ind_feature],trials=1
                        ,na.action=na.pass,control=C5.0Control(...) );
        
      }
      else{
        
        fit[[m]] =rpart(y~.,data=dat[xval,ind_feature],weights=w_n[xval]
                        ,method=method,x=FALSE,y=TRUE
                        ,na.action=na.action,control=rpart.control(cost=feat_cost,...));
      }
      
      f<-predict.type(fit[[m]],dat)
      d=dcal(f,y)
      alpha[m] = nu * coefs(d,w,m)
      fits<-fits+alpha[m]*f
      
      w=wfun(w,d,alpha[m])

      w_n=w/(sum(w)+10^-20)
      num_zero[m]=length(which(w_n==0));
   
      tab<-table(sign(fits),y)
      train.err[m]<-1-sum(diag(tab))/n
      train.kap[m]<-1-kapstat(tab)
      if(verbose)  print(c('err',round(train.err[m],8),'eta',round(alpha[m],8),'sum_w',round(sum(w),8),'num_zero',num_zero[m],'max_w*M',round(max(w)*n,4) ))
      
      indx<- setdiff(1:n,xval)
      btmp<-as.numeric(as.factor(sign(fits)[indx]))
      if(length(btmp)==1){
        oobm.mat[indx,btmp]<-oobm.mat[indx,btmp]+1
      }else{
        oobm.mat[indx,][btmp==1,1]<- oobm.mat[indx,][btmp==1,1]+1
        oobm.mat[indx,][btmp==2,2]<- oobm.mat[indx,][btmp==2,2]+1
        denom<-apply(oobm.mat,1,sum)
        vals<-denom>0
        if(sum(vals)==1){
          ytr<-c(-1,1)[which.max(oobm.mat[vals,])]
        }else{
          ytr<-c(-1,1)[apply(oobm.mat[vals,],1,which.max)]
        }
        oobm.err[m]<-sum(ytr!=y[vals])/length(vals)
      }
      if(is.null(test.y)){
        next
      }
      fit1<-predict.type(fit[[m]],test.x)
      fits.test<-fits.test +alpha[m]*fit1
    
      tab<-table(sign(fits.test),test.y)
      test.err[m]<- 1-sum(diag(tab))/test.n
      test.kap[m]<-1-kapstat(tab)

    }
   
    a1=(fits==0)
    if(sum(a1)>0)
      fits[a1]<-sample(c(-1,1),sum(a1),TRUE,c(.5,.5))
    errs<-cbind(train.err,train.kap)
    ans<-list()
    ans[[1]]=fits
    if(!is.null(test.y)){
      errs<-cbind(errs,test.err,test.kap)
      ans[[2]]=fits.test
    }
    obj=list(trees=fit,alpha=alpha,F=ans,errs=errs,oob.str=list(oobm.err=oobm.err,oobm.mat=oobm.mat),lw=w,lossObj=lossObj,num_zero=num_zero)
    return(obj)
  }