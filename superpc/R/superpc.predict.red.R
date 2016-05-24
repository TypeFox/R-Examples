superpc.predict.red <- function(fit, data, data.test, threshold, n.components=3, n.shrinkage=20, shrinkages=NULL, compute.lrtest=TRUE, sign.wt="both", prediction.type=c("continuous","discrete"), n.class=2){
# try reduced predictor on test set

  
  soft.thresh<- function(x,tt){ sign(x)*(abs(x)-tt)*(abs(x)>tt)}


  this.call<- match.call()
  
  prediction.type=match.arg(prediction.type)

  type=fit$type

  if(!is.null(shrinkages)){ n.shrinkage=length(shrinkages)}

  lrtest.reduced<- matrix(NA,ncol=n.components,nrow=n.shrinkage)
  
  cur.vall<- array(NA,c(n.shrinkage,ncol(data$x),n.components))
  cur.vall.test<- array(NA, c(n.shrinkage,ncol(data.test$x),n.components))
  cur.vall.test.groups<- array(NA, c(n.shrinkage,ncol(data.test$x),n.components))
  cur.vall.1df<- matrix(NA, nrow=n.shrinkage,ncol=ncol(data$x))
  cur.vall.test.1df<- matrix(NA, nrow=n.shrinkage,ncol=ncol(data.test$x))

  corr.with.full<-matrix(NA,nrow=n.shrinkage, ncol=n.components)
  which.features <- abs(fit$feature.scores) > threshold
  x.sml <- data$x[which.features, ]
  x.svd <- mysvd(x.sml, n.components=n.components)
scal=apply(scale(abs(x.svd$u),center=F,scale=x.svd$d),2,sum)

  cur.v <- scale(t(data$x[which.features, ]) %*%x.svd$u, center=FALSE,scale=scal*x.svd$d)

  # flip the sign of the latent factors, if a coef is neg

  junk<-superpc.fit.to.outcome(fit,data, cur.v,  print=FALSE)
  if(fit$type=="survival"){coef=junk$coeftable[,1]}
  if(fit$type=="regression"){coef=junk$coeftable[-1,1]}

  cur.v<-scale(cur.v, center=FALSE,scale= sign(coef))
  ##

  train.means=apply(data$x,1,mean)
  xcen=t(scale(t(data$x), center=train.means,scale=F))
  xtest.cen=t(scale(t(data.test$x), center=train.means,scale=F))
  
#  import<-cor(t(data$x), cur.v)

  sc= xcen%*%cur.v
sc[abs(fit$feature.scores)< threshold]=0

import=sc



                                        # don't shrink all of the way to zero 
  maxshrink=max(abs(sc))

  if(sign.wt=="positive"){ maxshrink=max(abs(sc[sc>0]))}
  
  if(sign.wt=="negative"){ maxshrink=max(abs(sc[sc<0]))}

  
  if(is.null(shrinkages)){
    shrinkages<- seq(0,maxshrink,length=n.shrinkage+1)
    shrinkages= shrinkages[-(n.shrinkage+1)]
  }

  num.features<-matrix(NA,nrow=n.shrinkage, ncol=n.components)
  
  feature.list<-vector("list", n.shrinkage)


  for(i in 1:n.shrinkage){
    cat(i)
    sc2<- soft.thresh(sc,shrinkages[i])
    if(sign.wt=="positive"){sc2[sc2<0]<-0}
    if(sign.wt=="negative"){sc2[sc2>0]<-0}
    nonzero<-sc2!=0
    

    num.features[i,]<- apply(nonzero,2,sum)
    
    junk=vector("list",n.components)
    for(ii in 1:n.components){
      junk[[ii]]<- (1:nrow(xcen))[nonzero[,ii]]
    }
    feature.list[[i]]=junk

    for(ii in 1:n.components){
      cur.vall[i,,ii]<-apply(t(scale(t(xcen[nonzero[,ii],,drop=FALSE]), center=FALSE,scale=1/sc2[nonzero[,ii],ii])),2,mean)
      cur.vall.test[i,,ii]<-apply(t(scale(t(xtest.cen[nonzero[,ii],,drop=FALSE]), center=FALSE,scale=1/sc2[nonzero[,ii],ii])),2,mean)

                                        # find quantile break points on training data; apply to test data
      if(prediction.type=="discrete") {
        if(sum(is.na(cur.vall.test[i,,ii]))==0){
          br=quantile(cur.vall[i,,ii], (0:n.class)/n.class)
          cur.vall.test.groups[i,,ii]<-cut(cur.vall.test[i,,ii],breaks=br,labels=FALSE)
          o=is.na(cur.vall.test.groups[i,,ii])
          cur.vall.test.groups[i,o,ii]<-1
        }
        else{ cur.vall.test.groups[i,,ii]<-1}
      }
    }   
    
    cur.vall.1df[i,]=apply( scale(cur.vall[i,,],center=F,scale=1/coef),1,mean)
    cur.vall.test.1df[i,]=apply( scale(cur.vall.test[i,,],center=F,scale=1/coef),1,mean)
    

    if(prediction.type=="discrete") {
      if(sum(is.na(cur.vall.1df[i,]))==0){
        br=quantile(cur.vall.1df[i,], (0:n.class)/n.class)
        cur.vall.test.1df[i,]<- cut( cur.vall.test.1df[i,], breaks=br,labels=FALSE)
        o=is.na(cur.vall.test.1df[i,])
        cur.vall.test.1df[i,o]<-1
      }
  else{cur.vall.test.1df[i,]<-1}
     
   
       }

  }


 if(prediction.type=="continuous"){
   dimnames(cur.vall.test)=list(NULL,NULL,rep("continuous", dim(cur.vall.test)[3]))
 }

 if(prediction.type=="discrete"){
   cur.vall.test<- cur.vall.test.groups
   dimnames(cur.vall.test)=list(NULL,NULL, rep("factor", dim(cur.vall.test)[3]))
 }

  cat("",fill=TRUE)


  if(compute.lrtest){ 
for(ii in 1:n.components){
    for(i in 1:n.shrinkage){
      if(type=="survival"){
        require(survival)
        
# with too much shrinkage,
# all predictors may be shrunk to zero and Cox model bombs. I check for this first

         if(sum(is.na(cur.vall.test[i, , 1:ii]))==0){
junk=superpc.fit.to.outcome(fit, data.test, cur.vall.test[i,,1:ii], print=FALSE)$results$loglik
        lrtest.reduced[i,ii]=2*(junk[2]-junk[1])
        }
      }
      else{
         if(sum(is.na(cur.vall.test[i, , 1:ii]))==0){
junk=superpc.fit.to.outcome(fit, data.test, cur.vall.test[i,,1:ii], print=FALSE)
        if(!is.null(junk$fstat)){lrtest.reduced[i,ii]<-junk$results$fstat[1]}
      }
      }
    }
  }
}
  

  return(list(shrinkages=shrinkages, lrtest.reduced=lrtest.reduced, 
num.features=num.features, feature.list=feature.list, import=import,  v.test=cur.vall.test, coef=coef, v.test.1df= cur.vall.test.1df,
              n.components=n.components, sign.wt=sign.wt, type=type,call=this.call))
  
}



