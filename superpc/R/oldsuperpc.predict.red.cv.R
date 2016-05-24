superpc.predict.red.cv <- function(fitred, fitcv, data, threshold, num.reduced.models=30, sign.wt="both"){

 # try reduced predictor on cv folds, via prevalidation

                           
  this.call=match.call()

  type=fitred$type

  n.components=fitred$n.components


  n.fold<-length(fitcv$folds)

  shrinkages<- fitred$shrinkages
  num.reduced.models<-length(shrinkages)
  cur.vall<- array(NA,c(num.reduced.models,ncol(data$x),n.components))

  for(j in 1:n.fold){
    cat(j,fill=TRUE)
    fit.temp<-list(feature.scores=fitcv$featurescores.fold[,j], type=type)
    ii<-fitcv$folds[[j]]
    
    data1<-list(x=data$x[,-ii],y=data$y[-ii],censoring.status=data$censoring.status[-ii])
    data2<-list(x=data$x[,ii],y=data$y[ii],censoring.status=data$censoring.status[ii])
    junk<- superpc.predict.red(fit.temp, data1,data2, threshold, num.reduced.models=num.reduced.models, n.components=n.components,compute.lrtest=FALSE, sign.wt=sign.wt)
    cur.vall[,ii,]<-junk$v.test
  }

  lrtest.reduced<-rep(NA,num.reduced.models)

  for(i in 1:num.reduced.models){
    if(type=="survival"){
      require(survival)
      junk<- coxph(Surv(data$y, data$censoring.status) ~cur.vall[i,,])$loglik
      lrtest.reduced[i]=2*(junk[2]-junk[1])
    }
    else{
      junk<- summary(lm(data$y~cur.vall[i,,]))
      if(!is.null(junk$fstat)){lrtest.reduced[i]<-junk$fstat[1]}
    }

  }

return(list(shrinkages=shrinkages, lrtest.reduced=lrtest.reduced,  n.components=n.components, num.features=fitred$num.features, v.preval.red=cur.vall, sign.wt=sign.wt, type=type,call=this.call))
}
