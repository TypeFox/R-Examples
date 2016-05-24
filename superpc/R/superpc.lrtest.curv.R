superpc.lrtest.curv <- function (object, data, newdata, n.components=1, threshold=NULL,  n.threshold=20)
{

  this.call <- match.call()

                                        # compute lrtest statistics based on fit "object", training data "data",
                                        # and test  data "newdata",
                                        # over a set of threshold values

  type=object$type

  if(!is.null(threshold)) {n.threshold=length(threshold)}
  if(is.null(threshold)){
    second.biggest<- -sort(-abs(object$feature.scores))[2]
    threshold<- seq(0,second.biggest, length=n.threshold)
  }

  n.pc <- n.components
  lrtest<-rep(NA, n.threshold)
  num.features<-rep(NA, n.threshold)


  cat("",fill=TRUE)
  for(ii in 1:n.threshold){
    cat(ii)

    object.temp<- superpc.predict(object, data, newdata,threshold=threshold[ii], n.components=n.pc)
    
    num.features[ii]<- sum(object.temp$which.features)
    
    v.pred<-object.temp$v.pred
    
    
    if(type=="survival"){
      require(survival)
      junk<- coxph(Surv(newdata$y, newdata$censoring.status) ~v.pred)$loglik
      lrtest[ii]<-2*(junk[2]-junk[1])
    }
    else{junk<- summary(lm(newdata$y~v.pred))
         lrtest[ii]<-junk$fstat[1]
       }

  }
  cat("",fill=TRUE)

  return(list(lrtest=lrtest,threshold=threshold,num.features=num.features, type=type, call=this.call))
}
