"predict.glmmNPML" <-
predict.glmmNPML<-function(object, newdata,  type="link", ...){

  if(missing(newdata)){ #Emp. Bayes Pred. (Aitkin, 96)
      if (type=="link") {
          return(round(object$ebp,digits=4))
      } else {
          #rebp <- object$fitted.values
          #rebp<- switch(object$family$link, "log"= exp(ebp),
          #                                  "logit"=exp(ebp)/(1+exp(ebp)),
          #                                  "identity"=ebp,
          #                                  "inverse"=1/ebp,
          #                                  "probit"= pnorm(ebp,0,1))
          return(round(object$fitted.values, digits=4))
      }

  } else  {
      m<-length(object$mass.points)
      k<-length(object$masses)
      Terms<-  delete.response(terms(object$formula))
      if (k==1){  # GLM
              X<-model.matrix(Terms, model.frame(Terms,newdata))
              dimnames(X)[[1]]<-dimnames(model.frame(Terms,newdata))[[1]]
              pred<-  as.vector(X%*%matrix(object$coef[1:dim(X)[2]]))
      }  else if (k==m){ #Overdispersion model
              X<-model.matrix(Terms, model.frame(Terms,newdata))[,-1,drop=FALSE]
              dimnames(X)[[1]]<-dimnames(model.frame(Terms,newdata))[[1]]
              if (dim(X)[2]!=0){
                    pred<- as.vector(X%*%matrix(object$coef[1:dim(X)[2]])+ sum(object$masses*object$mass.points))
              } else {
                    pred<- rep(0, dim(X)[1])+ sum(object$masses*object$mass.points)
              }
      } else { #Random coefficient models
              X<-model.matrix(Terms, model.frame(Terms,newdata))[,-1,drop=FALSE]
              dimnames(X)[[1]]<-dimnames(model.frame(Terms,newdata))[[1]]
              object$mass.points<- ifelse(is.na(object$mass.points),0,object$mass.points)
              #r<-  names(newdata) %in% gsub('~','',object$random)[2]
              r<-  names(newdata) %in% object$Misc$mform   #28/02/06
              if(is.factor(newdata[,r])){newdata[,r]<-as.numeric(newdata[,r])-1}  #28/02/06
              
              if (dim(X)[2]!=0){
                  pred<- as.vector(X%*%matrix(object$coef[1:dim(X)[2]])+ sum(object$masses*object$mass.points[1:k])+ newdata[,r]*sum(object$masses[1:(k-1)]*object$mass.points[(k+1):m])) #24-07-06
                  } else {
                  pred<- rep(0, dim(X)[1])+ sum(object$masses*object$mass.points[1:k])+ newdata[,r]*sum(object$masses[1:k]*object$mass.points[(k+1):m])
              }
      }
      if (type=="link"){
             rpred<-as.vector(pred)
      } else {
              rpred <- object$family$linkinv(pred)
      }
      names(rpred)<-dimnames(X)[[1]]
      return(round(rpred,digits=4))
  }

}

