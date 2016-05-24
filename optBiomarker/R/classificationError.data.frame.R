# $Id: errorest.R,v 1.25 2005/06/29 08:50:28 hothorn Exp $



## coding style greatfully accommodated from 'errorest' - ipred package
classificationError <- function(formula, data, ...)
{ 
UseMethod("classificationError", data)
}

classificationError.default <- function(formula, data, ...)
  stop(paste("Do not know how to handle objects of class", class(data)))
  

		   
		   
classificationError.data.frame <- function(formula, 
                                data, 
                                method=c("RF","SVM","LDA","KNN"), 
                                errorType = c("cv", "boot", "six32plus"),
				senSpec=TRUE,
				negLevLowest=TRUE,
				na.action=na.omit, 
				control=control.errorest(k=NROW(na.action(data)),nboot=100),
                                ...) 
{
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if(missing(formula)) 
       stop("'formula' missing")
  if(!identical(class(formula),"formula"))
       stop("incorrect formula")
  if (!identical(as.numeric(length(attr(terms(formula[-3], data = data), "term.labels"))),1))
       stop("response must be one dimensional")
  noPredictor <- (length(attr(terms(formula[-2], data = data), "term.labels")) < 1)
  if(noPredictor)
  stop("no predictor variable specified in 'formula'")
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
     m[[1]] <- as.name("model.frame")
     m$method<-NULL
     m$errorType<-NULL
     m$senSpec<-NULL
     m$negLevLowest<-NULL
     m$control<-NULL
     m$... <- NULL
     mf<-eval(m,parent.frame())
     response <- attr(attr(mf, "terms"), "response")
     attr(mf, "terms") <- NULL
     y <- mf[,response]
     data <- mf

## User supplied methods

  allmethod<-c("RF","SVM","LDA","KNN")
  alltype<-c("cv", "boot", "six32plus")
  
  
  if( any(c(length(method)<1, length(errorType)<1)))
    stop("'method' and 'errorType' should be of length >= 1") 

  usrmethod<-match.arg(method, allmethod, several.ok = TRUE)
  usrtype<-match.arg(errorType, alltype, several.ok = TRUE)
  nMethod<-length(usrmethod)
  nType<-length(usrtype)
  
  errorRate<-matrix(,nrow=nType, ncol=nMethod)
  rownames(errorRate)<-usrtype
  colnames(errorRate)<-usrmethod
  
  if(senSpec & ("cv" %in% usrtype)){
   rocData<-matrix(NA, ncol=nMethod,nrow=2)
   colnames(rocData)<-usrmethod
   rownames(rocData)<-c("sens","spec")
   } else rocData<-NULL
  
  for (i in usrtype){
      if(senSpec & identical(i,"cv")) control$predictions<-TRUE  else control$predictions<-FALSE
       for (j in usrmethod){
   est<-errorest(formula=formula, data=data, 
                            model=switch(j, RF=randomForest,SVM=svm,LDA=lda,KNN=ipredknn),
			    predict=switch(j, RF=predict,SVM=predict,LDA=classpredict.lda,KNN=classpredict.knn), 
                            estimator =switch(i, cv="cv",boot="boot",six32plus="632plus"), 
			    na.action=na.action,
                            est.para=control,...)
  errorRate[i,j]<-est$error
  if(senSpec & identical(i,"cv")){
    if(negLevLowest) {
     neg<-min(levels(y)) 
     pos<-max(levels(y)) 
     } else 
     {neg<-max(levels(y))
     pos<-min(levels(y))
     }
    rocData["sens",j]<-mean(est$predictions[y==pos] == pos)
    rocData["spec",j]<-mean(est$predictions[y==neg] == neg)
  }

}
}
results<-list(call=cl,errorRate=errorRate,rocData=rocData)
class(results)<-"classificationError"
results
}

