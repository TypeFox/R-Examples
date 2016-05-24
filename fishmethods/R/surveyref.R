surveyref<-function(x=NULL,refpt=25,compyear=NULL,reffix=FALSE,refrange=NULL,nboot=500,allboots=FALSE){
	 if(is.null(x)) stop ("x must be specified.")
	 if(is.null(x$ARIMA_output)) stop ("x is not a surveyfit output object.")
       if(is.null(refpt)) stop("refpt must be specified.")
       if(is.null(compyear)) stop("compyear must be specified.")
       if(reffix==TRUE & is.null(refrange)) stop("refrange must be specified if reffix is TRUE.")
       if(length(refpt)>1) refpt<-refpt[1]
       if(!all(compyear %in% x$index$year)) stop("compyear not found in time series.")

       if(reffix==TRUE){
          if(is.na(match(refrange[1],x$index$year))) stop("Can not find lower value of refrange survey fit object.")
          if(is.na(match(refrange[2],x$index$year))) stop("Can not find upper value of refrange survey fit object.")
         }
	 allruns<-NULL
	 newpred<-NULL
	 newq<-NULL
       if(reffix==TRUE){
     		  f_yr<-match(refrange[1],x$index$year)
       	  l_yr<-match(refrange[2],x$index$year)
	 	  qvalue<-quantile(x$index$fitted[f_yr:l_yr],refpt/100)
        }
       if(reffix==FALSE) qvalue<-quantile(x$index$fitted,refpt/100)
       nyr<-length(compyear)
	 n<-length(x$index$fitted)
	 for(b in 1:nboot){
   		newy<-x$index$fitted+sample(x$ARIMA_output$residuals,size=n,
                   replace=TRUE)
   		outboot<-arima(newy,order=c(0,1,1),fixed=x$ARIMA_output$coef)
    		for(i in 1:n){
         	 pis<-rep(NA,1,length(seq(i,n,1)))
         	 pis[1]<-1
     	 		if(i<n){
			for(j in 2:as.numeric(length(pis))){
      		   pis[j]<-(1-abs(outboot$coef))*(abs(outboot$coef)^(j-1-1))*-1
       		  }
      	      }
     	        newpred[i]<-newy[i]-abs(outboot$coef)*sum(outboot$residuals[i:n]*pis)     
            }   
         allruns<-cbind(allruns,newpred)  
        if(reffix==TRUE) newq<-rbind(newq,quantile(newpred[f_yr:l_yr],refpt/100))    
        if(reffix==FALSE) newq<-rbind(newq,quantile(newpred,refpt/100))    

}
runs<-as.data.frame(cbind(x$index$year,allruns))
testyear<-as.data.frame(stack(runs[runs[,1] %in% c(compyear),-1])$values)
testyear[,1]<-round(testyear[,1],2)
qdist<-as.data.frame(newq)[,1]
qdist<-round(qdist,2)

ref<-quantile(qdist,c(seq(0,1,0.05)))
probs<-NULL
loopcnt<-length(as.vector(ref))
  for(cl in 1:loopcnt){
   probs[cl]<-length(testyear[testyear[,1]<ref[cl],])/(nboot*length(compyear))
  }
  prob<-data.frame(confidence=c(sort(seq(0,100,5),decreasing=TRUE)),prob=probs)
  prob<-prob[order(prob[,1]),]
	ans<-NULL
      ans$comp_refpt<-as.data.frame(cbind(qvalue,mean(qdist)))
		names(ans$comp_refpt)<-c("orig_q","boot_mean_q")
      ans$comp_fitted<-as.data.frame(cbind(x$index$year,x$index$fitted,rowMeans(allruns)))
	    names(ans$comp_fitted)<-c("year","fitted","boot_mean")
      ans$emp_dist_index<-as.data.frame(table(testyear))
	    names(ans$emp_dist_index)<-c("value","freq")
  	ans$emp_dist_refpt<-as.data.frame(table(qdist))
	    names(ans$emp_dist_refpt)<-c("value","freq")
      ans$prob_index<-prob
      if(allboots==TRUE) ans$boot_runs<-runs
      return(ans)
} #end function

