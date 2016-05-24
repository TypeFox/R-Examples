surveyfit<-function(year=NULL,index=NULL,logtrans=TRUE,graph=TRUE){
	if(length(year)!=length(index)) stop ("Lengths of year and index are different.")
	if(any(is.na(index))) stop("Missing values are not allowed.")
	if(logtrans==TRUE) logn<-log(index) else logn<-index
	n<-length(logn)
	outacf<-acf(diff(ts(logn),lag=1),plot=FALSE)[1:3]
	out<-arima(logn,order=c(0,1,1))
	ests<-c(abs(out$coef),sqrt(out$var.coef))
	varc<-out$sigma2
	table1<-matrix(NA,1L,7L)	
	table1<-cbind(n,round(rbind(outacf$acf),2),round(ests[1],2),round(ests[2],2),
        round(varc,2))
	colnames(table1)<-c("n","r1","r2","r3","theta","SE","sigma2c")
      rownames(table1)<-"Results"
	pred<-NULL
	for(i in 1:n){
          pis<-rep(NA,1,length(seq(i,n,1)))
          pis[1]<-1
      	if(i<n){
			for(j in 2:as.numeric(length(pis))){
        		pis[j]<-(1-ests[1])*(ests[1]^(j-1-1))*-1
       		}
      	}
       pred[i]<-logn[i]-ests[1]*sum(out$residuals[i:n]*pis)
     }
         ans<-NULL  
         ans$results<-table1
	   ans$index<-as.data.frame(cbind(year,logn,pred))
		names(ans$index)<-c("year","log_index","fitted")
         ans$ARIMA_output<-out
    if(graph==TRUE){
       plot(logn~year,ylab="loge(index): (o) observed, (-) predicted",xlab="Year")
       lines(pred~year,col="black")

     }
         return(ans)
}#function end

