ypr<-function(age=NULL,wgt=NULL,partial=NULL,M=NULL,plus=FALSE,oldest=NULL,maxF=2,incrF=0.001,graph=TRUE){										          
	if(is.null(age)) 
         stop ("age vectoris missing") 
  	if(is.null(age)) 
         stop ("wgt vector is missing") 
     	if(is.null(partial)) 
         stop ("partial recruitment vector is missing")
  	if(is.null(M)) 
         stop ("M value or vector is missing")
      if(plus==TRUE & is.null(oldest)) stop("oldest must be specified for plus group calculation.")    
      if(any(length(age)!=c(length(age),length(wgt),length(partial))))
         stop("Length of vectors unequal")
      if(length(M)==1) M<-rep(M,length(age))
        data<-as.data.frame(cbind(age,wgt,partial,M))
	YPR<-as.data.frame(cbind(rep(NA,ceiling(maxF/incrF)+1),rep(NA,ceiling(maxF/incrF)+1)))
	names(YPR)<-c("F","YPR")
        if(plus==TRUE){
                 len<-oldest-min(data$age)+1
                 if(oldest>max(data$age)){
                 pdata<-data[rep(length(data$age),times=oldest-data$age[length(data$age)]), ] 
                 pdata$age<-seq(max(data$age)+1,oldest,1)
                 data<-rbind(data,pdata)}
           }
        if(plus==FALSE) len<-max(data$age)-min(data$age)+1
	  F<-0.0					      
	  for(i in 1:length(YPR$F)){
            data$S<-exp(-data$partial*F-data$M) 
            data$psb[1]<-1
  		for(y in 2:len){
      	   data$psb[y]<-data$S[y-1]
    		}
         	data$psb<-cumprod(data$psb)
   		data$YPR<-((data$partial*F)/(data$partial*F+data$M))*
               (1-exp(-data$partial*F-data$M))*data$psb*data$wgt
    		YPR$YPR[i]<-sum(data$YPR)
    		YPR$F[i]<-F
    		F<-F+incrF
  	 }
     	 Ymax<-max(YPR$YPR)
	 Fmax<-YPR$F[which(YPR$YPR==max(YPR$YPR))]
	 s10<-((YPR$YPR[2]-YPR$YPR[1])/(YPR$F[2]-YPR$F[1]))*0.10

	F10<-Fmax/2;df<-F10/2;ok<-0;fuzz<-0.0001                
	while (ok==0){
    		  data$s1<-exp(-data$partial*F10-data$M)
    		  data$s1p<-1
    		  for(y in 2:len){
       		data$s1p[y]<-data$s1[y-1]
      }
    	data$s1p<-cumprod(data$s1p)
   	d1<-sum(((data$partial*F10)/(data$partial*F10+data$M))*
             (1-exp(-data$partial*F10-data$M))*data$s1p*data$wgt)
    	data$s2<-exp(-data$partial*(F10+0.0001)-data$M)
    	data$s2p<-1
    	for(y in 2:len){
     	    data$s2p[y]<-data$s2[y-1]
     	}
    	data$s2p<-cumprod(data$s2p)
    	d2<-sum(((data$partial*(F10+0.0001))/(data$partial*(F10+0.0001)+data$M))*
            	(1-exp(-data$partial*(F10+0.0001)-data$M))*data$s2p*data$wgt)
   	slope<-(d2-d1)/((F10+0.0001)-F10)
   	if(abs(s10-slope)<=fuzz) ok<-1
      if(ok==0){
        if(slope>s10) F10<-F10+df
    	  if(slope<s10) F10<-F10-df
    	  df<-df/2 
    	}
    	Y10<-d1  
  }
         ans<-NULL  
         ans<-matrix(NA,2L,2L)
         ans<-rbind(cbind(F10,Y10),cbind(Fmax,Ymax))   
         dimnames(ans)<-list(cbind("F0.10","Fmax"),c("F","Yield_Per_Recruit"))  
     outpt<-list(ans,YPR);names(outpt)<-c("Reference_Points","F_vs_YPR")
    if(graph==TRUE) plot(YPR[,2]~YPR[,1],ylab="Yield-Per-Recruit",xlab="F",type="l")
    return(outpt)
}


