powertrend<-function(trend=1,A1=NULL,PSE=NULL,pserel=1,maxyrs=3,pR=100,
     step=5,alpha=0.05,tail=2,graph=TRUE){
    if(is.null(A1)) stop("Starting Abundance Required.")
    if(is.null(PSE)) stop("PSE Required.")
    if(pserel>3) stop("No such pserel.")
    if(trend>2|trend<1) stop("No such trend.")
    if(tail>2) stop("No such tailed test.")

   tlen<-length(seq(-100,pR,by=step))
   results<-data.frame(years=NA,trend=NA,pse=NA,alpha=NA,tail=NA,R=NA,power=NA)
   place<-0
    
   nyr<-maxyrs	
   R<--100
   for (j in 0:(tlen-1)){ 
   	  s2x<-((nyr+1)*(nyr-1))/12 
	  if(trend==1){
           r<-(R/100/(nyr-1))
	     b<-A1*r
           s2res<-ifelse(pserel==1,((PSE*A1)^2)*(1+(r/2)*(nyr-1)),
           ifelse(pserel==2,((PSE*A1)^2)*(1+r*(nyr-1)*(1+(r/6)*(2*nyr-1))),
            ((PSE*A1)^2)*(1+(3*r/2)*(nyr-1)*(1+(r/3)*(2*nyr-1)+(r^2/6)*nyr*(nyr-1)))))           
        }
        if (trend==2){
            r<-((R/100+1)^(1/(nyr-1)))-1
            b<-log(1+r)
            if(pserel==1){
		   sum<-0
               for(k in 1:maxyrs){
                  sum<-sum+log((PSE^2/((1+r)^(k-1)))+1)
               }
               s2res<-sum/maxyrs
             }
            if(pserel==2){
              s2res<-log(1+PSE^2)
             }
            if(pserel==3){
		  sum<-0
              for(k in 1:maxyrs){
                sum<-sum+log((PSE^2*(1+r)^(k-1))+1)
              }
              s2res<-sum/maxyrs
            }
       }
      sb<-sqrt(s2res/(nyr*s2x))
      delta<-abs(b/sb)
      tdist<-abs(qt(alpha/tail,nyr-2))
      powert<-round(1-pt(tdist,df=nyr-2,ncp=delta),2)
      place<-place+1
      results[place,1]<-maxyrs;results[place,2]<-trend;results[place,3]<-PSE;results[place,4]<-alpha
      results[place,5]<-tail;results[place,6]<-R;results[place,7]<-powert
      R<-R+step
   }
  if(graph==TRUE){
     labels<-paste("Trend: ",ifelse(trend==1,"Linear","Exponential"),
                   " PSE: ",PSE," alpha: ",alpha," tail: ",tail,
                   " pserel: ",pserel)
     plot(results$power~results$R,ylim=c(0,1),type="l",main=labels,
    xlim=c(min(results$R),max(results$R)),ylab="Power",
     xlab=paste("Percent change over ",maxyrs," years"),cex.main=1)

   }
   return(results)
} #end function

