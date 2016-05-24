slca<-function(x,type=1,fryr=NULL,Linf=NULL,K=NULL,t0=NULL,Lrange=NULL,Krange=NULL){
      if(!is.data.frame(x)) 
         stop ("dataframe does not exist") 
	if(is.null(fryr))
         stop("fractions of year values missing")
      if(type==2 & is.null(Linf))
         stop("Linf value absent")
      if(type==2 & is.null(K))
         stop("K value absent")
      if(type==2 & is.null(t0))
         stop("t0 value absent")
      if(type==1 & is.null(Lrange))   stop("Lrange values missing")
	if(type==1 & length(Lrange)!=3) stop ("Missing Lrange values.")
      if(type==1 & is.null(Krange)) stop("Krange values missing")
	if(type==1 & length(Krange)!=3) stop ("Missing Krange values.")
	seas<-as.numeric(fryr)
	rownames(x)<-c(seq(1,length(x[,1]),1))
	ndists<-ncol(x)-1;nlens<-length(x[,1])
	x<-rbind(x,c(x[nlens,1]+(x[2,1]-x[1,1])
        	,rep(0,ndists)))
	score<-numeric()
	ts<-data.frame(tmin=rep(0,nlens+1),tmax=rep(0,nlens+1))

	sumprod<-function(Ki,Li,tt){
   		RK<-1/Ki
  		L0<-Li*(1-exp(Ki*tt))
    		ts$tmin<<-ifelse(Li>x[,1],RK*log(ifelse(Li>x[,1],
         		(Li-L0)/(Li-x[,1]),NA)),0)
   		ts$tmax[seq(1,nlens,1)]<<-ifelse(x[seq(2,nlens+1,1),1]<Li,
                ts$tmin[seq(2,nlens+1,1)],0)
   		ts$delta<<-ifelse(ts$tmin==0|ts$tmax==0,0,ts$tmax-ts$tmin)
   		ts$tbar<<-ifelse(ts$delta!=0,(ts$tmax+ts$tmin)/2,0)
   		ts$diffn<<-ifelse(ts$delta!=0,sin(pi*ts$delta)/(pi*ts$delta),0)
   		discr<-data.frame(length=x[,1])
  		for(d in 1:ndists){
     		    discr[,ncol(discr)+1]<-ts$diffn*cos(2*pi*(ts$tbar-seas[d]))*
			sqrt(x[,d+1])
     		}
   		return(sum(rowSums(discr[,2:as.numeric(ndists+1)])))
	 }
  
     if(type==1){
    	   L1<-as.numeric(Lrange[1])
    	   L2<-as.numeric(Lrange[2])
         int_L<-as.numeric(Lrange[3])
         K1<-as.numeric(Krange[1])
         K2<-as.numeric(Krange[2])
         int_K<-as.numeric(Krange[3])
     	  Lcnt<-((L2-L1)/int_L)+1
     	  Kcnt<-((K2-K1)/int_K)+1
     	  LKtable<-data.frame(Li=-1,Ki=-1,t0i=-1,Smax=-1)
    	  for(i in 1:Lcnt){
       	if(i==1) Li<-L1 else Li<-Li+int_L
       	  for(j in 1:Kcnt){
         		if(j==1) Ki<-K1 else Ki<-Ki+int_K
         	    for(m in 1:2){
           		 if(m==1)tt<-0
            	 if(m==2)tt<-0.25
            	 score[m]<-sumprod(Ki,Li,tt)
            	 Smax<-round(sqrt(score[1]^2+score[2]^2),1)
           		 t0i<-round(ifelse((atan2(score[2],score[1])/(2*pi))<0,
                	 (atan2(score[2],score[1])/(2*pi))+1,
                 	 atan2(score[2],score[1])/(2*pi)),2)
                }
         	    LKtable<-rbind(LKtable,c(Li,Ki,t0i,Smax))
              }
         }
 	LKtable<-LKtable[LKtable$Ki!=-1,]
 	LKS<-tapply(LKtable$Smax,list(LKtable$Li,LKtable$Ki),sum)
 	LKt<-tapply(LKtable$t0i,list(LKtable$Li,LKtable$Ki),sum)
      output<-list(LKS,LKt);names(output)<-c("Smax by Linf and K",
            "t0 by Linf and K")
    }
	
    if(type==2){
     sumprod(Li=Linf,Ki=K,tt=t0)
     incage<-1-trunc(ts$tmin[1]-0.5)
     outs<-data.frame(d=-1,age=-1,num=-1)
     for(d in 1:ndists){
        ages<-x[-length(x[,1]),c(1,d+1)]
        for(l in 1:length(ages[,1])){
       	ages$agemin[l]<-ts$tmin[l]-seas[d]+0.5
       	ages$minage[l]<-trunc(ages$agemin[l])
       	ages$agemax[l]<-ts$tmax[l]-seas[d]+0.5
       	ages$maxage[l]<-trunc(ages$agemax[l])
       	ages$nage[l]<-ages$maxage[l]-ages$minage[l]+1
         	for(age in as.numeric(ages$minage[l]):as.numeric(ages$maxage[l])){
           		ages$fract[l]<-1
           		if(ages$nage[l]>1) ages$fract[l]<-1/(ages$agemax[l]-ages$agemin[l])
          	      if(ages$nage[l]>1 & age==ages$minage[l]) 
                     ages$fract[l]<-ages$fract[l]*(ages$minage[l]+1-ages$agemin[l])
                  if(ages$nage[l]>1 & age==ages$maxage[l]) 
                    ages$fract[l]<-ages$fract[l]*(ages$agemax[l]-ages$maxage[l])
           		outs<-rbind(outs,c(d,age,ages$fract[l]*ages[l,2]))
             }
        }  
     }
   ACout<-ifelse(is.na(tapply(outs[-1,3],list(outs[-1,2],outs[-1,1]), sum)
	),0,round(tapply(outs[-1,3],list(outs[-1,2],outs[-1,1]), sum),1))
   output<-list(ACout);names(output)<-c("Numbers-at-age by Seasonal Distribution")
  }
  return(output)
}
