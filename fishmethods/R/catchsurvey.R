catchsurvey<-function(year=NULL, catch=NULL, recr=NULL, post=NULL, M=NULL,T=NULL,
        phi=NULL, w=1, initial=c(NA,NA,NA),uprn=NA,graph=TRUE){
      if(is.null(year)) stop("Missing year vector")
  	if(is.null(catch)) stop("Missing catch vector")
	if(is.null(recr)) stop("Missing recruits vector")
	if(is.null(post)) stop("Missing post-recruits vector")
      if(is.null(M)) stop("M is missing")
      if(is.null(T)) stop("T is missing")
      if(is.null(phi)) stop("phi is missing")
      if(any(is.na(initial))) stop("initial values of recr and post missing")
      if(is.na(uprn)) stop("upper bounds for recr and post missing")
      if(any(c(length(catch),length(year),length(recr),length(post)) %in% length(catch)=="FALSE"))
         stop("vectors lengths are different")
      outs<-NULL
      yrs<-length(year)
      lower<-c(rep(1,yrs),0)	
      upper<-c(rep(uprn,yrs),1)
	Rs<-rep(initial[1],yrs-1)
	nest<-c(rep(0,yrs))
	Z<<-c(rep(0,yrs-1))
	parms<-c(Rs,initial[2],initial[3])			
      TPSS<-NULL;TRSS<-NULL;rest<-NULL
	cs<-function(x){						 			
   		rest<<-x[1:yrs-1]
 		for(i in 1:yrs){
   		   if(i==1) nest[1]<<-x[yrs]
   		   if(i>1){
    		     nest[i]<<-max(0.001,(nest[i-1]+rest[i-1]/phi)*
		      exp(-M)-x[yrs+1]*catch[i-1]*exp(-M*(1-T)))
    		   }
 	      }
           TPSS<-0;TRSS<-0
           TPSS<-sum((log(post[seq(1,yrs,1)])-log(nest[seq(1,yrs,1)]))^2,na.rm=T)
           TRSS<-sum(((log(recr[seq(1,yrs-1,1)])-log(rest[seq(1,yrs-1,1)]))*w)^2,na.rm=T)
           TPSS+TRSS
         } 
	outs<-optim(parms, cs, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"), 
          control=list(maxit=100000),hessian=T)
      cov<-solve(outs$hessian)

	R<-rest/(outs$par[yrs+1]*phi)
	N<-nest/outs$par[yrs+1]
	TA<-R+N[-c(yrs)]
	Z[seq(1,yrs-1,1)]<--log(N[seq(2,yrs,1)]/(R[seq(1,yrs-1,1)]+N[seq(1,yrs-1,1)]))
	Fmort<-Z-M
      q<-outs$par[yrs+1]
      out2<-cbind(rest, R, TA,Z,Fmort)
      out3<-cbind(nest,N)
      output<-list(q,out2,out3)
      names(output)<-c("q","Estimates 1","Estimates 2")
     if(graph==TRUE){
       par(mfrow=c(2,2))
       #Plot Observed versus Predicted Recruit Indices
       plot(x=year,y=recr,col="black",ylab="Survey Index", xlab="Year",
        main="Recruit (O=black, P=red)",ylim=c(0,max(recr,out2[,1])))
        lines(x=year[1:yrs-1],y=out2[,1],col="red")

		#Plot Observed versus Predicted Post-Recruit Indices
		plot(x=year,y=post,col="black",ylab="Survey Index", xlab="Year",
         main="Post-Recruit (O=black, P=red)",ylim=c(0,max(post,out3[,1])))
		lines(x=year,y=out3[,1],col="red")

		#Plot Total Abundance
		plot(x=year[1:yrs-1],y=out2[,3],type="l",col="black",ylab="Numbers", xlab="Year",
     	main="Stock Abundance",xlim=c(min(year),max(year)))

		#Plot Fishing Mortality
		plot(x=year[1:yrs-1],y=out2[,5],type="l",col="black",ylab="F", xlab="Year",
           main="Fishing Mortality",xlim=c(min(year),max(year[1:yrs-1])))
     }
    
return(output)
}
