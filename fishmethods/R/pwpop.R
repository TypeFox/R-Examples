pwpop<-function(abund=NULL,year=NULL,periods=NULL,Cs=NULL,
	startR=NULL,upperR=NULL,lowerR=NULL,graph=TRUE){   
     if(is.null(abund)) 
         stop ("abund vector does not exist.")
     if(is.null(year)) 
         stop ("year vector does not exist.")
      if(length(abund)!=length(year)) 
        stop("unequal abund/years vector lengths.")
 	 if(is.null(periods)) 
         stop ("periods does not exist.")
       if(periods>1){
        if(is.null(Cs)) stop ("CC does not exist.")
        if(length(Cs)!=as.numeric(periods-1)) stop(paste("There should be ", periods-1," CC values",sep=""))
        }
        if(is.null(startR)) stop ("startR vector does not exist.")
        if(length(startR)!=as.numeric(periods)) stop(paste("There should be ", periods," startR values",sep=""))
	   if(is.null(upperR)) stop ("upperR vector does not exist.")
        if(length(upperR)!=as.numeric(periods)) stop(paste("There should be ", periods," upperR values",sep=""))
	   if(is.null(lowerR)) stop ("lowerRb vector does not exist.")
        if(length(lowerR)!=as.numeric(periods)) stop(paste("There should be ", periods," lowerR values",sep=""))
      N<-abund
      t<-1:length(year)
      d<-as.data.frame(cbind(N,year,t))
      cc<-Cs-year[1]+1
      for(pp in 1:periods){
       if(pp==1){
          parms<-c(log(d$N[1]),startR[pp]) #a1,b1
          lower<-c(log(d$N[1])*-10,lowerR[pp])
          upper<-c(log(d$N[1])*10,upperR[pp]) 
       }
       if(pp>=2) {
        parms<-c(parms,cc[pp-1],startR[pp])#add C, b
        lower<-c(lower,1,lowerR[pp])
	   upper<-c(upper,t[length(t)],upperR[pp])
       }  
     }             
    fitmodel<-function(y){
         # Generate equations
          models<-NULL
         
          for(pp in 1:periods){
             if(pp==1) models[pp]<-paste("a1",sep="") 
             if(pp>1) models[pp]<-paste(models[pp-1],"+C",pp-1,"*(b",pp-1,"-b",pp,")",sep="")         
           }
          for(pp in 1:as.numeric(length(models))) models[pp]<-paste(models[pp],"+b",pp,"*d$t",sep="")
         #generate ifelse levels
          texter<-NULL
          if(periods==1) texter<-paste("d$lpred<<-a1+b1*d$t",sep="")
          if(periods>1){
            for(pp in 1:periods){
             if(pp==1) texter<-paste("d$lpred<<-ifelse(d$t<C",pp,",",models[pp],",",sep="") 
             if(pp>1 & pp<periods) texter<-paste(texter,"ifelse(d&t>=C",pp-1," & d&t<C",pp,",",models[pp],",",sep="")
             if(pp==periods) texter<-paste(texter,models[pp])
            }
            for(pp in 1:as.numeric(periods-1)) texter<-paste(texter,")",sep="") 
        }
        #### set parameters
         if(periods==1) eval(parse(text=paste("a1<-y[1];b1<-y[2]",sep="")))        
         if(periods>1){
            setvar<-NULL
            i<-2
            for(pp in 1:periods){
              if(pp==1) setvar<-paste("a1<-y[1];b1<-y[2]",sep="")
              if(pp>1){
                i<-i+1
                setvar<-paste(setvar,";C",pp-1,"<-y[",i,"]",sep="")
                i<-i+1
                setvar<-paste(setvar,";b",pp,"<-y[",i,"]",sep="")
              }       
            }
           eval(parse(text=setvar))
         }
         eval(parse(text=texter))              
         sum((log(d$N)-d$lpred)^2,na.rm=T)
   }
   results1<-try(optim(parms, fitmodel, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"), 
      control=list(maxit=100000),hessian=TRUE),TRUE)
       #generate labels
       for(pp in 1:periods){
         if(pp==1) labs<-c("logN0","logR1") #a1,b1
         if(pp>=2) labs<-c(labs,paste("C",pp-1,sep=""),paste("logR",pp,sep=""))
       }
       labs<-c(labs,"RSS","AIC","R2")
    if(class(results1)!="try-error"){
       cov1<-solve(results1$hessian)
         outpt<-data.frame(Parameters=labs,
           Estimate=c(round(results1$par,3),sum((log(d$N)-d$lpred)^2,na.rm=T),
                       length(N[!is.na(N)]) *log(sum((log(d$N)-d$lpred)^2,na.rm=T)/ length(N[!is.na(N)]))+(2*length(parms)+1),
                      1-(round(sum((log(d$N)-d$lpred)^2,na.rm=T),2)/round(sum((log(d$N)-mean(log(d$N),na.rm=T))^2,na.rm=T),2))),
                    SE=c(sqrt(diag(cov1)),NA,NA,NA))
       d$pred<-exp(d$lpred)
      if(graph==TRUE){
         plot(d$N~c(d$t+year[1]-1),type="b",col="black",pch=16,ylab="Abundance",xlab="Year")
         newt<-seq(1,t[length(t)],0.01) 
         models<-NULL
          for(pp in 1:periods){
             if(pp==1) models[pp]<-paste("a1",sep="") 
             if(pp>1) models[pp]<-paste(models[pp-1],"+C",pp-1,"*(b",pp-1,"-b",pp,")",sep="")         
           }
          for(pp in 1:as.numeric(length(models))) models[pp]<-paste(models[pp],"+b",pp,"*newt",sep="")
         #generate ifelse levels
          texter<-NULL
          if(periods==1) texter<-paste("lpred<-a1+b1*newt",sep="")
          if(periods>1){
            for(pp in 1:periods){
             if(pp==1) texter<-paste("lpred<-ifelse(newt<C",pp,",",models[pp],",",sep="") 
             if(pp>1 & pp<periods) texter<-paste(texter,"ifelse(newt>=C",pp-1," & newt<C",pp,",",models[pp],",",sep="")
             if(pp==periods) texter<-paste(texter,models[pp])
            }
            for(pp in 1:as.numeric(periods-1)) texter<-paste(texter,")",sep="") 
          }
         #### set parameters
         if(periods==1) eval(parse(text=paste("a1<-results1$par[1];b1<-results1$par[2]",sep="")))        
         if(periods>1){
            setvar<-NULL
            i<-2
            for(pp in 1:periods){
              if(pp==1) setvar<-paste("a1<-results1$par[1];b1<-results1$par[2]",sep="")
              if(pp>1){
                i<-i+1
                setvar<-paste(setvar,";C",pp-1,"<-results1$par[",i,"]",sep="")
                i<-i+1
                setvar<-paste(setvar,";b",pp,"<-results1$par[",i,"]",sep="")
              }       
            }
         
         }
          eval(parse(text=setvar))
          eval(parse(text=texter))
          lines(exp(lpred)~c(newt+year[1]-1),col="red")
     } #graph
    } #no error
    if(class(results1)=="try-error"){
        outpt<-data.frame(Parameters=labs,
           Estimate=c(rep(NA,length(parms)),NA,NA,NA),
                    SE=c(rep(NA,length(parms)),NA,NA,NA))
     }
      out<-list(outpt,d);names(out)<-c("Estimates","Data")
      return(out)
}     
