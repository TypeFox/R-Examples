#Objects of class 'avetest' and 'maxtest'
#representing results of AVE- and MAX-tests (Hung, 2000)
setClass("avetest",representation(p="numeric",
                                  stat="numeric",
                                  test="character",
                                  method="character",
                                  nboot="numeric",
                                  simerror="numeric",
                                  duration="numeric",
                                  call="call"))
setClass("maxtest",representation(p="numeric",
                                  stat="numeric",
                                  test="character",
                                  name="character",
                                  method="character",
                                  nboot="numeric",
                                  simerror="numeric",
                                  duration="numeric",
                                  call="call"))
#
#
#Calling the right method for the AVE-test and the respective type of
#data, selected test statistic and computation method
setGeneric("avetest",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(all(C@D==c(1,1,1,1,1,1))) warning("The global tests are just single min-tests for a 1x1 design.\n")
  if(method=="bootstrap"&&is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
  if(method=="hung"&&(!is.null(nboot)||!is.null(simerror))) warning("Arguments nboot or simerror ignored.")
  if(method!="hung"&&method!="bootstrap") stop(paste("Computation method '",method,"' unknown.",sep=""))
  if(is.null(test)) stop("Test statistic must be specified.")
  if(test!="ttest"&&test!="ztest") stop(paste("Test statistic '",test,"' unknown.",sep=""))
  if(is.null(simerror)) simerror<-9
  if(is.null(nboot)) nboot<-900
  if(!is.numeric(nboot) || !is.numeric(simerror)) stop("nboot and simerror must be numeric.")
  res<-standardGeneric("avetest")
  res@call=match.call()
  res
})
#
#
#'carpet' objects
setMethod("avetest","carpet",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(method=="bootstrap"){
    if(test=="ttest") result<-avestudent2Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-avebinomial2Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  if(method=="hung"){
    if(is.binary(C@data)) stop("No analytical approach implemented for binary data.")
    if(test=="ttest") result<-avestudent2(C=C,...)
    if(test=="ztest") stop("No analytical approach implemented for the Z-statistic.")
  }
  return(result)
})
#
#
#'cube' objects
setMethod("avetest","cube",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(all(C@D==c(1,1,1))) warning("The global tests are just single min-tests for a 1x1x1 design.\n")
  if(method=="hung") stop("No analytical approach implemented for k=3.")
  if(method=="bootstrap"){
    if(is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
    if(is.null(simerror)) simerror<-9
    if(is.null(nboot)) nboot<-900
    if(test=="ttest") result<-avestudent3Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-avebinomial3Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  return(result)
})
setMethod("avetest","ANY",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  stop("Argument must be an object of class 'carpet' or 'cube'.")
})
#                                        
#
#Calling the right method for the MAX-test and the respective type of
#data, selected test statistic and computation method
setGeneric("maxtest",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(all(C@D==c(1,1,1,1,1,1))) warning("The global tests are just single min-tests for a 1x1 design.\n")
  if(method=="bootstrap"&&is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
  if(method=="hung"&&(!is.null(nboot)||!is.null(simerror))) warning("Arguments nboot or simerror ignored.")
  if(method!="hung"&&method!="bootstrap") stop(paste("Computation method '",method,"' unknown.",sep=""))
  if(is.null(test)) stop("Test statistic must be specified.")
  if(test!="ttest"&&test!="ztest") stop(paste("Test statistic '",test,"' unknown.",sep=""))
  if(is.null(simerror)) simerror<-9
  if(is.null(nboot)) nboot<-900
  if(!is.numeric(nboot) || !is.numeric(simerror)) stop("nboot and simerror must be numeric.")
  res<-standardGeneric("maxtest")
  res@call=match.call()
  res
})
#
#
#'carpet' objects
setMethod("maxtest","carpet",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(method=="bootstrap"){
    if(test=="ttest") result<-maxstudent2Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-maxbinomial2Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  if(method=="hung"){
    if(is.binary(C@data)) stop("No analytical approach implemented for binary data.")
    if(test=="ttest") result<-maxstudent2(C=C,...)
    if(test=="ztest") stop("No analytical approach implemented for the Z-statistic.")
  }
  return(result)
})
#
#
#'cube' objects
setMethod("maxtest","cube",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(all(C@D==c(1,1,1))) warning("The global tests are just single 'min'-tests for a 1x1 design.\n")
  if(method=="hung") stop("No analytical approach implemented for k=3.")
  if(method=="bootstrap"){
    if(is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
    if(is.null(simerror)) simerror<-9
    if(is.null(nboot)) nboot<-900
    if(test=="ttest") result<-maxstudent3Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-maxbinomial3Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  return(result)
})
setMethod("maxtest","ANY",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  stop("Argument must be an object of class 'carpet' or 'cube'.")  
})
#
#
#S4 methods for the 'show' and 'summary' generic functions
#and objects of classes 'avetest' and 'maxtest'
setMethod("show","maxtest",function(object){
  cat("\nMAX-test for the existence of an efficacious combination\n")
  if(object@p<0.0001){
    cat("tmax=",round(object@stat,4),"; pmax<0.0001\n",sep="")
  }
  if(object@p>0.9999){
    cat("tmax=",round(object@stat,4),"; pmax>0.9999\n",sep="")
  }
  if(object@p>=0.0001 && object@p<=0.9999){
    cat("tmax=",round(object@stat,4),"; pmax=",object@p,"\n",sep="")
  }
  cat("Combination where maximum test statistic occurs:",object@name,"\n\n")
})
setMethod("show","avetest",function(object){
  cat("\nAVE-test for the existence of an efficacious combination\n")
  if(object@p<0.0001){
    cat("tave=",round(object@stat,4),"; pave<0.0001\n\n",sep="")
  }
  if(object@p>0.9999){
    cat("tave=",round(object@stat,4),"; pave>0.9999\n\n",sep="")
  }
  if(object@p>=0.0001 && object@p<=0.9999){
    cat("tave=",round(object@stat,4),"; pave=",object@p,"\n\n",sep="")
  }
})
setMethod("summary","maxtest",function(object){
  cat("\nMAX-test for the existence of an efficacious combination\n")
  cat("tmax=",round(object@stat,4),"\n",sep="")
  if(object@p<0.0001) cat("pmax<0.0001\n\n",sep="")
  if(object@p>0.9999) cat("pmax>0.9999\n\n",sep="")
  if(object@p>=0.0001 && object@p<=0.9999) cat("pmax=",object@p,"\n\n",sep="")
   cat("Method:",object@method,"\n")
  if(object@method=="Bootstrap"){
    cat("Combination where maximum test statistic occurs:",object@name,"\n")
    cat("Total number of simulations:",object@nboot,"\n")
  }
  cat("Simulation standard error:",round(max(object@simerror),6),"\n")
  cat("Total computation time:",(object@duration-object@duration%%3600)/3600,"hours,",
      ((object@duration-object@duration%%60)/60)-(object@duration-object@duration%%3600)/60,
      "minutes,",round(object@duration%%60,digits=0), "seconds\n\n")
})
setMethod("summary","avetest",function(object){
  cat("\nAVE-test for the existence of an efficacious combination\n")
  cat("tave=",round(object@stat,4),"\n",sep="")
  if(object@p<0.0001) cat("pave<0.0001\n\n",sep="")
  if(object@p>0.9999) cat("pave>0.9999\n\n",sep="")
  if(object@p>=0.0001 && object@p<=0.9999) cat("pave=",object@p,"\n\n",sep="")
  cat("Method:",object@method,"\n")
  if(object@method=="Bootstrap"){
    cat("Total number of simulations:",object@nboot,"\n")
    cat("Simulation standard error:",round(max(object@simerror),6),"\n")
  }
  cat("Total computation time:",(object@duration-object@duration%%3600)/3600,"hours,",
      ((object@duration-object@duration%%60)/60)-(object@duration-object@duration%%3600)/60,
      "minutes,",round(object@duration%%60,digits=0), "seconds\n\n")
})
