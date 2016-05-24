#Objects of class 'mintest' representing results for
#calculations of min-tests (Hung, 2000)
setClass("mintest",
         representation(p="numeric",
                        stat="numeric",
                        gnames="character",
                        test="character",
                        method="character",
                        nboot="numeric",
                        simerror="numeric",
                        duration="numeric",
                        call="call"))
#
#
#Calling the right method for the min-test and the respective type of
#data, selected test statistic and computation method
setGeneric("mintest",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(method=="bootstrap"&&is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
  if(method=="hung"&&(!is.null(nboot)||!is.null(simerror))) warning("Arguments nboot or simerror ignored.")
  if(method!="hung"&&method!="bootstrap") stop(paste("Computation method '",method,"' unknown.",sep=""))
  if(is.null(test)) stop("Test statistic must be specified.")
  if(test!="ttest"&&test!="ztest") stop(paste("Test statistic '",test,"' unknown.",sep=""))
  if(is.null(simerror)) simerror<-9
  if(is.null(nboot)) nboot<-900
  if(!is.numeric(nboot) || !is.numeric(simerror)) stop("nboot and simerror must be numeric.")
  res<-standardGeneric("mintest")
  res@call=match.call()
  res
})
#
#
#'carpet' objects
setMethod("mintest",signature(C="carpet"),function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(method=="bootstrap"){
    if(test=="ttest") result<-pvstudent2Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-pvbinomial2Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  if(method=="hung"){
    if(is.binary(C@data)) stop("No analytical approach implemented for binary data.")
    if(test=="ttest") result<-pvstudent2(C=C,...)
    if(test=="ztest") stop("No analytical approach implemented for the Z-statistic.")
  }
  return(result)
})
#
#
#'cube' objects
setMethod("mintest",signature(C="cube"),function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  if(method=="bootstrap"){
    if(test=="ttest") result<-pvstudent3Boot(C=C,nboot=nboot,simerror=simerror,...)
    if(test=="ztest") result<-pvbinomial3Boot(C=C,nboot=nboot,simerror=simerror,...)
  }
  if(method=="hung"){ stop("No analytical methods implemented for k=3.") }
  return(result)
})
setMethod("mintest",signature(C="ANY"),function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...){
  stop("Need an object of class 'carpet' or 'cube'.")
})
#
#
#S4 methods for the 'show', 'summary' and 'plot'
#generic functions and objects of class 'mintest'
setMethod("show","mintest",function(object){
  cat("\nCombination\t p-value\n")  
  for(i in 1:length(object@gnames)){
    if(object@p[i]<.0001){
      cat(fixdigit(object@gnames[i],7),"\t <0.0001\n")
    }
    if(object@p[i]>.9999){
      cat(fixdigit(object@gnames[i],7),"\t >0.9999\n")
    }
    if(object@p[i]>=.0001 && object@p[i]<=.9999){
      cat(fixdigit(object@gnames[i],7),"\t",fixdigit(object@p[i],6),"\n")
    }
  }
  cat("\n")
})
setMethod("summary","mintest",
function(object){
  cat("\nCombination\t p-value\t Statistic\n")
  for(i in 1:length(object@gnames)){
    if(object@p[i]<.0001){
      cat(fixdigit(object@gnames[i],7),"\t <0.0001\t",round(object@stat[i],4),"\n")
    }
    if(object@p[i]>.9999){
      cat(fixdigit(object@gnames[i],7),"\t >0.9999\t",round(object@stat[i],4),"\n")
    }
    if(object@p[i]>=.0001 && object@p[i]<=.9999){
      cat(fixdigit(object@gnames[i],7),"\t",fixdigit(object@p[i],6),"\t",round(object@stat[i],4),"\n")
    }
  }
  cat("\n")
  cat("Method:",object@method,"\n")
  if(object@method=="Bootstrap"){
    cat("Number of simulations: N=",object@nboot,"\n",sep="")
    cat("Maximum simulation standard error:",round(max(object@simerror),6),"\n",sep="")
  }
  cat("Total computation time:",(object@duration-object@duration%%3600)/3600,"hours,",
      ((object@duration-object@duration%%60)/60)-(object@duration-object@duration%%3600)/60,
      "minutes,",round(object@duration%%60,digits=0), "seconds\n\n")
})
setMethod("plot",signature("mintest","missing"),function(x,y){
  plot(c(1:length(x@gnames)),x@p,
       pch=3,
       col=rep("red",length(x@gnames)),
       xaxt="n",
       xlab="Combination groups",
       ylab="Significance level (multiplicity-adjusted)")
  lines(c(1,length(x@gnames)),c(.1,.1),lty=2,col="green3")
  lines(c(1,length(x@gnames)),c(.05,.05),lty=2,col="green4")
  lines(c(1,length(x@gnames)),c(.025,.025),lty=2)
  axis(1,
       at=c(1:length(x@gnames)),
       labels=x@gnames,
       par("las"=2),
       padj=rep(1,length(x@gnames)))
})
