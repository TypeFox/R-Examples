#Objects of class 'margint' representing results of
#confidence interval estimation on bifactorial designs
setClass("margint",
         representation(kiu="numeric",
                        kio="numeric",
                        alpha="numeric",
                        cnames="character",
                        test="character",
                        method="character",
                        nboot="numeric",
                        simerror="numeric",
                        duration="numeric",
                        call="call"))
#
#
#Calling the right method for the marginal intervals and the respective type of
#data, selected test statistic and computation method
setGeneric("margint",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,alpha=.05,...){
  if(method=="bootstrap"&&is.null(nboot)&&is.null(simerror)) stop("Either nboot or simerror must be specified.")
  if(method=="tdistr"&&(!is.null(nboot)||!is.null(simerror))) warning("Arguments nboot or simerror ignored.")
  if(method!="tdistr"&&method!="bootstrap") stop(paste("Computation method '",method,"' unknown.",sep=""))
  if(is.null(test)) stop("Test statistic must be specified.")
  if(test!="ttest"&&test!="ztest") stop(paste("Test statistic '",test,"' unknown.",sep=""))
  if(is.null(simerror)) simerror<-9
  if(is.null(nboot)) nboot<-900
  if(!is.numeric(nboot) || !is.numeric(simerror)) stop("nboot and simerror must be numeric.")
  res<-standardGeneric("margint")
  res@call=match.call()
  res
})
#
#
#'carpet' objects
setMethod("margint","carpet",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,alpha=.05,...){
  if(method=="bootstrap"){
    if(test=="ttest"){result<-intstudent2Boot(C,nboot=nboot,simerror=simerror,alpha=alpha,...)}
    if(test=="ztest"){result<-intbinomial2Boot(C,nboot=nboot,simerror=simerror,alpha=alpha,...)}
  }
  if(method=="tdistr"){
    if(is.binary(C@data)) stop("No analytical approach implemented for binary data.")
    if(test=="ttest") result<-intstudent2(C,alpha=alpha,...)
    if(test=="ztest") stop("No analytical approach implemented for the Z-statistic.")
  }
  return(result)
})
#
#
#'cube' objects
setMethod("margint","cube",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,alpha=.05,...){
  if(method=="bootstrap"){
    if(test=="ttest"){result<-intstudent3Boot(C,nboot=nboot,simerror=simerror,alpha=alpha,...)}
    if(test=="ztest"){result<-intbinomial3Boot(C,nboot=nboot,simerror=simerror,alpha=alpha,...)}
  }
  if(method=="tdistr"){ stop("No analytical methods implemented for k=3.") }
  return(result)
})
setMethod("margint","ANY",function(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,alpha=.05,...){
  stop("Need an object of class 'carpet' or 'cube'.")
})
#
#
#S4 methods for the 'show', 'summary' and 'plot'
#generic functions and objects of class 'margint'
setMethod("show","margint",
function(object){
  cat("\nContrast\tConfidence interval\n")
  for(l in 1:length(object@cnames)){
    cat(object@cnames[l],"\t[",round(object@kiu[l],digits=4),"; ",round(object@kio[l],digits=4),"]\n",sep="")
  }
  cat("\n")
})
setMethod("summary","margint",function(object){
  cat("\nContrast\t Confidence interval\n")
  for(l in 1:length(object@cnames)){
    cat(object@cnames[l],"\t [",round(object@kiu[l],digits=4),"; ",round(object@kio[l],digits=4),"]\n",sep="")
  }
  cat("\n")
  cat("Method:",object@method,"\n")
  if(object@method=="Bootstrap"){
    cat("Number of simulations: ",object@nboot,"\n",sep="")
    cat("Simulation standard error: ",round(object@simerror,6),"\n",sep="")
  }
  cat("Nominal level of the simultanous intervals: ",1-object@alpha,"\n",sep="")
  cat("Total computation time:",(object@duration-object@duration%%3600)/3600,"hours,",
      ((object@duration-object@duration%%60)/60)-(object@duration-object@duration%%3600)/60,
      "minutes,",round(object@duration%%60,digits=0), "seconds\n\n")
})
setMethod("plot",signature("margint","missing"),function(x,y){
  plot(c(1,1),c(x@kiu[1],x@kio[1]),type="l",xlim=c(0,length(x@cnames)+1),ylim=c(min(x@kiu)-1,max(x@kio)+1),xaxt="n",ylab="Value of contrast",xlab="")
  lines(c(.9,1.1),c(x@kiu[1],x@kiu[1]))
  lines(c(.9,1.1),c(x@kio[1],x@kio[1]))
  for(i in 2:length(x@cnames)){
    lines(c(i,i),c(x@kiu[i],x@kio[i]))
    lines(c(i-.1,i+.1),c(x@kiu[i],x@kiu[i]))
    lines(c(i-.1,i+.1),c(x@kio[i],x@kio[i]))
  }
  lines(c(0,length(x@cnames)+1),c(0,0),lty="dashed")
  axis(1,
       at=c(1:length(x@cnames)),
       labels=x@cnames,
       par("las"=2),
       padj=rep(1,length(x@cnames)))
})
