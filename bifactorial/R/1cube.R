setClass("cube",representation(data="list",D="numeric",n="numeric",call="call"))
setGeneric("cube",function(data,D,...){
  res<-standardGeneric("cube")
  res@call=match.call()
  res
})
setMethod("cube",signature(data="list",D="numeric"),function(data,D,...){
  n<-numeric(0)
  if(!is.numeric(D)) stop("Dimensions must be specified by integer values.")
  if(length(D)!=3) stop("Three dimensions must be specified for 'cube' objects.")
  for(i in 1:length(data)) n[i]<-length(data[[i]])
  new("cube",data=data,D=D,n=n)
})
setMethod("cube",signature(data="ANY",D="ANY"),function(data,D,...){
  stop("Need a list of numeric data and an integer Densions vector.")
})
setMethod("show",signature("cube"),function(object){
  cb<-function(a,b,c){(a*(object@D[2]+1)*(object@D[3]+1))+(b*(object@D[3]+1))+c+1}
  stdev<-numeric(0)
  cat("\n")
  cat("Cube size:",object@D[1],"x",object@D[2],"x",object@D[3],"\n\n")
  if(is.binary(object@data)){
    cat("Group\t\tn\t\trate\t\tstdev\n")
    for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
      p<-mean(object@data[[cb(a,b,c)]])
      stdev<-c(stdev,round(sqrt(object@n[cb(a,b,c)]*p*(1-p)),3))
    }}}
  }
  else{
    cat("Group\t\tn\t\tmean\t\tstdev\n")
    for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
      stdev<-c(stdev,round(sd(object@data[[cb(a,b,c)]]),2))
    }}}
  }
  for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
    cat("(",a,",",b,",",c,")\t\t",object@n[cb(a,b,c)],"\t\t",round(mean(object@data[[cb(a,b,c)]]),3),"\t\t",stdev[cb(a,b,c)],"\n",sep="")
  }}}
  cat("\n")
})
setMethod("summary",signature("cube"),function(object){
  cb<-function(a,b,c){(a*(object@D[2]+1)*(object@D[3]+1))+(b*(object@D[3]+1))+c+1}
  cat("\nCube dimensions:",object@D[1],"x",object@D[2],"x",object@D[3],"\n",sep="")
  cat("Total sample size:",sum(object@n),"\n\n")
  stdev<-numeric(0)
  if(is.binary(object@data)){
    cat("Group\trate (stdev)\n")
    for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
      p<-mean(object@data[[cb(a,b,c)]])
      stdev<-c(stdev,round(sqrt(object@n[cb(a,b,c)]*p*(1-p)),3))
    }}}
  }
  else{
    cat("Group\tmean (stdev)\n")
    for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
      stdev<-c(stdev,round(sd(object@data[[cb(a,b,c)]]),2))
    }}}
  }
  for(a in 0:object@D[1]){for(b in 0:object@D[2]){for(c in 0:object@D[3]){
    cat("(",a,",",b,",",c,")\t ",round(mean(object@data[[cb(a,b,c)]]),3)," (",stdev[cb(a,b,c)],")\n",sep="")
  }}}
  cat("\n")
})
setMethod("plot",signature(x="cube",y="missing"),function(x,y){
  stop("Plotting of factorial designs is available for k=2 only.")
})
