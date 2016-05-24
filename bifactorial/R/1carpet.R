setClass("carpet",representation(data="list",D="numeric",n="numeric",call="call"))
setGeneric("carpet",function(data,D,...){
  res<-standardGeneric("carpet")
  res@call=match.call()
  res
})
setMethod("carpet",signature(data="list",D="numeric"),function(data,D,...){
  n<-numeric(0)
  if(!is.numeric(D)) stop("Dimensions must be specified by integer values.")
  if(length(D)!=2) stop("Two Dimensions must be specified for 'carpet' objects.")
  for(i in 1:length(data)) n[i]<-length(data[[i]])
  res<-new("carpet",data=data,D=D,n=n)
  res
})
setMethod("carpet",signature(data="ANY",D="ANY"),function(data,D,...){
  stop("Data must be a list of numeric data vectors.")
})
setMethod("show","carpet",function(object){
  cb<-function(a,b) (a*(object@D[2]+1))+(b+1)
  cat("\n")
  cat("Carpet size:",object@D[1],"x",object@D[2],"\n\n")
  cat("Sample size allocation matrix:\n\t")
  for(b in 0:object@D[2]){cat(paste(b,"\t",sep=""))}
  cat("\n")
  for(a in 0:object@D[1]){
    cat(a,"\t")
    for(b in 0:object@D[2]){
      cat(object@n[cb(a,b)],"\t")
    }
    cat("\n")
  }
  if(is.binary(object@data[[cb(1,1)]])) cat("\nDescriptive statistics: Event rates\n\t")
  else cat("\nDescriptive statistics: Mean response values\n\t")
  for(b in 0:object@D[2]){cat(paste(b,"\t",sep=""))}
  cat("\n")
  for(a in 0:object@D[1]){
    cat(a,"\t")
    for(b in 0:object@D[2]){
      cat(round(mean(object@data[[cb(a,b)]]),3),"\t")
    }
    cat("\n")
  }
  cat("\nDescriptive statistics: Standard deviations\n\t")
  for(b in 0:object@D[2]){cat(paste(b,"\t",sep=""))}
  cat("\n")
  for(a in 0:object@D[1]){
    cat(a,"\t")
    for(b in 0:object@D[2]){
      if(is.binary(object@data[[cb(a,b)]])){
        p<-mean(object@data[[cb(a,b)]])
        cat(round(sqrt(object@n[cb(a,b)]*p*(1-p)),3),"\t")
      }
      else cat(round(sd(object@data[[cb(a,b)]]),2),"\t")
    }
    cat("\n")
  }
  cat("\n")
})
setMethod("summary",signature("carpet"),function(object){
  cb<-function(a,b) (a*(object@D[2]+1))+(b+1)
  cat("\nCarpet Densions:",object@D[1],"x",object@D[2],"\n")
  cat("Total sample size:",sum(object@n),"\n\n")
  if(is.binary(object@data[[cb(1,1)]])){
    cat("Event rates and standard deviations (in parentheses)\nin the treatment groups:\n\n\t")
  }
  else cat("Mean response values and standard deviations (in parentheses)\nin the treatment groups:\n\n\t")
  for(b in 0:object@D[2]){cat(paste(b,"\t\t",sep=""))}
  cat("\n")
  for(a in 0:object@D[1]){
    cat(a,"\t")
    for(b in 0:object@D[2]){
      if(is.binary(object@data[[cb(a,b)]])){
        p<-mean(object@data[[cb(a,b)]])
        cat(paste(round(p,3)," (",round(sqrt(object@n[[cb(a,b)]]*p*(1-p)),3),")\t",sep=""))
      }
      else{
        cat(paste(round(mean(object@data[[cb(a,b)]]),1)," (",round(sd(object@data[[cb(a,b)]]),2),")","\t",sep=""))
      }
    }
    cat("\n")
  }
  cat("\n")
})
setMethod("plot",signature(x="carpet",y="missing"),function(x,y){
  cb<-function(a,b) (a*(x@D[2]+1))+(b+1)
  xwert<-ywert<-zwert<-numeric(0)
  for(j in 0:x@D[2]){for(i in 0:x@D[1]){
    xwert<-c(xwert,i)
    ywert<-c(ywert,j)
    zwert<-c(zwert,mean(x@data[[cb(i,j)]]))
  }}
  wireframe(zwert~xwert*ywert,
            data=data.frame(zwert,xwert,ywert),
            scales=list(arrows=FALSE),
            xlab="Dose A",ylab="Dose B",zlab="Response",
            zlim=c(min(zwert)-(.25*(max(zwert)-min(zwert))),max(zwert)+(.25*(max(zwert)-min(zwert)))),
            drape=TRUE,colorkey=FALSE)#,bg="red")
})
