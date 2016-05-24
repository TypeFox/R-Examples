meanranks<-function(x,allranks=NULL,labels=NULL,rankx=FALSE) {
 if(missing(x))
  stop("x must be a vector, matrix or data frame of rankings")
 dimx<-dim(x)
 if(is.null(allranks)) {
  if(is.null(dimx)) allranks<-1:length(x)
  else allranks<-1:dimx[2]
 }
 if(is.null(dimx)) xx <- x
 else xx<-as.matrix(x)
 if(any(is.na(x))) xx<-muranks(x,allranks=allranks,rankx=rankx)
 if(is.null(dimx)) mranks <- mean(xx)
 else mranks<-apply(xx,2,mean)
 if(is.null(labels)) {
  if(is.null(colnames(xx))) labels<-allranks
  else labels<-format(colnames(xx))
 }
 mrx<-list(ranks=xx,labels=labels,mean.ranks=mranks)
 class(mrx)<-"meanranks"
 return(mrx)
}

print.meanranks<-function(x,...) {
 rankpos <- order(x$mean.ranks)
 cat("\n\tMean rank\n")
 for(i in 1:length(x$mean.ranks))
  cat(x$labels[rankpos[i]],"\t",round(x$mean.ranks[rankpos[i]],2),"\n",sep = "")
}
