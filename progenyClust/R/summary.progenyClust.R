summary.progenyClust <-
function(object,...){
  if(!is.null(object$call$score.invert) && (object$call$score.invert==T)){
    if(object$method=='gap'){
      ngap=which(object$mean.gap==min(object$mean.gap))
      nscore=NA
    }else if(object$method=='score'){
      nscore=which(object$mean.score==min(object$mean.score))
      ngap=NA
    }else{
      ngap=which(object$mean.gap==min(object$mean.gap))
      nscore=which(object$mean.score==min(object$mean.score))
    }
  }else{
    if(object$method=='gap'){
      ngap=which(object$mean.gap==max(object$mean.gap))
      nscore=NA
    }else if(object$method=='score'){
      nscore=which(object$mean.score==max(object$mean.score))
      ngap=NA
    }else{
      ngap=which(object$mean.gap==max(object$mean.gap))
      nscore=which(object$mean.score==max(object$mean.score))
    }
  }
  if(!is.na(ngap)){
    ngap=object$ncluster[ngap+1]
  }
  if(!is.na(nscore)){
    nscore=object$ncluster[nscore]
  }
  output<-list(call=object$call,n.gap=ngap,n.score=nscore)
  class(output)<-'summary.progenyClust'
  return(output)
}
