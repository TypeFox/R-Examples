##' @export
summary.confScoreSurv <- function(object,
                                  times,
                                  type=1,
                                  qScore=FALSE,
                                  ...){
  if (type==1){
    meanScore <- do.call("cbind",lapply(object$models,function(m){
      colMeans(m$score)
    }))}
  else{
    meanScore <- do.call("cbind",lapply(object$models,function(m){
      1-sqrt(colMeans((1-m$score)^2))
    }))}
  if (qScore==TRUE){
    qScore <- do.call("cbind",lapply(1:length(object$models),function(m){
      qq <- t(apply(object$models[[m]]$score,2,quantile,c(0.5,.25,.75,0,1)))
      colnames(qq) <- paste(names(object$models)[m],".",c("median","iqrLow","iqrUp","min","max"),"Score",sep="")
      qq
    }))
  }
  mm <- data.frame(times=object$times,
                   meanScore=meanScore)
  if (!missing(times)){
    mm <- rbind(0,mm)[1+prodlim::sindex(jump.times=object$times,eval.times=times),,drop=FALSE]
    if (qScore==TRUE)
      qScore <- rbind(0,qScore)[1+prodlim::sindex(jump.times=object$times,eval.times=times),,drop=FALSE]
  }
  if (qScore==TRUE)
    mm <- cbind(mm,qScore)
  ##   rownames(mm) <- rep("",NROW(mm))
  mm
}
