ConfInt.Cindex <- function(x,times,ref=1,level=.95,digits=3,...){
  lower <- (1-level)/2
  upper <- 1-lower
  # median <- .5
  oob=x$BootCvCindexList
  if (is.null(oob)) stop("Out of bag matrix missing. Set 'keep.matrix' to TRUE.")
  ttt=x$time
  mmm <- names(oob)
  at <- prodlim::sindex(jump.times=ttt,eval.times=times)
  meanOob <- do.call("cbind",x$PredCindex)
  out <- lapply(at,function(a){
    meanDiff <- meanOob[a,ref]-meanOob[a,-ref]
    aResult <- do.call("cbind",lapply(oob,function(x){x[a,]}))
    aref <- aResult[,ref]
    adiff <- data.frame(aref-aResult[,-ref])
    aCI <- do.call("rbind",lapply(adiff,function(x){
      # quantile(x,c(median,lower,upper),na.rm=TRUE)
      quantile(x,c(lower,upper),na.rm=TRUE)
    }))
    a.out <- cbind(meanDiff,aCI)
    colnames(a.out) <- c("diff",paste(c("lower","upper"),level*100,sep="."))
    rownames(a.out) <- paste(mmm[ref],mmm[-ref],sep=" vs ")
    a.out
  })
  names(out) <- paste("time:",times)
  lapply(1:length(out),function(i){
    cat("\n\n")
    cat(names(out)[i])
    cat("\n")
    print(out[[i]],digits=digits)})
  class(out) <- "CiCindex"
  invisible(out)
}
