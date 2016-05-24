predict.indicators<-function(object, newdata, ...) {
  C <- object$C
  nc<-nrow(C)
  taxnames <-colnames(C)
  if(sum(taxnames %in% colnames(newdata))<length(taxnames)) stop("Not all taxon names that form indicators could be found in compositional table.")
  postaxa <- numeric(length(taxnames))
  for(i in 1:length(taxnames)) postaxa[i] = which(colnames(newdata)==taxnames[i])
  stat<-matrix(0, nrow=nrow(newdata), ncol=nc)
  if(length(dim(object$A))==2) {
    lowerCI<-matrix(0, nrow=nrow(newdata), ncol=nc)
    upperCI<-matrix(0, nrow=nrow(newdata), ncol=nc)    
  }
  for(i in 1:nc) {
    ci <- as.numeric(C[i,])
    found<-rowSums(newdata[,postaxa[ci==1], drop=FALSE]>0)==sum(ci) #Is the species combination there? (for each site)
    if(length(dim(object$A))==2) {
      stat[,i]<- found*object$A$stat[i]
      lowerCI[,i]<- found*object$A$lowerCI[i]
      upperCI[,i]<- found*object$A$upperCI[i]
    }
    else {
      stat[,i]<- found*object$A[i]
    }
  }
  ## Returns the maximum probability among all combinations explored
  ## If confidence intervals are available
  if(length(dim(object$A))==2) {
    wm = apply(lowerCI,1,which.max)
    p = matrix(NA,nrow=nrow(newdata), ncol=3)
    for(i in 1:length(wm)) {
      p[i,] = c(stat[i,wm[i]],lowerCI[i,wm[i]], upperCI[i,wm[i]])
    }
    rownames(p)<-rownames(newdata)
    colnames(p)<-c("Prob.", "lowerCI", "upperCI")
  } else {
    p = apply(stat,1,max)
    names(p)<-rownames(newdata)
  }
  return(p)
}