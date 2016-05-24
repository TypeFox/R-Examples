jomo1con.MCMCchain<- function(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL,start.imp=NULL, nburn=100,output=1, out.iter=10) {
  if (is.null(X)) X=matrix(1,nrow(Y),1)
  if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
  if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
  if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
  for (i in 1:ncol(X)) {
    if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
  }
  stopifnot(nrow(Y)==nrow(X), nrow(beta.start)==ncol(X), ncol(beta.start)==ncol(Y),nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(Y), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start))
  betait=matrix(0,nrow(beta.start),ncol(beta.start))
  for (i in 1:nrow(beta.start)) {
    for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
  }
  covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  for (i in 1:nrow(l1cov.start)) {
    for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
  }   
  nimp=1
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  Y<-data.matrix(Y)
  storage.mode(Y) <- "numeric"    
  X<-data.matrix(X)
  storage.mode(X) <- "numeric"    
  if (output!=1) out.iter=nburn+2
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(X), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(0, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=1
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
  betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),nburn))
  omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),nburn))
  meanobs<-colMeans(Y,na.rm=TRUE)
  if (!is.null(start.imp)) {
    start.imp<-as.matrix(start.imp)
    if ((nrow(start.imp)!=nrow(Yimp))||(ncol(Yimp)!=ncol(start.imp))) {
      cat("start.imp dimensions incorrect. Not using start.imp as starting value for the imputed dataset.\n")
      start.imp=NULL
    } else {
      Yimp<-start.imp
    }
  }
  if (is.null(start.imp)) {
    for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
  } 
  .Call("MCMCjomo1con", Y, Yimp, Yimp2, X,betait,betapost,covit, omegapost, nburn, l1cov.prior,out.iter, PACKAGE = "jomo")
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  Yimp=Yimp2
  betapostmean<-apply(betapost, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(betapostmean)
    cat("The posterior covariance matrix is:\n")
    print(omegapostmean)
  }
  imp<-data.frame(imp)
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  colnames(imp)<-c(colnamy,colnamx,"Imputation","id")
  dimnames(betapost)[1] <- list(colnamx)
  dimnames(betapost)[2] <- list(colnamy)
  dimnames(omegapost)[1] <- list(colnamy)
  dimnames(omegapost)[2] <- list(colnamy)
  return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost))
}