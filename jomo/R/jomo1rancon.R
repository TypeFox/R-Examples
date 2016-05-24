jomo1rancon<- function(Y, X=NULL, Z=NULL, clus, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=100, nbetween=100, nimp=5, output=1, out.iter=10) {
  if (nimp<2) {
    nimp=2
    cat("Minimum number of imputations:2. For single imputation using function jomo1rancon.MCMCchain\n")
  }
  if (is.null(X)) X=matrix(1,nrow(Y),1)
  if (is.null(Z)) Z=matrix(1,nrow(Y),1)
  if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
  if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
  if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
  clus<-factor(unlist(clus))
  previous_levels_clus<-levels(clus)
  levels(clus)<-0:(nlevels(clus)-1)
  if (is.null(u.start)) u.start = matrix(0, nlevels(clus), ncol(Z)*ncol(Y))
  if (is.null(l2cov.start)) l2cov.start = diag(1, ncol(u.start))
  if (is.null(l2cov.prior)) l2cov.prior = diag(1, ncol(l2cov.start))
  for (i in 1:ncol(X)) {
    if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
  }
  for (i in 1:ncol(Z)) {
    if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
  }
  stopifnot(nrow(Y)==nrow(clus),nrow(Y)==nrow(X), nrow(beta.start)==ncol(X), ncol(beta.start)==ncol(Y),nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(Y), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start), nrow(Z)==nrow(Y), ncol(l2cov.start)==ncol(u.start), ncol(u.start)==ncol(Z)*ncol(Y))
  betait=matrix(0,nrow(beta.start),ncol(beta.start))
  for (i in 1:nrow(beta.start)) {
    for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
  }
  covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  for (i in 1:nrow(l1cov.start)) {
    for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
  }   
  uit=matrix(0,nrow(u.start),ncol(u.start))
  for (i in 1:nrow(u.start)) {
    for (j in 1:ncol(u.start)) uit[i,j]=u.start[i,j]
  }
  covuit=matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  for (i in 1:nrow(l2cov.start)) {
    for (j in 1:ncol(l2cov.start)) covuit[i,j]=l2cov.start[i,j]
  }   
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  colnamz<-colnames(Z)
  Y<-data.matrix(Y)
  storage.mode(Y) <- "numeric"    
  X<-data.matrix(X)
  storage.mode(X) <- "numeric"    
  Z<-data.matrix(Z)
  storage.mode(Z) <- "numeric"    
  clus <- matrix(as.integer(levels(clus))[clus], ncol=1)
  if (output!=1) out.iter=nburn+nbetween
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+ncol(Z)+3)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(Z), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[1:nrow(clus), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[1:nrow(X), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(0, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=1
  betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nimp-1)))
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  upost<-matrix(0,nrow(u.start),ncol(u.start))
  upostall<-array(0, dim=c(nrow(u.start),ncol(u.start),(nimp-1)))
  omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nimp-1)))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  covupost<- array(0, dim=c(nrow(l2cov.start),ncol(l2cov.start),(nimp-1)))
  cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  meanobs<-colMeans(Y,na.rm=TRUE)
  for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
  #for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=rnorm(1,mean=meanobs[j], sd=0.01)
  .Call("jomo1rancon", Y, Yimp, Yimp2, X, Z, clus, betait, uit, bpost, upost, covit, opost, covuit, cpost, nburn, l1cov.prior, l2cov.prior,out.iter, PACKAGE = "jomo")
  #betapost[,,1]=bpost
  #upostall[,,1]=upost
  #omegapost[,,1]=opost
  #covupost[,,1]=cpost
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  upost<-matrix(0,nrow(u.start),ncol(u.start))
  cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  Yimp=Yimp2
  if (output==1) cat("First imputation registered.", "\n")
  for (i in 2:nimp) {
    Yimp2=matrix(0, nrow(Y),ncol(Y))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[(i*nrow(clus)+1):((i+1)*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=i
    .Call("jomo1rancon", Y, Yimp, Yimp2, X, Z, clus, betait, uit, bpost, upost, covit,opost, covuit,cpost, nbetween, l1cov.prior, l2cov.prior,out.iter, PACKAGE = "jomo")
    betapost[,,(i-1)]=bpost
    upostall[,,(i-1)]=upost
    omegapost[,,(i-1)]=opost
    covupost[,,(i-1)]=cpost
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    upost<-matrix(0,nrow(u.start),ncol(u.start))
    cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
    imp[(i*nrow(Y)+1):((i+1)*nrow(Y)),1:ncol(Y)]=Yimp2
    Yimp=Yimp2
    if (output==1) cat("Imputation number ", i, "registered", "\n")
  }
  betapostmean<-apply(betapost, c(1,2), mean)
  upostmean<-apply(upostall, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  covupostmean<-apply(covupost, c(1,2), mean)
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(betapostmean)
    cat("The posterior mean of the random effects estimates is:\n")
    print(upostmean)
    cat("The posterior mean of the level 1 covariance matrices is:\n")
    print(omegapostmean)
    cat("The posterior mean of the level 2 covariance matrix is:\n")
    print(covupostmean)
  }
  imp<-data.frame(imp)
  imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)]<-factor(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])
  levels(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])<-previous_levels_clus
  clus<-factor(clus)
  levels(clus)<-previous_levels_clus
  for (j in 1:(ncol(Y)+ncol(X)+ncol(Z))) {
    imp[,j]=as.numeric(imp[,j])
  }
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  colnames(imp)<-c(colnamy,colnamx,colnamz,"clus","id","Imputation")
  return(imp)
}