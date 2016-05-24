jomo1cat <-
  function(Y.cat, Y.numcat, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, nburn=100, nbetween=100, nimp=5, output=1, out.iter=10) {
    if (nimp<2) {
      nimp=2
      cat("Minimum number of imputations: 2. For single imputation using function jomo1cat.MCMCchain\n")
    }
    if (is.null(X)) X=matrix(1,nrow(Y.cat),1)
    if (is.null(beta.start)) beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat))))
    if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
    if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
    previous_levels<-list()
    Y.cat<-data.frame(Y.cat)
    for (i in 1:ncol(Y.cat)) {
      Y.cat[,i]<-factor(Y.cat[,i])
      previous_levels[[i]]<-levels(Y.cat[,i])
      levels(Y.cat[,i])<-1:nlevels(Y.cat[,i])
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    stopifnot( nrow(beta.start)==ncol(X), ncol(beta.start)==((sum(Y.numcat)-length(Y.numcat))),nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(beta.start), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start))
    betait=matrix(0,nrow(beta.start),ncol(beta.start))
    for (i in 1:nrow(beta.start)) {
      for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
    }
    covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    for (i in 1:nrow(l1cov.start)) {
      for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
    }    
    colnamycat<-colnames(Y.cat)
    colnamx<-colnames(X)
    Y.cat<-data.matrix(Y.cat)
    storage.mode(Y.cat) <- "numeric"    
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"
    Y=cbind(Y.cat)
    Yi=cbind(matrix(0,nrow(Y.cat),(sum(Y.numcat)-length(Y.numcat))))
    h=1
    for (i in 1:length(Y.numcat)) {
      for (j in 1:nrow(Y)) {
        if (is.na(Y.cat[j,i])) {
          Yi[j,h:(h+Y.numcat[i]-2)]=NA
        }
      } 
      h=h+Y.numcat[i]-1
    }
    if (output!=1) out.iter=nburn+nbetween
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
    imp[1:nrow(Y),1:ncol(Y)]=Y
    imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[1:nrow(X), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=1
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nimp-1)))
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nimp-1)))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    .Call("jomo1mix", Y, Yimp, Yimp2, Y.cat, X,betait,bpost,covit,opost, nburn, l1cov.prior,Y.numcat, 0, out.iter, PACKAGE = "jomo")
    #betapost[,,1]=bpost
    #omegapost[,,1]=opost
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Y.cat
    if (output==1) cat("First imputation registered.", "\n")
    for (i in 2:nimp) {
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+1)]=i
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
      .Call("jomo1mix", Y, Yimp, Yimp2, Y.cat, X,betait,bpost,covit, opost, nbetween, l1cov.prior, Y.numcat, 0, out.iter, PACKAGE = "jomo") 
      
      betapost[,,(i-1)]=bpost
      omegapost[,,(i-1)]=opost
      bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
      opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))      
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y)]=Y.cat
      if (output==1) cat("Imputation number ", i, "registered", "\n")
    }

    betapostmean<-apply(betapost, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      cat("The posterior covariance matrix is:\n")
      print(omegapostmean)
    }
    imp<-data.frame(imp)
    for (i in 1:ncol(Y)) {
        imp[,i]<-as.factor(imp[,i]) 
        levels(imp[,i])<-previous_levels[[i]]
    }
    if (is.null(colnamycat)) colnamycat=paste("Y", 1:ncol(Y.cat), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycat,colnamx,"Imputation","id")
    return(imp)
  }
