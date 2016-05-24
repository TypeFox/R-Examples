jomo1rancat.MCMCchain <-
  function(Y.cat, Y.numcat, X=NULL, Z=NULL, clus, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, start.imp=NULL, nburn=100, output=1, out.iter=10) {
    if (is.null(X)) X=matrix(1,nrow(Y.cat),1)
    if (is.null(Z)) Z=matrix(1,nrow(Y.cat),1)
    if (is.null(beta.start)) beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat))))
    if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
    if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
    clus<-factor(unlist(clus))
    previous_levels_clus<-levels(clus)
    levels(clus)<-0:(nlevels(clus)-1)
    if (is.null(u.start)) u.start = matrix(0, nlevels(clus), ncol(Z) * ((sum(Y.numcat) - length(Y.numcat))))
    if (is.null(l2cov.start)) l2cov.start = diag(1, ncol(u.start))
    if (is.null(l2cov.prior)) l2cov.prior = diag(1, ncol(l2cov.start))
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
    for (i in 1:ncol(Z)) {
      if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
    }
    stopifnot( nrow(beta.start)==ncol(X), ncol(beta.start)==((sum(Y.numcat)-length(Y.numcat))),nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(beta.start), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start),nrow(Z)==nrow(Y.cat), ncol(l2cov.start)==ncol(u.start), ncol(u.start)==ncol(Z)*((sum(Y.numcat)-length(Y.numcat))))
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
    nimp=1
    colnamycat<-colnames(Y.cat)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    Y.cat<-data.matrix(Y.cat)
    storage.mode(Y.cat) <- "numeric"    
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"    
    Z<-data.matrix(Z)
    storage.mode(Z) <- "numeric"    
    clus <- matrix(as.integer(levels(clus))[clus], ncol=1)
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
    if (output!=1) out.iter=nburn+2
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+ncol(Z)+3)
    imp[1:nrow(Y),1:ncol(Y)]=Y
    imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[1:nrow(Z), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[1:nrow(clus), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[1:nrow(X), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=1  
    betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),nburn))
    omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),nburn))
    upostall<-array(0, dim=c(nrow(u.start),ncol(u.start),nburn))
    covupost<- array(0, dim=c(nrow(l2cov.start),ncol(l2cov.start),nburn))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    if (!is.null(start.imp)) {
      start.imp<-as.matrix(start.imp)
      if ((nrow(start.imp)!=nrow(Yimp2))||(ncol(Yimp2)!=ncol(start.imp))) {
        cat("start.imp dimensions incorrect. Not using start.imp as starting value for the imputed dataset.\n")
        start.imp=NULL
      } else {
        Yimp2<-start.imp
      }
    }
    if (is.null(start.imp)) {
      for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    } 
    .Call("MCMCjomo1ranmix", Y, Yimp, Yimp2, Y.cat, X, Z, clus,betait,uit,betapost,upostall,covit,omegapost, covuit, covupost, nburn, l1cov.prior,l2cov.prior,Y.numcat, 0,out.iter, PACKAGE = "jomo")
    
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Y.cat
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
    for (i in 1:ncol(Y)) {
      imp[,i]<-as.factor(imp[,i]) 
      levels(imp[,i])<-previous_levels[[i]]
    }
    imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)]<-factor(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])
    levels(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])<-previous_levels_clus
    clus<-factor(clus)
    levels(clus)<-previous_levels_clus
    for (j in (ncol(Y.cat)+1):(ncol(Y.cat)+ncol(X)+ncol(Z))) {
      imp[,j]=as.numeric(imp[,j])
    }
    if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y.cat), sep = "")
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycat,colnamx,colnamz,"clus","id","Imputation")
    cnycatcomp<-rep(NA,(sum(Y.numcat)-length(Y.numcat)))
    count=0
    for ( j in 1:ncol(Y.cat)) {
      for (k in 1:(Y.numcat[j]-1)) {
        cnycatcomp[count+k]<-paste(colnamycat[j],k,sep=".")
      }
      count=count+Y.numcat[j]-1
    }
    dimnames(betapost)[1] <- list(colnamx)
    dimnames(betapost)[2] <- list(cnycatcomp)
    dimnames(omegapost)[1] <- list(cnycatcomp)
    dimnames(omegapost)[2] <- list(cnycatcomp)
    colnamcovu<-paste(cnycatcomp,rep(colnamz,each=ncol(omegapost)),sep="*")
    dimnames(covupost)[1] <- list(colnamcovu)
    dimnames(covupost)[2] <- list(colnamcovu)
    dimnames(upostall)[1]<-list(levels(clus))
    dimnames(upostall)[2]<-list(colnamcovu)
    dimnames(Yimp2)[2] <- list(cnycatcomp)
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost, "finimp.latnorm" = Yimp2))
  }
