jomo1ran <-
  function(Y, X=NULL, Z=NULL,clus, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=100, nbetween=100, nimp=5, a=NULL, meth="common",output=1, out.iter=10) {
    stopifnot(meth=="common"|meth=="fixed"|meth=="random")
    ncon=0
    ncat=0
    Y.con=NULL
    Y.cat=NULL
    Y.numcat=NULL
    for (i in 1:ncol(Y)) {
      if (is.numeric(Y[,i])) {
        ncon=ncon+1
        if (is.null(Y.con)) {
          Y.con<-data.frame(Y[,i])
        } else {
          Y.con<-data.frame(Y.con,Y[,i])
        }
        colnames(Y.con)[ncon]<-colnames(Y)[i]
      }
      else {
        if (is.factor(Y[,i])) {
        ncat=ncat+1
        if (is.null(Y.cat)) {
          Y.cat<-data.frame(Y[,i])
        } else {
          Y.cat<-data.frame(Y.cat,Y[,i])
        }
        colnames(Y.cat)[ncat]<-colnames(Y)[i]
        Y.numcat<-cbind(Y.numcat,nlevels(Y[,i]))
        }
      }
    }
    if (is.null(X)) X=matrix(1,nrow(Y),1)
    if (is.null(Z)) Z=matrix(1,nrow(Y),1)
    if (meth=="common") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1con.", "\n")
        imp<-jomo1rancon(Y=Y.con, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancat.", "\n")
        imp<-jomo1rancat(Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmix.", "\n")
        imp<-jomo1ranmix(Y.con, Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output,out.iter=out.iter)
      }
    }
    if (meth=="fixed") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1ranconhr(Y=Y.con, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=0, meth="fixed", output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1rancathr(Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=0, meth="fixed",output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1ranmixhr(Y.con, Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=0, meth="fixed", output=output,out.iter=out.iter)
      }
    }
    if (meth=="random") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1ranconhr(Y=Y.con, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth="random", output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1rancathr(Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth="random",output=output,out.iter=out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        imp<-jomo1ranmixhr(Y.con, Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, Z=Z, clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth="random", output=output,out.iter=out.iter)
      }
    }
    return(imp)
  }
