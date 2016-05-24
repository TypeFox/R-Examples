jomo2 <-
  function(Y, Y2, X=NULL, X2=NULL, Z=NULL,clus, beta.start=NULL,l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=100, nbetween=100, nimp=5, a=NULL, meth="common",output=1, out.iter=10) {
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
    ncon2=0
    ncat2=0
    Y2.con=NULL
    Y2.cat=NULL
    Y2.numcat=NULL
    for (i in 1:ncol(Y2)) {
      if (is.numeric(Y2[,i])) {
        ncon2=ncon2+1
        if (is.null(Y2.con)) {
          Y2.con<-data.frame(Y2[,i])
        } else {
          Y2.con<-data.frame(Y2.con,Y2[,i])
        }
        colnames(Y2.con)[ncon2]<-colnames(Y2)[i]
      }
      else {
        if (is.factor(Y2[,i])) {
          ncat2=ncat2+1
          if (is.null(Y2.cat)) {
            Y2.cat<-data.frame(Y2[,i])
          } else {
            Y2.cat<-data.frame(Y2.cat,Y2[,i])
          }
          colnames(Y2.cat)[ncat2]<-colnames(Y2)[i]
          Y2.numcat<-cbind(Y2.numcat,nlevels(Y2[,i]))
        }
      }
    }
    if (is.null(X2)) X2=matrix(1,nrow(Y2),1)
    if (meth=="common") {
        cat("Found ", ncon, "level 1 continuous and ", ncat, "level 1 categorical outcomes, ",ncon2," level 2 continuous and ",ncat2," level 2 categorical outcomes. Using function jomo2com, assuming common covariance matrix across clusters", "\n")
        imp<-jomo2com(Y.con=Y.con, Y.cat=Y.cat, Y.numcat=Y.numcat, Y2.con=Y2.con, Y2.cat=Y2.cat, Y2.numcat=Y2.numcat, X=X, X2=X2, Z=Z, clus=clus, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter)
    }
    if (meth=="fixed") {
      cat("Found ", ncon, "level 1 continuous and ", ncat, "level 1 categorical outcomes, ",ncon2," level 2 continuous and ",ncat2," level 2 categorical outcomes. Using function jomo2hr with fixed cluster-specific covariance matrices.", "\n")
      imp<-jomo2hr(Y.con=Y.con, Y.cat=Y.cat, Y.numcat=Y.numcat, Y2.con=Y2.con, Y2.cat=Y2.cat, Y2.numcat=Y2.numcat, X=X, X2=X2, Z=Z, clus=clus, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp,a=a, meth="fixed", output=output, out.iter=out.iter)
    }
    if (meth=="random") {
      cat("Found ", ncon, "level 1 continuous and ", ncat, "level 1 categorical outcomes, ",ncon2," level 2 continuous and ",ncat2," level 2 categorical outcomes. Using function jomo2hr with random cluster-specific covariance matrices.", "\n")
      imp<-jomo2hr(Y.con=Y.con, Y.cat=Y.cat, Y.numcat=Y.numcat, Y2.con=Y2.con, Y2.cat=Y2.cat, Y2.numcat=Y2.numcat, X=X, X2=X2, Z=Z, clus=clus, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth="random", output=output, out.iter=out.iter)
    }
    return(imp)
  }
