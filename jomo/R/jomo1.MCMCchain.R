jomo1.MCMCchain <-
  function(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL,start.imp=NULL, nburn=100, output=1, out.iter=10) {
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
    if (ncat==0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1con.", "\n")
      imp<-jomo1con.MCMCchain(Y=Y.con, X=X, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, start.imp=start.imp,nburn=nburn, output=output,out.iter=out.iter)
      attr(imp, "function") = "jomo1con.MCMCchain"
    }
    if (ncat>0 & ncon==0) {
      cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1cat.", "\n")
      imp<-jomo1cat.MCMCchain(Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, start.imp=start.imp,nburn=nburn, output=output,out.iter=out.iter)
      attr(imp, "function") = "jomo1cat.MCMCchain"
    }
    if (ncat>0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1mix.", "\n")
      imp<-jomo1mix.MCMCchain(Y.con=Y.con,Y.cat=Y.cat,Y.numcat=Y.numcat, X=X, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, start.imp=start.imp,nburn=nburn, output=output,out.iter=out.iter)
      attr(imp, "function") = "jomo1mix.MCMCchain"
    }
    return(imp)
  }
