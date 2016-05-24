jomo <-
  function(Y, Y2=NULL, X=NULL, X2=NULL, Z=NULL,clus=NULL, beta.start=NULL,l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=100, nbetween=100, nimp=5, a=NULL, meth="common",output=1, out.iter=10) {
    if (is.null(Y2)) {
      if (is.null(clus)) {
        cat("No clustering, using functions for single level imputation.\n")
        imp<-jomo1(Y=Y, X=X, beta.start=beta.start, l1cov.start=l1cov.start, l1cov.prior=l1cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, output=output, out.iter=out.iter)
      }
      if (!is.null(clus)) {
        cat("Clustered data, using functions for two-level imputation.\n")
        imp<-jomo1ran(Y=Y, X=X, Z=Z,clus=clus, beta.start=beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth=meth, output=output, out.iter=out.iter) 
      }
    } else {
      cat("2-level data, using functions for two-level imputation.\n")
      imp<-jomo2(Y=Y, Y2=Y2, X=X, X2=X2, Z=Z,clus=clus, beta.start=beta.start, l2.beta.start=l2.beta.start, u.start=u.start, l1cov.start=l1cov.start, l2cov.start=l2cov.start, l1cov.prior=l1cov.prior, l2cov.prior=l2cov.prior, nburn=nburn, nbetween=nbetween, nimp=nimp, a=a, meth=meth, output=output, out.iter=out.iter) 
    }
    return(imp)
  }
