clusterVCV <-
function(data, fm, cluster1, cluster2=NULL) {
  
  # Calculation shared by covariance estimates
  est.fun = estfun(fm)
  inc.obs = complete.cases(data[,all.vars(formula(fm))])
  
  # Shared data for degrees-of-freedom corrections
  N  = dim(fm$model)[1]
  NROW = NROW(est.fun)
  K  = dim(vcov(fm))[1]
  
  # Calculate the covariance matrix estimate for the first cluster.
  cluster1 = data[inc.obs,cluster1]
  cov1 = covc(cluster1, estfun=est.fun, N1=N, K1=K, NROW1=NROW, fm1=fm)
  
  if(is.null(cluster2)) {
    # If only one cluster supplied, return single cluster
    # results
    return(cov1)
  } else {
    # Otherwise do the calculations for the second cluster
    # and the "intersection" cluster.
    cluster2 = data[inc.obs,cluster2]
    cluster12 = paste(cluster1,cluster2, sep="")
    
    # Calculate the covariance matrices for cluster2, the "intersection"
    # cluster, then then put all the pieces together.
    cov2   = covc(cluster2, estfun=est.fun, N1=N, K1=K, NROW1=NROW, fm1=fm)
    cov12  = covc(cluster12, estfun=est.fun, N1=N, K1=K, NROW1=NROW, fm1=fm)
    covMCL = (cov1 + cov2 - cov12)
    
    # Return the output of coeftest using two-way cluster-robust
    # standard errors.
    return(covMCL)
  }
}
