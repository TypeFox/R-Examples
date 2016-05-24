### kClass: Generates k-class estimators (point est, se, p-value, and confint)
###         Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        k, a vector (or scalar) of k values for estimation
###        beta0, a vector (or scalar) of null values for testing and p-values
###        alpha, significance level for confidence intervals
###        heteroSE, use heteroscedastic robust standard errors?
###        clusterID, if cluster-robust standard errors are desired, please provide a vector of
###                   length that's identical to the sample size.
###                   For example, if n = 6 and clusterID = c(1,1,1,2,2,2),
###                   there would be two clusters where the first cluster
###                   is formed by the first three observations and the
###                   second cluster is formed by the last three observations
###                   clusterID can be numeric, character, or factor.
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
KClass = function(ivmodel,
                  beta0=0,alpha=0.05,k=c(0,1),
                  heteroSE=FALSE,clusterID=NULL) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("You must supply an ivmodel class. Run ivmodel() and see ivmodel() function for details")
    return(NULL)
  }
  if(missing(k)) {
    print("You must specify a value for k")
    return(NULL)
  }
 
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj =ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR
  degF = ivmodel$n - ivmodel$p - 1
  
  # Compute k-class estimator
  denom = (sum(Dadj^2) - k * sum(qr.resid(ZadjQR,Dadj)^2))
  kPointEst = (sum(Yadj * Dadj) - k * sum(Yadj * qr.resid(ZadjQR,Dadj))) / denom 
  
  # Compute std error, testStat, and confidence intervals
  kVarPointEst = rep(0,length(k))
  kCI = matrix(0,length(k),2)
  for(i in 1:length(k)) {
    if(denom[i] <= 0) { #k exceeds the maximum possible value
	  kVarPointEst[i] = NA; kCI[i,] = c(NA,NA)
	} else {
      if(heteroSE) {
	    kVarPointEst[i] = ivmodel$n/degF * sum( (Yadj - Dadj * kPointEst[i])^2 * (Dadj - k[i]*qr.resid(ZadjQR,Dadj))^2 ) / (denom[i])^2
	  } else if(!is.null(clusterID)){
	  	if(length(clusterID) != ivmodel$n) {
	  		### Missing problem here needs to be taken care of ###
	  		print("Cluster ID vector is not the same length as the sample size")
	  		return(NULL)
	  	}
	  	if(!is.character(clusterID) && !is.factor(clusterID) && !is.numeric(clusterID)) {
	  		print("Cluster ID must be either a character vector, a factor vector, or a numeric vector")
	  		return(NULL)
	  	}
	  	clusterID = as.factor(clusterID)
	  	clusterID = as.numeric(clusterID)
	  	nCluster <- length(unique(clusterID))
        dfc <- nCluster/(nCluster-1)*(ivmodel$n-1)/degF
        residFit = Yadj - Dadj * kPointEst[i]
        umat = tapply(residFit*(Dadj - k[i]*qr.resid(ZadjQR,Dadj)),clusterID, sum)
        kVarPointEst[i] = dfc * t(umat) %*% umat / denom[i]^2
	  }
	  else {
        kVarPointEst[i] = 1/(degF) *sum((Yadj - Dadj * kPointEst[i])^2) / denom[i]
	  }
	  # Note that R's AER package uses the z-scores instead of the t distribution. They are basically the same as n gets large.
	  kCI[i,] =  c(kPointEst[i] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]),
	               kPointEst[i] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]))
    }
  }
  
  # Compute test statistics
  # This operation takes a vector of k estimates (pointEst, k*1 dim) and 
  # another vector of null values (beta0,1 * length(beta0))
  # and subtracts pointEst from a vector of null values.
  # The last operation of dividing by sqrt() is done column-wise 
  # E.g. [1,2,3 ; 2,3,4] / (2,3,4) --> [0.5,2/3,0.75;1,1,1] 
  kTestStat = outer(kPointEst,beta0,FUN="-") / sqrt(kVarPointEst)
  
  # Compute p-value
  kPValue = 2*(1 - pt(abs(kTestStat),degF))
  
  # Package output 
  kPointEst = matrix(kPointEst,length(k),1)
  kVarPointEst = matrix(sqrt(kVarPointEst),length(k),1)
  rownames(kPointEst) = rownames(kVarPointEst) = rownames(kCI) = rownames(kPValue) = rownames(kTestStat) = k
  colnames(kPValue) = colnames(kTestStat) = beta0
  colnames(kPointEst) = "Estimate"
  colnames(kVarPointEst) = "Std. Error"
  colnames(kCI) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(point.est = kPointEst,std.err = kVarPointEst,test.stat = kTestStat,p.value = kPValue,ci = kCI))
}