# condPvalue: This is the conditional p-value integral described in Andrews, Moreira, and Stock (2007, Journal of Econometrics)
#             P(m,qT) = 1 - Pr(LR > m |QT = qT)
#          m: test stat value
#         qT:  conditional value of Qt
#          k: number of moment conditions (i.e. number of instruments)
#          
#        eps: epsilon (for numerical stability, default set based on Andrews, Moreira, and Stock's 2007 recommendation (0.02)
# OUTPUT: a p-value
# NOTE: mimics condivreg.ado's new_try function in Mikusheva and Poi (2006)'s STATA function
condPvalue <- function(m,qT,k,df2,eps = 0.02) {
  if(k == 1) { # Eqivalent to AR #
    return(1 - pf(m,k,df2))
  }
  K = gamma(k/2) / (sqrt(pi) * gamma((k-1)/2))  
  if(k == 2) {
    return(1 - 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * sin(x)^2/m),k)},lower=0,upper=pi/2)$value) 
  } else if(k == 3) {
    return(1 - 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k)},lower=0,upper=1)$value)
  } else if(k == 4) {
    nonapproxRegion = 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k) * (1 - x^2)^((k-3)/2)},lower=0,upper=1-eps)$value
    approxRegion =  2 * K * pchisq( (qT + m) / (1 + qT *(1 - eps/2)^2 /m),k) * (1/2 * (pi/2 - asin(1 - eps)) - (1 - eps)/2 * sqrt( 1- (1 - eps)^2)) 
    return(1 - nonapproxRegion - approxRegion)
  } else {
    return(1 - 2*K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k) * (1 - x^2)^((k-3)/2)},lower=0,upper=1)$value)
  }
}

### CLR: Generates the CLR confidnece interval
###         Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        alpha, significance level for confidence intervals
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
CLR = function(ivmodel,beta0=0,alpha=0.05) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("CLR: You must supply an ivmodel class. See ivmodel function for details")
	return(NULL)
  }
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; ZadjQR = ivmodel$ZadjQR
  Y = ivmodel$Y; D = ivmodel$D
  
  # Qs, Qst, Qst, Qt #
  YadjANDDadj = cbind(Yadj,Dadj) #You can also use the original Y and D
  PZYadjANDDadj = qr.fitted(ZadjQR,YadjANDDadj); 
  RZYadjANDDadj = qr.resid(ZadjQR,YadjANDDadj); 
  
  sigmaHat = (t(RZYadjANDDadj) %*% RZYadjANDDadj) / (ivmodel$n - ivmodel$p - ivmodel$L)
  sigmaHatInv = invTwobyTwoSymMatrix(sigmaHat) 
  
  a0 = c(beta0,1); b0 = c(1,-beta0)
  denomS = t(b0) %*% sigmaHat %*% b0; denomT = t(a0) %*% sigmaHatInv %*% a0
  numT = PZYadjANDDadj %*% sigmaHatInv 
  QS = sum( (PZYadjANDDadj %*% b0)^2) / denomS
  QT = sum( (PZYadjANDDadj %*% sigmaHatInv %*% a0)^2) / denomT
  QTS = sum( (PZYadjANDDadj %*% b0) * (PZYadjANDDadj %*% sigmaHatInv %*% a0)) / (sqrt(denomS) * sqrt(denomT))
  
  LRtest = 1/2 * (QS - QT + sqrt((QS + QT)^2 - 4*(QS *QT - QTS^2)))
  
  test.stat = matrix(LRtest,1,1)
  p.value = tryCatch({condPvalue(LRtest,QT,ivmodel$L,ivmodel$n - ivmodel$p- ivmodel$L)},error=function(e){0})
  p.value = matrix(p.value,1,1) 
  maxEigen = max(quadSolver(a=1,b=-1*(QS + QT),c=QS *QT - QTS^2)) #of Q matrix
  C = tryCatch({uniroot(function(y){condPvalue(m=maxEigen - y,qT = y,k = ivmodel$L,df2=ivmodel$n - ivmodel$p- ivmodel$L) - alpha},
  	            lower=10*.Machine$double.eps,upper=maxEigen - 5*.Machine$double.eps,maxiter=5000)$root},error=function(e){0})
  quadMatrix.CLR = sigmaHatInv %*% t(YadjANDDadj) %*% PZYadjANDDadj %*% sigmaHatInv  - C * sigmaHatInv
  ci.CLR = quadSolver(a=quadMatrix.CLR[1,1],b=2*quadMatrix.CLR[1,2],c=quadMatrix.CLR[2,2])
  if(quadMatrix.CLR[1,1] > 0) {
    if(is.na(ci.CLR[1])) {
      ci.CLR = matrix(c(-Inf,Inf),1,2)
	  info = c("Whole Real Line")
    } else {
	    ci.CLR = matrix(c(-Inf,ci.CLR[1],ci.CLR[2],Inf),2,2,byrow=TRUE)
		info = c(paste("(-Infinity, ",ci.CLR[1,2],"] union [",ci.CLR[2,1],", Infinity)",sep=""))
	  }
  } else {
    ci.CLR = matrix(c(ci.CLR[1],ci.CLR[2]),1,2)
	info = paste("[",ci.CLR[1,1],", ",ci.CLR[1,2],"]",sep="")
  }
  
  # Package output 
  colnames(p.value) = colnames(test.stat) = beta0
  colnames(ci.CLR) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(test.stat = test.stat,p.value = p.value,ci = ci.CLR,ci.info = info))
}
