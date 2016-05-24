call.netGSA <-
function(
  D1,  		            #Influence matrix for condition 1
  D2,		        	#Influence matrix for condition 2
  x, 			      	#the p x n data matrix
  y, 			      	#vector of class indicators of length n
  B, 		    	  	#indicator matix for pathways (npath x p)
  varEstCntrl           #parameters to pass for variance estimation
){
  
  p = dim(x)[1] #No. of genes
  n = length(y) #No. of samples in total
  npath = dim(B)[1] #No. of gene sets to consider
  
  y = as.integer(as.factor(y))
  n1 = sum(y == 1)
  n2 = sum(y == 2)
  
  X1 = x[, (y == 1)]
  X2 = x[, (y == 2)]
  
  ##Matrices Needed 
  Ip = diag(rep(1, p))
  D1Inv = solve(D1)
  D2Inv = solve(D2)
  
  D1D1 = D1 %*% t(D1)
  D2D2 = D2 %*% t(D2)
  
  D1D1Inv = solve(t(D1) %*% D1)
  D2D2Inv = solve(t(D2) %*% D2)
  
  ##Building the "contrast" matrix L, see Result in the paper
  L1 = (B %*% D1) * B
  L2 = (B %*% D2) * B
  LN = cbind(-L1, L2)
  
  ##Initialzing test matrices, each test is explained where it's calculated
  teststat = matrix(0, npath, 1)
  num.tstat = matrix(0, npath, 1)
  
  ##Initializing the vector for degrees of freedom for the test statistics.
  df = matrix(0, npath, 1)
  
  #Initialzing pvalues of each test
  pvals = matrix(0, npath, 1)
  
  ##------------------
  ##ESTIMATION OF BETA
  ##------------------
  beta1 = (1/n1) * D1Inv %*% apply(X1, 1, sum)
  beta2 = (1/n2) * D2Inv %*% apply(X2, 1, sum)
  
  ##-----------------
  ##ESTIMATION OF VAR COMPONENTS
  ##-----------------
  #First calculate residuals: 
  resid = matrix(0, p, n)
  for (i in 1:n1) {
    resid[, i] = X1[, i] - D1 %*% beta1
  }
  for (i in 1:n2) {
    resid[, n1 + i] = X2[, i] - D2 %*% beta2
  }
  
  #initial estimates
  s2g = (mean(apply(resid, 2, sd)))^2
  s2e = (sd(diag(resid)))^2
  
  ##Estimation of variance components is based on profile likelihood with Newton's method.
  ##It can be done using either ML or REML.
  ##However, the method uses the approximate method for datasets with
  ##more than 5K genes to prevent computational issues and memory overflow.    
  approxEst = approxVarEst(s2e, s2g, D1, D2, resid, n1, n2, control = varEstCntrl) 
  se0 = approxEst$s2e
  sg0 = approxEst$s2g
  if (p<5000){
	  S = profileVarEst(se0, sg0, D1, D2, resid, n1, n2, control = varEstCntrl)
	  s2e = S$s2e
	  s2g = S$s2g  	
  } else {
	  s2e = se0
	  s2g = sg0  	
  } 
 
  ##-----------------
  ##CALCULATING DEGREES OF FREEDOM & TEST STATS
  ##-----------------
  #matrices needed in calculatoin of degrees of freedom
  mc = chol2inv(chol(s2g * D1D1 + s2e * Ip))
  mca = mc %*% D1D1
  mt = chol2inv(chol(s2g * D2D2 + s2e * Ip))
  mta = mt %*% D2D2
  
  #trace of matrices needed for H (and therefore KK)
  t11c = sum(diag(mca %*% mca))
  t11t = sum(diag(mta %*% mta))
  t22c = sum(diag(mc %*% mc))
  t22t = sum(diag(mt %*% mt))
  t12c = sum(diag(mca %*% mc))
  t12t = sum(diag(mta %*% mt))
  
  #These are the elements of the empirical Information matrix
  EH11 = (1/2) * (t11c + t11t)
  EH22 = (1/2) * (t22c + t22t)
  EH12 = (1/2) * (t12c + t12t) #(-1/2) * (tt1 + tt2) + (1/2) * (tt3 + tt4)
  
  #In this version of the code, K matrix is calculated directly!
  ## Kmat will be the expected information matrix. Here we need its inverse 
  ## in calculating the degrees of freedom.
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  KmatInv = solve(Kmat)
  
  #These matrices are needed in the calculation of test statistics
  mctildi = s2e * D1D1Inv + s2g * Ip
  mttildi = s2e * D2D2Inv + s2g * Ip
  
  ##-----------------
  for (rr in 1:npath) {
    Lrow = t(as.matrix(LN[rr, ])) #single row of L3
    
    Lrow1 = t(as.matrix(Lrow[, 1:p]))
    Lrow2 = t(as.matrix(Lrow[, (p + 1):(2 * p)]))
    
    LC11Lprime = (1/n2) * Lrow2 %*% mttildi %*% t(Lrow2) + (1/n1) * Lrow1 %*% mctildi %*% t(Lrow1)
    
    g1 = (1/n2) * (Lrow2 %*% t(Lrow2)) + (1/n1) * (Lrow1 %*% t(Lrow1))
    g2 = (1/n2) * Lrow2 %*% D2D2Inv %*% t(Lrow2) + (1/n1) * Lrow1 %*% D1D1Inv %*% t(Lrow1)
    g = matrix(c(g1, g2), 2, 1)
    
    #test statistic
    num.tstat[rr] = Lrow2 %*% beta2 + Lrow1 %*% beta1
    teststat[rr] = num.tstat[rr]/sqrt(LC11Lprime)	
    
    #calculating df based on the Satterthwaite approximation method 
    #using the formula nu=(2*LCL'^2)/g'Kg with K being the empirical covariance matrix.
    #NOTE: If df2<2, it is set to 2
    
    df[rr] = 2 * (LC11Lprime)^2/(t(g) %*% KmatInv %*% g)    
    if (df[rr] < 2) 
      df[rr] = 2
    
  }
  ##-----------------
  
  #p-values for each test
  pvals = 1 - pt(abs(teststat), df) + pt(-abs(teststat), df)
  
  output = list(beta = c(beta1, beta2), teststat = teststat, df = df, p.value = pvals, 
                s2.epsilon = s2e, s2.gamma = s2g)
  
  return(output)
}
