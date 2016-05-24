PenalizedLMM_opt <-
function(y, X, subject, lambda_pop, lambda_random, max_iter=100, tol = 1e-5){
  # estimate penalised linear mixed model - has been optimized (very slightly)
  #
  # INPUT:
  #      - y: response vector (centered to have mean zero)
  #      - X: fixed effect design matrix (for all random effects - ie should be rbind of all design matrices)
  #      - lambda_pop: population (fixed) effects regularisation parameters
  #      - lambda_random: regularisation parameter for variance of random effects
  #      - max_iter, tol: max number of iterations & convergence tolerance
  #
  #
  #
  
  # prepare some variables:
  p = ncol(X) # number of predictors
  n = length(y) # number of observations
  N = length(unique(subject)) # number of subjects
  nsub = n/N # number of observations per subject (assumed equal across all subjects)
  beta = ginv(t(X)%*%X)%*%t(X)%*%y # initial estimate fixed effects - will be p dimensional vector
  betaOld = beta # to determine convergence
  sig = rep(1, p) # initial estimate variances of random effects - also a p dimensional vector
  sigma = sqrt( sum((y - X%*%beta)**2)/n ) # initial estimate of variance
  D = diag(sig)
  NLL = rep(0, max_iter+1) # to track NLL over iterations
  
  # build D matrix:
  Dtilde = rep(1, p) # initial random effects estimates
  # initial estimate of random effects. Each row a random effects vector for the corresponding subject. bvec is the vectorized version
  bvec = rep(0, p*N)
  for (i in 1:N){
    ID = seq(1+(i-1)*p,(i*p))
    ID2 = seq(1+(i-1)*nsub,(i*nsub))
    bvec[ID] = c(solve(t(diag(Dtilde))%*% t(X[ID2,]) %*% X[ID2,] %*% diag(Dtilde) + diag(p)) %*% t(diag(Dtilde))%*%t(X[ID2,])%*%(y[ID2] -  X[ID2,] %*%beta))
  }
  b = matrix(bvec, ncol=p, byrow = TRUE)
  
  # prepare penalty coefficients
  penW = c(rep(1, p), (lambda_random/lambda_pop)*rep(1,p))
  
  # some additional definitions:
  iter = 0 # iteration count
  conv = FALSE # convergence track
  PhiOld = c(beta, Dtilde, sigma) # used to determine convergence
  
  ## begin EM algorithm:
  while((iter < max_iter)&(!conv)){
    
    ### perform M step:
    fixedDesign = X # this is very inefficient, will tidy up later - keeping it this way to help with interpretability for now
    randomDesign = BuildRandomDesign(fixedDes=X, blups=bvec, N=N)
    #randomDesignOld = Z %*% diag(bvec) %*% kronecker(rep(1,N), diag(p))
    
    # solve our lasso problem:
    coefs = glmnet(x=cbind(fixedDesign, randomDesign), y=y, family="gaussian",  alpha=1, lambda=lambda_pop*as.numeric(sigma**2), penalty.factor=penW, intercept = FALSE, lower.limits = c(rep(-Inf, p), rep(0, p)), standardize = FALSE)
    # note that we have ensured that the estimated variances are positive!
    beta = as.numeric(coefs$beta[1:p])
    #sig = as.numeric(abs(coefs$beta[(p+1):(2*p)]))
    Dtilde =  as.numeric((coefs$beta[(p+1):(2*p)]))
    
    # update white noise estiamte:
    sigma = 0 # we update in loop below #sqrt(t(y-X%*%beta) %*% solve(Z%*%Dtilde %*% Dtilde %*% t(Z)+diag(n)) %*% (y-X%*%beta)/n)
    #Dtilde = sig # this is a bit silly but leave it for now
    
    ### perform E step:
    #bvec = c(solve(t(Dtilde)%*% t(Z) %*% Z %*% Dtilde + diag(N*p)) %*% t(Dtilde)%*%t(Z)%*%(y -  X %*%beta))
    for (i in 1:N){
      ID = seq(1+(i-1)*p,(i*p))
      ID2 = seq(1+(i-1)*nsub,(i*nsub))
      bvec[ID] = c(solve(t(diag(Dtilde))%*% t(X[ID2,]) %*% X[ID2,] %*% diag(Dtilde) + diag(p)) %*% t(diag(Dtilde))%*%t(X[ID2,])%*%(y[ID2] -  X[ID2,] %*%beta))
      sigma = sigma + t(y[ID2]- X[ID2, ]%*%beta) %*% solve(X[ID2,]%*%diag(Dtilde) %*% diag(Dtilde) %*% t(X[ID2,])+diag(nsub)) %*% (y[ID2]-X[ID2,]%*%beta)
      ####sigma  = sigma + t(y[ID2]- X[ID2, ]%*%beta - X[ID2,] %*% bvec[ID]) %*% (y[ID2]- X[ID2, ]%*%beta - X[ID2,] %*% bvec[ID]) + sum(diag(t(X[ID2,]) %*% X[ID2,] %*% ( as.numeric(sigma**(-2)) * t(X[ID2,]) %*% X[ID2,] + diag())))
    }
    b = matrix(bvec, ncol=p, byrow = TRUE)
    sigma = sqrt(sigma/n)
    
    # get new NLL - should be decreasing!
    for (i in 1:N){
      ID = seq(1+(i-1)*p,(i*p))
      ID2 = seq(1+(i-1)*nsub,(i*nsub))
      NLL[(iter+1)] = NLL[(iter+1)] + calculateNLL(resp = y[ID2], des = X[ID2, ], fixEf = beta, reVar = Dtilde, resVar = as.numeric(sigma)) #+ lambda_pop*sum(abs(beta)) + lambda_random * sum(abs(Dtilde))
      #cat(NLL[(iter+1)], "\n")
    }
    #NLL[(iter+1)] = -1*NLL[(iter+1)] + lambda_pop*sum(abs(beta)) + lambda_random * sum(abs(Dtilde))
    
    # define vector of parameters to be estimated:
    Phi = c(beta, Dtilde, sigma)
    Res = sum(sapply(Phi-PhiOld, FUN=function(x){x**2}))
    Denom = sum(sapply(PhiOld, FUN=function(x){x**2}))
    # check convergence:
    if (sum(abs(Res/Denom))< tol){
      conv = TRUE
    } else {
      iter = iter + 1
      betaOld = beta
      PhiOld = Phi
    }
  }
  return(list(beta=beta, sigmaRE=Dtilde, sigmaNoise=sigma, BLUPs=b, it=iter, NLL=NLL[1:(iter+1)]))
}
