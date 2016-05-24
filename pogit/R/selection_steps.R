#### (invisible function)
#### sub-function "draw_indicators()" for select_logit() and select_pois()
#### last version: 2016/03/24
#### This function is not meant to be called directly by the user. 

#### Sample the vector of indicators for regression effects (for BVS)

# ... ARGUMENTS ...
# y             response vector of latent variables 
#               (=representation as Gaussian regression model in auxiliary variables)
# X             design matrix in the regression model 
# delta         vector of indicators for selection of regression effects
#               (with delta_j=1 if alpha_j is allocated to the slab)
# gamma         indicator variable for selection of random intercept variance parameter
#               (with gamma=1 if theta is allocated to the slab)
# omega         (mixture) weight of slab component for regression effects
# pi            (mixture) weight of slab component for random intercept parameter
# model         model (list) object
# prior         prior (list) object
# invA0         inverse prior variance of regression effects 


draw_indicators <- function(y, X, delta, gamma, omega, pi, model, prior, invA0){
  
  nDelta <- model$d - sum(model$deltafix)
  nGamma <- model$ri - sum(model$gammafix)
  
  ## update indicators for regression effects
  if (nDelta > 0){
    iDel <- which(model$deltafix==0)  
    ranOrdDelta <- sample(nDelta)
    pdelta <- matrix(NA, 1, model$d)
    pdelta[which(model$deltafix==1)] <- 1
    
    for (i in 1:nDelta){
      j         <- iDel[ranOrdDelta[i]]
      delta.new <- delta
      lp <- matrix(0, 2, 1)
      
      for (ii in 0:1){
        delta.new[j] <- ii
        lprior       <- ii*log(omega) + (1 - ii)*log(1 - omega) 
        if (model$ri == 0){
          llik <- lmarglik(y, X, delta.new, prior$a0, invA0)
        } else {
          llik <- lmarglik(y, X, c(delta.new, gamma), prior$a0, invA0)
        }
        lp[ii+1] <- llik + lprior
      }
      maxL  <- max(lp)
      l     <- exp(lp-maxL)
      lprob <- l/sum(l) 
      
      deltaj <- runif(1) > lprob[1]
      if (deltaj != delta[j]) delta[j] <- deltaj
      pdelta[j] <- lprob[2]
    }    
  } else {
    pdelta <- NULL
  }
  
  ## update indicator for random intercept variance parameter
  if (model$ri == 1){
    if (model$gammafix == 0){
      gamma.new <- gamma
      lpg <- matrix(0, 2, 1)
      
      for (jj in 0:1){
        gamma.new <- jj
        lprior    <- jj*log(pi) + (1 - jj)*log(1 - pi)
        llik      <- lmarglik(y, X, c(delta, gamma.new), prior$a0, invA0)
        lpg[jj+1] <- llik + lprior
      }
      maxL  <- max(lpg)
      l     <- exp(lpg - maxL)
      lprob <- l/sum(l)
      
      gammaj <- runif(1)>lprob[1]
      if (gammaj != gamma) gamma <- gammaj
      pgamma <- lprob[2]
    } else {
      pgamma <- 1
    }
  } else {
    pgamma <- NULL
  }
  
  return(list(deltanew = delta, pdeltanew = pdelta, gammanew = gamma, 
              pgammanew = pgamma))
}



draw_indicators_nb <- function(y, X, delta, omega, model, prior, invA0){
  
  nDelta <- model$d - sum(model$deltafix)
  
  ## update indicators for regression effects
  if (nDelta > 0){
    iDel <- which(model$deltafix == 0)  
    ranOrdDelta <- sample(nDelta)
    pdelta <- matrix(NA, 1, model$d)
    pdelta[which(model$deltafix==1)] <- 1
    
    for (i in 1:nDelta){
      j         <- iDel[ranOrdDelta[i]]
      delta.new <- delta
      lp <- matrix(0, 2, 1)
      
      for (ii in 0:1){
        delta.new[j] <- ii
        lprior       <- ii*log(omega) + (1 - ii)*log(1 - omega) 
        llik <- lmarglik(y, X, delta.new, prior$a0, invA0)
        lp[ii+1] <- llik + lprior
      }
      maxL  <- max(lp)
      l     <- exp(lp-maxL)
      lprob <- l/sum(l) 
      
      deltaj <- runif(1) > lprob[1]
      if (deltaj != delta[j]) delta[j] <- deltaj
      pdelta[j] <- lprob[2]
    }    
  } else {
    pdelta <- NULL
  }
  
  return(list(deltanew = delta, pdeltanew = pdelta))
}




#### (invisible function)
#### sub-function "lmarglik()" for select_logit() and select_poisson()
#### last version: 2016/03/24
#### This function is not meant to be called directly by the user. 

#### Compute marginal likelihood of a Gaussion regression model

#### ... ARGUMENTS ...  
# y   		response vector
# X      	design matrix
# ind			index vector for (selected) regression effects
# a0prior	prior mean of regression effects
# iA0			inverse prior variance of regression effects

lmarglik <- function(y, X, ind, a0prior, iA0){
  
  index <- c(1, which(ind==1) + 1)	     
  invA0 <- iA0[index, index, drop = FALSE] 
  a0    <- a0prior[index,]			
  Xsel  <- X[, index, drop = FALSE] 		 
  
  Apost <- solve(invA0 + t(Xsel)%*%Xsel)  # Apost=(A0^-1 + Z*'Sigma^-1Z*)	
  apost <- invA0%*%a0 + t(Xsel)%*%y       # apost=Apost*(A0^-1*a0 + Z*'Sigma^-1*y) 
  
  # conditional (log) marginal likelihood 
  #h <- log(det(Apost))-(-log(det(invA0)))
  h <- 2*log(det(chol(Apost))) - (-log(det(invA0)))
  Q <- t(y)%*%y - t(apost)%*%Apost%*%apost + t(a0)%*%invA0%*%a0 
  lml <- 0.5*(h - Q)
  
  return(lml)
}


#### (invisible function)
#### sub-function "draw_psi()" for select_logit() and select_pois()
#### last version: 2015/03/11
#### This function is not meant to be called directly by the user. 

## Sample the scale parameters psi of the slab component (for BVS)
## Arguments: 
## ... alpha (=regression effects), 
## ... index (=indices of included effects), 
## ... prior (=prior (list) object)


draw_psi <- function(alpha, index, prior){
  d <- length(alpha)
  psiv <- matrix(0, d, 1) 
  
  if (prior$slab == "Student"){
    psiv <- 1/rgamma(d, shape = prior$psi.nu + t(index)/2, rate = (prior$psi.Q + 0.5*alpha^2))
  } else if (prior$slab == "Normal"){
    psiv <- prior$psi.Q*matrix(1, d, 1)
  }
  
  psi <- t(psiv)
  return(psi)
}



