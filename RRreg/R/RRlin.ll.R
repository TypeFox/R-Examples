# loglikelihood up to an unknown constant
# phi=c(beta,sigma,pi) : vector with parameters
# nbeta : number of regression parameters
# y : response variable (vector)
# u : nonRR variables (matrix)
# gidx : index for individual group membership (vector)
# U : number of nonRR variables included? (logical)
# intercept : include? (logical)
# n.star : aggregated frequencies. rows: response patterns; cols: groups
# PWarray: missclassification matrices as array[response,true state, group]
# PWpp : individual probabilities taken from PWarray-rows [individual,state]
# M : number of RR variables
RRlin.ll <- function(phi, nbeta, y, u, gidx, U, intercept, n.star, PWarray, PWpp, M){
  
  # INPUT : different names for parameters
  beta <- phi[1:nbeta]
  sigma <- phi[nbeta+1]
  pi.true <- c(phi[(nbeta+2):length(phi)] , 0)
  # sum(pi) has to be smaller than 1 !
#   pi.true[pi.true<=0] <- 1e-6
#   pi.true[pi.true>=1] <- 1-1e-6
#   if (sum(pi.true)!=)
  
  pi.true[length(pi.true)] <- 1-sum(pi.true)

  G <- ncol(n.star) # number of groups
  J <- nrow(n.star)  # number of response patterns
  K <- ncol(PWarray[,,1])   # number of true states (for Kukrep>1)
  N <- length(y)
  if(M==1){
    patterns <- dimnames(PWarray)$true   # only 1 RR variable
  }else{
    patterns <- unlist(strsplit(dimnames(PWarray)$true,":"))  # more RR variables
  }
  patterns <- matrix(as.numeric(patterns),nrow=K,byrow=T)
  
  # loop across all group-combinations
  first <- rep(0,G)
  sec <- rep(0,G)
  for (g in 1:G){
    # first part of sum
    #regression part
    Ng <- colSums(n.star)[g]
    t1 <- matrix(NA,Ng,K)
    for (j in 1:K){
      xx <- matrix(patterns[j,], Ng, M, T)
      if (U>0){
        xx <- cbind(xx, u[gidx == g,,drop=FALSE])
      }
      if (intercept){
        xx <- cbind(1,xx)
      }
      t1[,j] <- exp(-.5 * ( ( y[gidx == g] - (xx %*% beta) ) /sigma)^2 )
    }
    
    # missclassification part: N x K matrix
    t2 <- PWpp[gidx == g,] * matrix(pi.true, nrow=Ng, ncol=K, byrow=T)
    # constants for all persons:
    
    temp <- rowSums(t1 * t2) / (sigma*sqrt(2*pi) * PWpp[gidx == g,]%*%pi.true )
    first[g] <- sum(log( temp ) )
    
    # second part of loglik
    sec[g] <- t(n.star[,g]) %*% log( PWarray[,,g] %*% pi.true )
  }
  ll <- sum (first,sec)
  if(ll == -Inf)
    ll <- -9e300
  ll
}