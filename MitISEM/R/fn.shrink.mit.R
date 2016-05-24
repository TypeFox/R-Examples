# Function to shrink a mixture of multivariate t densities (mit)
# Shrinkage from an H comp. mixture to an H'<= H component mixture
#
# inputs:  
#   mit     : [list] mit density of H student-t components
#   tol.pr  : [double in [0,1)] minimum probability required to keep mix. components
#             (default: 0, component probability does not alter mixture)
# outputs:
#   [list] containing
#     mit.s  : [list] 'shrank' mit density with H'<=H components
#     crash  : [logical] 'TRUE' if any component is removed, 'FALSE' else
#
# note: components are removed if their scale matrices are not pds
#       (optional) components with probability p < tol.pr are removed
#
# author : Nalan Basturk
# date   : 20120912

fn.shrink.mit <- function(mit,tol.pr=0){
   # test PDS invertible covariance Sigma
   fn.testSigma <- function(vSigma){
      k <- sqrt(length(vSigma))
      r <- fn.isPDS(matrix(vSigma,k,k))
      return(r)
   }
   # test for mixture components with too small probability
   fn.testProbs <- function(probs,tol.pr){
      (probs<=tol.pr)
   }
   # indicator for zero variance for each component
   crash.Sigma <- apply(mit$Sigma,1,FUN=fn.testSigma) # indicator for crashing variance
   # indicator for too small mixture probability
   crash.Probs <- fn.testProbs(mit$p,tol.pr)
   
   crash       <- (crash.Sigma|crash.Probs)     
   if(all(crash))
      stop("mixture density shrinks to 'NULL', try different probability tolerance 'tol.pr' or increase 'N'");
   if(any(crash)){
     # 'shrank' mixture density
     mit.s      <- mapply(function(x)as.matrix(as.matrix(x)[!crash,]),mit,SIMPLIFY=FALSE)
     # correct mixture definition
     mit.s$p    = as.vector(mit.s$p) / sum(mit.s$p) 
     mit.s$df   = as.vector(mit.s$df)
	 if(ncol(mit.s$mu)==1) mit.s$mu = matrix(mit.s$mu,nrow=1)
     if(ncol(mit.s$Sigma)==1) mit.s$Sigma = matrix(mit.s$Sigma,nrow=1)     
   }else{
     mit.s      <- mit
   }
   return(list(mit=mit.s,crash=any(crash)))  
}
