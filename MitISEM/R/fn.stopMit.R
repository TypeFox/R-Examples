# Function to assess convergence of the MitISEM algorithm
# depending on the CV convergence or AR acceptance rate convergence
#
# inputs:
#   method  : [string] 'CV' or 'AR' defining convergence criterion
#   w       : [vector size N] IS weights given N parameter draws
#   CV_last : [double], last CV value to compare new CV
#   CVtol   : [double >0] convergence tolerance for CV
#   AR_last : [double], last AR acceptance prob. to compare new CV
#   ARtol   : [double >0] convergence tolerance for AR prob.
# outputs:
#   out     : [vector size 3], with new CV value, new AR prob value, 
#              stop (logical for convergence achieved)
# author : Nalan Basturk
# date   : 20120912


# Compute new CV from IS weights
# and indicator to finalize the number of mixture comps in MitISEM
fn.CVstop <- function(w,CV_last,CVtol){
  if(sd(w)==0)
    stop ("IS weights 'w' are constant, try increasing number of draws 'N'")
  CV_new  <- sd(w)/mean(w)
  hstop   <- abs((CV_new-CV_last)/CV_last) <= CVtol
  return(c(CV_new,hstop))
}
# Calculate (approximate expected) acceptance rate given IS weights from draws
# and indicator to finalize the number of mixture comps in MitISEM
fn.ARstop <- function(w,AR_last,ARtol){
  AR_new  <- mean(w) / max(w)
  hstop   <- abs((AR_new-AR_last)/AR_last) <= ARtol
  return(c(AR_new,hstop))
}

# Function indicating end of adding mixture comps in MitISEM candidate
fn.stopMit <- function(method='CV',w,CV_last,CVtol,AR_last,ARtol){
  if(all(method!=c('CV','AR')))
    stop("Stopping method for MitISEM should be 'CV' or 'AR'")
  
  r1 <- fn.CVstop(w,CV_last,CVtol)
  r2 <- fn.ARstop(w,AR_last,ARtol)
  
  r      <- rbind(r1,r2)
  tmp    <- which(method==c('CV','AR')) 
  hstop  <- r[tmp,2]    
  
  out        <- c(r1[1],r2[1],hstop)
  names(out) <- c('CV','AR','stop')
  return(out)
}