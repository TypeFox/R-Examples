#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       credible.region.R                   ####
####                                                 ####
#### FUNCTIONS:  credible.region                     ####
#########################################################

### ======================================
### credible.region
### ======================================
##
## Function to compute a simultaneous credible region for a vector parameter
## from the MCMC sample
##
## Reference: Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995)
##            Bayesian computation and stochastic systems (with Discussion)
##            Statistical Science, vol. 10, 3 - 66, page 30
## and        Held, L. (2004)
##            Simultaneous inference in risk assessment; a Bayesian perspective
##            In: COMPSTAT 2004, Proceedings in Computational Statistics (J. Antoch, Ed.)
##            213 - 222, page 214
##
## \item{sample}{a data frame or matrix with sampled values (one column = one parameter)}
## \item{probs}{probabilities for which the credible regions are to be computed}
##
credible.region <- function(sample, probs=c(0.90, 0.975))
{
  if (!length(dim(sample))) stop("Incorrect sample parameter supplied")
  if (sum(probs < 0.50)) stop("probs must be each at least 0.50")
  if (sum(probs >= 1)) stop("probs must be strictly lower than 1")
  
  nn <- dim(sample)[1]
  p <- dim(sample)[2]
  k <- floor(nn*probs)

  ## Get ranks separately for each parameter
  ranks <- apply(sample, 2, rank, ties.method="first")

  ## Compute sets S_j = {n + 1 - min_i r_i^{(j)}, max_i r_i^{(j)}}, j=1,...,n
  ## and the set S
  S.left <- nn + 1 - apply(ranks, 1, min)
  S.right <- apply(ranks, 1, max)
  SS <- apply(cbind(S.left, S.right), 1, max)
  SS <- SS[order(SS)]

  ## Get jstar and the credible region
  jstar <- SS[k]
  result <- list()
  for (kk in 1:length(jstar)){
    up <- sample[ranks == jstar[kk]]
    low <- sample[ranks == nn + 1 - jstar[kk]]    
    result[[kk]] <- rbind(low, up)
    rownames(result[[kk]]) <- c("Lower", "Upper")
    colnames(result[[kk]]) <- colnames(sample)
  }    
  names(result) <- paste(probs)

  return(result)
}


