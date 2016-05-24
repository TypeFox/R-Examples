estimate_lambda <-
function(n.sim, J, prior.z, overlap, pi0){

n.zones <- length(prior.z)

# q_0 = q_1 = 1, since there is no concept of overlap for j=0, 1
q <- c(1, 1, rep(0, J-1))
var.q <- c(NA, NA, rep(0, J-1))


#-------------------------------------------------------------------------------
# All k's, sample configurations of k single zones, and see if the k single
# zones overlap or not
#-------------------------------------------------------------------------------
for(k in 2:J){
  unit.zones <- sample(1:n.zones, k*n.sim, replace=TRUE, prob=prior.z)
  unit.zones <- matrix(unit.zones, ncol=k)
  no.overlap <- check_overlap(unit.zones, overlap)
  
  denom <- rep(factorial(k), n.sim)
  
  # Compute multinomial coefficient watching out for when certain single zones
  # are sampled more than once
  duplicates <- apply(unit.zones, 1, function(x){length(unique(x))})
  duplicates <- which(duplicates != k)
  for(i in duplicates) {
    denom[i] <-  factorial(k)/prod(factorial(table(unit.zones[i,])))
  }
  
  q[k+1] <- mean(no.overlap/denom)
  var.q[k+1] <- var(no.overlap/denom)/n.sim
}


#-------------------------------------------------------------------------------
# Using q compute estimate of lambda and prior probability of j clusters
#-------------------------------------------------------------------------------
lambda <- (1-pi0)/( (1-pi0)*J + pi0*sum(q[-1]) )
lambda <- rep(lambda, J)
lambda <- c(1-sum(lambda), lambda)
prior.j <- normalize(lambda*q)

return(list(lambda=lambda, prior.j=prior.j))
}
