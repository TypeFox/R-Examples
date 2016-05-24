nminmax <-
function(l, u, theta, sigma, n1.lower, n2.lower, t.vec, type1, type2, gamma=rep(-4,2), 
                    binding=FALSE, n1.upper=ceiling(2*n1.lower), n2.upper=ceiling(2*n2.lower),
                    n.sim=1e4, seed=NULL) 
{
  n1L <- n1.lower; n1U <- n1.upper
  n2L <- n2.lower; n2U <- n2.upper
  K <- length(t.vec)
  
  if(binding==FALSE){
    continue <- TRUE
    prev.n <- n1L
    while(continue){
       n1M <- ceiling((n1L+n1U)/2)
       n2M <- ceiling((n2L+n2U)/2)
       cat("*")
       bound3 <- nonbinding(l,u, theta, sigma, n1M, n2M, t.vec, type1, type2, gamma=gamma, 
                            plot=FALSE, force = FALSE, n.sim=n.sim, seed=seed)   
       if(bound3$futilL[K]-bound3$equivL[K] > 0.001){n1U <- n1M; n2U <- n2M}
       else if(bound3$futilL[K]-bound3$equivL[K] < -0.001){n1L <- n1M; n2L <- n2M}
       else continue <- FALSE
       if(n1M == prev.n) {
         continue <-FALSE
         bound3$futilL[K]<- (bound3$futilL[K]+bound3$equivL[K])/2
         bound3$equivL[K]<-  bound3$futilL[K]
         bound3$futilU[K]<- -bound3$futilL[K]
         bound3$equivU[K]<- -bound3$equivL[K]
       }
       prev.n <- n1M 
    }
    if(continue) print("no solution found; increase the values of n1.upper and n2.upper")
  }
  else{
    continue <- TRUE
    prev.n <- n1L
    while(continue){
      n1M <- ceiling((n1L+n1U)/2)
      n2M <- ceiling((n2L+n2U)/2)
      bound3 <- binding(l,u, theta, sigma, n1M, n2M, t.vec, type1, type2, gamma=gamma, 
                        plot=FALSE, force = FALSE, n.sim=n.sim, seed=seed)    
      if(bound3$futilL[K]-bound3$equivL[K] > 0.001){n1U <- n1M; n2U <- n2M}
      else if(bound3$futilL[K]-bound3$equivL[K] < -0.001){n1L <- n1M; n2L <- n2M}
      else continue <- FALSE
      if(n1M == prev.n) {
        continue <-FALSE
        bound3$futilL[K]<- (bound3$futilL[K]+bound3$equivL[K])/2
        bound3$equivL[K]<-  bound3$futilL[K]
        bound3$futilU[K]<- -bound3$futilL[K]
        bound3$equivU[K]<- -bound3$equivL[K]
      }
      prev.n <- n1M 
    }
    if(continue) print("no solution found; increase the values of n1.upper and n2.upper")
  }
  return(list(n1minmax = n1M, n2minmax = n2M, 
              typeI= bound3$typeI, typeII= bound3$typeII, 
              equivL= bound3$equivL, equivU= bound3$equivU, 
              futilL= bound3$futilL, futilU= bound3$futilU))
}
