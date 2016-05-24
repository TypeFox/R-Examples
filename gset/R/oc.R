oc <-
function(l, u, theta, sigma, K, n1, n2, boundaries, futility=TRUE, binding=FALSE, n.sim=1e4)
{
  reject.t2 <- rep(0, n.sim)
  accept.t2 <- rep(0, n.sim)

  t <- 1/K*(1:K) 
  nn1<- ceiling(n1*t); nn1[K]<- n1      
  nn2<- ceiling(n2*t); nn2[K]<- n2   
 
  for(i in 1:n.sim) {
    x1 <- rnorm(n1, 0, sigma)
    x2 <- rnorm(n2, theta, sigma)
    for(k in 1:K) {
      x1.stage   <- x1[1:nn1[k]]   
      mean1.temp <- mean(x1.stage)
      v1.temp    <- var(x1.stage)
      
      x2.stage   <- x2[1:nn2[k]]   
      mean2.temp <- mean(x2.stage)
      v2.temp    <- var(x2.stage)  
      
      v.pool <- (v1.temp*(nn1[k]-1)+ v2.temp*(nn2[k]-1))/(nn1[k]+ nn2[k]-2)
      se.pool <- sqrt(v.pool*(1/nn1[k]+1/nn2[k]))
      diff <- mean2.temp -  mean1.temp
      
      t.L <-  (diff - l) / se.pool
      t.U <-  (diff - u) / se.pool  
      
      if(futility){
        if(binding){   
      		reject <- (t.L > boundaries$equivL[k]) & (t.U < boundaries$equivU[k])
      		accept <- (t.L <= boundaries$futilL[k]) | (t.U>= boundaries$futilU[k])
      		if(reject.t2[i] == 0 & accept.t2[i] == 0 & reject) reject.t2[i] <- nn1[k]
      		if(reject.t2[i] == 0 & accept.t2[i] == 0 & accept) accept.t2[i] <- nn1[k]
      	}	
        else{
          reject <- (t.L > boundaries$equivL[k]) & (t.U < boundaries$equivU[k])
      		accept <- (t.L <= boundaries$futilL[k]) | (t.U >= boundaries$futilU[k])
      		if(reject.t2[i] == 0 & reject) reject.t2[i] <- nn1[k]
      		if(reject.t2[i] == 0 & accept.t2[i] == 0 & accept) accept.t2[i] <- nn1[k]
        }
  	  }
  	  else {   		      		
    		reject <- (t.L > boundaries$equivL[k] & t.U < boundaries$equivU[k])
    		if(reject.t2[i] == 0 & reject) reject.t2[i] <- nn1[k]
      }
    }
  }
  if(futility) {
    R <-	sum(reject.t2 > 0)/n.sim
   	A <-  sum(accept.t2 > 0)/n.sim	
    stopn1 <- reject.t2 + accept.t2
    stopn2 <- stopn1
    for(i in 1:n.sim) { 
      if(stopn1[i]==0) {stopn1[i] <- nn1[K]; stopn2[i] <- nn2[K]}}
    En1 <- sum(stopn1)/n.sim
    En2 <- sum(stopn2)/n.sim
    
    stopp <- rep(0,K)
    for(k in 1:(K-1)) stopp[k]<- mean(reject.t2==nn1[k]| accept.t2==nn1[k])   
    stopp[K] <- 1-sum(stopp[1:(K-1)])
    
    stoppE <- rep(0,K)
    for(k in 1:K) stoppE[k]<- mean(reject.t2==nn1[k])   
    
    stoppF <- rep(0,K)
    for(k in 1:K) stoppF[k]<- mean(accept.t2==nn1[k])   
    
    return(list(reject.rate=R, accept.rate=A, En1=round(En1,1), En2=round(En2,1), 
                prob.stop=stopp, prob.stopE=stoppE, prob.stopF=stoppF)) 
  }
  else {
    R  <-	sum(reject.t2 > 0)/n.sim
    stopn1 <- reject.t2
    stopn2 <- stopn1
    for(i in 1:n.sim) { 
      if(stopn1[i]==0) {stopn1[i] <- nn1[K]; stopn2[i] <- nn2[K]}}
    En1 <- sum(stopn1)/n.sim
    En2 <- sum(stopn2)/n.sim
    
    stopp <- rep(0,K)
    for(k in 1:K)  stopp[k]<- mean(reject.t2==nn1[k])
    stopp[K] <- 1-sum(stopp[1:(K-1)])
    
    stoppE <- rep(0,K)
    for(k in 1:K) stoppE[k]<- mean(reject.t2==nn1[k])   
    
    stoppF <- rep(0,K)
    for(k in 1:K) stoppF[k]<- mean(accept.t2==nn1[k])   
    
    return(list(reject.rate=R,  En1=round(En1,1),  En2=round(En2,1),
                prob.stop=stopp, prob.stopE=stoppE, prob.stopF=stoppF)) 
  }     		                        
}
