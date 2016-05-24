probm <- function(eta, margin, only.pr = TRUE, bc = FALSE){ # bc stands for binary continuous case
 
  epsilon <- 0.0000001 
  max.p   <- 0.9999999 
  
  derp1.dereta1 <- der2p1.dereta1eta1 <- d.n <- der2p.dereta <- NULL
 
if( margin == "probit" ){
 
  pr  <- pnorm(eta)
  
  
  if(only.pr == FALSE){
  
  d.n <- dnorm(eta)
  d.n <- ifelse(d.n < epsilon, epsilon, d.n )
  der2p.dereta <- -(eta * dnorm(eta))          # second deriv
  
  }
  
  if(bc == TRUE){ # RECALL that this is deriv of Y = 0 wrt eta1 and not the usual case of Y = 1!!
  derp1.dereta1      <- -dnorm(-eta)         # First derivative of 1-prob(eta1) respect to eta1.
  der2p1.dereta1eta1 <- eta * dnorm(-eta)  ## This is the second derivative of 1 - p1 respect to eta1  
  }  
    
    
    
}


if( margin == "logit" ){
 
  pr  <- plogis(eta)
  
  if(only.pr == FALSE){
  
  d.n <- dlogis(eta)
  d.n <- ifelse(d.n < epsilon, epsilon, d.n )
  der2p.dereta <- -((1 - 2 * (exp(-eta)/(1 + exp(-eta)))) * exp(-eta)/(1 + exp(-eta))^2)
  
  }
  
  if(bc == TRUE){
  derp1.dereta1    <- -((1 - exp(-eta)/(1 + exp(-eta))) * exp(-eta)/(1 + exp(-eta))) # First derivative of 1-prob(eta1) respect to eta1.
      
  der2p1.dereta1eta1 <- (1 - (3 - 2 * (exp(-eta)/(1 + exp(-eta)))) * exp(-eta)/(1 + exp(-eta))) * 
                           exp(-eta)/(1 + exp(-eta))
 ## This is the second derivative of 1 - p1 respect to eta1  
  }
  
}



#if( margin == "loglog" ){
#
#  pr  <- exp(-exp(eta)) 
#  d.n <- -(exp(-exp(eta)) * exp(eta))
#  d.n <- ifelse(d.n < epsilon, epsilon, d.n )
#  der2p.dereta <- -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta)) 
#    
#} #not in glm

if( margin == "cloglog" ){
 
  pr  <- 1-exp(-exp(eta))
  
  
  if(only.pr == FALSE){
  
  d.n <- exp(-exp(eta)) * exp(eta)
  d.n <- ifelse(d.n < epsilon, epsilon, d.n )
  der2p.dereta <- (1 - exp(eta)) * exp(-exp(eta)) * exp(eta) 
  
  
  }
  
  
  if(bc == TRUE){
  
  derp1.dereta1    <-  -(exp(-exp(eta)) * exp(eta))
                   
  der2p1.dereta1eta1 <- -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta))
 
                }
 
}

if( margin == "cauchit" ){
 
  pr  <- 1 / pi * atan(eta) + 0.5
  
  
  if(only.pr == FALSE){
  
  d.n <- 1 / (pi * (1 + eta^2))
  d.n <- ifelse(d.n < epsilon, epsilon, d.n )
  der2p.dereta <- -(2 * (eta/(pi * (1 + eta^2)^2)))
  
  }
  
  
  if(bc == TRUE){
  
  derp1.dereta1    <- -(1/(pi * (1 + eta^2)))
         
  der2p1.dereta1eta1 <- 2 * (eta/(pi * (1 + eta^2)^2))
 ## This is the second derivative of 1 - p1 respect to eta1 

  }
  
  
}


  pr <- mm(pr) 
  
    
    
 list(pr = pr, d.n = d.n, der2p.dereta = der2p.dereta, 
      derp1.dereta1 = derp1.dereta1, der2p1.dereta1eta1 = der2p1.dereta1eta1 )  
 
}    