"epi.bohning" <- function(obs, exp, alpha = 0.05){
  
   J <- length(obs)
   smr <- obs / exp
   smr.bar <- sum(smr) / J

   # Bohning's test:   
   top <- (1 / (J - 1)) * sum(((obs - (smr.bar * exp))^2) / exp) - smr.bar
   bottom <- sqrt((2 * smr.bar) / (J - 1))
   bohning <- top / bottom
   p <- 1 - pnorm(q = bohning, mean = 0, sd = 1)
               
   # Results:
   rval <- as.data.frame(cbind(test.statistic = bohning, p.value = p)) 
   return(rval)
}