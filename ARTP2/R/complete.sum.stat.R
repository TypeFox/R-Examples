
complete.sum.stat <- function(sum.stat, options){
  
  msg <- paste("Calculating P and SE if not provided:", date())
  if(options$print) message(msg)
  
  nf <- length(sum.stat$stat)
  lambda <- sum.stat$lambda
  for(i in 1:nf){
    
    st <- sum.stat$stat[[i]]
    id.no.SE <- which(is.na(st$SE))
    id.no.P <- which(is.na(st$P))
    
    if(length(id.no.SE) > 0){
      z2 <- qchisq(st$P[id.no.SE], df = 1, lower.tail = FALSE)
      st$SE[id.no.SE] <- abs(st$BETA[id.no.SE]/sqrt(z2))
    }
    
    if(length(id.no.P) > 0){
      st$P[id.no.P] <- pchisq((st$BETA[id.no.P]/st$SE[id.no.P])^2, df = 1, lower.tail = FALSE)
    }
    
    st$SE <- sqrt(lambda[i]) * st$SE
    st$P <- pchisq((st$BETA/st$SE)^2, df = 1, lower.tail = FALSE)
    
    sum.stat$stat[[i]] <- st
    rm(st)
    gc()
  }
  
  sum.stat
  
}

