emc <-
function(n1m, n0m, n2f, n1f, n0f, nm, nf, rho, dv, itertime) {

  sum.all <- nm + 2 * nf
  
  pc <- (n1m + n2f * 2 + n1f)/sum.all    # the initial value of pc
  
  iter <- 0      
    
  nonstop <- TRUE 
    
  while(nonstop & (iter <= itertime))
      
  {
      
    temp <- rho * pc * (1 - pc)

    pc.sq <- pc^2

    qc.sq <- (1 - pc)^2
      
    deno1 <- pc.sq + temp
      
    deno2 <- qc.sq + temp
      
    e1 <- n2f * pc.sq / deno1
      
    e2 <- n2f * temp / deno1
      
    e3 <- n0f * qc.sq / deno2
      
    e4 <- n0f * temp / deno2

    sume2e4 <- e4 + e2
      
    trho <- sume2e4 / (sume2e4 + n1f)    
      
    tpc <- (2 * e1 + n1f + sume2e4 + n1m) / sum.all

    nonstop <- (abs(trho - rho) > dv | abs(tpc - pc)> dv)

    if(is.nan(nonstop)|is.na(nonstop))  break
      
    if(nonstop & (iter <= itertime)) 
        
    {
        
      rho <- trho
        
      pc <- tpc    
        
      iter <- iter + 1   
        
     }

  }
  
  rho.last <- rho
  
  pc.last <- pc
  
  list(rho.last = rho.last, pc.last = pc.last)
  
}
