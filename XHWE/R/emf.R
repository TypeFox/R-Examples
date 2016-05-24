emf <-
function(n1m, n0m, n2f, n1f, n0f, nm, nf, rho, dv, itertime)
  
{


  iter <- 0 

  nfT2 <- 2 * nf

  pf <- (n2f * 2 + n1f) / nfT2    # the initial value of pf

  nonstop <- TRUE 
    
  while(nonstop & (iter <= itertime))
      
  {
      
    temp <- rho * pf * (1 - pf)

    pf.sq <- pf^2

    qf.sq <- (1 - pf)^2
      
    deno1 <- pf.sq + temp
      
    deno2 <- qf.sq + temp
      
    e1 <- n2f * pf.sq / deno1
      
    e2 <- n2f * temp / deno1
      
    e3 <- n0f * qf.sq / deno2
      
    e4 <- n0f * temp / deno2

    sume2e4 <- e4 + e2
      
    trho <- sume2e4 / (sume2e4 + n1f)    
      
    tpf <- (2 * e1 + n1f + sume2e4) / nfT2
 
    nonstop <- (abs(trho - rho) > dv | abs(tpf - pf)> dv)

    if(is.nan(nonstop)|is.na(nonstop))  break
      
    if(nonstop & (iter <= itertime))  
        
    {
        
      rho <- trho
        
      pf <- tpf    
        
      iter <- iter + 1   
        
    }
      
  }

  rho.last <- rho
  
  pf.last <- pf

  if (nm == 0){ 

  list(rho.last = rho.last, pf.last = pf.last)

  }

  else{  

  pm.last <- n1m / nm
 
  list(rho.last = rho.last,pm.last = pm.last,pf.last = pf.last)

  }

}
