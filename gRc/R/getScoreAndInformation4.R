getScore   <- function(m, K, scale='original'){
  UseMethod('getScore')
}

getScore.rcon <- function(m, K, scale='original'){ ### OK !!!

  ir    <- m$intRep  
  S     <- m$dataRep$S
  n     <- m$dataRep$n
  
  Sigma <- solve.default(K)

  vccTerms <- ir$vccI
  eccTerms <- ir$eccI
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc

  score    <- rep(NA, nparm)
  J        <- matrix(0,nrow=nparm, ncol=nparm)
  
  f  <- n-1;
  f2 <- f/2

  DSigma <- Sigma-S

  ## Find Score vector
  ##
  for (u in 1:lvcc){
    Ku <- vccTerms[[u]]
    ##val <-  f2 * trAW(Ku, DSigma)
    val <-  f2 * .Call("trAW", Ku, DSigma, PACKAGE="gRc")
    if (scale=='free')
      score[u] <- val * K[Ku[[1]],Ku[[1]]]
    else
      score[u] <- val 
  }    

  if (lecc>0){
    for (u in 1:lecc){
      Ku  <- eccTerms[[u]]
      ##val <-  f2 * trAW(Ku, DSigma)
      val <-  f2 * .Call("trAW", Ku, DSigma, PACKAGE="gRc")
      score[u+lvcc] <- val
    }    
  }  

  ## Find Fisher information
  ##
  xxxx <- .Call("trAWBWlist", vccTerms, Sigma, vccTerms, 1, PACKAGE="gRc")
  kk <- 1
  for (u in 1:lvcc){    ##print("V,V")
    Ku  <- vccTerms[[u]]
    for (v in u:lvcc){
      Kv  <- vccTerms[[v]]
      val <- f2* xxxx[kk]
      kk  <- kk +1
      #val2 <- f2* xxxx[(u-1)+lvcc*(v-1)+1] 
      #val <- f2*.Call("trAWBW", Ku, Sigma, Kv, PACKAGE="gRc")
      #print(c(u, v, val, val2))
      
      if (scale=='free')
        J[v,u] <- J[u,v] <- val * (K[Ku[[1]],Ku[[1]]]*K[Kv[[1]],Kv[[1]]])
      else
        J[v,u] <- J[u,v] <- val 
    }    
  }

  
  if (lecc>0){    ##print("V,E")    
    
    xxxx <- .Call("trAWBWlist", vccTerms, Sigma, eccTerms, 0, PACKAGE="gRc")
    
    for (u in 1:lvcc){
      Ku  <- vccTerms[[u]]
      for (v in 1:lecc){
        Kv  <- eccTerms[[v]]
        ##val <- f2* trAWBW(Ku, Sigma, Kv)
        #val <- f2* .Call("trAWBW", Ku, Sigma, Kv, PACKAGE="gRc")
        val <- f2* xxxx[(u-1)+lvcc*(v-1)+1] 
        #print(c(u, v, val,val2))
        idx = Ku[[1]]
        if (scale=='free'){
          J[u,v+lvcc] <- J[v+lvcc,u] <- val * K[idx,idx] 
        } else {
          J[u,v+lvcc] <- J[v+lvcc,u] <- val
        }
      }    
    }

    xxxx <- .Call("trAWBWlist", eccTerms, Sigma, eccTerms, 1, PACKAGE="gRc")

    kk = 1
    for (u in 1:lecc){    ##print("E,E")    
      Ku  <- eccTerms[[u]]
      for (v in u:lecc){
        Kv  <- eccTerms[[v]]
        #val <- f2* .Call("trAWBW", Ku, Sigma, Kv, PACKAGE="gRc")
                                        #val <- f2* xxxx[(u-1)+lecc*(v-1)+1]
        val <- f2* xxxx[kk]; kk=kk+1
        #print(c(u, v, val,val2))
        J[v+lvcc,u+lvcc] <-J[u+lvcc,v+lvcc] <- val #f2* trAWBW(Ku, Sigma, Kv)      
      }    
    }

    
  }
  
  return(list(score=score, J=J))
}

getScore.rcor <- function(m,K, scale='original'){
  
  md    <- m$dataRep
  ir    <- m$intRep

  S     <- md$S;
  f     <- md$n - 1

  vccTerms <- ir$vccI
  eccTerms <- ir$eccI
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc

  score    <- rep(NA, nparm)
  J        <- matrix(0,nrow=nparm, ncol=nparm)

  C        <- cov2cor(K); 
  Cinv     <- solve.default(C)             ## Brug IKKE cholSolve(C) - numerisk ustabil

  a        <- sqrt(diag(K))        ## a indeholder eta'erne
  A        <- diag(a)              

  ASA      <- a* t(a*S)
  CASA     <- C %*% ASA
  ACAS     <- (a * C) %*% (a * S)
  CinvASA  <- Cinv-ASA

  ## Find score vector
  ##
  ##cat("Score - VCC\n")
  for (u in 1:lvcc){    
    ## Score for VCC 
    Ku  <-  vccTerms[[u]]
    #val <-  f*(trA(Ku) - trAW(Ku, ACAS))
    val <-  f*(trA(Ku) - .Call("trAW", Ku, ACAS, PACKAGE="gRc")) 

    if (scale=='free')
      score[u] <- val
    else
      score[u] <- val / (A[Ku[[1]],Ku[[1]]]) ## OK, dec 07
  }

  if (lecc>0){
    ##cat("Score - ECC\n")
    ## Score for ECC:
    for (u in 1:lecc){
      Ku <- term.u <- eccTerms[[u]]
      #score[u+lvcc] <- (f/2) * (trAW(Ku, CinvASA))
      score[u+lvcc] <- (f/2) * .Call("trAW", Ku, CinvASA, PACKAGE="gRc")
    }    
  }


  ## Fisher information matrix
  ##

##   Cinv<<-Cinv
##   CC <<- C
##   vcc <<- vccTerms
##   ecc <<- eccTerms
##   break
  
  ##cat("Information - VCC x VCC\n")
  for (u in 1:lvcc){
    Ku  <-  vccTerms[[u]]
    ##cat("u:", u, "\n")
    for (v in 1:lvcc){
      ##cat("v:", v, "\n")
      Kv  <- vccTerms[[v]]      

      if (u==v)
        ##val <- 2*f* trAWBV(Ku, Cinv, Kv, C) ## OK, sept 07
        val <- 2*f* .Call("trAWBV", Ku, Cinv, Kv, C, PACKAGE="gRc") ## OK, sept 07
      else
        ##val <- f * trAWBV(Ku, Cinv, Kv, C) ## OK, sept 07
        val <- f * .Call("trAWBV", Ku, Cinv, Kv, C, PACKAGE="gRc") ## OK, sept 07
      
      if (scale=='original')
        val <- val / (A[Ku[[1]],Ku[[1]]]*A[Kv[[1]],Kv[[1]]])

      J[u,v] <- J[v,u] <-val
    }    
  }

  if (lecc>0){

    ##xxxx <- .Call("trAWBlist", vccTerms, Cinv, eccTerms, 0, PACKAGE="gRc")
    ##cat("Information - VCC x ECC\n")
    for (u in 1:lvcc){
      #cat("u:", u, "\n")
      ## Score for VCC x ECC
      Ku <-  vccTerms[[u]]
      for (v in 1:lecc){
        #cat("v:", v, "\n")
        Kv <- eccTerms[[v]]
        ##val <- f * trAWB(Ku, Cinv, Kv) ## OK, sept 07
        val <- f * .Call("trAWB", Ku, Cinv, Kv, PACKAGE="gRc") ## OK, sept 07
        ##val <- f* xxxx[(u-1)+lvcc*(v-1)+1] 
        #print(c(val,val2))
        if (scale=='free')
          J[u,v+lvcc] <- J[v+lvcc,u] <- val
        else
          J[u,v+lvcc] <- J[v+lvcc,u] <- val /(A[Ku[[1]],Ku[[1]]]);
      }    
    }


    #xxxx <- .Call("trAWBWlist", eccTerms, Cinv, eccTerms, 1, PACKAGE="gRc")
    ##cat("Information - ECC x ECC\n")
    kk = 1
    for (u in 1:lecc){
      #cat("u:", u, "\n")
      Ku <- eccTerms[[u]]
      for (v in u:lecc){
       # cat("v:", v, "\n")
        Kv  <- eccTerms[[v]]
        ##val  <- (f/2)* trAWBW(Ku, Cinv, Kv) ## OK, sept 07
        val <- (f/2) * .Call("trAWBW", Ku, Cinv, Kv, PACKAGE="gRc")
        #val <- (f/2)* xxxx[kk]; kk=kk+1
        ##print(c(val, val2))
        J[u+lvcc,v+lvcc] <- J[v+lvcc,u+lvcc] <- val 
      }    
    }

  }
  ##print(J)
  return(list(score=score, J=J))
}









