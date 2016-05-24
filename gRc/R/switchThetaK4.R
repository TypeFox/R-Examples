## New !!
# theta2K2 <- function(theta, vccTerms, eccTerms, nrK, type, scale='original'){
#   if (type=="rcon")
#     theta2K2.rcon(theta, vccTerms, eccTerms, nrK, scale)
# }

# K2theta2 <- function(K, vccTerms, eccTerms, type, scale='original'){
#   if (type=="rcon")
#     K2theta2.rcon(K, vccTerms, eccTerms, scale)

# }

# theta2K2.rcon <- function(theta, vccTerms, eccTerms, nrK, scale='original'){
  
#   K        <- matrix(0, nrow=nrK, ncol=nrK)
#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc
  
#   for (u in 1:lvcc){
#     term.u <- vccTerms[[u]]
#     val    <- theta[u]
#     for (j in 1:nrow(term.u)){
#       term.uj <- term.u[j,]
#       K[term.uj,term.uj] <- val 
#     }  
#   }

#   if (scale=='free')
#     diag(K) <- exp(diag(K))
  
#   if (lecc>0){
#     for (u in 1:lecc){
#       term.u <- eccTerms[[u]]
#       val <- theta[u+lvcc]
#       for (j in 1:nrow(term.u)){
#         term.uj <- term.u[j,]
#         term.uj <- rep(term.uj,2)[1:2]
#         K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val
#       }  
#     }
#   }  
#   return(K)
# }

# K2theta2.rcon <- function(m, K, scale='original'){

#   vccTerms <- m$intRep$vccI
#   eccTerms <- m$intRep$eccI

#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc
#   theta <- rep(NA, nparm)

#   for (u in 1:length(vccTerms)){
#     term.u <- vccTerms[[u]]
#     #print(term.u)
#     term.u <- term.u[[1]]
#     #print(term.u)
#     term.u <- rep(term.u,2)[1:2]
#     if (scale=='free')
#       val <- log(K[term.u[1],term.u[2]]) # lambda
#     else
#       val <- K[term.u[1],term.u[2]]      # eta
#     theta[u] <- val
#   }
#   if (length(eccTerms)>0){
#     for (u in 1:length(eccTerms)){
#       term.u <- eccTerms[[u]]
#       term.u <- term.u[1,]
#       term.u <- rep(term.u,2)[1:2]
#       val <- K[term.u[1],term.u[2]] 
#       theta[u+lvcc] <- val 
#     }
#   }
#   return(theta)
# }

# ## !!



theta2K <- function(m, theta, scale='original'){
  UseMethod("theta2K")
}

K2theta <- function(m, K, scale='original'){
  UseMethod("K2theta")
}


## Matrices
theta2K.rcon <- function(m, theta, scale='original'){

  p        <- nrow(m$dataRep$S)
  vccTerms <- m$intRep$vccI
  eccTerms <- m$intRep$eccI
    
  K        <- matrix(0, nrow=p, ncol=p)
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc
 
  for (u in 1:lvcc){
    term.u <- vccTerms[[u]]
    val    <- theta[u]
    for (j in 1:nrow(term.u)){
      term.uj <- term.u[j,]
      ## term.uj <- rep(term.uj,2)[1:2]
      ##K[term.uj[1],term.uj[2]] <- val
      K[term.uj,term.uj] <- val 
    }  
  }

  if (scale=='free')
    diag(K) <- exp(diag(K))
  
  if (lecc>0){
    for (u in 1:lecc){
      term.u <- eccTerms[[u]]
      val <- theta[u+lvcc]
      #print(val)
      for (j in 1:nrow(term.u)){
        ##for (j in 1:length(term.u)){
        ##term.uj <- term.u[[j]]
        term.uj <- term.u[j,]
        term.uj <- rep(term.uj,2)[1:2]
        #print(term.uj)
        K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val
      }  
    }
  }
  return(K)
}



theta2K.rcor <- function(m, theta, scale='original'){

  p        <- nrow(m$dataRep$S)
  vccTerms <- m$intRep$vccI
  eccTerms <- m$intRep$eccI

  K        <- matrix(0, nrow=p, ncol=p)
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc

  ##print(theta); print(scale)
  for (u in 1:lvcc){
    term.u <- vccTerms[[u]]
    #print(term.u)
    for (j in 1:length(term.u)){
      term.uj <- term.u[[j]]
      term.uj <- rep(term.uj,2)[1:2]
      if (scale=='free')
        val <- exp(theta[u]*2)
      else
        val <- theta[u]^2

      K[term.uj[1],term.uj[2]] <- val 
    }  
  }
  ##print(K)

  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]];       
      val    <- theta[u+lvcc]
      for (j in 1:nrow(term.u)){
        term.uj <- term.u[j,]
        term.uj <- rep(term.uj,2)[1:2]
        tuj1 <- term.uj[1]
        tuj2 <- term.uj[2]
        #print(c(tuj1,tuj2))
        #print(K[tuj1,tuj1]); print(K[tuj2,tuj2])
        K[tuj1,tuj2] <- K[tuj2,tuj1] <- val*(sqrt(K[tuj1,tuj1] * K[tuj2,tuj2]))        
        #print("OK here")
      }
    }
  }
  return(K)
}




K2theta.rcon <- function(m, K, scale='original'){

  vccTerms <- m$intRep$vccI
  eccTerms <- m$intRep$eccI

  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc
  theta <- rep(NA, nparm)

  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    #print(term.u)
    term.u <- term.u[[1]]
    #print(term.u)
    term.u <- rep(term.u,2)[1:2]
    if (scale=='free')
      val <- log(K[term.u[1],term.u[2]]) # lambda
    else
      val <- K[term.u[1],term.u[2]]      # eta
    theta[u] <- val
  }
  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      term.u <- term.u[1,]
      term.u <- rep(term.u,2)[1:2]
      val <- K[term.u[1],term.u[2]] 
      theta[u+lvcc] <- val 
    }
  }
  return(theta)
}


## For diagonal elements K[i,i] of K, the corresponding parameters are
## log(sqrt(K[i,i])) : scale = 'free'
## sqrt(K[i,i])      : scale = 'original'
## For off-diagonal elements K[i,j] the parameters are those values
## scaled with the diagonals

K2theta.rcor <- function(m, K, scale='original'){
    
  vccTerms <- m$intRep$vccI
  eccTerms <- m$intRep$eccI
  
  lvcc     <- length(vccTerms)
  lecc     <- length(eccTerms)
  nparm    <- lvcc + lecc
  theta <- rep(NA, nparm)
  
  C <- cov2cor(K)

  for (u in 1:lvcc){
    term.u <- vccTerms[[u]]
    #print(term.u)
    term.u <- term.u[[1]]
    #print(term.u)
    term.u <- rep(term.u,2)[1:2]
    if (scale=='free')
      val <- log(sqrt(K[term.u[1],term.u[2]])) # lambda
    else
      val <- sqrt(K[term.u[1],term.u[2]])      # eta
    theta[u] <- val
  }
  if (lecc>0){
    for (u in 1:lecc){
      term.u <- eccTerms[[u]]
      #print(term.u)
      ##term.u <- term.u[[1]]
      term.u <- term.u[1,]
      #print(term.u)
      term.u <- rep(term.u,2)[1:2]

      ##val <- K[term.u[1],term.u[2]] / sqrt(K[term.u[1],term.u[1]]*K[term.u[2],term.u[2]])
      val <- C[term.u[1],term.u[2]] 
      theta[u+lvcc] <- val 
    }
  }
  return(theta)
}
















# theta2K.rcon <- function(m, theta, scale='original'){
#   md       <- getSlot(m, 'dataRep')
#   ir       <- getSlot(m, 'intRep')  
#   S        <- md$S;
#   p        <- nrow(S)
#   K        <- matrix(0, nrow=p, ncol=p)
#   vccTerms <- ir$vccI
#   eccTerms <- ir$eccI
#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc
 
#   for (u in 1:lvcc){
#     term.u <- vccTerms[[u]]
#     for (j in 1:length(term.u)){
#       term.uj <- term.u[[j]]
#       term.uj <- rep(term.uj,2)[1:2]
#       val <- theta[u]
#       K[term.uj[1],term.uj[2]] <- val 
#     }  
#   }

#   if (scale=='free')
#     diag(K) <- exp(diag(K))
  
#   if (lecc>0){
#     for (u in 1:lecc){
#       term.u <- eccTerms[[u]];       
#       val <- theta[u+lvcc]
#       for (j in 1:length(term.u)){
#         term.uj <- term.u[[j]]
#         term.uj <- rep(term.uj,2)[1:2]
#         K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val
#       }  
#     }
#   }
#   return(K)
# }





# theta2K.rcor <- function(m, theta, scale='original'){

#   md       <- getSlot(m, 'dataRep')
#   ir       <- getSlot(m, 'intRep')  
#   S        <- md$S;
#   p        <- nrow(S)
#   K        <- matrix(0, nrow=p, ncol=p)
#   vccTerms <- ir$vccI
#   eccTerms <- ir$eccI
#   ## print(eccTerms)
#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc

#   ## print(theta)
#   for (u in 1:lvcc){
#     term.u <- vccTerms[[u]]
#     for (j in 1:length(term.u)){
#       term.uj <- term.u[[j]]
#       term.uj <- rep(term.uj,2)[1:2]
#       if (scale=='free')
#         val <- exp(theta[u]*2)
#       else
#         val <- theta[u]^2
#       K[term.uj[1],term.uj[2]] <- val 
#     }  
#   }

#   if (length(eccTerms)>0){
#     for (u in 1:length(eccTerms)){
#       term.u <- eccTerms[[u]];       
#       val    <- theta[u+lvcc]
#       for (j in 1:length(term.u)){
#         term.uj <- term.u[[j]]
#         term.uj <- rep(term.uj,2)[1:2]
#         tuj1 <- term.uj[1]
#         tuj2 <- term.uj[2]
#         K[tuj1,tuj2] <- K[tuj2,tuj1] <- val*(sqrt(K[tuj1,tuj1] * K[tuj2,tuj2]))        
#       }  
#     }
#   }
#   return(K)
# }

# K2theta.rcon <- function(m, K, scale='original'){

#   md       <- getSlot(m, 'dataRep')
#   ir       <- getSlot(m, 'intRep')  
#   S        <- md$S;  
#   vccTerms <- ir$vccI
#   eccTerms <- ir$eccI
#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc
#   theta <- rep(NA, nparm)

#   for (u in 1:length(vccTerms)){
#     term.u <- vccTerms[[u]]
#     term.u <- term.u[[1]]
#     term.u <- rep(term.u,2)[1:2]
#     if (scale=='free')
#       val <- log(K[term.u[1],term.u[2]]) # lambda
#     else
#       val <- K[term.u[1],term.u[2]]      # eta
#     theta[u] <- val
#   }
#   if (length(eccTerms)>0){
#     for (u in 1:length(eccTerms)){
#       term.u <- eccTerms[[u]]
#       term.u <- term.u[[1]]
#       term.u <- rep(term.u,2)[1:2]
#       val <- K[term.u[1],term.u[2]] 
#       theta[u+lvcc] <- val 
#     }
#   }
#   return(theta)
# }





# ## For diagonal elements K[i,i] of K, the corresponding parameters are
# ## log(sqrt(K[i,i])) : scale = 'free'
# ## sqrt(K[i,i])      : scale = 'original'
# ## For off-diagonal elements K[i,j] the parameters are those values
# ## scaled with the diagonals

# K2theta.rcor <- function(m, K, scale='original'){

#   md       <- getSlot(m, 'dataRep')
#   ir       <- getSlot(m, 'intRep')  
#   S        <- md$S;  
#   vccTerms <- ir$vccI
#   eccTerms <- ir$eccI
#   lvcc     <- length(vccTerms)
#   lecc     <- length(eccTerms)
#   nparm    <- lvcc + lecc
#   theta <- rep(NA, nparm)
  
#   C <- cov2cor(K)

#   for (u in 1:lvcc){
#     term.u <- vccTerms[[u]]
#     term.u <- term.u[[1]]
#     term.u <- rep(term.u,2)[1:2]
#     if (scale=='free')
#       val <- log(sqrt(K[term.u[1],term.u[2]])) # lambda
#     else
#       val <- sqrt(K[term.u[1],term.u[2]])      # eta
#     theta[u] <- val
#   }
#   if (lecc>0){
#     for (u in 1:lecc){
#       term.u <- eccTerms[[u]]
#       term.u <- term.u[[1]]
#       term.u <- rep(term.u,2)[1:2]
#       #val <- K[term.u[1],term.u[2]] / sqrt(K[term.u[1],term.u[1]]*K[term.u[2],term.u[2]])
#       val <- C[term.u[1],term.u[2]] 
#       theta[u+lvcc] <- val 
#     }
#   }
#   return(theta)
# }














# ##getThetaFromK2.rcor <- function(m,K=m$fit$K,scale='original'){
# K2theta.rcor <- function(m, K, scale='original'){
# ##cat("getThetaFromK.rcor\n")
#   S<-m$S
#   dimnames(K)<- dimnames(S)
#   v <- .modelInfo(m,'stdrepI')
#   vccTerms <- v$vcc
#   eccTerms <- v$ecc

#   lvcc <- length(vccTerms)
#   lecc <- length(eccTerms)

#   theta <- rep(NA, lvcc+lecc)

#   C <- cov2cor(K)

#   for (u in 1:lvcc){
#     term.u <- vccTerms[[u]]
#     term.u <- term.u[[1]]
#     term.u <- rep(term.u,2)[1:2]
#     if (scale=='free')
#       val <- log(sqrt(K[term.u[1],term.u[2]])) # lambda
#     else
#       val <- sqrt(K[term.u[1],term.u[2]])      # eta
#     theta[u] <- val
#   }
#   if (lecc>0){
#     for (u in 1:lecc){
#       term.u <- eccTerms[[u]]
#       term.u <- term.u[[1]]
#       term.u <- rep(term.u,2)[1:2]
#       #val <- K[term.u[1],term.u[2]] / sqrt(K[term.u[1],term.u[1]]*K[term.u[2],term.u[2]])
#       val <- C[term.u[1],term.u[2]] 
#       theta[u+lvcc] <- val 
#     }
#   }
#   return(theta)
# }



# ##getKFromTheta2.rcor <- function(m,theta=m$theta,scale='original'){
# theta2K.rcor <- function(m, theta, scale='original'){
  
# ##cat("getKFromTheta.rcor\n")
  
#   S <- m$S
#   v <- .modelInfo(m,'stdrepI')
#   vccTerms <- v$vcc
#   eccTerms <- v$ecc
#   lvcc <- length(vccTerms)
#   K <- S; K[,] <- 0;   

  
#   for (u in 1:length(vccTerms)){
#     term.u <- vccTerms[[u]]
#     for (j in 1:length(term.u)){
#       term.uj <- term.u[[j]]
#       term.uj <- rep(term.uj,2)[1:2]
#       if (scale=='free')
#         val <- exp(theta[u]*2)
#       else
#         val <- theta[u]^2
#       K[term.uj[1],term.uj[2]] <- val 
#     }  
#   }

#   if (length(eccTerms)>0){
#     for (u in 1:length(eccTerms)){
#       term.u <- eccTerms[[u]];       
#       val <- theta[u+lvcc]
#       for (j in 1:length(term.u)){
#         term.uj <- term.u[[j]]
#         term.uj <- rep(term.uj,2)[1:2]
#         K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val*(sqrt(K[term.uj[1],term.uj[1]] * K[term.uj[2],term.uj[2]]))        
#       }  
#     }
#   }
#   return(K)
# }




#####################################################################################




##getThetaFromK2.rcon <- function(m,K=m$fit$K,scale='original'){ ## OK !!!
  ##cat("getThetaFromK.rcon\n")


##getKFromTheta2 <- function(m,theta=m$theta,scale='original') UseMethod('getKFromTheta2')

##getKFromTheta2.rcon <- function(m,theta=m$theta,scale='original'){
##cat("getKFromTheta.rcon\n")


#getThetaFromK2   <- function(m,K=m$fit$K,scale='original') UseMethod('getThetaFromK2')

# getThetaFromK2.rcon <- function(m,K=m$fit$K,scale='original'){ ## OK !!!

  
#   ##cat("getThetaFromK.rcon\n")
  
#   S<-m$S
#   dimnames(K)<- dimnames(S)

#   v <- .modelInfo(m,'stdrepI')
#   vccTerms <- v$vcc
#   eccTerms <- v$ecc
#   lvcc <- length(vccTerms)
#   theta <- rep(NA, length(vccTerms)+length(eccTerms))
  
#   for (u in 1:length(vccTerms)){
#     term.u <- vccTerms[[u]]
#     term.u <- term.u[[1]]
#     term.u <- rep(term.u,2)[1:2]
#     if (scale=='free')
#       val <- log(K[term.u[1],term.u[2]]) # lambda
#     else
#       val <- K[term.u[1],term.u[2]]      # eta
#     theta[u] <- val
#   }
#   if (length(eccTerms)>0){
#     for (u in 1:length(eccTerms)){
#       term.u <- eccTerms[[u]]
#       term.u <- term.u[[1]]
#       term.u <- rep(term.u,2)[1:2]
#       val <- K[term.u[1],term.u[2]] 
#       theta[u+lvcc] <- val 
#     }
#   }
#   return(theta)
# }







