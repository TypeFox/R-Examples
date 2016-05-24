
rconScoreMatch <- function(m, control=m$control, trace=0){

  tstart <- proc.time()
  if (trace>=2)
    cat("..Fitting with score matching\n")

  theta <- rconScoreTheta(m)
  allcc <- list(vcc=m$vcc, ecc=m$ecc)
  vn           <- unlist(lapply(allcc,names))
  names(theta) <- vn
  
  S <- m$dataRep$S
  n <- m$dataRep$n

  K            <- theta2K(m, theta, scale='original')
  dimnames(K)  <- dimnames(S)

  diagK <- diag(K)
  ## print("After score matching"); print(ellK(K,S,n-1))
  ## Ensure diagonals are positive
  ## 
  if (min(diagK)<0){
    vccI <- intRep(m)$vccI
    dS <- diag(S)
    for (i in 1:length(vccI)){
      cc <- as.numeric(vccI[[i]])
      dS[cc] <- mean(dS[cc])
    }
    diag(K) <- 1/dS
  }
  
  if (min(eigen(K)$values)<0){
    K     <- regularizeK(K)
  }

  logL  <- ellK(K,S,n-1) #print("After redoing diagonals"); print(logL)
  
  if (trace>=3)
    cat("..Score matching, logL:", logL, "Time:", proc.time()-tstart, "\n")

  ans <- list(K=K, logL=logL, coef=theta)
  return(ans)
}




rconScoreTheta <- function(m){

#  cat("C-implementation of scorematching\n")
  vcc <- m$intRep$vccI
  ecc <- m$intRep$eccI
  S <- m$dataRep$S
  n <- m$dataRep$n
  logL <- -99999

  ans <- .Call("scorematching", 
               S     = S, 
               nobs  = n, 
               vcc   = vcc, 
               ecc   = ecc,
               logL  = logL, 
               debug=5, PACKAGE="gRc")
  return(ans)
}


## This one is now obsolete
##
.rconScoreTheta <- function(m){

  vccTerms <- m$intRep$vccI
  eccTerms <- m$intRep$eccI

  lvcc <- length(vccTerms)
  lecc <- length(eccTerms)
  ltot <- lvcc+lecc
  
  S <- m$dataRep$S
  n <- m$dataRep$n
  
  A <- matrix(0, ncol=ltot, nrow=ltot)
  B <- rep(0,ltot)

  
  for (u in 1:lvcc){
    Ku     <- vccTerms[[u]]
    #bu     <- trA(Ku)

    #bu  <- nrow(Ku) ## This is tr(A)
    #B[u]   <- bu
    B[u] <- nrow(Ku)
    ##A[u,u] <- trAWB(Ku, S, Ku)
    ##A[u,u] <- .Call("trAWB", Ku, S, Ku, PACKAGE="gRc")
  }

  xxxx <- .Call("trAWBlist", vccTerms, S, vccTerms, 1, PACKAGE="gRc") # OK
  kk <- 1
  for (u in 1:lvcc){
    A[u,u] <- xxxx[kk]
    kk <- kk + lvcc - (u-1)
  }
  
  if (lecc>0){
    
    xxxx <- .Call("trAWBlist", vccTerms, S, eccTerms, 0, PACKAGE="gRc") #OK
    zzzz <- .Call("trAWBlist", eccTerms, S, eccTerms, 1, PACKAGE="gRc") #OK
                                        #print(xxxx)
    for (u in 1:lvcc){
      #Ku      <- vccTerms[[u]]
      for (v in 1:lecc){
        #Kv     <- eccTerms[[v]]
        auv <- xxxx[(u-1)+lvcc*(v-1)+1]
        #auv2 <- .Call("trAWB", Ku, S, Kv, PACKAGE="gRc")
        #auv3 <- trAWB(Ku, S, Kv)
        #print(c(u, v, auv, auv2,auv3))        
        A[u,v+lvcc] <- A[v+lvcc,u] <- auv

      }
    }


    kk <- 1
    for (u in 1:lecc){
      #Ku      <- eccTerms[[u]]
      for (v in u:lecc){
        #Kv     <- eccTerms[[v]]
        #auv <- xxxx[(u-1)+lecc*(v-1)+1]
        auv  <- zzzz[kk]; kk <- kk+1
        #auv2 <- .Call("trAWB", Ku, S, Kv, PACKAGE="gRc")
        #auv3 <- trAWB(Ku, S, Kv)
        #print(c(auv, auv2,auv3))
        A[u+lvcc,v+lvcc] <- A[v+lvcc,u+lvcc] <- auv
      }
    }

  } # if (lecc>0)

  
  theta <- solve.default(A,B)
##   print(A)
##   print(B)
##   print(theta)
                                        #theta <- qr.solve(A,B)  
  return(theta)
}

#     for (u in 1:lvcc){
#       Ku      <- vccTerms[[u]]
#       for (v in 1:lecc){
#         Kv     <- eccTerms[[v]]
#         auv    <- trAWB(Ku, S, Kv)
#         #auv <- .Call("trAWB", Ku, S, Kv, PACKAGE="gRc")
#         A[u,v+lvcc] <- A[v+lvcc,u] <- auv
#       }
#     }

#     for (u in 1:lecc){
#       Ku      <- eccTerms[[u]]
#       for (v in u:lecc){
#         Kv     <- eccTerms[[v]]
#         auv    <- trAWB(Ku, S, Kv)
#         #auv <- .Call("trAWB", Ku, S, Kv, PACKAGE="gRc")
#         A[u+lvcc,v+lvcc] <- A[v+lvcc,u+lvcc] <- auv      
#       }
#     }






rcorScoreMatch <- function(m, control=m$control, trace=0){
  tstart <- proc.time()
  if (trace>=2)
    cat("..Fitting with score matching (rcor)\n")
  #theta     <- rcorScoreTheta(m)
  theta        <- rconScoreTheta(m)
  vn           <- unlist(lapply(getcc(m),names))
  names(theta) <- vn

  S        <- getSlot(m, "dataRep")$S
  ##Sorig    <- getSlot(m, "dataRep")$Sorig ## BRIS
  n        <- getSlot(m, "dataRep")$n
  
  oclass       <- class(m)
  class(m)     <- c("rcon","rcox")
  K            <- theta2K(m, theta, scale='original')
  ## print("As RCON"); print(K)

  ## Ensure diagonals are positive
  if (min(diag(K))<0){  
    vccI <- intRep(m)$vccI
    dS <- diag(S)
    for (i in 1:length(vccI)){
      cc <- as.numeric(vccI[[i]])
      dS[cc] <- mean(dS[cc])
    }
    diag(K) <- 1/dS
  }
  ##print("After redoing diagonals");  print(K)

  if (min(eigen(K)$values)<0)  
    K     <- regularizeK(K)
  ##print("After regularizing diagonals"); print(K)

  dimnames(K)  <- dimnames(S)
  ##print(ellK(K,S,n-1))
  
  class(m) <- oclass
  K <- findKinModel(m,K,type="rcor")

  #print(K)  
  #print(ellK(K,S,n-1))
  #print(1/sqrt(diag(S)))
  
  ## K            <- theta2K(m, theta, scale="original")
  ## dimnames(K)  <- dimnames(S)

  logL   <- ellK(K, S, n-1)

  ##cat(">>logL (scale):", logL, "logL (orig):", ellK(K,Sorig,n-1),"\n")
  if (trace>=3)
    cat("..Score matching, logL:", logL, "Time:", proc.time()-tstart, "\n")


  return(list(K=K, logL=logL, coef=theta))
}

rcorScoreTheta <- function(m){


  theta        <- rconScoreTheta(m)
  vn           <- unlist(lapply(getcc(m),names))
  names(theta) <- vn

  S        <- getSlot(m, "dataRep")$S
  Sorig        <- getSlot(m, "dataRep")$Sorig
  n        <- getSlot(m, "dataRep")$n

  oclass <- class(m)
  class(m) <- c("rcon","rcox")
  K            <- theta2K(m, theta, scale='original')

  dimnames(K)  <- dimnames(S)
  ##print(ellK(K,S,n-1))

  ## Ensure diagonals are positive
  vccI <- intRep(m)$vccI
  dS <- diag(S)
  for (i in 1:length(vccI)){
    cc <- as.numeric(vccI[[i]])
    dS[cc] <- mean(dS[cc])
  }
  diag(K) <- 1/dS
  
  K     <- regularizeK(K) 
  logL  <- ellK(K,S,n-1)
  ##print(ellK(K,S,n-1))
  
  class(m) <- oclass
  theta    <- K2theta(m, K, scale="original")

  vn       <- unlist(lapply(getcc(m),names))
  names(theta) <- vn
  return(theta)
}


  

#   ctrl   <- m$control
#   ctrl$vcov <- NULL
#   class(m) <- c("rcon","rcox")
#   KS       <- rconScoreMatch(m, control=ctrl)$K
#   class(m) <- oclass
#   K.curr   <- findKinModel(m, KS=KS, type='rcor', regularize=TRUE)

#   theta    <- K2theta(m, K.curr, scale="original")
#   vn       <- unlist(lapply(getcc(m),names))
#   names(theta) <- vn
#   return(theta)

#   if (!is.null(control$vcov)){
#     V <- calculateVCOV(m, K=K, vcov=control$vcov, nboot=control$nboot)
#   }




#   print(K)  
#   K <<- K
#   S <<- S
#   print(S)
#   K <<- K; S<<-S

#  print(ellK(K,S,n-1))
  
#  if (min(eigen(K)$values)<0){
#    cat("Score matching gave n.d. K; regularizing to make K p.d.\n")
#    K <- regularizeK(K)
#  }   



## iR <- getSlot(m, "intRep")
## vccTerms <- iR$vcc
## eccTerms <- iR$ecc
## allTerms <- c(vccTerms, eccTerms)

  #ctrl      <- control
  ##ctrl$vcov <- NULL ## We don't want vcov calculated until at the end...
  #m$control <- ctrl


# rconScoreThetaDiag <- function(m){
#   iR <- getSlot(m, "intRep")
#   vccTerms <- iR$vcc
#   S    <- getSlot(m, "dataRep")$S
#   n    <- getSlot(m, "dataRep")$n

#   a <- B <- rep(0,length(vccTerms))
#    for (u in 1:length(vccTerms)){
#      Ku <- vccTerms[[u]]
#      bu     <- trA(Ku)
#      B[u]   <- bu
#      a[u]   <- trAWB(Ku, S, Ku)
#    }
#   theta <- B/a

#   theta <- c(theta, rep(0,length(iR$ecc)))
  
#   K            <- theta2K(m, theta, scale='original')
#   dimnames(K)  <- dimnames(S)
#   return(K)

# }
