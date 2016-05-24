
calculateVCOV <- function(m, K, vcov="boot", nboot=250){
  vn          <- unlist(lapply(getcc(m),names))
  switch(vcov,
         "boot"={
           #cat(".bootstrapVCOV", nboot, "\n")
           V <- .bootstrapVCOV(m,K,nboot)
         },
         "inf" ={
           #cat("scoreVar", nboot, "\n")
           V <- solve(getScore(m,K=K)$J)
         }
         )
  dimnames(V) <- list(vn,vn)  
  V
}



calculateVCOV <- function(m, K, vcov="boot", nboot=250){
  vn          <- unlist(lapply(getcc(m),names))
  V <- solve(getScore(m,K=K)$J)
  
#   switch(vcov,
#          "boot"={
#            #cat(".bootstrapVCOV", nboot, "\n")
#            V <- .bootstrapVCOV(m,K,nboot)
#          },
#          "inf" ={
#            #cat("scoreVar", nboot, "\n")
#            V <- solve(getScore(m,K=K)$J)
#          }
#          )
  dimnames(V) <- list(vn,vn)  
  V
}



.bootstrapVCOV <- function(m, K, nboot=250){

  ##cat("bootstrapVar", nboot, "\n")
  res  <- matrix(NA, ncol=length(c(getSlot(m,"vcc"),getSlot(m,"ecc"))), nrow=nboot)
  Sf   <- solve(K)
  n    <- dataRep(m, "n")
  for(i in 1:nboot){
    d<- mvrnorm(n, mu=rep(0, nrow(Sf)),
                Sigma=Sf, tol = 1e-6, empirical = FALSE)
    m$dataRep$S <- cov(d)
    if (class(m)[1]=="rcon")
      th <- rconScoreTheta(m)
    else
      th <- rcorScoreTheta(m)
    res[i,] <- th
  }  
  V           <- cov(res)  
} 
