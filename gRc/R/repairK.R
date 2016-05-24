## Based on a start value KS an estimate K in the model is returned.
##
## Entries in K which should be zero according to the model are so.
## Non-zero elements in K are obtained by averaging relevant entries in KS
##

findKinModel <- function(m, KS, type="rcon", regularize=TRUE)
  UseMethod("findKinModel")


## MATRIX VERSION
findKinModel.rcox <- function(m, KS, type="rcon", regularize=TRUE){
  if (is.null(KS))
    return(NULL)
  
  switch(type,
         'rcor'={
           ##C <- cov2cor(KS)
           ##print(C)
           C <- findKinModel(m,cov2cor(KS),   type='rcon')
           A <- findKinModel(m,diag(diag(KS)),type='rcon')
           a <- sqrt(diag(A))
           KKmod <- a * t(a*C) ## Short for KKmod <- A%*%C%*%A
           KKmod           
         },
         
         'rcon'={
           KS2 <- KS;
           KS2[,] <- 0; diag(KS2) <- diag(KS)
           VCC <- intRep(m,"vccI")
           ECC <- intRep(m,"eccI")

           for (i in seq(along=ECC)){
             x      <- ECC[[i]]
             x2    <- x[,2:1,drop=FALSE]
             KS2[x] <- KS2[x2] <- mean(KS[x])
           }
           
           for (i in seq(along=VCC)){
             x           <-unlist(VCC[[i]])
             if (nrow(x)>1){
               x           <- as.numeric(x)
               diag(KS2)[x] <- mean(diag(KS)[x])
             }
           }
           
           if (regularize && (min(eigen(KS2)$values) < 0)){             
             KKmod <- regularizeK(KS2)
           } else {
             KKmod <- KS2
           }
         }
         )
  return(KKmod)
}



## K is made positive definite
##

regularizeK <- function(K){
  Kdiag <- abs(diag(diag(K)))           
  Krest <- K-Kdiag
  alpha <- 0.9
  repeat{
    Kmod  <- Kdiag + alpha*Krest    
    if (min(eigen(Kmod)$values) > 0){
      Kmod <- Kdiag + 0.95*alpha*Krest ## Be on the safe side...
      ##cat("The end", alpha,"\n")
      break()
    }
    alpha <- alpha-0.1    
  }
  Kmod
}


