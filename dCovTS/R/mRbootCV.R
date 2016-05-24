mRbootCV <- function(n,q,MaxLag,b,parallel=FALSE){   
 x <- replicate(q,rnorm(n))
 if (missing(MaxLag) || MaxLag < 0) 
      stop("'MaxLag' must be greater than 0")
 A0 <- crossDist(x,0)$A
 Atilde0 <- lapply(1:q, FUN=function(j) ATilde(A0[[j]]))
  rstar <- function(k){
   A <- crossDist(x,k)$A
   B <- crossDist(x,k)$B
   Atilde <- lapply(1:q, FUN=function(j) ATilde(A[[j]]))
   Btilde <- lapply(1:q, FUN=function(j) ATilde(B[[j]]))
   cv <- vector()
   Rarray <- function(Atilde,Btilde,Atilde0,k){
    Wtstar <- rbind(rnorm(n-k))
    Wt <- rbind(rnorm(n))
    Rm <- matrix(NA,q,q)
     for (i in 1:q){
     for (j in 1:q){
       dcov <- sqrt((Wtstar%*%(Atilde[[i]]*Btilde[[j]])%*%t(Wtstar))/((n-k)^2))
       dvarx <- sqrt(mean((Atilde0[[i]]*Atilde0[[i]]))*mean((Atilde0[[j]]*Atilde0[[j]])))
        Rm[i,j] <- dcov/sqrt(dvarx)
      }
     }
    return(Rm)
   }
   result <- replicate(b,Rarray(Atilde,Btilde,Atilde0,k))
   s <- 1
   for(i in 1:q){
    for (j in 1:q){
     quant <- quantile(result[i,j,],0.95)
     pv <- mean(result[i,j,]>=quant)
     pvadj <- p.adjust(pv,method="fdr")
     cv[s] <- quantile(result[i,j,],1-pvadj)
     s <- s+1
    }
   }
   return(cv)
  }
 if(parallel==TRUE){
  closeAllConnections()
  #cl <- makeCluster(detectCores())
  cl <- makeCluster(2)
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- 1:MaxLag
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="rbind",.export=c("crossDist","ATilde")) ))
  fe <- eval(fe_call)
  Rstar <- fe %dopar% rstar(i)
  stopCluster(cl)
  mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[,i])),nrow=q,ncol=q,byrow=T)
 }
 else {
 res <- lapply(1:MaxLag,FUN=function(i) rstar(i))
 Rstar <- t(sapply(1:2,FUN=function(i) unlist(res[[i]])))
 mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[[i]])),nrow=q,ncol=q,byrow=T) 
 }
 return(mcv)
}
