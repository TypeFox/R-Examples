"check.pmatrix" <-
function(X, minvi = .03){

   compute.Ppp <- function(X,P1,N,J,m){
      label <- as.vector(t(outer(paste("P(X",1:J,">=",sep=""),paste(1:(m-1),")",sep=""), paste, sep="")))
      Ppp <- matrix(0,J*(m-1),J*(m-1))
      i <- 0
      j <- 0
      for(i in 1:(J-1)) for(j in (i+1):J)
       Ppp[((i-1)*(m-1)+1):((i-1)*(m-1)+(m-1)),((j-1)*(m-1)+1):((j-1)*(m-1)+(m-1))] <-
       t(outer(X[,i],0:(m-2),">")) %*% outer(X[,j],0:(m-2),">")/N
      Ppp <- Ppp + t(Ppp) + kronecker(diag(J),matrix(-1,m-1,m-1))
      Ppp[Ppp < -.5] <- NA
      dimnames(Ppp) <- list(label,label)
      return(Ppp[order(P1),order(P1)])
   }

   compute.Pmm <- function(X,P1,N,J,m){
      label <- as.vector(t(outer(paste("P(X",1:J,"<",sep=""),paste(1:(m-1),")",sep=""), paste, sep="")))
      Pmm <- matrix(0,J*(m-1),J*(m-1))
      i <- 0
      j <- 0
      for(i in 1:(J-1)) for(j in (i+1):J)
       Pmm[((i-1)*(m-1)+1):((i-1)*(m-1)+(m-1)),((j-1)*(m-1)+1):((j-1)*(m-1)+(m-1))] <-
       t(outer(X[,i],0:(m-2),"<=")) %*% outer(X[,j],0:(m-2),"<=")/N
      Pmm <- Pmm + t(Pmm) + kronecker(diag(J),matrix(-1,m-1,m-1))
      Pmm[Pmm < -.5] <- NA
      dimnames(Pmm) <- list(label,label)
      return(Pmm[order(P1),order(P1)])
   }

   Scores2Steps <- function(X,P1 = unit.order){
      N <- nrow(X) 
      J <- ncol(X) 
      maxx <- max(X)
      Y <- matrix(t(X),1,N*J)
      Y <- matrix(rep(Y,maxx),maxx,N*J,TRUE)
      Y <- ifelse(Y < row(Y),0,1)
      Y <- matrix(as.vector(Y),N,maxx*J,TRUE)
      unit.order = 1:ncol(Y)
      Y <- Y[,order(P1)]
      dimnames(Y)[[2]] <- as.vector(t(outer(paste("X",1:J,">=",sep=""), paste(1:maxx, sep = ""), paste, sep = "")))[order(P1)]
      return(Y)
   }   

  X <- check.data(X)
  J <- ncol(X)
  N <- nrow(X)
  ncat <- max(X) + 1
  maxx <- ncat - 1
  P1 <- matrix(t(apply(outer(as.matrix(X), 1:maxx, ">=") * 1, c(2,3), mean)), nrow = maxx * J)
  X.ISRF <- Scores2Steps(X,P1) 
  I.item <- rep(1:J,each=maxx)[order(P1)]
  I.step <- dimnames(X.ISRF)[[2]]
  I.labels <- dimnames(X)[[2]]
  if(length(I.labels) == 0) I.labels <- paste("C", 1:ncol(X))

  sorted.P1 <- sort(P1)
  Ppp <- compute.Ppp(X,P1,N,J,ncat)
  Pmm <- compute.Pmm(X,P1,N,J,ncat)
  # NIEUW
  K <- length(I.step)
  res <- list()
  res$Ppp <- Ppp
  res$Pmm <- Pmm
  res$ac <- (J-1) * (J-2) * maxx^3
  res$vi <- list()
  res$n.vi <- list()
  res$max.vi <- list()
  res$sum.vi <- list()
  res$z <- list()
  res$n.z <- list()
  res$max.z <- list()

  nvi.Ppp <- nvi.Pmm <- NULL
  for (i in 1:K){
    res$vi$Pmm[[i]] <- res$vi$Ppp[[i]] <- list()
    res$z$Pmm[[i]] <- res$z$Ppp[[i]] <- list()
    names(res$vi$Ppp)[[i]] <- names(res$z$Ppp)[[i]] <- paste("Subgroup ",I.step[i],sep="")
    names(res$vi$Pmm)[[i]] <- names(res$z$Pmm)[[i]] <- paste("Subgroup ",dimnames(Pmm)[[1]][i],sep="")

    TFpp  <- TFmm <- matrix(FALSE,nrow(Ppp),ncol(Ppp))
    Dpp <- outer(Ppp[i,],Ppp[,i],"-")
    Dpp[is.na(Dpp)] <- 0
    TFpp[upper.tri(Dpp) & Dpp > minvi] <- TRUE
    TFpp[lower.tri(Dpp) & Dpp < -minvi] <- TRUE
    Dmm <- outer(Pmm[i,],Pmm[,i],"-")
    TFmm[upper.tri(Dmm) & Dmm < -minvi] <- TRUE
    TFmm[lower.tri(Dmm) & Dmm > minvi] <- TRUE
    Dmm[is.na(Dmm)] <- 0
    
    
    res$vi$Ppp[[i]] <- abs(TFpp * Dpp)
    res$vi$Pmm[[i]] <- abs(TFmm * Dmm)
    nvi.Ppp <- rbind(nvi.Ppp,res$vi$Ppp[[i]])
    nvi.Pmm <- rbind(nvi.Pmm,res$vi$Pmm[[i]])
    # res$n.vi$Ppp[i,] <- apply(TFpp,1,sum)
    # res$n.vi$Pmm[i,] <- apply(TFmm,1,sum)
  }
  res$n.vi$total <- apply(sign(rbind(nvi.Ppp,nvi.Pmm)),2,sum)
  res$n.vi$Ppp <- apply(sign(nvi.Ppp),2,sum)
  res$n.vi$Pmm <- apply(sign(nvi.Pmm),2,sum)
  res$max.vi$total <- apply(rbind(nvi.Ppp,nvi.Pmm),2,max)
  res$max.vi$Ppp <- apply(nvi.Ppp,2,max)
  res$max.vi$Pmm <- apply(nvi.Pmm,2,max)
  res$sum.vi$total <- apply(rbind(nvi.Ppp,nvi.Pmm),2,sum)
  res$sum.vi$Ppp <- apply(nvi.Ppp,2,sum)
  res$sum.vi$Pmm <- apply(nvi.Pmm,2,sum)

  nz.Ppp <- nz.Pmm <- NULL
  for (i in 1:K){
    sample11 <- X.ISRF[,i]==1
    if (sum(sample11) < 2) Z <- matrix(0,K,K) else {
       Np11 <- compute.Ppp(X[sample11,],P1,N,J,ncat) * N
       Np1  <- matrix(rep(apply(X.ISRF[sample11,],2,sum),K),K,K)
       Np01 <- round(Np1 - Np11)
       Np10 <- t(Np01) 
       k <- pmin(Np01,Np10,na.rm=TRUE)
       n <- Np01 + Np10
       n[n < 1] <- .5
       B <- ((2 * k + 1 - n)^2 - 10 * n) / (12 * n)
       Z <- abs(sqrt(2 * k + 2 + B) - sqrt(2 * n - 2 * k + B))
       Z <- Z * sign(res$vi$Ppp[[i]])
    }   
    Z[Z < qnorm(.95)] <- NA
    res$z$Ppp[[i]] <- Z
    nz.Ppp <- rbind(nz.Ppp,Z)

    sample00 <- X.ISRF[,i]==0
    if (sum(sample00) < 2) Z <- matrix(0,K,K) else {
       Nm00 <- compute.Pmm(X[sample00,],P1,N,J,ncat) * N
       Nm1  <- matrix(rep(apply(X.ISRF[sample00,],2,sum),K),K,K)
       Nm0  <- nrow(X[sample00,]) - Nm1
       Nm01 <- round(t(Nm0) - Nm00)
       Nm10 <- t(Nm01) 
       k <- pmin(Nm01,Nm10,na.rm=TRUE)
       n <- Nm01 + Nm10
       n[n < 1] <- .5
       B <- ((2 * k + 1 - n)^2 - 10 * n) / (12 * n)
       Z <- abs(sqrt(2 * k + 2 + B) - sqrt(2 * n - 2 * k + B))
       Z <- Z * sign(res$vi$Pmm[[i]])
    }  
    Z[Z < qnorm(.95)] <- NA
    res$z$Pmm[[i]] <- Z
    nz.Pmm <- rbind(nz.Pmm,Z)
  }
  nz.Ppp[is.na(nz.Ppp)] <- 0
  nz.Ppp[is.nan(nz.Ppp)] <- 0
  nz.Pmm[is.na(nz.Pmm)] <- 0
  nz.Pmm[is.nan(nz.Pmm)] <- 0
  res$n.z$Ppp <- apply(sign(nz.Ppp),2,sum)
  res$n.z$Pmm <- apply(sign(nz.Pmm),2,sum)
  res$n.z$total <- res$n.z$Ppp + res$n.z$Pmm
  res$max.z$Ppp <- apply(nz.Ppp,2,max)
  res$max.z$Pmm <- apply(nz.Pmm,2,max)
  res$max.z$total <- pmax(res$max.z$Ppp,res$max.z$Pmm)
  
  Hi <- coefH(X,FALSE)$Hi
  pmatrix.list <- list(results=res, I.item=I.item, I.step=I.step, I.labels=I.labels, Hi=Hi, minvi=minvi, ncat=ncat, N=N)
  class(pmatrix.list) <- "pmatrix.class"  
  return(pmatrix.list)
}
