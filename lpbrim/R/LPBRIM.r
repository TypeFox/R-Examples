#' @title Bipartite label propagation (LP)
#' @description Performs LP
#' @export
#' @param x Bipartite network
#' @param as.adjacency TRUE if x is a matrix, FALSE if it its an object with three columns
bLP <- function (x,as.adjacency=TRUE) {
   if(as.adjacency) x[x>0] <- 1
   OrderVec <- c(rownames(x),colnames(x))
   x <- x[sample(c(1:NROW(x))),sample(c(1:NCOL(x)))]
   # HTL labels
   lT <- c(1:NROW(x))
   names(lT) <- rownames(x)
   # LTL labels
   lB <- rep(0,NCOL(x))
   names(lB) <- colnames(x)
   # Seeding the initial modularity
   oldQ <- 0
   newQ <- 1e-10
   Nsteps <- 0
   # Labels propagation loop
   while(oldQ<newQ)
   {
      Nsteps <- Nsteps + 1
      oldQ <- newQ
      # Step 1 : update lB
      for(lsp in 1:NCOL(x))
      {
         Nei <- rownames(x)[x[,lsp]>0]
         NeiLab <- lT[Nei]
         if(as.adjacency)
         {
            lB[lsp] <- NeiLab[mostFrequent(NeiLab, NA)]
         } else {
            lB[lsp] <- NeiLab[mostFrequent(NeiLab, x[Nei,lsp])]
         }
      }
      names(lB) <- colnames(x)
      # Step 2 : update lT
      for(tsp in 1:NROW(x))
      {
         Nei <- colnames(x)[x[tsp,]>0]
         NeiLab <- lB[Nei]
         if(as.adjacency)
         {
            lT[tsp] <- NeiLab[mostFrequent(NeiLab, NA)]
         } else {
            lT[tsp] <- NeiLab[mostFrequent(NeiLab, x[tsp,Nei])]
         }
      }
      names(lT) <- rownames(x)
      # Shaping the vectors
      Modules <- c(lT,lB)
      Comms <- unique(Modules)
      NComm <- length(Comms)
      Smat <- matrix(0,ncol=NComm,nrow=sum(dim(x)))
      colnames(Smat) <- Comms
      rownames(Smat) <- c(rownames(x),colnames(x))
      for(i in 1:length(Modules)) Smat[names(Modules)[i],as.character(Modules[i])]<-1
      newQ <- Qbip(x,Smat)
   }
   Modules <- as.numeric(as.factor(Modules))
   names(Modules) <- c(names(lT),names(lB))
   return(Modules[OrderVec])
}

#' @title Barber's modularity
#' @description Returns Q
#' @export
#' @param x Bipartite adjacency matrix
#' @param s Community partition
Qbip <- function(x,s)
{
   Q <- 0
   x[x>0] <- 1
   p <- NROW(x)
   h <- NCOL(x)
   m <- sum(x)
   nc <- NCOL(s)
   A <- x
   P <- matrix(kronecker(colSums(x),rowSums(x)),nrow=NROW(x),ncol=NCOL(x))/m
   B <- A-P
   Rm <- s[c(1:p),]
   # If the network is not modular
   if(is.null(dim(Rm))) return(0)
   Tm <- s[c(p+(1:h)),]
   # Induce Qr from T
   BT <- B%*%Tm
   Isum <- NULL
   for(i in 1:p)
   {
      Ksum <- NULL
      for(k in 1:nc)
      {
         Ksum[k] <- Rm[i,k]*(BT[i,k])
      }
      Isum[i] <- sum(Ksum)
   }
   Q = (1/m)*sum(Isum)
   return(Q)
}

#' @title BRIM
#' @description Performs the BRIM algorithm
#' @export
#' @param x Bipartite adjacency matrix
bBRIM = function(x)
{
   Nsteps <- 0
   CommDiv <- bLP(x)
   x[x>0] <- 1
   Comms <- unique(CommDiv)
   NComm <- length(Comms)
   Smat = matrix(0,ncol=NComm,nrow=sum(dim(x)))
   colnames(Smat) <- Comms
   rownames(Smat) <- c(rownames(x),colnames(x))
   # Fill the S matrix
   for(i in 1:length(CommDiv)) Smat[names(CommDiv)[i],as.character(CommDiv[i])]<-1
   # Initial modularity
   FromR <- TRUE
   cBM <- -10
   preBM <- 10
   # Some important values
   p <- NROW(x)
   h <- NCOL(x)
   m <- sum(x)
   nc <- NCOL(Smat)
   A <- x
   P <- matrix(kronecker(colSums(x),rowSums(x)),nrow=NROW(x),ncol=NCOL(x))/m
   B <- A-P
   # Optimization loop
   while((Nsteps<3)|(preBM!=cBM))
   {
      Nsteps <- Nsteps + 1
      preBM <- Qbip(x,Smat)
      # Modularity matrix
      Rm <- as.matrix(Smat[c(1:p),])
      Tm <- as.matrix(Smat[c(p+(1:h)),])
      # Product matrix for T & R
      rBT <- B%*%Tm
      tBT <- t(B)%*%Rm
      # Optimization
      if(FromR)
      {
         Rm[,] <- 0
         Rm[cbind(1:NROW(x),apply(rBT,1,which.max))] <- 1
      } else {
         Tm[,] <- 0
         Tm[cbind(1:NCOL(x),apply(tBT,1,which.max))] <- 1
      }
      Smat[c(1:NROW(x)),] <- Rm
      Smat[(NROW(x)+c(1:NCOL(x))),] <- Tm
      # New bipartition
      cBM <- Qbip(x,Smat)
      FromR <- !FromR
   }
   return(list(S=Smat,M=x,Q=cBM,c=NCOL(Smat)))
}
