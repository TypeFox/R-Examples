# $Id: fdh.R 140 2015-05-15 21:48:02Z B002961 $

# FDH efficiency uden brug af LP.
# Der er ingen kontrol af argumenter, den taenkes at blive kaldt fra dea.R
# der udfoerer alle kontroller


fdh <- function(X,Y, ORIENTATION="in", XREF=NULL, YREF=NULL, 
                FRONT.IDX=NULL, DIRECT=NULL, TRANSPOSE=FALSE, oKr=0)  {

   if ( ORIENTATION=="graph" )
      stop("ORIENTATION==\"graph\" does not work for fdh,",
            "use dea( ...,RTS=\"fdh\", ... ) ")

   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }

   if ( TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }
   orgKr <- dim(XREF)

   if ( length(FRONT.IDX) > 0 && oKr==dim(XREF)[1] )  {
      # Brug kun FRONT.IDX hvis den ikke allerede er brugt og saa vil
      # oKr vaare forskellig fra dim(XREF)[1]
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in 'dea'")
      XREF <- XREF[FRONT.IDX,, drop=FALSE]
      YREF <- YREF[FRONT.IDX,, drop=FALSE]
   }
   rNames <- rownames(XREF)
   if ( is.null(rNames) & !is.null(colnames(YREF)) )
      rNames <- rownames(YREF)

   K <- dim(X)[1]
   m <- dim(X)[2]
   n <- dim(Y)[2]
   Kr <- dim(XREF)[1]
   # Saat kun oKr hvis den ikke allerede er sat i kaldet af fdh
   if (oKr==0) oKr <- orgKr[1]
   eff <- rep(NA, K)
   peer <- rep(NA, K)


# Directional efficiency
if ( !is.null(DIRECT) )  {

   if ( class(DIRECT)=="matrix" && dim(DIRECT)[1] > 1 ) {
      if ( ORIENTATION=="in" )  {
         dirX <- DIRECT  # matrix(DIRECT, nrow=K, ncol=m)
         # dirY <- matrix(.Machine$double.xmin, nrow=K, ncol=n)
         dirY <- matrix(NA, nrow=K, ncol=n)
      } else if ( ORIENTATION=="out" )  {
         # dirX <- matrix(.Machine$double.xmin, nrow=K, ncol=m)
         dirX <- matrix(NA, nrow=K, ncol=m)
         dirY <- DIRECT  # matrix(DIRECT, nrow=K, ncol=n)
      } else if ( ORIENTATION=="in-out" )  {
         dirX <- DIRECT[,1:m,drop=FALSE]   
                 # matrix(DIRECT[,1:m], nrow=K, ncol=m)
         dirY <- DIRECT[,(m+1):(m+n), drop=FALSE]  
                 # matrix(DIRECT[,(m+1):(m+n)], nrow=K, ncol=n)
      }
   } else {
      # Her er DIRECT en vektor og derfor ens for alle firms, dvs.
      # alle raekker skal vaere ens
      if ( ORIENTATION=="in" )  {
         dirX <- matrix(DIRECT, nrow=K, ncol=m, byrow=T)
         # dirY <- matrix(.Machine$double.xmin, nrow=K, ncol=n)
         dirY <- matrix(NA, nrow=K, ncol=n)
      } else if ( ORIENTATION=="out" )  {
         # dirX <- matrix(.Machine$double.xmin, nrow=K, ncol=m)
         dirX <- matrix(NA, nrow=K, ncol=m)
         dirY <- matrix(DIRECT, nrow=K, ncol=n, byrow=T)
      } else if ( ORIENTATION=="in-out" )  {
         dirX <- matrix(DIRECT[1:m], nrow=K, ncol=m, byrow=T)
         dirY <- matrix(DIRECT[(m+1):(m+n)], nrow=K, ncol=n, byrow=T)
      }
   }

   eps <- 1e-6

   for ( k in 1:K )  {
      # For each firm find max(XREF/X) over inputs (rows)
      # Se kun for de firmaer der dominerer firm k i den retning der 
      # ikke aendres
      # xk <-  NULL  # matrix(X[k,,drop=FALSE], nrow=Kr, ncol=m, byrow=TRUE)
      # yk <-  NULL  # matrix(Y[k,,drop=FALSE], nrow=Kr, ncol=n, byrow=TRUE)
      xk <- matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
      yk <- matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
      if ( ORIENTATION=="in" )  {
         idx <- 
          rowSums( (xk >= XREF-eps & dirX[k,]==0) | abs(dirX[k,])>0+eps) & 
                rowSums(yk <= YREF+eps) ==  n
      } else if ( ORIENTATION=="out" )  {
         idx <- rowSums(xk >= XREF) == m & 
           rowSums( (yk <= YREF+eps & dirY[k,]==0) | abs(dirY[k,])>0+eps)
      } else if ( ORIENTATION=="in-out" )  {
         idx <- 
           rowSums( (xk >= XREF-eps & dirX[k,]==0) | abs(dirX[k,])>0+eps) & 
           rowSums( (yk <= YREF+eps & dirY[k,]==0) | abs(dirY[k,])>0+eps)
      }
      if ( is.null(xk) ) 
         xk <-  matrix(X[k,,drop=FALSE], nrow=sum(idx,na.rm=TRUE), 
            ncol=m, byrow=TRUE)
      else
         xk <-  xk[idx,,drop=FALSE]

      if ( is.null(yk) )  
         yk <-  matrix(Y[k,,drop=FALSE], nrow=sum(idx,na.rm=TRUE), 
            ncol=n, byrow=TRUE)
      else 
         yk <-  yk[idx,,drop=FALSE]

      allDir <- cbind( (xk-XREF[idx,])/dirX[rep(k,sum(idx)),], 
                       (YREF[idx,]-yk)/dirY[rep(k,sum(idx)),])
      minDir <- apply(allDir,1,min, na.rm=TRUE)
      # eff[k] <- max(minDir)
	tryCatch ( eff[k] <- max(minDir, na.rm=TRUE), warning = function(w) NULL,
          finaly = (eff[k] <- Inf)  )
      # der er kun gjort plads til een peer per firm
      if ( !is.na(eff[k]) && abs(eff[k]) < Inf  )  { 
         peer[k] <- (1:Kr)[idx][which.max(minDir)]
      } else { 
         peer[k] <- NA 
      }
      # peer[k] <- (1:Kr)[idx][which.max(minDir)]
   }
   # we only need lambda to be able to call peers() to get peers.
   lam <- matrix(0, nrow=K, ncol=Kr)
   rownames(lam) <- rownames(X)
   colnames(lam) <- rNames
   for (k in 1:K)  {
       lam[k, peer[k]] <- 1
   }
   if ( is.null(rownames(lam)) )  {
      if ( length(FRONT.IDX)>0 )  {
         colnames(lam) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
      } else {
         colnames(lam) <- paste("L",1:Kr,sep="")
      }
   } else {
       colnames(lam) <- paste("L",rNames,sep="_")
   }

   e <- list(eff=eff, objval=eff, peers=peer, lambda=lam, RTS="fdh",
             direct=DIRECT, ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
   class(e) <- "Farrell"
   return(e)
}  # if !is.null(DIRECT)




# Find efficiency when compared to each dominating firm; the actual 
# efficiency is then then min/max over the dominating firms.
if ( ORIENTATION=="in" )  {
   for ( k in 1:K )  {
      # For each firm find max(XREF/X) over inputs (rows)
      yk <- matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
      if ( dim(YREF)[1] > 0 )
          idx <- rowSums(yk <= YREF) == n
      else 
          idx <- 0
      if ( sum(idx) == 0 )  {
         # Der er ingen loesning, eff er NA
         eff[k] <- Inf
         peer[k] <- NA
         next
      } 
      maxIn <- apply( XREF[idx,,drop=FALSE] / 
        matrix(X[k,], nrow=sum(idx), ncol=m, byrow=TRUE)
                   , 1, max )
      eff[k] <- min(maxIn)
      # der er kun gjort plads til een peer
      peer[k] <- (1:Kr)[idx][which.min(maxIn)]
   }
} else {
   for ( k in 1:K )  {
      xk <- matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
      idx <- rowSums(xk >= XREF) == m
      if ( sum(idx) == 0 )  {
         # Der er ingen loesning, eff er NA
         eff[k] <- -Inf
         peer[k] <- NA
         next
      } 
      minOut <- apply( YREF[idx,,drop=FALSE] / 
        matrix(Y[k,,drop=FALSE], nrow=sum(idx), ncol=n, byrow=TRUE)
                      , 1, min)
      eff[k] <- max(minOut)
      peer[k] <- (1:Kr)[idx][which.max(minOut)]
   }
}

# we only need lambda to be able to call peers() to get peers.
lam <- matrix(0, nrow=K, ncol=Kr)
rownames(lam) <- rownames(X)
colnames(lam) <- rNames

for (k in 1:K)  {
    # print(peer[k])
    if ( !is.na(peer[k]) )
       lam[k, peer[k]] <- 1
}
if ( is.null(rownames(lam)) )  {
   if ( length(FRONT.IDX)>0 )  {
      colnames(lam) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
   } else {
      colnames(lam) <- paste("L",1:Kr,sep="")
   }
} else {
    colnames(lam) <- paste("L",rNames,sep="_")
}

# print(lam)


e <- list(eff=eff, objval=eff, peers=peer, lambda=lam, RTS="fdh", 
          ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
class(e) <- "Farrell"

return(e)
}  # fdh function

