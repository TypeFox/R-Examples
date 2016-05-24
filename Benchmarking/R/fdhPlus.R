# $Id: fdhPlus.R 149 2015-07-03 14:53:18Z b002961 $


# FDH+.  Beregn foerst CRS efficiency med kun en mulig peer; se om
# lambda ligger inden for graenderne; hvis er den goer er alt ok,
# ellers saet lambda svarende til den overtradte graense og beregn
# efficiency svarende.

dea.fdhPlus <- function(X, Y, ORIENTATION="in",
          XREF=NULL, YREF=NULL, FRONT.IDX=NULL, 
          DIRECT=NULL, param=0.15, TRANSPOSE=FALSE, oKr=0)
{

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop(paste("Unknown value for ORIENTATION:",ORIENTATION))
   }

   if ( ORIENTATION=="graph" )
      stop("ORIENTATION==\"graph\" does not work for fdh+")

   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }

   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in 'dea'")
      XREF <- XREF[,FRONT.IDX, drop=FALSE]
      YREF <- YREF[,FRONT.IDX, drop=FALSE]
   }
   rNames <- colnames(XREF)
   if ( is.null(rNames) & !is.null(colnames(YREF)) )
      rNames <- colnames(YREF)


   # Saet parametrene low og high
   if ( is.null(param) )  {
      param <- .15
   }
   if ( length(param) == 1 )  {
      low <- 1-param
      high <- 1+param
   } else {
      low <- param[1]
      high <- param[2]
   }


   m <- dim(X)[2]  # number of inputs
   n <- dim(Y)[2]  # number of outputs
   K <- dim(X)[1]  # number of units, firms, DMUs
   Kr <- dim(XREF)[1] # number of units,firms in the reference technology


# # Find dominating reference firms for each of the Kr reference firms
# Dom <- list(NA,Kr)  # list of dominating reference firms
# for ( i in 1:Kr )  {
#  # Stoerre for
#  Dom[[i]] <- which(
#   rowSums(matrix(XREF[i,,drop=FALSE],nrow=Kr,ncol=m,byrow=TRUE) > 
#           XREF) == m
#   &
#   rowSums(matrix(YREF[i,,drop=FALSE],nrow=Kr,ncol=n,byrow=TRUE) < 
#           YREF) == n
#   )
# }


E <- rep(NA,K)
peer <- rep(NA, K)
lambda <- matrix(0, nrow=K, ncol=Kr)
lamR <- matrix(NA, nrow=K, ncol=Kr)




# Directional efficiency
if ( !is.null(DIRECT) )  {

stop("Directional efficiency does not yet work for RTS='fdh+'")

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

   for ( k in 1:K )  {
      # For each firm find max(XREF/X) over inputs (rows)
      # Se kun for de firmaer der dominerer firm k i den retning der 
      # ikke aendres
      xk <-  NULL  # matrix(X[k,,drop=FALSE], nrow=Kr, ncol=m, byrow=TRUE)
      yk <-  NULL  # matrix(Y[k,,drop=FALSE], nrow=Kr, ncol=n, byrow=TRUE)
      if ( ORIENTATION=="in" )  {
         yk <- matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
         idx <- rowSums(yk <= YREF) == n
      } else if ( ORIENTATION=="out" )  {
         xk <-  matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
         idx <- rowSums(xk >= XREF) == m
      } else if ( ORIENTATION=="in-out" )  {
         xk <-  matrix(X[k,], nrow=Kr, ncol=m, byrow=TRUE)
         yk <-  matrix(Y[k,], nrow=Kr, ncol=n, byrow=TRUE)
         idx <- rowSums(xk >= XREF) == m & rowSums(yk <= YREF) == n
      }
      if ( is.null(xk) ) 
         xk <-  matrix(X[k,,drop=FALSE], nrow=sum(idx), ncol=m, byrow=TRUE)
      else
         xk <-  xk[idx,,drop=FALSE]

      if ( is.null(yk) )  
         yk <-  matrix(Y[k,,drop=FALSE], nrow=sum(idx), ncol=n, byrow=TRUE)
      else 
         yk <-  yk[idx,,drop=FALSE]

      allDir <- cbind( (xk-XREF[idx,])/dirX[rep(k,sum(idx)),], 
                       (YREF[idx,]-yk)/dirY[rep(k,sum(idx)),])
      minDir <- apply(allDir,1,min, na.rm=TRUE)
      eff[k] <- max(minDir)
      # der er kun gjort plads til een peer per firm
      peer[k] <- (1:Kr)[idx][which.max(minDir)]
   }
   # we only need lambda to be able to call peers() to get peers.
   lam <- matrix(0, nrow=K, ncol=Kr)
   for (k in 1:K)  {
       lam[k, peer[k]] <- 1
   }
   if ( !is.null(dimnames(X)[[1]]) )  {
	   names(eff) <- dimnames(X)[[1]]
	}
   e <- list(eff=eff, objval=eff, peers=peer, lambda=lam, RTS="fdh",
             direct=DIRECT, ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
   class(e) <- "Farrell"
   return(e)
}  # if !is.null(DIRECT)




for ( k in 1:K )  {
   # Loeb hver firm der skal beregnes for, igennem
   # xk <- X[k,,drop=FALSE]
   # yk <- Y[k,,drop=FALSE]
   xk <- X[k,]
   yk <- Y[k,]
   ek <- rep(NA,Kr)
   # For alle mulige reference firms
   for ( r in 1:Kr )  {
     xref <- XREF[r,,drop=FALSE]
     yref <- YREF[r,,drop=FALSE] 
     # Drop this firm as reference if it does not dominate the firm at
     # question at the ends of the possible lambda interval.
     # Tager hensyn saa det ogsaa virker med super efficiens
     if ( ORIENTATION=="in" && yk > high*yref )  next
     else if ( ORIENTATION=="out" && xk < low*xref )  next
     etemp <- dea.csrOne(xk,yk, xref, yref, ORIENTATION)  
     #etemp <- dea(xk,yk, RTS="crs", ORIENTATION, XREF=xref, YREF=yref)
     lamR[k,r] <- etemp$lambda
     ek[r] <- etemp$eff
     if ( ORIENTATION=="in" && etemp$lambda < low )  {
           ek[r] <- max( low*xref/xk )
           lamR[k,r] <- low
     } else if ( ORIENTATION=="out" && etemp$lambda > high )  {
           ek[r] <- min( high*yref/yk )
           lamR[k,r] <- high
     }
   }  # for r
   # Lad vaere med at printe warning hvis min/max er Inf/-Inf, 
   # men saet vaerdien direkte i en finaly clause
   if ( ORIENTATION=="in" )  {
      tryCatch ( E[k] <- min(ek, na.rm=TRUE), warning = function(w) NULL, 
          finaly = (E[k] <- Inf) )
   } else {
	tryCatch ( E[k] <- max(ek, na.rm=TRUE), warning = function(w) NULL,
          finaly = (E[k] <- -Inf)  )
   }
   # print(ek)
   if ( !is.na(E[k]) && abs(E[k]) < Inf  )  { 
      peer[k] <- which.min(ek) 
   } else { 
      peer[k] <- NA 
   }
} # for k

for (k in 1:K)  { 
   if (!is.na(peer[k])) { lambda[k,peer[k]] <- lamR[k,peer[k]] }
}

if ( !is.null(dimnames(X)[[1]]) )  {
	names(E) <- dimnames(X)[[1]]
}

obj <- list(eff=E, lambda=lambda,RTS="fdh+",  lamR=lamR, peers=peer, 
            ORIENTATION=ORIENTATION, TRANSPOSE=FALSE)
class(obj) <- "Farrell"

return(obj)

} # function




# dea.csrOne only works for one firm when there is ONE peer for that firm
dea.csrOne <- function(X,Y, XREF, YREF, ORIENTATION="in")  {
   # X and Y are vectors of input and outout for ONE firm
   m <- length(X)
   n <- length(Y)
   Kr <- dim(XREF)[1]

   # alle kombinationer af x og y sae vi kan finde all partielle
   # produktiviteter. Index ix og iy sae vi blot kan bruge y[iy]/x[ix]

   ix <- gl(m,1,m*n)
   iy <- gl(n,m)

   # Produktivity for each input-output combination
   yx0 <- Y[iy]/X[ix]
   ek <- rep(NA,Kr)
   for ( k in 1:Kr )  {
      # Find max relative productivity compared to each potential peer
      yxk <- YREF[k,iy] / XREF[k,ix]
      ek[k] <- max(yx0/yxk)
   }
   if ( ORIENTATION == "in" )  {
      peer <- which.min(ek)
      E <- ek[peer]
      lam <- max(Y/YREF[peer,])
   } else {
      peer <- which.min(1/ek)
      E <- 1/ek[peer]
      lam <- min(X/XREF[peer,])
   }
   # print(paste(h,":: E =",E[h],"; peer =", peer[h],"; lambda =",lam[h]),quote=FALSE)
   
   obj <- list(eff=E, lambda=lam, peer=peer)
   return(obj)
} # dea.csrOne

