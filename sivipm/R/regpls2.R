###################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
# regression pls2
# @title regression pls multivarie
# @regpls2
# @param Y outputs data.frame (which may contain missing values)
# @param dataX.exp inputs data.frame (expanded polynomial) (which may contain missing values)
# @param nc number of components
# @param fast (if no missing only): fast procedure
# @param output option to control what is returned.
# @return all scores loadings cross validation weights VIP coefficients
# Q2, Q2cum, PRESS: returned when there are no missing, only 

regpls2 <- function(Y, dataX.exp, nc = 2,
                    na.miss=NA, fast = FALSE, output = FALSE) {
    X <- as.matrix(dataX.exp)
    n <- nrow(X)
    p <- ncol(X)
    
    YY <- as.matrix(Y)
    q <- ncol(YY)

   
    # with/without missing values
    if (is.na(na.miss)) {
    if (any(is.na(X)) || any(is.na(YY))) 
        na.miss <- TRUE else na.miss <- FALSE
  }
    
    # centrage reduction
    X.old <- scale(X)
    YY.old <- scale(YY)
    Xx <- X.old
    YYy <- YY.old 
    
    T <- matrix(nrow = nc, ncol = n)
    W <- matrix(nrow = nc, ncol = p)
    
    P <- matrix(nrow = nc, ncol = p)
    C <- matrix(nrow = nc, ncol = q)
    U <- matrix(nrow = nc, ncol = n)
    
    h <- 1

    # Algo without missing values
    if (!na.miss) {
      if (fast) {
         ret <- fastregpls2nomissing(X.old, YY.old, C, P, T, U, W,  h, n, p, q, nc)
       } else {
      ret <- regpls2nomissing(X.old, YY.old, C, P, T, U, W,  h, n, p, q, nc)
    }
    } # fin na.miss

    # Algo with missing values
    if (na.miss) {
      ret <- regpls2missing(X.old, YY.old,C, P, T, U, W, h, n, p, q, nc)
    } # end missing


    
    if (output == TRUE) 
        {
            x.scores <- as.data.frame(ret$T)
            y.scores <- as.data.frame(ret$U)
            x.loads <- as.data.frame(ret$P)
            y.loads <- as.data.frame(ret$C)
            weights <- as.data.frame(ret$W)
            
            
            
            cor.tx <- matrix(nrow = nc, ncol = p)
            for (j in 1:p) {
                i.exist <- which(complete.cases(X[, j]))
                #  1st solution 
                for (k in 1:nc) {
                  cor.tx[k, j] <- cor(ret$T[k, i.exist], X[i.exist, j])
                }

                 # 2nd solution, not more quick
#                cor.tx[,j] <- apply(ret$T, 1, function(T,X,ind) cor(T[ind], X), X[i.exist, j], i.exist)
                
            }
            R2x <- cor.tx^2
            Rdx <- rowMeans(R2x)
        }  # end output == TRUE
    

    cor.ty <- matrix(nrow = nc, ncol = q)
    for (j in 1:q) {
        i.exist <- which(complete.cases(YY[, j]))
        # 1st solution 
        for (k in 1:nc) {
            cor.ty[k, j] <- cor(ret$T[k, i.exist], YY[i.exist, j])
        }

        # 2nd solution, not more quick
#        cor.ty[,j] <- apply(ret$T, 1, function(X,YY,ind) cor(X[ind], YY), YY[i.exist, j], i.exist)
    }

 
                       
    R2y <- cor.ty^2
    Rdy <- rowMeans(R2y)
    
    
    
    Rd.mat <- matrix(0, nc, nc)
    for (j in 1:nc) {
        Rd.mat[1:j, j] <- Rdy[1:j]
    }
    VIP <- sqrt((t(ret$W)^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), nc, nc))
    dimnames(VIP) <- list(colnames(X), paste("c", 1:nc, sep=""))
    
    if (output == TRUE) 
        {
            VIPind <- matrix(nrow = ncol(YY), ncol = ncol(X))
            dimnames(VIPind) <- list(colnames(YY), colnames(X))

            for (k in 1:ncol(YY)) {
                for (j in 1:nc) {
                  Rd.mat[1:j, j] <- R2y[1:j, k]
                }
                VIPind[k, ] <- sqrt((t(ret$W)^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), nc, 
                  nc))[, nc]
            }
            
            EV <- rbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy))
            
            # predictions
            
            mu.x <- attributes(Xx)$"scaled:center"
            sd.x <- attributes(Xx)$"scaled:scale"
            mu.y <- attributes(YYy)$"scaled:center"
            sd.y <- attributes(YYy)$"scaled:scale"

            X.hat <- t(ret$T) %*% ret$P %*% diag(sd.x, p, p) + matrix(rep(mu.x, each = n), 
                n, p)
            # Not used: Dx <- sqrt(rowSums((X - X.hat)^2)) 
            
            if (is.null(ret$YY.hat)) {
              ret$YY.hat <- t(ret$T) %*% ret$C %*% diag(sd.y, q, q) + matrix(rep(mu.y, each = n), 
                n, q)
            }
            # Not used: Dy <- sqrt(rowSums((YY - YY.hat)^2)) 
            

            
            # coeffs beta
            if (is.null(ret$beta)) {
                 #p is the number of monomials
                ret$Ws <- matrix(nrow = nc, ncol = p)
                ret$Ws[1, ] <- ret$W[1, ]

                if (nc >= 2 ) {
                    for (h in 2:nc) {
                        wt <- diag(1, p)
                        for (i in 1:(h - 1)) {
                            wt <- wt %*% (diag(1, p) - ret$W[i, ] %o% ret$P[i, ])
                        }
                        ret$Ws[h, ] <- wt %*% ret$W[h, ]
                    } # fin h
                } # fin (nc >2 )
                ret$beta <- t(ret$Ws) %*% ret$C
            }    # fin (is.null(ret$beta))
            
            betaNat <- calcbetaNat(ret$beta, YY, X, mu.x, sd.x, sd.y)
            # Q2cumule
            if (!na.miss) 
                {
                  Q2cum <- rep(0, nc)
                  Q2ckh <- matrix(NA, nrow=nc, ncol= ncol(YY))
                  for (h in 1:nc) {
                    a <- matrix(ret$PRESS[1:h,]/ret$RSS[1:h,], ncol= ncol(YY))
                    Q2ckh[h,] <- 1 - apply(a, 2, prod)

                    Q2cum[h] <- 1 - prod(rowSums(ret$PRESS)[1:h]/rowSums(ret$RSS)[1:h])
                    if ( h > 1) {
                      Q2cum[h] <- max(Q2cum[h], Q2cum[h - 1])
                    }
                    
                 }
                  
                 }  # end (!na.miss)
            # noms
            
            # The dimnames
            if (is.null(colnames(YY))) 
                colnames(YY) <- paste("Y", 1:ncol(YY), sep = "")
            if (is.null(colnames(X))) 
                colnames(X) <- paste("X", 1:ncol(X), sep = "")
            
            
            dimnames(weights) <- list(paste("w", 1:nc, sep = ""), colnames(X))
            dimnames(ret$Ws) <- list(paste("w*", 1:nc, sep = ""), colnames(X))
            dimnames(x.scores) <- list(paste("t", 1:nc, sep = ""), 1:n)
            dimnames(x.loads) <- list(paste("p", 1:nc, sep = ""), colnames(X))
            dimnames(y.scores) <- list(paste("u", 1:nc, sep = ""), 1:n)
            dimnames(y.loads) <- list(paste("c", 1:nc, sep = ""), colnames(YY))
            dimnames(cor.tx) <- list(paste("t", 1:nc, sep = ""), colnames(X))
            dimnames(cor.ty) <- list(paste("t", 1:nc, sep = ""), colnames(YY))
            colnames(EV) <- seq(1, ncol(EV))
            dimnames(X.hat) <- list(1:n, colnames(X))
            dimnames(ret$YY.hat) <- list(1:n, colnames(YY))
            dimnames(ret$beta) <- list(colnames(X), colnames(YY))
            dimnames(betaNat$betaNat) <-  dimnames(ret$beta)
            names(betaNat$betaNat0) <- colnames(YY)
            
            retour <- list(betaCR = ret$beta, betaNat= betaNat$betaNat,
                           betaNat0 = betaNat$betaNat0,
                           mweights = as.data.frame(ret$Ws), x.scores = x.scores, 
                x.loadings = x.loads, y.scores = y.scores, y.loadings = y.loads, 
                weights = weights, cor.tx = cor.tx, cor.ty = cor.ty, expvar = EV, 
                VIP = VIP, VIPind = VIPind, x.hat = X.hat, y.hat = ret$YY.hat)
            
            if (!na.miss) {
	        labnc <-  paste("c", 1:nc, sep="")
                  # RSS has nc+1 rows
		dimnames(ret$RSS) <- list(
		  paste("c", 1:nrow(ret$RSS), sep=""), 
		  colnames(YY))
                retour$RSS <- ret$RSS
                retour$PRESS <- ret$PRESS
		dimnames(retour$PRESS) <- list(labnc,  colnames(YY))
                 colnames(ret$Q2) <- colnames(YY)
                retour$Q2 <- ret$Q2
                # Q2 less than zero are set to zero
                retour$Q2[retour$Q2<0] <- 0
                rownames(retour$Q2) <- labnc
                 retour$Q2cum <- cbind(Q2ckh,Q2cum)
                colQ2 <- paste("Q2-", colnames(retour$PRESS), sep="")
                 colnames(retour$Q2cum) <- c(colQ2, "total-Q2cum")
                  rownames(retour$Q2cum) <- labnc
		colnames(retour$expvar) <- labnc
            }
            return(retour)
        }  # end output == TRUE
 else {
   ret <- list(VIP = VIP)
        return(ret)
    }
    
    
    
}  # end regpls2 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Regression PLS without missing values
# Internal function
# Return: C,P,T,U, W, RSS,  PRESS, Q2

regpls2nomissing <- function(X.old, YY.old,C, P, T, U, W, h, n, p, q, nc)
  {      
##        RSS <- matrix(nrow = nc + 1, ncol = q)
##        RSS[1, ] <- rep(n - 1, q)
        # For version C, structures are transposed
        RSS <- matrix(nrow = q, ncol = nc + 1)
        RSS[, 1] <- rep(n - 1, q)

        # Allocation of the working structures for the C version
         w.old <- rep(NA, p)
        w.new <- rep(NA, p)
        w.dif <- rep(NA, p)
        wh.si  <- rep(NA, p)
        wh.siold  <- rep(NA, p)
        wsi.dif  <- rep(NA, p)
         th.si <- rep(NA, n-1)
        uh.si <- rep(NA, n-1)
        ch.si  <- rep(NA,q)
        YYhat.si  <- rep(NA,q)
         pression <- matrix(NA, nrow=n, ncol=q)
        C <- matrix(NA, nrow=q, ncol=nc)
        PRESS <- matrix(NA, nrow=q, ncol=nc)
        Q2 <- matrix(NA, nrow=q, ncol=nc)
        
        
        a=.C("boucle", as.double(X.old), as.double(YY.old),
          as.integer(n), as.integer(p), as.integer(q),
          as.integer(nc),
	    as.double(w.old), as.double(w.new), as.double(w.dif),
          as.double(wh.si), as.double(wh.siold), as.double(wsi.dif),
            as.double(th.si), as.double(uh.si),
	    as.double(ch.si),  as.double(YYhat.si),  
	    pression=as.double(pression),
	    RSS=as.double(RSS), C=as.double(C),
            PRESS=as.double(PRESS), 
	    Q2=as.double(Q2), P=as.double(P),
          W=as.double(W), T=as.double(T), 
	    U=as.double(U), NAOK=TRUE)

        RSS  <-  matrix(a$RSS, nrow=q)
        PRESS <-  matrix(a$PRESS, nrow=q)
        pression <- matrix(a$pression, nrow=n)
        Q2 <- matrix(a$Q2, ncol=nc)
        P <- matrix(a$P, ncol=nc)
        C <- matrix(a$C, ncol=nc)
        U <- matrix(a$U, ncol=nc)
        T <- matrix(a$T, ncol=nc)
        W <- matrix(a$W, ncol=nc)
        
        # All the results are transposed compared to the R version,
        # except 'pression'
        

# VERSION R        
##         
##         PRESS <- matrix(NA, nc, q)
##         Q2 <- matrix(NA, nc, q)
##         
##         repeat {
##             # NIPALS PLS2
##             u.new <- YY.old[, 1]
##             w.old <- rep(1, p)
## 
##             repeat {
##                  w.new <- t(X.old) %*% u.new/sum(u.new^2)
##                 w.new <- w.new/sqrt(sum(w.new^2))
##                 t.new <- X.old %*% w.new
##                 c.new <- t(YY.old) %*% t.new/sum(t.new^2)
##                 u.new <- YY.old %*% c.new/sum(c.new^2)
##                 
##                 w.dif <- w.new - w.old
##                 w.old <- w.new
##                 if (sum(w.dif^2) < 1e-12) 
##                   # if ( iter==100)
##                 break
##                 # iter <- iter + 1
##               } # fin repeat
##             
##             p.new <- t(X.old) %*% t.new/sum(t.new^2)
##             c.new <- t(YY.old) %*% t.new/sum(t.new^2)
##             
##             RSS[h + 1, ] <- colSums((YY.old - t.new %*% t(c.new))^2)
##             
##             pression <- matrix(0, n, q)
## 
##             for (i in 1:n) {
##                 uh.si <- YY.old[-i, 1]
##                 wh.siold <- rep(1, p)
##                 itcv <- 1
##                 repeat {
##                    wh.si <- t(X.old[-i, ]) %*% uh.si/sum(uh.si^2)
##                    wh.si <- wh.si/sqrt(sum(wh.si^2))
##                   th.si <- X.old[-i, ] %*% wh.si
##                   ch.si <- t(YY.old[-i, ]) %*% th.si/sum(th.si^2)
##                  uh.si <- YY.old[-i, ] %*% ch.si/sum(ch.si^2)
##                    wsi.dif <- wh.si - wh.siold
##                  wh.siold <- wh.si
##                  if (sum(wsi.dif^2) < 1e-12) 
##                     break
##                   itcv <- itcv + 1  # not used
##                 }
##                 YYhat.si <- (X.old[i, ] %*% wh.si) %*% t(ch.si)
##                 pression[i, ] <- (YY.old[i, ] - YYhat.si)^2
##             }  # end i
##             
##             PRESS[h, ] <- colSums(pression)
##             Q2[h, ] <- 1 - PRESS[h, ]/RSS[h, ]
##             
## 
##             
##             Xres <- X.old - t.new %*% t(p.new)
##             YYres <- YY.old - t.new %*% t(c.new)
##             X.old <- Xres
##             YY.old <- YYres
##             P[h, ] <- p.new
##             C[h, ] <- c.new
##             U[h, ] <- u.new
##             T[h, ] <- t.new
##             W[h, ] <- w.new
##             
##             if (h == nc) 
##                 break
##             h <- h + 1
##           } # fin repeat

        ret <- list( C=t(C),P=t(P),T=t(T),U=t(U), W=t(W), RSS=t(RSS), PRESS=t(PRESS), Q2=t(Q2))
            
return(ret)
} # end reglps2nomissing
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Regression PLS with missing values
# Internal function
# Return: C,P,T,U, W

regpls2missing <- function(X.old, YY.old, C, P, T, U, W,h, n, p, q, nc)
  {

        
        repeat {
            w.old <- rep(1, p)
            w.new <- rep(0, p)
            t.new <- rep(0, n)
            u.new <- YY.old[, 1]
            c.new <- rep(0, q)
            p.new <- rep(0, p)
            
            repeat {
                for (j in 1:p) {
                  k.exist <- which(complete.cases(u.new))
                  i.exist <- which(complete.cases(X.old[, j]))
                  i.exist <- intersect(k.exist, i.exist)
                  w.new[j] <- sum(X.old[i.exist, j] * u.new[i.exist])/sum(u.new[i.exist]^2)
                }

                w.new <- w.new/sqrt(sum(w.new^2))
                
                for (i in 1:n) {
                  j.exist <- which(complete.cases(t(X.old)[, i]))
                  t.new[i] <- sum(X.old[i, j.exist] * w.new[j.exist])/sum(w.new[j.exist]^2)
                }
                
                
                for (j in 1:q) {
                  i.exist <- which(complete.cases(YY.old[, j]))
                  c.new[j] <- sum(YY.old[i.exist, j] * t.new[i.exist])/sum(t.new[i.exist]^2)
                }
                
                for (i in 1:n) {
                  j.exist <- which(complete.cases(t(YY.old)[, i]))
                  u.new[i] <- sum(YY.old[i, j.exist] * c.new[j.exist])/sum(c.new[j.exist]^2)
                }
                
                
                w.dif <- w.new - w.old
                w.old <- w.new
                
                if (sum(w.dif^2) < 1e-12) {
                  break
                }
            }
            
            for (j in 1:p) {
                i.exist <- which(complete.cases(X.old[, j]))
                p.new[j] <- sum(X.old[i.exist, j] * t.new[i.exist])/sum(t.new[i.exist]^2)
            }
            
            Xres <- X.old - t.new %*% t(p.new)
            YYres <- YY.old - t.new %*% t(c.new)
            X.old <- Xres
            YY.old <- YYres
            P[h, ] <- p.new
            C[h, ] <- c.new
            U[h, ] <- u.new
            
            W[h, ] <- w.new
            T[h, ] <- t.new
            
            
            if (h == nc) 
                break
            h <- h + 1
        }
        ret <- list( C=C,P=P,T=T,U=U, W=W)
return(ret)
} # end reglps2missing
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# calcbetaNat compute the natural beta
# Input
# beta: matrix nmonomes X nreponses
# YY:   matrix nobs X nreponses
# X:    matrix nobs X nmonomes
# Return
# betaNat: matrix nmonomes X nreponses
# betaNat0: vector nrep
# Internal function
calcbetaNat <- function(beta, YY, X, mu.x, sd.x, sd.y) {
  nrep <- ncol(YY)
  nmono <- nrow(beta)
  betaNat0 <- rep(NA, nrep)
  betaNat <- matrix(nrow=nmono, ncol=nrep)

  for (irep in 1:nrep) {
    som <- 0
    for (imono in 1:nmono) {
      fact <- sd.y[irep]/sd.x[imono]
      betaNat[imono, irep] <- beta[imono, irep] * fact
       som <- som + betaNat[imono, irep] * mu.x[imono]
    }
    betaNat0[irep] <- mean(YY[, irep], na.rm=TRUE) - som
  }
  return(list(betaNat=betaNat, betaNat0=betaNat0))
} # end calcbetaNat
