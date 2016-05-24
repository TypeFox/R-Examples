###################################################################
# plspolychaos R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78350 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/plspolychaos
#
# This file is part of plspolychaos R package.
# This is free software: you can redistribute it and/or modify
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
# This file is a modified version of the one of the sivipm package
###################################################################
# regression pls2
# @title regression pls multivarie without missing
# @regpls2
# @param Y outputs data.frame 
# @param dataX.exp inputs data.frame (expanded polynomial) 
# @param nc number of components
# @return ret, R2y, Q2, Q2cum, PRESS
################################################################
regpls2 <- function(Y, dataX.exp, nc = 2) {
    X <- as.matrix(dataX.exp)
    n <- nrow(X)
    p <- ncol(X)
    
    YY <- as.matrix(Y)
    q <- ncol(YY)
    
    
    # centrage reduction
    X.old <- scale(X)
    YY.old <- scale(YY)
    Xx <- X.old
    YYy <- YY.old
    
    T <- matrix(0, nrow = nc, ncol = n)
    W <- matrix(0, nrow = nc, ncol = p)
    
    P <- matrix(0, nrow = nc, ncol = p)
    C <- matrix(0, nrow = nc, ncol = q)
    U <- matrix(0, nrow = nc, ncol = n)
    
    
    # Algo without missing values
        ret <- fastregpls2nomissing(X.old, YY.old, C, P, T, U, W,  n, p, q, nc)
    

      # fonction pour apply
    onek <- function( T, X, i.exist) {
      return(cor(T[i.exist], X[i.exist]))
    }
    
    onej <- function(X, nc, T) {
      # fonction pour apply
      i.exist <- which(complete.cases(X))
      cor.tx <- apply(T, 1, onek, X, i.exist)
      return(cor.tx)
    } # fin onej
    cor.tx <- apply(X, 2, onej, nc, ret$T)

##   cor.tx <- matrix(nrow = nc, ncol = p)
##     for (j in 1:p) {
##         i.exist <- which(complete.cases(X[, j]))
##         # 1st solution
##         for (k in 1:nc) {
##             cor.tx[k, j] <- cor(ret$T[k, i.exist], X[i.exist, j])
##         }
##     }
    
    
    R2x <- cor.tx^2
    Rdx <- rowMeans(R2x)
    
    cor.ty <- apply(YY, 2, onej, nc, ret$T)
##   cor.ty <- matrix(nrow = nc, ncol = q)
##     for (j in 1:q) {
##         i.exist <- which(complete.cases(YY[, j]))
##         # 1st solution
##         for (k in 1:nc) {
##             cor.ty[k, j] <- cor(ret$T[k, i.exist], YY[i.exist, j])
##         }
##     }
    
    
    R2y <- cor.ty^2
    Rdy <- rowMeans(R2y)
    
    
    Rd.mat <- matrix(0, nc, nc)
    for (j in 1:nc) {
        Rd.mat[1:j, j] <- Rdy[1:j]
    }
    
    # predictions
    
    mu.x <- attributes(Xx)$"scaled:center"
    sd.x <- attributes(Xx)$"scaled:scale"
    mu.y <- attributes(YYy)$"scaled:center"
    sd.y <- attributes(YYy)$"scaled:scale"
    
# Not used:    X.hat <- t(ret$T) %*% ret$P %*% diag(sd.x, p, p) + matrix(rep(mu.x, each = n),         n, p)
    # Not used: Dx <- sqrt(rowSums((X - X.hat)^2))
    
    if (is.null(ret$YY.hat)) {
        ret$YY.hat <- t(ret$T) %*% ret$C %*% diag(sd.y, q, q) + matrix(rep(mu.y, 
            each = n), n, q)
    }
    # Not used: Dy <- sqrt(rowSums((YY - YY.hat)^2))
    
    
    # Q2cumule
    
    Q2cum <- rep(0, nc)
    Q2ckh <- matrix(NA, nrow = nc, ncol = ncol(YY))
    for (h in 1:nc) {
        a <- matrix(ret$PRESS[1:h, ]/ret$RSS[1:h, ], ncol = ncol(YY))
        Q2ckh[h, ] <- 1 - apply(a, 2, prod)
        
        Q2cum[h] <- 1 - prod(rowSums(ret$PRESS)[1:h]/rowSums(ret$RSS)[1:h])
        if (h > 1) {
            Q2cum[h] <- max(Q2cum[h], Q2cum[h - 1])
        }
        
    }
    
    
    
    # The dimnames
    labnc <- paste("c", 1:nc, sep = "")
    if (is.null(colnames(YY))) 
        colnames(YY) <- paste("Y", 1:ncol(YY), sep = "")
    if (is.null(colnames(X))) 
        colnames(X) <- paste("X", 1:ncol(X), sep = "")
    
    
    dimnames(ret$Ws) <- list(paste("w*", 1:nc, sep = ""), colnames(X))
    dimnames(cor.tx) <- list(paste("t", 1:nc, sep = ""), colnames(X))
    dimnames(cor.ty) <- list(paste("t", 1:nc, sep = ""), colnames(YY))
# Not used   dimnames(X.hat) <- list(1:n, colnames(X))
    dimnames(ret$YY.hat) <- list(1:n, colnames(YY))
    rownames(R2y) <- labnc
    
    retour <- list(ret = ret, mweights = as.data.frame(ret$Ws), cor.tx = cor.tx, 
        cor.ty = cor.ty, mu.x = mu.x, sd.x = sd.x, mu.y = mu.y, sd.y = sd.y,
                   #  Not used x.hat = X.hat, 
        y.hat = ret$YY.hat, R2x = R2x, R2y = R2y, PRESS = ret$PRESS)
    
    
    # RSS has nc+1 rows
    dimnames(ret$RSS) <- list(paste("c", 1:nrow(ret$RSS), sep = ""), colnames(YY))
    retour$RSS <- ret$RSS
    retour$PRESS <- ret$PRESS
    dimnames(retour$PRESS) <- list(labnc, colnames(YY))
    colnames(ret$Q2) <- colnames(YY)
    retour$Q2 <- ret$Q2
    # Q2 less than zero are set to zero
    retour$Q2[retour$Q2 < 0] <- 0
    rownames(retour$Q2) <- labnc
    retour$Q2cum <- cbind(Q2ckh, Q2cum)
    colQ2 <- paste("Q2-", colnames(retour$PRESS), sep = "")
    colnames(retour$Q2cum) <- c(colQ2, "total-Q2cum")
    rownames(retour$Q2cum) <- labnc
    
    return(retour)
    
    
    
}  # end regpls2 
####################################################
####################################################
## calcbeta compute the beta
calcbeta <- function(W, P, C, nc, p) {
    
    # p is the number of monomials
    Ws <- matrix(nrow = nc, ncol = p)
    Ws[1, ] <- W[1, ]
    
    if (nc >= 2) 
        {
            for (h in 2:nc) {
                wt <- diag(1, p)
                for (i in 1:(h - 1)) {
                  wt <- wt %*% (diag(1, p) - W[i, ] %o% P[i, ])
                }
                Ws[h, ] <- wt %*% W[h, ]
            }  # fin h
        }  # fin (nc >2 )
    beta <- t(Ws) %*% C
    return(beta)
}
####################################################
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
    betaNat <- matrix(nrow = nmono, ncol = nrep)
    
    for (irep in 1:nrep) {
        som <- 0
        for (imono in 1:nmono) {
            fact <- sd.y[irep]/sd.x[imono]
            betaNat[imono, irep] <- beta[imono, irep] * fact
            som <- som + betaNat[imono, irep] * mu.x[imono]
        }
        betaNat0[irep] <- mean(YY[, irep], na.rm = TRUE) - som
    }
    return(list(betaNat = betaNat, betaNat0 = betaNat0))
}  # end calcbetaNat
###################################################################
  # Fast regression PLS without missing values
# Internal function
# h: first current component indice
# n: number of X and Y rows
# p: number of monomials (ncol(X))
# q: number of response variables (ncol(Y))
# nc: required number of components
# Return: C,P,T,U, W, RSS,  PRESS, Q2, beta


fastregpls2nomissing <- function(X.old, YY.old, C, P, T, U, W,  n, p, q, nc) {
    seuil <- 1e-12
    RSS <- matrix(0, nrow = nc + 1, ncol = q)
    RSS[1, ] <- rep(n - 1, q)
    PRESS <- matrix(NA, nc, q)
    Q2 <- matrix(NA, nc, q)
    Ws <- matrix(nrow = p, ncol = nc)
    ## save the genuine data and their mean and std
    mu.y <- attributes(YY.old)$"scaled:center"
    sd.y <- attributes(YY.old)$"scaled:scale"
    YY <- YY.old
    XX <- X.old
    YY.hat <- matrix(NA, nrow = n, ncol = q)
    beta <- matrix(NA, nrow = p, ncol = nc)

    # Init 
    w.dif <- rep(0, p)
    RSS <- matrix(0, nrow = nc + 1, ncol = q)
    RSS[1, ] <- rep(n - 1, q)
    u.new <- YY.old[, 1]
    somunew <- sum(u.new^2)
    w.old <- rep(1, p)
    leRSS <- 0
    # fin init 

    
    for (hcur in 1:nc) {
        retC <- .C("boucleregpls", as.integer(n),
                   as.integer(p),
                   Xold=as.double(X.old),
                   YYold=as.double(YY.old),
                   W=as.double(W[hcur, ]),
                   as.double(w.old),
                   C=as.double(C[hcur,]),
                   U=as.double(u.new),
                   T=as.double(T[hcur, ]),
                   P=as.double(P[hcur, ]),
                   as.double(w.dif),
                   as.double(somunew),
                   RSS=as.double(leRSS))
        RSS[hcur + 1, ] <- retC$RSS

        ## Maj des resultats
        P[hcur, ] <- retC$P[1:p]
        C[hcur, ] <- retC$C[1:q]
        U[hcur, ] <- retC$U[1:n]
        T[hcur, ] <- retC$T[1:n]
        W[hcur, ] <-retC$W[1:p]
## Le C a modifie X.old et YY.old
 #      X.old <- X.old - t.new %*% t(p.new)
        X.old <- matrix(retC$Xold, n, p)
#      YY.old <- YY.old - t.new %*% t(c.new)
        YY.old <- matrix(retC$YYold, n, q)
        
      
      A <- T[1:hcur, , drop = FALSE] %*% t(T[1:hcur, , drop = FALSE])
        ## A is now a matrix (hcur X hcur)
        A <- solve(A)
        ## A is now a matrix (hcur,hcur)
        
#         Hii <- rep(0, n)
#         for (i in 1:n) {
#             Hii[i] = sum(A[i, ] * T[1:hcur, i])
#         }
        Hii <- rep(0, n)
        retC <- .C("calcHii", as.integer(nc),
           as.integer(hcur), as.integer(n),
           as.double(T), as.double(A),
           Hii=as.double(Hii))
        
        YY.hat <- t(T[1:hcur, , drop = FALSE]) %*% C[1:hcur, , drop = FALSE] %*% 
            diag(sd.y, q, q) + matrix(rep(mu.y, each = n), n, q)
        
        
        
        if (hcur == 1) {
            Ws[, 1] <- W[1, , drop = FALSE]
        } else {
            wt <- diag(1, p)
            for (i in 1:(hcur - 1)) {
                wt <- wt %*% (diag(1, p) - W[i, ] %o% P[i, ])
            }
            Ws[, hcur] <- wt %*% W[hcur, ]
        }
        beta[, hcur] <- Ws[, 1:hcur, drop = FALSE] %*% C[1:hcur, , drop = FALSE]
        
        # beta is a matrix (nmono,1)
        YYz <- XX %*% beta[, hcur, drop = FALSE]
        
        A <- (YY - YYz)
        A <- (A * A)/(1 - retC$Hii)^2
        ## A is now a matrix (n, ncol(Y))
        
        PRESS[hcur, ] <- colSums(A)
        Q2[hcur, ] <- 1 - PRESS[hcur, ]/RSS[hcur, ]
        
    
        

    }  # fin hcur
    
    ret <- list(C = C, P = P, T = T, U = U, W = W, RSS = RSS, PRESS = PRESS, Q2 = Q2, 
        Ws = t(Ws), YY.hat = YY.hat, beta = beta)
    return(ret)
}  # end fastreglps2nomissing  
