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
  # Fast regression PLS without missing values
# Internal function
# h: first current component indice
# n: number of X and Y rows
# p: number of monomials (ncol(X))
# q: number of response variables (ncol(Y))
# nc: required number of components
# Return: C,P,T,U, W, RSS,  PRESS, Q2

fastregpls2nomissing <- function(X.old, YY.old,C, P, T, U, W, h, n, p, q, nc)
  {      
  seuil <- 1.0e-12
  RSS <- matrix(nrow = nc + 1, ncol = q)
  RSS[1, ] <- rep(n - 1, q)
  PRESS <- matrix(NA, nc, q)
  Q2 <- matrix(NA, nc, q)
  Ws <- matrix(nrow = p, ncol=nc)
  ## save the genuine data and their mean and std
  mu.y <- attributes(YY.old)$"scaled:center"
  sd.y <- attributes(YY.old)$"scaled:scale"
  YY <- YY.old
  XX <- X.old
  YY.hat <- matrix(NA, nrow=n, ncol=q)

  
  for (hcur in h:nc) {
            u.new <- YY.old[, 1]
            w.old <- rep(1, p)
            repeat {
                w.new <- t(X.old) %*% u.new/sum(u.new^2)
                w.new <- w.new/sqrt(sum(w.new^2))
                t.new <- X.old %*% w.new
                c.new <- t(YY.old) %*% t.new/sum(t.new^2)
                u.new <- YY.old %*% c.new/sum(c.new^2)
                
                w.dif <- w.new - w.old
                w.old <- w.new
                if (sum(w.dif^2) < seuil) 
                  break
            }
            
            p.new <- t(X.old) %*% t.new/sum(t.new^2)
            c.new <- t(YY.old) %*% t.new/sum(t.new^2)
            RSS[hcur + 1, ] <- colSums((YY.old - t.new %*% t(c.new))^2)

            P[hcur, ] <- p.new
            C[hcur, ] <- c.new
            U[hcur, ] <- u.new
            T[hcur, ] <- t.new
            W[hcur, ] <- w.new
            A <- T[1:hcur, , drop=FALSE] %*% t(T[1:hcur, , drop=FALSE])
            ## A is now a matrix (hcur X hcur) 

            A <- solve(A)            
            ## A  is now a matrix (hcur,hcur)

            A <- t(T[1:hcur, , drop=FALSE]) %*% A
            ## A is now a matrix (n,hcur)
            Hii <- rep(0,n)
            for (i in 1:n) {
               Hii[i] = sum(A[i,]*T[1:hcur, i])
            }
             # Note. Use of 'sapply' is longer than the preceding loop
            # zz <- sapply(1:n,
#                         FUN=function(i,x,y,hcur){sum(x[i,]*y[1:hcur,i])}, A, T, hcur)

            
             YY.hat <- t(T[1:hcur, , drop=FALSE]) %*% C[1:hcur,, drop=FALSE]  %*% diag(sd.y, q, q) + matrix(rep(mu.y, each = n), 
                n, q)
            

            
            if (hcur == 1) {
              Ws[, 1] <- W[1, , drop=FALSE]
            } else {
               wt <- diag(1, p)
                for (i in 1:(hcur - 1)) {
                  wt <- wt %*% (diag(1, p) - W[i, ] %o% P[i, ])
                }
               Ws[, hcur ] <- wt %*% W[ hcur, ]
            }
              beta <- Ws[, 1:hcur, drop=FALSE ] %*% C[1:hcur,, drop=FALSE]
              # beta is a matrix (nmono,1)
              YYz <- XX %*% beta
            
            A <- (YY - YYz)
            A <- (A*A) / (1-Hii)^2
            ## A  is now a matrix  (n, ncol(Y))
            
            PRESS[hcur, ] <- colSums(A)
            Q2[hcur, ] <- 1 - PRESS[hcur, ]/RSS[hcur, ]
            
            X.old <- X.old - t.new %*% t(p.new)
            YY.old <- YY.old - t.new %*% t(c.new) 
        } # fin hcur
  
        ret <- list( C=C,P=P,T=T,U=U, W=W, RSS=RSS, PRESS=PRESS, Q2=Q2, Ws=t(Ws), YY.hat=YY.hat, beta = beta)
return(ret)
} # end fastreglps2nomissing
