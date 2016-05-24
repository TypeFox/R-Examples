# Helper functions for GWRR

# Do bi-section section search for phi using LARS
gwl.bw.cv <- function(lb, ub, eps, X, y, S, kernel){
   a <- lb
   b <- ub
   c <- (a+b)/2
   diff <- b - a
   N <- dim(X)[1]
   
   while (diff > eps){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      CV <- gwl.cv.err(a.c, X, y, S, N, kernel)
      RMSE.a.c <- CV$err
      CV <- gwl.cv.err(c.b, X, y, S, N, kernel)
      RMSE.c.b <- CV$err

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
      }

      if (RMSE.a.c > RMSE.c.b){
         a <- a.c
         RMSE.a <- RMSE.a.c
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }

   CV <- gwl.cv.err(lb, X, y, S, N, kernel)
   RMSE.lb <- CV$err
   sol.lb <- CV$sol
   frac.lb <- CV$frac
   bool.lb <- CV$bool
   CV <- gwl.cv.err(ub, X, y, S, N, kernel)
   RMSE.ub <- CV$err
   sol.ub <- CV$sol
   frac.ub <- CV$frac
   bool.ub <- CV$bool
   CV <- gwl.cv.err(c, X, y, S, N, kernel)
   RMSE.c <- CV$err
   sol.c <- CV$sol
   frac.c <- CV$frac
   bool.c <- CV$bool

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- lb
      sol.c <- sol.lb
      RMSE.c <- RMSE.lb
      frac.c <- frac.lb
      bool.c <- bool.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- ub
      sol.c <- sol.ub
      RMSE.c <- RMSE.ub
      frac.c <- frac.ub
      bool.c <- bool.ub
   }
   params <- list(c, sol.c, RMSE.c, frac.c, bool.c)
   names(params) <- c("phi", "sol", "RMSPE", "frac", "bool")
   params
}

# Find the best LARS solution at each location by minimizing CV RMSPE using LARS
gwl.cv.err <- function(phi, X, y, S, N, kernel){
   LARGE <- 1000000
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)
   W.sqrt <- sqrt(W)
   p <- dim(X)[2]
   sol.i <- array(0, N)
   err.i <- array(0, N)
   frac.i <- array(0, N)
   bool.i <- array(0, dim=c(N,p))

   for (i in 1:N){
      W.i <- diag(W.sqrt[i,])
      W.i <- W.i[-i,-i]
      X.W <- W.i %*% X[-i,]
      y.W <- W.i %*% y[-i]
      lars.obj <- lars(X.W,y.W)   # LARS centers x and y and scales x to have length = 1 for each x-square
      betas <- lars.obj$beta

      # Find best LARS solution w.r.t. RMSPE
      sol.n <- dim(betas)[1]
      # 1st solution is all 0's
      err.min <- LARGE
      sol <- 0
      frac <- 1
      bool <- rep(1,p)   # Use OLS/WLS solution as default
      for (k in 2:sol.n){
         yhat <- X[i,] %*% betas[k,]
         err <- (y[i] - yhat)^2
         f <- sum(abs(betas[k,])) / sum(abs(betas[sol.n,]))
         if (err < err.min){
            err.min <- err
            sol <- k
            frac <- f
            bool <- ifelse(betas[k,] == 0, 0, 1)
         }
      }
      sol.i[i] <- sol
      err.i[i] <- err.min
      frac.i[i] <- frac
      bool.i[i,] <- bool
   }
   RMSPE <- sqrt(sum(err.i) / N)
   params <- list(RMSPE, sol.i, frac.i, bool.i)
   names(params) <- c("err", "sol", "frac", "bool")
   params
}

# Estimate GWR beta coefficients given phi and LARS solution
gwl.beta <- function(phi, sol, frac, bool, X, y, S, N, kernel){
   p <- dim(X)[2]   
   betas.i <- array(0, dim=c(p,N))
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)
   W.sqrt <- sqrt(W)

   # Estimate regression coefficients
   for (i in 1:N){
      # Weight X and y by spatial weights
      W.i <- diag(W.sqrt[i,])
      X.W <- W.i %*% X
      y.W <- W.i %*% y
      lars.obj <- lars(X.W, y.W)   # LARS centers x and y and scales x to have length = 1 for each x-square
      betas <- lars.obj$beta
      sol.n <- dim(betas)[1]
      # Match estimates to CV solution as closely as possible
      if (frac[i] == 1){
         # Use last solution
         betas.i[,i] <- betas[sol.n,]
      }else{
         # Find the number of solutions that match CV bool
         match <- 0
         k.match <- sol.n   # Use WLS solution by default
         for (k in 2:sol.n){
            b <- ifelse(betas[k,] == 0, 0, 1)
            if (isTRUE(all.equal(b, bool[i,]))){
               match <- match + 1
               k.match <- k
            }
         }
         if (match == 1){
            betas.i[,i] <- betas[k.match,]
         }else{
            # Match as close to CV fraction as possible
            k.use <- sol.n   # Use WLS solution by default
            diff.min <- 1
            for (k in 2:sol.n){
               f <- sum(abs(betas[k,])) / sum(abs(betas[sol.n,]))
               diff <- abs(frac[i] - f)
               if (diff < diff.min){
                  diff.min <- diff
                  k.use <- k
               }
            }
            betas.i[,i] <- betas[k.use,]
         }
      }
   }
   betas.i
}

# Do bi-section section search for phi and ridge parameter in GWRR
gwrr.bw.rd.cv <- function(phi.lb, phi.ub, eps, lam.lb, lam.ub, lam.eps, X, y, S, kernel){
   a <- phi.lb
   b <- phi.ub
   c <- (a+b)/2
   diff <- b - a
   N <- dim(X)[1]

   while (diff > eps){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      CV <- gwrr.rd.cv(a.c, lam.lb, lam.ub, lam.eps, X, y, S, N, kernel)   # Was Estimate.lambda.bi.LARS.scale
      RMSE.a.c <- CV$RMSPE
      CV <- gwrr.rd.cv(c.b, lam.lb, lam.ub, lam.eps, X, y, S, N, kernel)
      RMSE.c.b <- CV$RMSPE

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
      }else{
         a <- a.c
         RMSE.a <- RMSE.a.c
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }
   CV <- gwrr.rd.cv(phi.lb, lam.lb, lam.ub, lam.eps, X, y, S, N, kernel)
   RMSE.lb <- CV$RMSPE
   lambda.lb <- CV$lambda
   CV <- gwrr.rd.cv(phi.ub, lam.lb, lam.ub, lam.eps, X, y, S, N, kernel)
   RMSE.ub <- CV$RMSPE
   lambda.ub <- CV$lambda
   CV <- gwrr.rd.cv(c, lam.lb, lam.ub, lam.eps, X, y, S, N, kernel)
   RMSE.c <- CV$RMSPE
   lambda <- CV$lambda

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- phi.lb
      RMSE.c <- RMSE.lb
      lambda <- lambda.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- phi.ub
      RMSE.c <- RMSE.ub
      lambda <- lambda.ub
   }
   params <- list(c, lambda, RMSE.c,  RMSE.c * N)
   names(params) <- c("phi", "lambda", "RMSPE", "cv.score")
   params
}

# Do bi-section section search for ridge shrinkage parameter in GWRR
gwrr.rd.cv <- function(phi, lb, ub, eps, X, y, S, N, kernel){
   a <- lb
   b <- ub
   c <- (a+b)/2
   diff <- b - a

   while (diff > eps){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      RMSE.a.c <- gwrr.cv.err(phi, a.c, X, y, S, N, kernel)   # Was CV.err.ridge.LARS.scale
      RMSE.c.b <- gwrr.cv.err(phi, c.b, X, y, S, N, kernel)  

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
      }

      if (RMSE.a.c > RMSE.c.b){
         a <- a.c
         RMSE.a <- RMSE.a.c
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }

   RMSE.lb <- gwrr.cv.err(phi, lb, X, y, S, N, kernel)
   RMSE.ub <- gwrr.cv.err(phi, ub, X, y, S, N, kernel)
   RMSE.c <- gwrr.cv.err(phi, c, X, y, S, N, kernel)

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- lb
      RMSE.c <- RMSE.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- ub
      RMSE.c <- RMSE.ub
   }
   params <- list(c, RMSE.c)
   names(params) <- c("lambda", "RMSPE")
   params
}

# Calculate CV RMSPE for all observations using GWRR with LARS-type scaling
gwrr.cv.err <- function(phi, lambda, X, y, S, N, kernel){
   yhat <- array(0,N)
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)
   W.sqrt <- sqrt(W)
   one <- rep(1, N)
   p <- dim(X)[2]
   I <- diag(1, p)

   for (i in 1:N){
      W.i <- diag(W.sqrt[i,])
      W.i[i,i] <- 0
      X.W <- W.i %*% X
      y.W <- W.i %*% y
      # Scale X and center X and y as is done in LARS
      meanx <- drop(one %*% X.W)/N
      X.W <- scale(X.W, meanx, FALSE)
      normx <- sqrt(drop(one %*% (X.W^2)))
      X.W <- scale(X.W, FALSE, normx)
      y.mu <- mean(y.W)
      y.W <- drop(y.W - y.mu)
      betas <- t(chol2inv(chol(t(X.W) %*% X.W + I * lambda)) %*% t(X.W) %*% y.W)
      betas <- scale(betas, FALSE, normx)[1:p]  # Scale betas by std(x)
      yhat[i] <- X[i,] %*% betas
   }
   sqrt(sum((y - as.vector(yhat))^2) / N)
}

# Estimate GWR beta coefficients given phi using local LARS scaling
gwrr.beta <- function(phi, lambda, X, y, S, N, kernel){
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)
   p <- dim(X)[2]
   
   W.sqrt <- sqrt(W)
   one <- rep(1, N)
   betas.i <- array(0, dim=c(p,N))
   I <- diag(1, p)

   for (i in 1:N){
      W.i <- diag(W.sqrt[i,])
      X.W <- W.i %*% X
      y.W <- W.i %*% y
      # Scale X and center X and y as is done in LARS
      meanx <- drop(one %*% X.W)/N
      X.W <- scale(X.W, meanx, FALSE)
      normx <- sqrt(drop(one %*% (X.W^2)))
      X.W <- scale(X.W, FALSE, normx)
      y.mu <- mean(y.W)
      y.W <- drop(y.W - y.mu)
      betas.i[,i] <- chol2inv(chol(t(X.W) %*% X.W + I * lambda)) %*% t(X.W) %*% y.W
      betas.i[,i] <- t(scale(t(betas.i[,i]), FALSE, normx)[1:p])  # Scale betas by std(x)
   }
   betas.i
}
   
# Do bi-section section search for phi
gwr.bw.cv <- function(lb, ub, eps, X, y, S, kernel){
   a <- lb
   b <- ub
   c <- (a+b)/2
   diff <- b - a
   N <- dim(X)[1]

   while (diff > eps){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      RMSE.a.c <- gwr.cv.err(a.c, X, y, S, N, kernel)
      RMSE.c.b <- gwr.cv.err(c.b, X, y, S, N, kernel)

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
         print(paste("Bandwidth: ", format(b,digits=4), " RMSPE :", format(RMSE.b,digits=4)))
      }

      if (RMSE.a.c > RMSE.c.b){
         a <- a.c
         RMSE.a <- RMSE.a.c
         print(paste("Bandwidth: ", format(a,digits=4), " RMSPE :", format(RMSE.a,digits=4)))
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }

   RMSE.lb <- gwr.cv.err(lb, X, y, S, N, kernel)
   RMSE.ub <- gwr.cv.err(ub, X, y, S, N, kernel)
   RMSE.c <- gwr.cv.err(c, X, y, S, N, kernel)

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- lb
      RMSE.c <- RMSE.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- ub
      RMSE.c <- RMSE.ub
   }
   print(paste("Bandwidth: ", format(c,digits=4), " RMSPE :", format(RMSE.c,digits=4)))
   params <- list(c, RMSE.c, RMSE.c * N)
   names(params) <- c("phi", "RMSPE", "cv.score")
   params
}

# Calculate CV RMSPE for all observations
gwr.cv.err <- function(phi, X, y, S, N, kernel){
   yhat <- array(0,N)
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)

   for (i in 1:N){
      W.i <- diag(W[i,])
      W.i[i,i] <- 0
      beta <- chol2inv(chol(t(X) %*% W.i %*% X)) %*% t(X) %*% W.i %*% y
      yhat[i] <- X[i,] %*% beta
   }
   sqrt(sum((y - as.vector(yhat))^2) / N)
}

# Function to calculate spatial correlation
w.exp <- function(phi, S) exp(-(S/phi))

# Function to calculate spatial correlation
w.gauss <- function(phi, S) exp(-(S/phi)^2)

# Estimate GWR beta coefficients given bandwidth
gwr.beta <- function(phi, X, y, S, N, kernel){
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)

   p <- dim(X)[2]
   betas <- array(0, dim=c(p,N))

   # Estimate regression coefficients
   for (i in 1:N){
      W.i <- diag(W[i,])
      betas[,i] <- t(chol2inv(chol(t(X) %*% W.i %*% X)) %*% t(X) %*% W.i %*% y)
   }
   betas
}

gwr.yhat <- function(betas, X){
   N <- dim(X)[1]
   yhat <- array(0, N)
   for (i in 1:N) yhat[i] <- X[i,] %*% betas[,i]
   yhat
}

gwr.rmse <- function(y, yhat){
   N <- length(y)
   sqrt(sum((y - yhat)^2) / N)
}

gwr.rsquare <- function(y, yhat){
   rs1 <- (cor(y, yhat))^2
   SS.err <- sum((y - yhat)^2)
   ybar <- mean(y)
   SS.tot <- sum((y - ybar)^2)
   rs2 <- 1 - (SS.err/SS.tot)
   rs2
}

# Calculate the variance decomposition proportions and scaled condition indexes at one location
local.vdp <- function(W.X, N){
   # Scale input matrix W.X by norm first
   one <- rep(1,N)
   normx <- sqrt(drop(one %*% (W.X^2)))
   W.X <- scale(W.X, FALSE, normx)

   # svd returns: P is vector of singular values,
   # U is matrix of left singular vectors, V is matrix of right singular vectors
   svd.list <- svd(W.X)
   P <- svd.list$d
   #U <- svd.list$u
   V <- svd.list$v
   p <- length(P)
   phi.ij <- array(0, dim = c(p,p))
   phi.i <- array(0,p)
   pi.ij <- array(0, dim = c(p,p))
   k.X <- array(0, dim=c(p,1))

   #Use the X'X = VS^2V' result
   for (j in 1:p){
      phi.ij[,j] <- V[,j]^2 / P[j]^2
   }

   phi.i <- apply(phi.ij, 1, sum)

   for (i in 1:p){
      pi.ij[,i] <- phi.ij[i,] / phi.i[i]
   }

   for (i in 1:p){
      k.X[i] <- P[1] / P[i]
   }

   params <- list(k.X, pi.ij)
   names(params) <- c("k.X","pi.ij")
   params
}
