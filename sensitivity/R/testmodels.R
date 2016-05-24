# Test models for sensitivity analysis 
#
# Gilles Pujol 2006
# Bertrand Iooss (2016)


# The non-monotonic Sobol g-function (Saltelli 2000)

sobol.fun <- function(X) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j]) / (1 + a[j])
  }
  y
}


# The non-monotonic Ishigami function (Saltelli 2000)

ishigami.fun <- function(X) {
  A <- 7
  B <- 0.1
  sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
}


# The non-monotonic function of Morris (Saltelli 2000)

morris.fun <- function(X) {
  w <- 2 * (X - 0.5)
  w[, c(3, 5, 7)] <- 2 * (1.1 * X[, c(3, 5, 7)] /
                          (X[, c(3, 5, 7)] + .1) - 0.5)
  y <- b0
  for (i in 1 : 20) {
    y <- y + b1[i] * w[, i]
  }
  for (i in 1 : 19) {
    for (j in (i + 1) : 20) {
      y <- y + b2[i, j] * w[, i] * w[, j]
    }
  }
  for (i in 1 : 18) {
    for (j in (i + 1) : 19) {
      for (k in (j + 1) : 20) {
        y <- y + b3[i, j, k] * w[, i] * w[, j] * w[, k]
      }
    }
  }
  for (i in 1 : 17) {
    for (j in (i + 1) : 18) {
      for (k in (j + 1) : 19) {
        for (l in (k+1):20) {
          y <- y + b4[i, j, k, l] * w[, i] * w[, j] * w[, k] * w[, l]
        }
      }
    }
  }
  y
}

b0 <- rnorm(1)
b1 <- c(rep(20, 10), rnorm(10))
b2 <- rbind(cbind(matrix(-15, 6, 6),
                  matrix(rnorm(6 * 14), 6, 14)),
            matrix(rnorm(14 * 20), 14, 20))
b3 <- array(0, c(20, 20, 20))
b3[1 : 5, 1 : 5, 1 : 5] <- - 10
b4 <- array(0, c(20, 20, 20, 20))
b4[1 : 4, 1 : 4, 1 : 4, 1 : 4] <- 5

environment(morris.fun) <- new.env()
assign("b0", b0, envir = environment(morris.fun))
assign("b1", b1, envir = environment(morris.fun))
assign("b2", b2, envir = environment(morris.fun))
assign("b3", b3, envir = environment(morris.fun))
assign("b4", b4, envir = environment(morris.fun))

# Functional toy function: Arctangent temporal function (Auder, 2011)
# B. Auder, Classification et modelisation de sorties fonctionnelles de codes de calcul, 
#   These de l'Universite Paris VI, 2011.
# X: input matrix (in [-7,7]^2)
# q: number of discretization steps of [0,2pi] interval
# output: vector of q values

atantemp.fun <- function(X, q = 100){
  
  n <- dim(X)[[1]]
  t <- (0:(q-1)) * (2*pi) / (q-1)
  
  res <- matrix(0,ncol=q,nrow=n)
  for (i in 1:n) res[i,] <- atan(X[i,1]) * cos(t) + atan(X[i,2]) * sin(t)
  
  return(res)
  
}

# Functional toy function 
# Inspired from Campbell K, McKay M,Williams B. 2006. 
#      Sensitivity analysis when model outputs are functions. 
#      Reliability Engineering and System Safety 91: 1468-1472.
# See also  A. Marrel, B. Iooss, M. Jullien, B. Laurent and E. Volkova. 
#     Global sensitivity analysis for models with spatially dependent outputs. 
#     Environmetrics, 22:383-397, 2011

campbell1D.fun <- function(X,theta=-90:90){
  
  K1=60
  K2=0.002
  
  return(10+X[,1]*exp(-(theta-10*X[,2])^2/(K1*X[,1]^2+X[,3]^2))
         +(X[,2]+X[,4])*exp(K2*X[,1]*(theta)))
  
}


