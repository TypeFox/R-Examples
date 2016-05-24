# 13) Compute IF of RIV.

# Inputs:
# (L,V) is the location/scatter multivariate estimator of a data matrix z.
# IFL is the IF of the location of the n data points, a (rxn) matrix
# IFV  is the IF of the scatter of the n data points, a (rxrxn) array

# Remarks:
# k=number of endogenous variables
# r=number of variables in the dataset

# The intercept=F is needed for the M/S estimator when dummies are present

IF.RIV <- function(L, V, IFL, IFV, kend, k, n, intercept) {
    r <- length(L)
    p <- r - k - 1      # number of endogenous and exogenous covariates
    
    # extract elements needed for IF.riv
    Vm <- matrix(V[(kend + 1):(nrow(V) - 1), -((kend + 1):(kend + k))],
                 nrow = nrow(V) - kend - 1)
    Swx <- Vm[, -ncol(Vm)]
    Sxw <- t(Swx)
    Sww <- matrix(V[(kend + 1):(nrow(V) - 1), (kend + 1):(nrow(V) - 1)],
                  nrow = nrow(V) - kend - 1)
    Swy <- Vm[, ncol(Vm)]
    Lm <- L[-((kend + 1):(kend + k))]
    Mx <- Lm[1:(length(Lm) - 1)]
    My <- Lm[length(Lm)]

    part1 <- Sxw %*% solve(Sww) %*% Swx
    b1 <- solve(part1) %*% Sxw %*% solve(Sww) %*% Swy
    
    IFL <- IFL[-((k + 1):(kend + k)), ]
    IF.Mx <- IFL[1:(length(Lm) - 1), ]
    IF.My <- IFL[length(Lm), ]
    
    # Function updated by Beat for the non-exactly identified case
    Swwinv <- solve(Sww)
    A <- Sxw %*% Swwinv %*% Swx
    B <- Sxw %*% Swwinv %*% Swy
    
    IF.Swx <- IFV[(kend + 1):(r - 1), -c((kend + 1):(kend + k), r), ]
    IF.Sxw <- IFV[-c((kend + 1):(kend + k), r), (kend + 1):(r - 1), ]
    if (p + k > 2) {
      IF.Sww <- IFV[(kend + 1):(r - 1), (kend + 1):(r - 1), ]
      IF.Swy <- IFV[(kend + 1):(r - 1), r, ]
      IF.prov <- array(NA, dim = c(p - kend + k, 1, n))
      for (i in 1:n) IF.prov[, , i] <- t(IF.Swy[, i])
      IF.Swy <- IF.prov
      
      if (p == 1) {
        IF.Xprov <- array(NA, dim = c(p, k, n))
        IF.Wprov <- array(NA, dim = c(k, p, n))
        for (i in 1:n) {
          IF.Xprov[, , i] <- t(IF.Sxw[, i])
          IF.Wprov[, , i] <- IF.Swx[, i]
        }
        IF.Swx <- IF.Wprov
        IF.Sxw <- IF.Xprov
      }
      
      IF.A1 <- multarray(IF.Sxw, solve(Sww) %*% Swx)
      IF.A2a <- multarray(IF.Sww, solve(Sww) %*% Swx)
      IF.A2 <- -tmultarray(IF.A2a, Sxw %*% solve(Sww))
      
      IF.A3 <- IF.A1
      for (i in 1:n) IF.A3[, , i] <- t(IF.A1[, , i])
      
      IF.A <- IF.A1 + IF.A2 + IF.A3
      IF.solveA <- multarray(IF.A, solve(A) %*% B)
      IF.solveA <- -tmultarray(IF.solveA, solve(A))
      
      IF.B1 <- multarray(IF.Sxw, solve(Sww) %*% Swy)
      IF.B2a <- multarray(IF.Sww, solve(Sww) %*% Swy)
      IF.B2 <- -tmultarray(IF.B2a, Sxw %*% solve(Sww))
      IF.B3 <- tmultarray(IF.Swy, Sxw %*% solve(Sww))
      
      IF.B <- IF.B1 + IF.B2 + IF.B3
      IF.B <- tmultarray(IF.B, solve(A))
      
      IF.b1.array <- IF.solveA + IF.B
      
      IF.b1 <- t(matrix(IF.b1.array, ncol = p, byrow = TRUE))
      if (intercept) {
        term1 <- tmultarray(IF.b1.array, t(Mx))
        term1 <- c(term1[1, 1, ])

        term2 <- as.vector(t(b1) %*% IF.Mx)
        IF.b0 <- matrix(IF.My - term1 - term2, nrow=1)
      }
    } else {
      IF.Sww <- IFV[(kend + 1):(kend + k), (kend + 1):(kend + k), ]
      IF.Swy <- IFV[(kend + 1):(kend + k), r, ]
      
      IF.A1 <- IF.Sxw * solve(Sww) * Swx
      IF.A2 <- -Sxw * solve(Sww) * IF.Sww * solve(Sww) * Swx
      IF.A3 <- Sxw * solve(Sww) * IF.Swx
      IF.A <- IF.A1 + IF.A2 + IF.A3
      IF.A <- -solve(A) * IF.A * solve(A) * B
      
      IF.B1 <- IF.Sxw * solve(Sww) * Swy
      IF.B2 <- -Sxw * solve(Sww) * IF.Sww * solve(Sww) * Swy
      IF.B3 <- Sxw * solve(Sww) * IF.Swy
      
      IF.B <- IF.B1 + IF.B2 + IF.B3
      IF.B <- solve(A) * IF.B
      
      IF.b1 <- IF.A + IF.B
      if (intercept) {
        IF.b0 <- IF.My - IF.b1 * Mx - b1 * IF.Mx
      }
    }

    if (intercept) {
      IFriv <- rbind(IF.b0, IF.b1)
      rownames(IFriv) <- c('Intercept', names(Mx))
    } else {
      IFriv <- rbind(IF.b1)
      rownames(IFriv) <- names(Mx)
    }

    IFriv
}
