#' Computation Of The Sufficient Statistics 
#' 
#' @description Computation of U and V, the two sufficient statistics of the likelihood of the mixed SDE
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma a(X_j(t)) dW_j(t)}.
#' @param X matrix of the M trajectories.
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param fixed fixed effects in the drift: value of the fixed effect when there is only one random effect, 0 otherwise.
#' @param times times vector of observation times.
#' @return
#' \item{U}{vector of the M statistics U(Tend)}
#' \item{V}{list of the M matrices V(Tend)}
#' @references See Bidimensional random effect estimation in mixed stochastic differential model, C. Dion and V. Genon-Catalot, \emph{Stochastic Inference for Stochastic Processes 2015, Springer Netherlands} \bold{1--28}



#' @details 
#' Computation of U and V, the two sufficient statistics of the likelihood of the mixed SDE
#' \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma a(X_j(t)) dW_j(t) = (\alpha_j, \beta_j)b(X_j(t))dt + \sigma a(X_j(t)) dW_j(t)} with \eqn{b(x)=(1,-x)^t}:
#' 
#' U : \eqn{U(Tend) = \int_0^{Tend} b(X(s))/a^2(X(s))dX(s) }
#'  
#' V : \eqn{V(Tend) = \int_0^{Tend} b(X(s))^2/a^2(X(s))ds }





UV <- function(X, model, random, fixed, times) {
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    Xm <- X[, -K]  
    delta <- diff(times)  
    Tend <- times[K]
    
    if (sum(random) > 2) {
        U <- matrix(0, 2, M)
        
        V <- as.list(1:M)
        b <- as.list(1:M)
        
        Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))  #Int1 <- apply(Xm * delta, 1, sum)
        
        if (model == "OU") {
            
            Int2 <- rowSums(Xm^2 * matrix(delta, M, length(delta)))
            
            
            for (j in 1:M) {
                b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, fixed, random), 2, K)  # 2xK  matrix
                
                U[, j] <- rowSums((b[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 1)]), 2, K - 1, byrow = TRUE)))
                
                V[[j]] <- matrix(c(Tend, -Int1[j], -Int1[j], Int2[j]), 2, 2)
            }
            
        }
        
        if (model == "CIR") {
            
            Int3 <- rowSums(1/Xm * matrix(delta, M, length(delta)))
            
            b <- as.list(1:M)
            bsig <- as.list(1:M)
            
            for (j in 1:M) {
                b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, fixed, random), 2, K)
                
                bsig[[j]] <- matrix(-1, 2, K)
                bsig[[j]][1, ] <- 1/X[j, ]
                
                U[, j] = rowSums((bsig[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 1)]), 2, K - 1, byrow = TRUE)))

                V[[j]] <- matrix(c(Int3[j], -Tend, -Tend, Int1[j]), 2, 2)
                
            }
        }
    }
    if (sum(random) == 2) {
        
        U <- numeric(M)
        V <- numeric(M)
        
        Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))
          
        if (model == "OU") {
            bxj <- apply(X, c(1, 2), bx, fixed, random)
            
             U <- rowSums(bxj[, 1:(K - 1)] * (X[ , 2:K] - X[ , 1:(K - 1)]))   + fixed * Int1
             V <- bxj[ , 1:(K - 1)]^2 %*% delta
        }
        
        if (model == "CIR") {
            
            U <- -(X[, K] - X[, 1]) + fixed * Tend
            V <- Int1
        }
        
        
    }
    
    if (sum(random) == 1) {
        U <- numeric(M)
        V <- numeric(M)
        
        
        if (model == "OU") {
            Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))
           
            U <- (X[, K] - X[, 1]) + fixed * Int1
            
            V <- rep(Tend, M)
        }
        
        if (model == "CIR") {
            
            Int2 <- rowSums((1/Xm) * matrix(delta, M, length(delta)))

            for (j in 1:M) {
                U[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)]) * (1/X[j, 2:K])) + fixed * Tend
            }
            
            V <- Int2
        }
    }
    return(list(U = U, V = V))
} 
