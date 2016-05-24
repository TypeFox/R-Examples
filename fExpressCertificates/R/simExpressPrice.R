# Simuliere Kurse (S(t_1),...,S(t_n)) für 
# alle Bewertungstage (t_1,...,t_n) ausgehend von S(t0) und m(t0) über 
# die multivariate Normalverteilung der Returns
#
# Falls wir für einen Zeitpunkt t0 > 0 simulieren wollen (z.B. 6M nach Emission),
# dann brauchen wir für pfadabhängige Zertifikate auch den bisherigen Verlauf/Minimum bis t0 (m(t0))
# @param m = m(t0) bisheriges Minimum, falls wir ausgehend von einem t0 > 0 simulieren
simExpressPriceMVN <- function(S, m=Inf, X, K, B, T, 
  r, r_d, sigma, mc.loops = 100000, payoffFunction, ...) {
  
    # simulate returns for the T valuation times
    n <- length(T)
    # 
    N <- mc.loops
    
    mu   <- r - r_d
    mean <- (mu - sigma^2/2) * T
    Sigma <- outer(T, T, pmin) * sigma^2
    R <- rmvnorm(n = mc.loops, mean = mean, sigma = Sigma)  
    
    # Matrix (N x n) of simulated prices (N = number of MC iterations and n number of valuation times)
    E <- S * exp(R)
    
    # if necessary also sample from running minimum m_t
    period_min <- matrix(NA, N, n)
    for (i in 1:n) {
      if (i == 1) {
         period_min[, i] <- rBrownianBridgeMinimum(N, t0 = 0, T = T[i], a = 0, b = R[, i], sigma = sigma)
      }
      else {
         period_min[, i] <- rBrownianBridgeMinimum(N, t0 = T[i - 1], T = T[i], a = R[, i - 1], b = R[, i], sigma = sigma)
      }
    }
    #running_min <- t(apply(t(apply(period_min, 1, cummin)), 1, pmin, m)) # (N x n)
    running_min <- matrix(NA, N, n)
    running_min[,1] <- pmin(period_min[,1], m)
    for (i in 2:n) {
      running_min[,i] <- pmin(running_min[,i-1], period_min[,i], m)
    }

    # transform running minimum returns to running minimum prices
    running_min <- S * exp(running_min)  # matrix (N x n) of running mininum m(t_k) = min(S_t, t < t_k)
    period_min <- S * exp(period_min)    # matrix (N x n) of period minimum m(t_{i-1},t_i) zwischen t_{i-1} und t_i
    
    # get redemption times t_i (N x 1) and payoff matrix (N x n)
    if (class(payoffFunction) == "character") {
      res <- eval(call(payoffFunction, S=E, m=running_min, X=X, K=K, B=B))
    } else if (class(payoffFunction) == "function") {
      res <- payoffFunction(S=E, m=running_min, X=X, K=K, B=B)
    } else {
      stop("payoffFunction must be either a function name or a function object")
    }
       
    # redemption times index i for redemption at t_i (N x 1)  
    redemptionIndex <- res$redemptionIndex
            
    # payoff matrix (N x n)
    payoffMatrix <- res$payoffMatrix
    
    # discounted payoffs
    prices <- rowSums(sweep(payoffMatrix, 2, FUN="*", exp(-r * T)))
        
    l = list(payoffs = payoffMatrix,
        underlying = E,
        underlying_period_minimum = period_min,
        underlying_running_minimum = running_min,
        prices = prices, 
        price = mean(prices), 
        redemptionTimes = redemptionIndex, 
        n = n, S = S, X = X, K = K, T = T)
    class(l) = "express.certificate"
    return(l)
}


# Monte Carlo simulation via truncated multivariate normal (TMVN)
#
# Simuliere nur die Kurse am LZE S(t_n) | S(t_1) < X(t_1),...,S(t_{n-1}) < X(t_{n-1}) 
# über die gestutzte multivariate Normalverteilung der konditionierten Returns (Ret(t_1),...,Ret(t_n))
#
# Die Methode funktioniert nur, wenn die erwarteten vorzeitigen Fälligkeiten direkt berechnet werden können:
# K(t_i) * 1_{S(t_i) > X(t_i)} * 1_{ forall j < i : S(t_j) < X(t_j)}   für t_i < t_n
#
simExpressPriceTMVN <- function (S, m=Inf, X, K, B, T, 
  r, r_d, sigma, mc.loops=100000, payoffFunction, ...) 
{
    n <- length(T)
    # 
    N <- mc.loops
    
    # Calculate the stop probabilities p(t_i) = P(0, t_i) = P(S(t_i) > X(t_i) | S(t_1) < X(t_1),...,S(t_{i-1}) < X(t_{i-1})
    p_0i <- calcRedemptionProbabilities(S = S, X = X, T = T, r = r, 
        r_d = r_d, sigma = sigma)$stop.probs
        
    Sigma <- sigma^2 * outer(T, T, pmin)
    mean <- ((r - r_d) - sigma^2/2) * T
    a <- rep(-Inf, n)
    b <- c(log(X[1:(n-1)]) - log(S), Inf)
    
    # Simuliere (Ret(t_1),...,Ret(t_n)) | S(t_1) < X(t_1),...,S(t_{n-1}) < X(t_{n-1})
    R <- rtmvnorm(n = mc.loops, mean = mean, sigma = Sigma, lower=a, upper=b)  
    
    # Matrix (N x n) of simulated prices (N = number of MC iterations and n number of valuation times)
    E <- S * exp(R)
    
    # if necessary also sample from running minimum m_t
    period_min <- matrix(NA, N, n)
    for (i in 1:n) {
      if (i == 1) {
         period_min[, i] <- rBrownianBridgeMinimum(N, t0 = 0, T = T[i], a = 0, b = R[, i], sigma = sigma)
      }
      else {
         period_min[, i] <- rBrownianBridgeMinimum(N, t0 = T[i - 1], T = T[i], a = R[, i - 1], b = R[, i], sigma = sigma)
      }
    }
    #running_min <- t(apply(t(apply(period_min, 1, cummin)), 1, pmin, m)) # (N x n)
    running_min <- matrix(NA, N, n)
    running_min[,1] <- pmin(period_min[,1], m)
    for (i in 2:n) {
      running_min[,i] <- pmin(running_min[,i-1], period_min[,i], m)
    }

    # transform running minimum returns to running minimum prices
    running_min <- S * exp(running_min)  # matrix (N x n) of running mininum m(t_k) = min(S_t, t < t_k)
    period_min <- S * exp(period_min)    # matrix (N x n) of period minimum m(t_{i-1},t_i) zwischen t_{i-1} und t_i
    
    # get redemption times and payoff matrix
    if (class(payoffFunction) == "character") {
      res <- eval(call(payoffFunction, S=E, m=running_min, X=X, K=K, B=B))
    } else if (class(payoffFunction) == "function") {
      res <- payoffFunction(S=E, m=running_min, X=X, K=K, B=B)
    } else {
      stop("payoffFunction must be either a function name or a function object")
    }

    # redemption times   
    redemptionIndex <- res$redemptionIndex           # 1000x an t3 per Konstruktion, da wir nur am LZE simulieren ;-)
        
    # payoffs at maturity h(S)
    payoffs <- res$payoffMatrix[,n]                  # N x 1

    # discounted payoffs
    prices <- sum(exp(-r * T[1:(n-1)]) * p_0i[1:(n-1)] * K[1:(n-1)]) + p_0i[n] * payoffs * exp(-r * T[n])
     
    l <- list(
        method="MC_TMVN",
        redemption_probabilities=p_0i,
        payoffsMaturity = payoffs,
        underlying = E,
        underlying_period_minimum = period_min,
        underlying_running_minimum = running_min,
        prices = prices, 
        price = mean(prices),
        n = n, S = S, X = X, K = K, T = T)
    class(l) = "express.certificate"
    return(l)
}