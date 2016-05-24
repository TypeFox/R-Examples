# The class "market model" and associated methods.

marketmodel <- function(name, n, gamma, diffmatrix, diag = FALSE) {
  # Create a market model object
  #
  # Args:
  #   name: name of the model.
  #   n: number of stocks in the market.
  #   gamma: growth rates.
  #   diffmatrix: diffusion matrix.
  #
  # Returns:
  #   A market model object.
  
  model <- list()
  model$name <- name
  model$n <- n
  model$gamma <- gamma
  model$diffmatrix <- diffmatrix
  model$diag <- diag
  
  class(model) <- "marketmodel"
  return(model)
}


print.marketmodel <- function(x, ...) {
  # Print a market model object.
  #
  # Args:
  #   x: a marketmodel object.  
  
  cat("\nName of model: ", x$name, "\n", sep = "")
  cat("Number of assets: ", x$n, "\n", sep = "")
  
  not <- ""
  if (x$diag == FALSE) {
    not <- "not "
  }
  cat("The diffusion matrix is ", not, "diagonal.\n\n", sep = "")
  cat("Other information (use $):\n")
  cat(names(x), "\n\n")
}



SimMarketModel <- function(model, n.years = 10,
                           frequency = 12,
                           initial.weight = rep(1/model$n, model$n),
                           sub.freq = 1) {
  # Simulate a market model
  #
  # Args:
  #   n.years: number of years.
  #   frequency: number of periods in each year (default is 12, i.e., 
  #              monthly).
  #   initial.weight: initial market weights. The default is equal-weight.
  #   sub.freq: number of sub-periods. This allows for more accurate
  #             simulation.
  #
  # Reurn:
  #   A toymkt object.
  
  n.asset <- model$n
  n.period <- n.years*frequency*sub.freq
  dt <-  1/(frequency*sub.freq)
  sd.BM <- sqrt(dt)
    
  # store prices
  price <- matrix(0, nrow = n.period + 1, ncol = n.asset)
  price[1, ] <- initial.weight
  
  # Simulation
  # Use simple Euler approximation, with time step = dt
  
  # Case 1: Constant coefficients
  if ((class(model$gamma) != "function") &
      (class(model$diffmatrix) != "function")) {
    gamma <- as.numeric(model$gamma)  # it is a vector
    Sigma <- model$diffmatrix  # it is a vector or a matrix
    
    if (model$diag == TRUE) {
      # Case 1.1: Diffusion matrix is diagonal
      Sigma <- as.numeric(Sigma)
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        dlogX <- gamma*dt + Sigma*dW
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    } else {
      # Case 1.2: Diffusion matrix is a matrix
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        dlogX <- as.numeric(gamma*dt + Sigma%*%dW)
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    }    
  }
  
  # Case 2: gamma is a function, Sigma is constant
  if ((class(model$gamma) == "function") &
      (class(model$diffmatrix) != "function")) {
    gamma <- model$gamma
    
    if (model$diag == TRUE) {
      # Case 2.1: Diffusion matrix is diagonal
      Sigma <- as.numeric(model$diffmatrix)
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        gamma.t <- gamma(price[i, ])  # compute gamma(t)
        dlogX <- gamma.t*dt + Sigma*dW
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }      
    } else {
      # Case 2.2: Diffusion matrix is a matrix
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        gamma.t <- gamma(price[i, ])  # compute gamma(t)
        dlogX <- as.numeric(gamma.t*dt + Sigma%*%dW)
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      } 
    }
  }
  
  # Case 3: gamma is constant, Sigma is a function
  if ((class(model$gamma) != "function") &
        (class(model$diffmatrix) == "function")) {
    gamma <- as.numeric(model$gamma)
    Sigma <- model$diffmatrix
    if (model$diag == TRUE) {
      # Case 3.1: Sigma is diagonal
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        Sigma.t <- Sigma(price[i, ])
        dlogX <- gamma*dt + Sigma.t*dW
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    } else {
      # Case 3.2: Sigma is a matrix
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        Sigma.t <- Sigma(price[i, ])
        dlogX <- as.numeric(gamma*dt + Sigma.t%*%dW)
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    }
  }
    
  # Case 4: gamma and Sigma are functions
  if ((class(model$gamma) == "function") &
      (class(model$diffmatrix) == "function")) {
    gamma <- model$gamma
    Sigma <- model$diffmatrix
    
    if (model$diag == TRUE) {
      # Case 4.1: Sigma is diagonal
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        gamma.t <- gamma(price[i, ])
        Sigma.t <- Sigma(price[i, ])
        dlogX <- gamma.t*dt + Sigma.t*dW
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    } else {
      # Case 4.2: Sigma is a matrix
      for (i in 1:n.period) {
        dW <- rnorm(n = n.asset, mean = 0, sd = sd.BM)
        gamma.t <- gamma(price[i, ])
        Sigma.t <- Sigma(price[i, ])
        dlogX <- as.numeric(gamma.t*dt + Sigma.t%*%dW)
        price[i + 1, ] <- price[i, ] * exp(dlogX)
      }
    }
  }
    
  # sample the result at the correct frequency
  price <- price[seq(1, n.period, by = sub.freq), ]
  colnames(price) <- paste("Asset", 1:n.asset, sep = "")
  price <- zoo(price,
               order.by = seq(0, by = 1/frequency, length.out = dim(price)[1]))  
  output <- toymkt(price = price, initial.weight = initial.weight,
                   buy.and.hold = TRUE)
  return(output)
}


AtlasModel <- function(n, g, sigma) {
  # Atlas model
  # Example 5.3.3 of Fernholz (2002)
  #
  # Args:
  #   n: number of stocks
  #   g: growth rate parameter
  #   sigma: common volatility
  #
  # Returns:
  #   A market model object
  
  gamma <- function(X) {
    gamma.t <- rep(0, n)
    gamma.t[X == min(X)] <- n*g
    return(gamma.t)
  }
  diffmatrix <- rep(sigma, n)
  model <- marketmodel(name = "Atlas model", n = n, gamma = gamma,
                       diffmatrix = diffmatrix, diag = TRUE)
  return(model)
}


VolStabModel <- function(n, alpha, sigma) {
  # Volatility-stabalized model
  # Section 12 of Fernholz and Karatzas (2009)
  #
  # Args:
  #   n: number of stocks
  #   alpha: growth rate parameter
  #   sigma: volatility parameter
  #
  # Returns:
  #   A market model object
  
  gamma <- function(X) {
             mu <- X/sum(X)
             return(0.5 * alpha / mu)
           }
  diffmatrix <- function(X) {
                  mu <- X/sum(X)
                  return(sigma / sqrt(mu))
                }
  model <- marketmodel(name = "Volatility-stabilized model", n = n,
                       gamma = gamma, diffmatrix = diffmatrix, diag = TRUE)
}