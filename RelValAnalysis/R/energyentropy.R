GetGroupWeight <- function(pi, group.index) {
  # Calculate group weights and conditional weights
  #
  # Args:
  #   pi: a portfolio vector
  #   group.index: a list whose components from a partition of 1:n, where
  #                n is the length of pi.
  #
  # Returns:
  #   A list:
  #     group.weight: the "first level" weights
  #     cond.weight: the second level weight (another list)

  pi <- as.numeric(pi)
  
  # Compute group weights
  m <- length(group.index)  
  lambda <- numeric(m)
  for (k in 1:m) {
    lambda[k] <- sum(pi[group.index[[k]]])
  }
  
  # Compute conditional weights
  pi.sub <- list()
  for (k in 1:m) {
    if (lambda[k] > 0) {
      pi.sub[[k]] <- pi[group.index[[k]]] / lambda[k]
    } else {
      n.group.k <- length(group.index[[k]])
      pi.sub[[k]] <- rep(1/n.group.k, n.group.k)
    }
  }
  

  return(list(group.weight = lambda, cond.weight = pi.sub))
}

RelativeEntropy <- function(p, q, group.index = NULL) {
  # Relative entropy H(pi|mu)
  #
  # Args:
  #   pi: probability vector
  #   mu: probability vector
  #   group.index: if given, the chain rule will be performed
  #
  # Returns:
  #   A non-negative number or Inf (if group.index = NULL)
  #   A numeric vector (otherwise)
  
  if (is.null(group.index)) {
    # Case 1: group.index not given
    p <- as.numeric(p)
    q <- as.numeric(q)
    return(sum(p*log(p), na.rm = TRUE) - sum(p*log(q), na.rm = TRUE))
  } else {
    # Case 2: Chain rule
    p.decomp <- GetGroupWeight(p, group.index)
    q.decomp <- GetGroupWeight(q, group.index)
    H.among <- RelativeEntropy(p.decomp$group.weight, q.decomp$group.weight)
    K <- length(group.index)
    H.within <- numeric(K)
    for (k in 1:K) {
      H.within[k] <- RelativeEntropy(p.decomp$cond.weight[[k]],
                                     q.decomp$cond.weight[[k]])
    }
    H.total <- RelativeEntropy(p, q)
    return(c(H.total, H.among, p.decomp$group.weight*H.within))
  }
}

FreeEnergy <- function(pi, R, group.index = NULL) {
  # Free energy
  #
  # Args:
  #  pi: portfolio weights
  #  R: simple returns of assets
  #
  # Returns:
  #  A non-negative number.
  
  pi <- as.numeric(pi)  
  R <- as.numeric(R)  
  if (is.null(group.index)) {
    r <- log(1 + R) 
    r.pi <- log(1 + sum(pi*R, na.rm = TRUE))  # portfolio log return
    r.WA <- sum(pi*r, na.rm = TRUE)  # weighted average log return
    return(r.pi - r.WA)
  } else {
    pi.decomp <- GetGroupWeight(pi, group.index)
    energy.total <- FreeEnergy(pi, R)
    K <- length(group.index)
    R.sub <- numeric(K)
    # Calculate returns of second level portfolios
    for (k in 1:K) {
      R.sub[k] <- sum((pi.decomp$cond.weight[[k]])*(R[group.index[[k]]]),
                      na.rm = TRUE)
    }
    energy.among <- FreeEnergy(pi.decomp$group.weight, R.sub)
    energy.within <- numeric(K)
    for (k in 1:K) {
      energy.within[k] <- FreeEnergy(pi.decomp$cond.weight[[k]], R[group.index[[k]]])
    }
    return(c(energy.total, energy.among, pi.decomp$group.weight*energy.within))
  }
}


EEControl <- function(pi.current, pi.next, nu.next, nu.implied = nu.next) {
  # Control term of the energy-entropy decomposition
  # See Pal and Wong (2013)
  #
  # Args:
  #   pi.current: current portfolio weights
  #   pi.next: portfolio weights for the next period
  #   nu.next: benchmark weights for the next period
  #   nu.implied: implied benchmark weights for the next period. The default is
  #               nu.next
  
  return(RelativeEntropy(pi.next, nu.next) -
         RelativeEntropy(pi.current, nu.implied))
}

EnergyEntropyDecomp <- function(market, weight, grouping = NULL, plot = TRUE) {
  # Energy-entropy decomposition
  # See Pal and Wong (2013)
  #
  # Args:
  #   market: a toymkt object
  #   weight: a zoo/matrix/dataframe representing the portfolio weights
  #   grouping: a matrix representing the grouping
  #   plot: TRUE or FALSE
  #
  # Returns:
  #   A list
  
  if (market$buy.and.hold == FALSE) {
    warning("The market is not buy-and-hold. Although the decomposition still holds, it is less informative.")
  }
  
  n.period <- dim(market$benchmark.weight)[1] - 1
  benchmark.weight <- as.matrix(data.frame(market$benchmark.weight))
  weight <- as.matrix(data.frame(weight))
  R <- as.matrix(data.frame(market$R))
  R[is.na(R)] <- 0
  
  
  if (is.null(grouping)) {
    # Case 1: single stage decomposition
    
    decomposition <- matrix(0, nrow = n.period, ncol = 5)

    for (i in 1:n.period) {
      pi.current <- weight[i, ]
      pi.next <- weight[i + 1, ]
      nu.current <- benchmark.weight[i, ]
      nu.next <- benchmark.weight[i + 1, ]
      R.current <- R[i, ]
      
      R.pi <- sum(pi.current*R.current, na.rm = TRUE)  # portfolio return
      R.nu <- sum(nu.current*R.current, na.rm = TRUE)  # benchmark return
      DlogV <- log((1 + R.pi)/(1 + R.nu))  # excess log return
      
      Denergy <- FreeEnergy(pi.current, R.current)  # free energy
      
      nu.t <- nu.current
      nu.t[is.na(nu.t)] <- 0  # precaution
      nu.implied <- nu.t*(1 + R.current)/sum(nu.t*(1 + R.current))  # implied weight

      # Control
      Dcontrol <- EEControl(pi.current, pi.next, nu.next, nu.implied)
      Ddrift <- Denergy + Dcontrol
      Drelative.entropy <- RelativeEntropy(pi.next, nu.next) -
                           RelativeEntropy(pi.current, nu.current)
      
      # Store results
      decomposition[i, ] <- c(DlogV, Denergy, -Drelative.entropy, Dcontrol, Ddrift)
    }
    
    colnames(decomposition) <- c("Excess log return",
                                 "Free energy",
                                 "Relative entropy",
                                 "Control",
                                 "Drift")
    time.index <- index(market$benchmark.weight)[1:n.period]
    decomposition <- zoo(decomposition, order.by = time.index)
    
    
    # if plot == TRUE, plot the decomposition
    if (plot == TRUE) {
      decomposition.cumsum <- rbind(0, apply(decomposition, 2, cumsum))
      decomposition.cumsum <- zoo(decomposition.cumsum,
                                  order.by = index(market$benchmark.weight))
      # extra margin
      extra <- 0.15*(max(decomposition.cumsum) - min(decomposition.cumsum))
      plot(decomposition.cumsum, plot.type = "single",
           xlab = "", ylab = "",
           ylim = c(min(decomposition.cumsum), max(decomposition.cumsum) + extra),
           main = "Energy-entropy decomposition",
           lwd = c(2, 2, 2, 2, 2),
           col = c("blue", "black", "orange", "green", "purple"))
      legend(x = "topleft",
             legend = c("log V", "Free energy",
                        "Relative entropy", "Control", "Drift"),
             lwd = c(2, 2, 2, 2, 2), cex = 0.7,
             col = c("blue", "black", "orange", "green", "purple"))
    }
    return(decomposition)
  } else {
    # Case 2: 2-stage decomposition
    
    # Error handling
    no.repeat <- as.numeric(names(table(grouping)))
    no.repeat <- sort(no.repeat)
    if (!all(no.repeat == 1:max(grouping))) {
      stop("The input must be a vector of positive integers between 1 and m, the number of groups. Each of 1, ..., m must appear at least once.")
    }
    n.assets <- length(grouping)
    if (n.assets != market$n) {
      stop("The length of grouping must be the same as the number of assets in market.")
    }    
    
    # Define group indices
    group.index <- list()
    for (k in 1:max(grouping)) {
      group.index[[k]] <- (1:n.assets)[grouping == k]
    }
    # number of groups
    n.group <- max(grouping)

    # Define storage variables
    # relative log return
    dlogV <- numeric(n.period)
    # portfolio group weights
    lambda <- matrix(0, ncol = n.group, nrow = n.period) 
    # benchmark group weights
    alpha <- matrix(0, ncol = n.group, nrow = n.period)
    # free energies
    free.energies <- matrix(0, ncol = 2 + n.group, nrow = n.period)
    # relative.entropies
    relative.entropies <- matrix(0, ncol = 2 + n.group, nrow = n.period)
    # controls
    controls <- matrix(0, ncol = 2 + n.group, nrow = n.period)
    
    # Colnames
    colnames(free.energies) <- c("FE.total", "FE.among",
                                 paste("FE.group", 1:n.group, sep = ""))
    colnames(relative.entropies) <- c("RE.total", "RE.among",
                                 paste("RE.group", 1:n.group, sep = ""))
    colnames(controls) <- c("C.total", "C.among",
                            paste("C.group", 1:n.group, sep = ""))
    

    # Energy-entropy decomposition
    for (i in 1:n.period) {
      pi.current <- weight[i, ]
      pi.next <- weight[i + 1, ]
      nu.current <- benchmark.weight[i, ]
      nu.next <- benchmark.weight[i + 1, ]      
      R.current <- R[i, ]
      
      # Excess log return
      R.pi <- sum(pi.current*R.current, na.rm = TRUE)  # portfolio return
      R.nu <- sum(nu.current*R.current, na.rm = TRUE)  # benchmark return
      dlogV[i] <- log((1 + R.pi)/(1 + R.nu))  # excess log return
      
      # Compute group weights for portfolio and benchmark
      pi.current.decomp <- GetGroupWeight(pi.current, group.index)
      pi.next.decomp <- GetGroupWeight(pi.next, group.index)
      nu.current.decomp <- GetGroupWeight(nu.current, group.index)
      nu.next.decomp <- GetGroupWeight(nu.next, group.index)      
      lambda[i, ] <- pi.current.decomp$group.weight
      alpha[i, ] <- nu.current.decomp$group.weight
      nu.implied <- nu.current*(1 + R.current)/sum(nu.current*(1 + R.current),
                                                   na.rm = TRUE)
      nu.implied.decomp <- GetGroupWeight(nu.implied, group.index)
    
      
      # Free energy chain rule
      free.energies[i, ] <- FreeEnergy(pi.current, R.current, group.index)

      # Relative entropy chain rule
      entropy.among <- RelativeEntropy(pi.current.decomp$group.weight,
                                       nu.current.decomp$group.weight) -
                       RelativeEntropy(pi.next.decomp$group.weight,
                                       nu.next.decomp$group.weight)
      entropy.within <- numeric(n.group)
      for (k in 1:n.group) {
        entropy.within[k] <- RelativeEntropy(pi.current.decomp$cond.weight[[k]],
                                             nu.current.decomp$cond.weight[[k]]) -
                             RelativeEntropy(pi.next.decomp$cond.weight[[k]],
                                             nu.next.decomp$cond.weight[[k]])
      }
      entropy.total <- entropy.among +
                       sum(pi.current.decomp$group.weight*entropy.within)
      relative.entropies[i, ] <- c(entropy.total, entropy.among,
                                   pi.current.decomp$group.weight*entropy.within)
      
      # Control chain rule
      control.among <- RelativeEntropy(pi.next.decomp$group.weight,
                                       nu.next.decomp$group.weight) -
                       RelativeEntropy(pi.current.decomp$group.weight,
                                       nu.implied.decomp$group.weight)
      control.within <- numeric(n.group)
      for (k in 1:n.group) {
        control.within[k] <- RelativeEntropy(pi.next.decomp$cond.weight[[k]],
                                             nu.next.decomp$cond.weight[[k]]) -
                             RelativeEntropy(pi.current.decomp$cond.weight[[k]],
                                             nu.implied.decomp$cond.weight[[k]])
      }
      control.total <- control.among +
                       sum(pi.current.decomp$group.weight*control.within)  
      controls[i, ] <- c(control.total, control.among,
                    pi.current.decomp$group.weight*control.within)
    }
    
    time.index <- index(market$R)
    dlogV <- zoo(dlogV, order.by = time.index)
    free.energies <- zoo(free.energies, order.by = time.index)
    relative.entropies <- zoo(relative.entropies, order.by = time.index)
    controls <- zoo(controls, order.by = time.index)
    
    decomposition <- list(dlogV = dlogV,
                          free.energies = free.energies,
                          relative.entropies = relative.entropies,
                          controls = controls)
    
    if (plot == TRUE) {
      energy.among <- free.energies[, 2]
      energy.within <- free.energies[, 1] - free.energies[, 2]
      entropy.among <- relative.entropies[, 2]
      entropy.within <- relative.entropies[, 1] - relative.entropies[, 2]
      control.among <- controls[, 2]
      control.within <- controls[, 1] - controls[, 2]
      combined <- cbind(dlogV,
                        energy.among, energy.within,
                        entropy.among, entropy.within,
                        control.among, control.within)
      combined <- rbind(0, apply(combined, 2, cumsum))
      combined <- zoo(combined, order.by = index(market$benchmark.weight))
      extra <- 0.2*(max(combined) - min(combined))
      plot(combined, plot.type = "single",
           xlab = "", ylab = "",
           ylim = c(min(combined), max(combined) + extra),
           main = "Energy-entropy decomposition",
           lwd = rep(2, 7),
           col = c("blue",
                   "black", "black",
                   "orange", "orange",
                   "green", "green"),
           lty = c(1, 1, 2, 1, 2, 1, 2))
      legend(x = "topleft",
             legend = c("Relative log return",
                        "Free energy: top level",
                        "Free energy: second level",
                        "Relative entropy: top level",
                        "Relative entropy: second level",
                        "Control: top level",
                        "Control: second level"),
             lwd = rep(2, 7), cex = 0.6,
             lty = c(1, 1, 2, 1, 2, 1, 2),
             col = c("blue",
                     "black", "black",
                     "orange", "orange",
                     "green", "green"))
    }
    
    return(decomposition)
  }
}


GetNewLambdaWeight <- function(pi.current, mu.next,
                                    energy, lambda = 0.5) {
  # New weights for the lambda-strategy
  # Section 6 of Pal and Wong (2013)
  # Use SDE approximation, good only for lambda small
  #
  # Args:
  #   pi.current: current portfolio weights
  #   mu.next: market weights for the next period (with full support)
  #   energy: free energy of the previous period
  #   lambda: parameter
  #
  # Returns:
  #   A portfolio weight vector
  
  pi.current <- as.numeric(pi.current)
  mu.next <- as.numeric(mu.next)
  
  # J-divergence
  J <- 0.5*(RelativeEntropy(pi.current, mu.next) +
            RelativeEntropy(mu.next, pi.current))
  
  # multiplier 
  multiplier <- min(lambda*energy/(2*J), 1)
  return(pi.current + multiplier*(mu.next - pi.current))
}
  

GetLambdaWeight <- function(market,
                            initial.weight = market$benchmark.weight[1, ],
                            lambda) {
  # Compute the weights of the lambda-strategy
  weight <- matrix(0, nrow = dim(market$benchmark.weight)[1],
                   ncol = dim(market$benchmark.weight)[2])
  
  weight[1, ] <- as.numeric(initial.weight)  
  for (i in 1:(dim(weight)[1] - 1)) {
    energy <- FreeEnergy(weight[i, ], market$R[i, ])
    weight[i + 1, ] <- GetNewLambdaWeight(weight[i, ],
                                          market$benchmark.weight[i + 1, ],
                                          energy, lambda)
  }
  
  return(weight)
}