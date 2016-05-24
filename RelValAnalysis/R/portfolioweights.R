Invest <- function(market, weight, plot = TRUE) {
  # Calculate portfolio value and returns together with those for the benchmark
  #
  # Args:
  #   market: a toymkt object
  #   weight: a vector/matrix/dataframe/zoo object containing the portfolio
  #           weights. If weight is a vector, it is assumed that the portfolio
  #           is constant-weighted throughout.
  #
  # Returns:
  #   A list containing the returns of the portfolio and the benchmark.
  
  # number of periods
  n.period <- dim(market$R)[1]  
  
  
  R <- as.matrix(data.frame(market$R[1:n.period, ]))
  benchmark.weight <- as.matrix(data.frame(market$benchmark.weight[1:n.period, ]))
  R.nu <- rowSums(R * benchmark.weight, na.rm = TRUE)
  r.nu <- log(1 + R.nu)
  Z.nu <- exp(c(0, cumsum(r.nu)))
  
  # Simulate portfolios
  # If weight is a numeric vector, is.vector(A)
  if (is.vector(weight)) {
    warning("Since only one weight vector is supplied, the portfolio is assumed to be constant-weighted.")
    weight.matrix <- matrix(0, nrow = n.period, ncol = market$n)
    for (i in 1:n.period) {
      weight.matrix[i, ] <- weight
    }
  } else {
    weight.matrix <- as.matrix(data.frame(weight[1:n.period, ]))
  }
  R.pi <- rowSums(R * weight.matrix, na.rm = TRUE)
  r.pi <- log(1 + R.pi)  
  Z.pi <- exp(c(0, cumsum(r.pi)))
  
  
  # Gather results
  growth <- cbind(Z.pi, Z.nu)
  colnames(growth) <- c("portfolio", "Benchmark")
  growth <- zoo(growth, order.by = index(market$growth))
  
  R.combined <- cbind(R.pi, R.nu)
  colnames(R.combined) <- c("portfolio", "Benchmark")
  R.combined <- zoo(R.combined, order.by = index(market$growth))
  
  r.combined <- cbind(r.pi, r.nu)
  colnames(r.combined) <- c("portfolio", "Benchmark")
  r.combined <- zoo(r.combined, order.by = index(market$growth))
  
  output <- list(growth = growth, R = R.combined, r = r.combined)
  
  # plot
  if (plot == TRUE) {
    par(mfrow = c(1, 2))
    extra <- 0.5*(max(growth) - min(growth))
    plot(growth, plot.type = "single",
         xlab = "", ylab = "",
         ylim = c(min(growth), max(growth) + extra),
         lwd = 2, col = c("blue", "darkgrey"),
         main = "Growth of $1")
    legend(x = "topleft", legend = c("portfolio", "Benchmark"),
           lwd = c(2, 2), col = c("blue", "darkgrey"), cex = 0.7)
    V <- log(growth[, 1] / growth[, 2])
    plot(V, xlab = "", ylab = "",
         lwd = 2, col = "blue",
         main = "Log relative value")
    par(mfrow = c(1, 1))
  }
  
  return(output)
}


GetWeight <- function(market, weight.function, ...) {
  # Compute portfolio weights from a weight function
  #
  # Args: 
  #   market: a toymkt object
  #   weight.function: a weight function
  #
  # Returns:
  #   A matrix of weights.
  
  if ("toymkt" %in% class(market)) {
    benchmark <- market$benchmark.weight
  } else {
    benchmark <- market
  }
  
  weight.matrix <- matrix(0, nrow = dim(benchmark)[1], ncol = dim(benchmark)[2])
  
  for (i in 1:dim(benchmark)[1]) {
    weight.matrix[i, ] <- weight.function(benchmark[i, ], ...)
  }
  
  colnames(weight.matrix) <- colnames(market$benchmark.weight)
  
  return(weight.matrix)
}
