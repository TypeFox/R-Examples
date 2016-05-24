# The class toymkt and associated methods.

toymkt <- function(price = NULL, R = NULL,
                   benchmark.weight = NULL,
                   initial.weight = NULL,
                   buy.and.hold = TRUE) {
  # Create a toymkt object from price or return data.
  #
  # Args:
  #   price: A zoo object containing the price or market captitalizations of
  #          the assets. Matrix and dataframes are also accepted but the
  #          dates will only be shown properly if a zoo object is supplied.
  #   R: Simple returns of the assets (same format as price).
  #   benchmark.weight: Portfolio weights of the benchmark (same format as
  #                     price).
  #   initial.weight: Initial market weights (normalization).
  #   buy.and.hold: If TRUE, the benchmark is a buy-and-hold portfolio. If not,
  #                 not. Default is TRUE.
  #
  # Returns:
  #   A toymkt object containing the data and additional information.  
  
  # Error handling
  # (At least one of price and R must be supplied.)
  if (is.null(price) & is.null(R)) {
    stop("At least one of price and R must be given.")
  }
  if (is.null(initial.weight) == FALSE) {
    if (all(initial.weight >= 0) == FALSE) {
      stop("The entries of initial.weight must be non-negative.")
    }    
  }

  # price, R, and benchmark.weight are not zoo objects, they will be converted
  # to zoo objects (but dates may be lost).
  if ((is.null(price) == FALSE) & (("zoo" %in% class(price)) == FALSE)) {
    price <- zoo(price)
  }
  if ((is.null(R) == FALSE) & (("zoo" %in% class(R)) == FALSE)) {
    R <- zoo(R)
  } 
  if ((is.null(benchmark.weight) == FALSE) & 
      (class(benchmark.weight) != "zoo")) {
    benchmark.weight <- zoo(benchmark.weight)
  } 
  
  market <- list()
  # Components:
  #   growth: growth of $1 for each asset
  #   R: simple return
  #   r: log return
  #   benchmark.weight: benchmark weights
  #   n: number of assets
  #   buy.and.hold: TRUE if market is buy-and-hold
  
 
  # Case 1: price is given.
  if (!is.null(price)) {
    n <- dim(price)[2]  # number of stocks
    
    # If price is given, market is assumed to be buy-and-hold even if
    # buy.and.hold = FALSE.
    if (buy.and.hold == FALSE) {
      warning("Since price is given, the market is assumed to be buy-and-hold.")
    }
    
    # All entries of price must be strictly positive
    if (all(price > 0) == FALSE) {
      stop("All entries of price must be strictly positive.")
    }
    
    # If names of assets are not given, they are set to Asset1, Asset2, ...,
    # Assetn.
    if (is.null(colnames(price))) {
      colnames(price) <- paste("Asset", 1:dim(price)[2], sep = "")
    }
    
    #  Compute growth of $1
    growth <- price
    for (j in 1:dim(growth)[2]) {  # normalize at $1 at beginning
      growth[, j] <- growth[, j] / as.numeric(growth[1, j])
    }
    market$growth <- growth
    
    # Normalize by initial.weight if it is given. Otherwise, initial.weight is
    # assumed to be equal-weighted.
    if (is.null(initial.weight) == FALSE) {  # normalize price
      initial.weight <- as.numeric(initial.weight)
    } else {
      warning("Since initial.weight is not given, the benchmark is assumed to be equal-weighted initially.")
      initial.weight <- rep(1/n, n)
    }
    for (j in 1:dim(price)[2]) {
      price[, j] <- price[, j] * (initial.weight[j] / as.numeric(price[1, j]))
    }

    # Compute log and simple returns
    # The convention is that R[t] is the return between date t and date t + 1.
    log.price <- log(price)
    dates <- index(price)[1:(dim(price)[1] - 1)]
    r <- diff(log.price)
    r <- zoo(r, order.by = dates)
    R <- exp(r) - 1
    market$R <- R
    market$r <- r
    
    # Compute benchmark.weight
    market$benchmark.weight <- price / rowSums(price, na.rm = TRUE)
    
    # Other attributes
    market$n <- n
    market$buy.and.hold <- TRUE
  } else {
    # Case 2: price is not given.
    # Then R must be given.
    
    # The number of periods must be at least 2.
    n.periods <- dim(R)[1]
    if (n.periods <= 1) {
      stop("The number of rows of R must be at least 2.")
    }
    
    n <- dim(R)[2]  # number of assets
    
    # If names of assets are not given, they are set to Asset1, Asset2, ...,
    # Assetn.
    if (is.null(colnames(R))) {
      colnames(R) <- paste("Asset", 1:n, sep = "")
    }    
    
    # Compute growth of $1
    R.matrix <- (as.matrix(data.frame(R)))
    growth <- apply(rbind(rep(1, n), as.matrix(1 + R.matrix)), 2, cumprod)
    # Define the final date
    new.date <- index(R)[n.periods] +
      as.numeric(index(R)[n.periods] - index(R)[n.periods - 1])
    growth <- zoo(growth, order.by = c(index(R), new.date))
    
    if (buy.and.hold == TRUE) {
      # If buy.and.hold = TRUE, then things are computed from initial.weight
      # and R.

      # If initial.weight and benchmark.weight are not specified, the initial
      # weights will be assumed to be 1/n.
      if (is.null(initial.weight)) {
        if (is.null(benchmark.weight)) {
          warning("Since initial.weight is not given, the benchmark is assumed to be equal-weighted initially.")
          initial.weight <- rep(1/n, n)
        } else {
          initial.weight <- as.numeric(benchmark.weight[1, ])
        }
      }
      
      price <- growth
      for (j in 1:dim(price)[2]) {
        price[, j] <- price[, j] * (initial.weight[j] / as.numeric(price[1, j]))
      }
      
      market$growth <- growth
      market$R <- R
      market$r <- log(1 + R)
      market$benchmark.weight <- price / rowSums(price, na.rm = TRUE)
      market$n <- n
      market$buy.and.hold <- TRUE
    } else {
      # Buy-and-hold = FALSE
      if (is.null(benchmark.weight)) {
        stop("benchmark.weight must be supplied when buy.and.hold = FAlSE.")
      }
      n <- dim(R)[2]  # number of assets
      market$growth <- growth
      market$R <- R
      market$r <- log(1 + R)
      market$benchmark.weight <- benchmark.weight
      market$n <- n
      market$buy.and.hold <- FALSE
    }
  }
  
  class(market) <- "toymkt"
  return(market)
}

print.toymkt <- function(x, ...) {
  # Print a toymkt object. (Names of the assets and pointers to stored
  # information.)
  #
  # Args:
  #   x: a toymkt object.  

  cat("\nThere are ", x$n, " assets:\n", sep = "")
  cat(colnames(x$growth), "\n\n")
  
  not <- ""
  if (x$buy.and.hold == FALSE) {
    not <- "not "
  }
  cat("The benchmark is ", not, "buy-and-hold.\n\n", sep = "")
  cat("Other information (use $):\n")
  cat(names(x), "\n\n")
}

plot.toymkt <- function(x, ...) {
  # Plot a toymkt object. 
  #
  # Args:
  #   x: a toymkt object.
  #
  # Returns:
  #   NULL (a few plots are shown)
  
  n <- x$n
  asset.names <- colnames(x$growth)
  op <- par(ask = TRUE)
  par(mfrow = c(1, 1))
  # Plot time series of growth of $1
  extra <- 0.15*(max(x$growth) - min(x$growth))
  plot(x$growth, plot.type = "single",
       xlab = "", ylab = "", col = 1:n, lwd = rep(2, n),
       ylim = c(min(x$growth), max(x$growth) + extra),
       main = "Growth of $1 of the assets")
  legend(x = "topleft", legend = asset.names,
           col = 1:n, lwd = rep(2, n))

  # plot time series of log returns
  extra <- 0.15*(max(x$r) - min(x$r))
  plot(x$r, plot.type = "single",
       ylim = c(min(x$r), max(x$R) + extra),
       xlab = "", ylab = "", col = 1:n, lwd = rep(2, n),
       main = "Time series of log returns")
  legend(x = "topleft", legend = asset.names,
         col = 1:n, lwd = rep(2, n))
  
  # plot time series of benchmark weights
  extra <- 0.15*(max(x$benchmark.weight) - min(x$benchmark.weight))
  plot(x$benchmark.weight, plot.type = "single",
       xlab = "", ylab = "", col = 1:n, lwd = rep(2, n),
       ylim = c(min(x$benchmark.weight),
                max(x$benchmark.weight) + extra),
       main = "Time series of benchmark weights")
  legend(x = "topleft", legend = asset.names,
         col = 1:n, lwd = rep(2, n))
  
  # plot time series of market entropy
  entropy.ts <- apply(x$benchmark.weight, 1, ShannonEntropy)
  entropy.ts <- zoo(entropy.ts, order.by = index(x$benchmark.weight))
  plot(entropy.ts, xlab = "", ylab = "",
       main = "Time series of Shannon entropy",
       lwd = 2, col = "blue")
  par(op)
}