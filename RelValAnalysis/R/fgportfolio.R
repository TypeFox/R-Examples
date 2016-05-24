# The class "fgp" and associated methods.

fgp <- function(name = NULL, gen.function, weight.function) {
  # Create a fgp object.
  #
  # Args:
  #   name: Name of the portfolio (a string).
  #   Phi: The generating function.
  #   weight.function: 
  #
  # Returns:
  #   A fgp object.
  
  portfolio <- list()
  portfolio$name <- name
  portfolio$gen.function <- gen.function
  portfolio$weight.function <- weight.function
  class(portfolio) <- "fgp"
  
  return(portfolio)
}

print.fgp <- function(x, ...) {
  # just print the name for now
  cat("A functionally generated portfolio: ", x$name, sep = "")
  cat("\n\nComponents of the fgp objects (use $):\n")
  cat(names(x))
}


DiversityPortfolio <- function(p = 0.5) {
  # Create an fgp object represented the
  # diversity-weighted portfolio.
  #
  # Args:
  #   p: interpolation paramter
  #
  # Returns:
  #   A fgp object
  name <- paste("Diversity-weighted portfolio with p = ", p, sep = "")
  
  Phi <- function(x) {
    x <- as.numeric(x)
    return((sum(x^p))^(1/p))
  }
  
  weight.function <- function(x) {
    x <- as.numeric(x)
    return(x^p / sum(x^p))
  }
  
  portfolio <- fgp(name = name, gen.function = Phi, weight.function = weight.function)
  
  return(portfolio)
}

ConstantPortfolio <- function(weight) {
  # Create an fgp object representing the constant-weighted portfolio
  # Note: The number of stocks is equal to the length of weight
  #
  # Args:
  #   weight: weights of the constant-weighted portfolio
  #
  # Returns:
  #   A fgp object
  n <- length(weight)
  
  name <- paste("Constant-weighted portfolio (", n, " stocks)", sep = "")
  
  Phi <- function(x) {
           x <- as.numeric(x)
           return(prod(x^weight))
         }
  
  weight.function <- function(x) {
                       return(weight)
                     }

  portfolio <- fgp(name = name, gen.function = Phi,
                   weight.function = weight.function)

  return(portfolio)
}


# The entropy-weighted portfolio
EntropyPortfolio <- fgp(name = "Entropy-weighted portfolio",
                        gen.function = function(x) {
                                         x <- as.numeric(x)
                                         return(return(sum(-x * log(x), na.rm = TRUE)))
                                       },
                        weight.function = function(x) {
                                            x <- as.numeric(x)
                                            y <- log(x)
                                            y[y == -Inf] <- 0
                                            return((x*y)/sum(x*y))
                                          })



FernholzDecomp <- function(market, portfolio, plot = TRUE) {
  # Fernholz decomposition
  #
  # Args:
  #   market: a toymkt object.
  #   portfolio: a fgp object.
  #   plot: if TRUE, the result will be plotted.
  #
  # Returns:
  #   A list containing 5 components.

  
  # The market must be buy-and-hold
  if (market$buy.and.hold == FALSE) {
    stop("The decomposition does not make sense if the market is not buy-and-hold.")
  }
  
  # Simulate portfolios
  weight.matrix <- GetWeight(market, portfolio$weight.function)
  growth <- Invest(market, weight.matrix, plot = FALSE)$growth
  
  # Comput terms in the decomposition
  Z.pi <- growth[, 1]  # portfolio
  Z.mu <- growth[, 2]  # market portfolio

  relative.return <- log(Z.pi / Z.mu)  # relative log return
  
  log.Phi <- log(apply(market$benchmark.weight, 1, portfolio$gen.function))
  log.Phi <- log.Phi - log.Phi[1]  # fluctuation of generating function
  
  output <- list()
  output$portfolio <- Z.pi
  output$benchmark <- Z.mu
  output$relative.return <- relative.return
  output$log.gen.funtion <- log.Phi
  output$drift <- relative.return - log.Phi
  

  # If plot == TRUE, plot the decomposition.
  combined <- cbind(output$relative.return, output$log.gen.funtion, output$drift)
  if (plot == TRUE) {
    extra <- 0.15*(max(combined) - min(combined))
    plot(combined, plot.type = "single", ylab = "", xlab = "",
         ylim = c(min(combined), max(combined) + extra),
         col = c("blue", "orange", "purple"), lwd = c(2, 2, 2),
         main = "Fernholz decomposition")
    legend(x = "topleft",
           legend = c("relative log return",
                      "generating function",
                      "drift process"),
           col = c("blue", "orange", "purple"),
           lwd = c(2, 2, 2), cex = 0.7)
  }
  
  return(output)
}

