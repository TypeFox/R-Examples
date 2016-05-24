# Class "capdist" - capital distribution

capdist <- function(x) {
  # if sum of weights is not close enough to 1, output warning message
  if (!all(x >= 0)) {
    stop("x must be a non-negative vector.")
  } else if (abs(sum(x) - 1) > 0.001) {
    warning("The sum of the weights is not close to 1.")
  }
  class(x) <- "capdist"
  return(x)
}

plot.capdist <- function(x, draw.line = TRUE, cut.end = 0.1, ...) {
  # plot cap weights against rank in log-log scale 
  x <- sort(x, decreasing = TRUE)
  log.x <- log(as.numeric(x))
  log.rank <- log(1:length(log.x))
  plot(log.rank, log.x,
       xlab = "log rank", ylab = "log cap weight",
       main = "Capital distribution curve",
       col = "blue", type = "l", lwd = 2)
  
  # if line == TRUE, add a regression line
  if (draw.line == TRUE) {
    n <- round((1 - cut.end) * length(log.x))
    abline(lm(log.x[1:n] ~ log.rank[1:n]), lwd = 2, lty = 2)
  }
}

print.capdist <- function(x, m = 5L, cut.end = 0.1, ...) {
  # summary for capdist objects
  n <- length(x)
  cat("Number of assets: ", n, "\n\n")
  if (n <= 2 * m) {
    print(x)
  } else {
    cat("Largest ", m, ":\n", sep = "")
    print(x[1:m])
    cat("\nSmallest ", m, ":\n", sep = "")
    print(x[(n - (m - 1)):n])
    cat("\nSlope of capital distribution curve:\n")
    cat(CapDistSlope(x, cut.end = cut.end))
  }
}

CapDistSlope <- function(x, cut.end = 0.1) {
  # slope of the capital distribution curve
  x <- sort(x, decreasing = TRUE)
  log.x <- log(as.numeric(x))
  log.rank <- log(1:length(log.x))
  n <- round((1 - cut.end) * length(log.x))
  return(as.numeric(coef(lm(log.x[1:n] ~ log.rank[1:n]))[2]))
}

ParetoCapDist <- function(n, index = 1) {
  # generate a Pareto capital distribution with a
  # user-defined slope for the log-log curve
  x <- (1:n)^(-index)
  x <- x/sum(x)
  return(capdist(x))
}

