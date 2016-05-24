# Measures of diversity

Diversity <- function(x, p = 0.5) {
  # Diversity as defined by Fernholz
  # Formally the same as the Lp norm.
  # 
  # Args:
  #   x: a numeric vector of non-negative numbers
  #   p: a non-negative parameter
  
  x <- as.numeric(x)
  
  if (!all(x >= 0)) {
    stop("The entries must be non-negative.")
  }
  
  return((sum(x^p))^(1/p))
}


ShannonEntropy <- function(x) {
  # Shannon entropy
  #
  # Args:
  #   x: a numeric vector of non-negative numbers summing to 1
  
  x <- as.numeric(x)
  
  if (!all(x >= 0) | (abs(sum(x) - 1) > 0.01)) {
    stop("The input is not a probability vector.")
  }
  
  return(sum(-x * log(x), na.rm = TRUE))
}

GeometricMean <- function(x, weight = rep(1/length(x), length(x))) {
  # Geometric mean
  #
  # Args:
  #   x: a numeric vector of non-negative numbers
  #   w: a vector of probability weights 
  
  x <- as.numeric(x)
  
  if (!all(x >= 0)) {
    stop("The entries must be non-negative.")
  }
  
  return(prod(x^weight))
}


RenyiEntropy <- function(x, p = 0.5) {
  # Renyi Entropy
  # Example 3.4.8 of Fernholz (2002)
  #
  # Args:
  #   x: a probability vector
  #   p: any number not equal to 1. The default is 0.5.
  
  if (p == 1) stop("p cannot be 1")
  if (p > 1) warning("for this value of p it is not a measure of diversity")
  
  return(log(sum(x^p) / (1 - p)))
}
