# File mipfp/R/binary.R
# by Johan Barthelemy and Thomas Suesse
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ------------------------------------------------------------------------------
# This files provides the following functions focusing on multivariate Bernoulli
# (binary) variables:
#   * Odds2PairProbs: Converts odds ratio between K binary variables to pairwise
#                     probability.
#   * Corr2PairProbs: Converts correlation between K binary variables to 
#                     pairwise probability.
#   * Odds2Corr: Converts odds ratios between K binary variables to correlation.
#   * Corr2Odds: Converts correlation between K binary variables to odds ratios.
#   * ObtainMultBinaryDist: Determines a multivariate binary joint-distribution.
#   * RMultBinary: Simulate multivariate binary data using Ipfp.
# ------------------------------------------------------------------------------

Odds2PairProbs <- function(odds, marg.probs) {
  # Converts odds ratio between K Bernoulli variables to pairwise probability.
  # 
  # This function transforms the odds ratio measure of association O_ij between
  # K binary (Bernoulli) random variables X_i and X_j to the pairwise 
  # probability P(X_i = 1, X_j = 1) where O_ij is defined as
  # O_ij = P(X_i = 1, X_j = 1) * P(X_i = 0, X_j = 0) 
  #        / P(X_i = 1, X_j = 0) * P(X_i = 0, X_j = 1).
  #  
  # Author: T. Suesse
  #  
  # Args:
  #   odds: A K x K matrix where the i-th row and the j-th colum represents the
  #         odds ratio O_ij.
  #   marg.probs: A vector with K elements of marginal probabilities with the
  #               i-th entry refering to P(X_i = 1).
  #
  # Returns: A matrix of the same dimension as odds containing the pairwise
  #          probabilities.
  
  # checking if input argument odds is valid
  if (is.null(odds) == TRUE | min(odds) < 0.0) {
    stop('Error: odds ratio odds misspecified!')
  }
  
  # checking if input argument marg.probs is specified
  if (is.null(marg.probs) == TRUE | min(marg.probs < 0.0)) {
    stop('Error: marg.probs missspecified!')
  }
  
  # determining some initial variables  
  K <- length(marg.probs)
  pair.proba <- matrix(0.0, K, K)
  
  # construction the upper triangular part of the final matrix
  # ... loop over the rows
  for (i in 1:(K - 1)) {
    
    # ... loop over the columns
    for (j in (i + 1):K) {
      
      # find pairwise probabilities by solving a quadratic polynomial      
      coeff <- c(odds[i, j] * marg.probs[i] * marg.probs[j], 
                 - (1 + (odds[i, j] - 1) * (marg.probs[i] + marg.probs[j])),                 
                 odds[i,j] - 1)
      roots <- Re(polyroot(coeff))
      
      # determining the solution wrt to the constraints sets by marg.probs
      s1 <- min(marg.probs[i], marg.probs[j])
      s2 <- max(0.0, marg.probs[i] + marg.probs[j] - 1.0)      
      cond <- (roots <= s1 & roots >= s2)      
      
      # ... no suitable solution
      if (max(cond) == FALSE) {
        stop("ERROR for variable ", i, " and variable ", j,
             ": Pairwise probability exceed limits given by marg.probs!\n")
      }
      
      # ... first suitable solution
      if (cond[1] == TRUE) {
        pair.proba[i, j] <- roots[1]
      } else {
        pair.proba[i, j] <- roots[2]
      }
      
    }
    
  }
  
  # adding the lower triangular part of the final matrix
  pair.proba <- pair.proba + t(pair.proba)
  
  # adding the diagonal of the final matrix
  diag(pair.proba) <- marg.probs
  
  # adding the dimension names
  dimnames(pair.proba) <- dimnames(odds)
  
  # returning the result
  return(pair.proba)
  
}

Corr2PairProbs <- function(corr, marg.probs) {
  # Converts correlation between K Bernoulli variables to pairwise probability.
  # 
  # This function transforms the correlation measure of association C_ij between
  # K binary (Bernoulli) random variables X_i and X_j to the pairwise
  # probability P(X_i = 1, X_j = 1) where C_ij is defined as
  # C_ij = cov(X_i, X_j) / sqrt(var(X_i) * var(X_j)).
  #  
  # Author: T. Suesse
  #  
  # Args:
  #   corr: A K x K matrix where the i-th row and the j-th colum represents the
  #         correlation C_ij.
  #   marg.probs: A vector with K elements of marginal probabilities with the
  #               i-th entry refering to P(X_i = 1).
  #
  # Returns: A matrix of the same dimension as corr containing the pairwise
  #          probabilities.
  
  # checking if input argument corr is specified
  if (is.null(corr) == TRUE ) {
    stop('Error: correlations corr not specified!')
  }
  
  # checking if input argument marg.probs is specified
  if (is.null(marg.probs) == TRUE | min(marg.probs < 0.0)) {
    stop('Error: marg.probs missspecified!')
  }
  
  K <- length(marg.probs)
  
  # determing the pairwise probabilities
  pair.proba <- corr * sqrt((marg.probs * (1 - marg.probs)) %*% 
                              t((marg.probs * (1 - marg.probs)))) + 
    marg.probs %*% t(marg.probs)
  
  # determining the constrains
  p.col <- matrix(marg.probs, K, K, byrow = TRUE)
  p.row <- t(p.col)
  s1 <- pmin(p.col, p.row)
  s2 <- pmax(matrix(0, K, K), p.col + p.row - 1)
  
  # assessing constraints
  cond <- matrix(as.logical((pair.proba <= s1) & (pair.proba >= s2)), K, K)
  diag(cond) <- TRUE
  
  if (min(cond) == FALSE) {
    warning("Correlation exceeds constrains set by marg.probs, i.e. ",
            "pair.proba[i, j] <= marg.probs[i]\n")
    cat("Problematic pairs:\n")
    print(which(cond == min(cond), arr.ind = TRUE))
  }
  
  # adding the dimension names
  dimnames(pair.proba) <- dimnames(corr)
  
  # returning the pairwise probabilities matrix
  return(pair.proba)
  
}

Odds2Corr <- function(odds, marg.probs) {
  # Converts odds ratios between K Bernoulli variables to correlation.
  # 
  # Assuming that there is K binary (Bernoulli) random variables X_1, ..., X_K, 
  # this function transforms the odds ratios measure of association O_ij of  
  # every pair (X_i, X_j) to the correlation C_ij where C_ij is defined as 
  # C_ij = cov(X_i, X_j) / sqrt(var(X_i) * var(X_j)) and
  # O_ij = P(X_i = 1, X_j = 1) * P(X_i = 0, X_j = 0) /
  #            P(X_i = 1, X_j = 0) * P(X_i = 0, X_j = 1).
  #  
  # Author: T. Suesse
  #  
  # Args:
  #   odds: A K x K matrix where the i-th row and the j-th colum represents the
  #         odds ratio O_ij.
  #   marg.probs: A vector with K elements of marginal probabilities with the
  #               i-th entry refering to P(X_i = 1).
  #
  # Returns: A list whose elements are
  #   corr: A matrix of the same dimension as odds containing the correlation
  #   pair.proba: A matrix of the same dimension as odds containing the pairwise
  #               probabilities.
  
  # checking if input argument odds is valid
  if (is.null(odds) == TRUE | min(odds) < 0.0) {
    stop('Error: odds ratio odds misspecified!')
  }
  
  # checking if input argument marg.probs is specified
  if (is.null(marg.probs) == TRUE | min(marg.probs < 0.0)) {
    stop('Error: marg.probs missspecified!')
  }
  
  # computing the pairwise probabilities and correlation matrix
  pair.proba <- Odds2PairProbs(odds, marg.probs)
  corr <- (pair.proba - marg.probs %*% t(marg.probs)) / 
    sqrt((marg.probs * (1 - marg.probs)) %*% 
           t((marg.probs * (1 - marg.probs))))  
  
  # adding the dimension names
  dimnames(corr) <- dimnames(odds)
  
  # returning the correlation and pairwise probabilities matrices
  return(list("corr" = corr, "pair.proba" = pair.proba)) 
  
}

Corr2Odds <- function(corr, marg.probs) {
  # Converts correlation between K Bernoulli variables to odds ratio.
  # 
  # This function transforms the correlation measure of association C_ij between
  # K binary (Bernoulli) random variables X_i and X_j to the odds ratios O_ij
  # where C_ij is defined as C_ij = cov(X_i, X_j) / sqrt(var(X_i) * var(X_j))
  # and O_ij = P(X_i = 1, Y_j = 1) * P(X_i = 0, Y_j = 0) /
  #            P(X_i = 1, Y_j = 0) * P(X_i = 0, Y_j = 1).
  #  
  # Author: T. Suesse
  #  
  # Args:
  #   corr: a K x K matrix where the i-th row and the j-th colum represents the
  #         correlation C_ij.
  #   marg.probs: a vector with K elements of marginal probabilities with the
  #               i-th entry refering to P(X_i = 1).
  #
  # Returns: A list whose elements are
  #   corr: A matrix of the same dimension as odds containing the correlation.
  #   pair.proba: A matrix of the same dimension as odds containing the pairwise
  #               probabilities.
  
  # checking if input argument corr is specified
  if (is.null(corr) == TRUE ) {
    stop('Error: correlations corr not specified!')
  }
  
  # checking if input argument marg.probs is specified
  if (is.null(marg.probs) == TRUE | min(marg.probs < 0.0)) {
    stop('Error: marg.probs missspecified!')
  }
  
  K <- length(marg.probs)
  
  # computing the pairwise probabilitie
  pair.proba <- Corr2PairProbs(corr, marg.probs)
  
  # computing the odds ratios
  mu.x <- matrix(marg.probs, K, K, byrow = TRUE)
  mu.y <- matrix(marg.probs, K, K, byrow = FALSE)
  odds <- pair.proba * (1 - mu.x - mu.y + pair.proba) / 
    ((mu.x - pair.proba) * (mu.y - pair.proba))
  diag(odds) <- Inf
  
  # adding the dimension names
  dimnames(odds) <- dimnames(corr)
  
  # returning the matrices of correlation and pairwise probabilities
  return(list("odds" = odds, "pair.proba" = pair.proba))
  
}

ObtainMultBinaryDist <- function(odds = NULL, corr = NULL, marg.probs, ...) {
  # Generates a multivariate binary joint-distribution.
  # 
  # Applies the IPFP procedure to obtain a joint distribution of K multivariate 
  # binary variables. It requires as input the odds ratio or alternatively the 
  # correlation as a measure of association between all the binary variables and
  # a vector of marginal probabilities. This function is useful when one wants 
  # to simulate from a multivariate binary distribution when only first order
  # (marginal probabilities) and second order moments (correlation OR 
  # odds ratio) are available.
  #
  # Author: T. Suesse
  #  
  # Args:
  #   odds: A K x K matrix where the i-th row and the j-th colum represents the
  #         odds ratio O_ij.
  #   corr: A K x K matrix where the i-th row and the j-th colum represents the
  #         correlation C_ij.
  #   marg.probs: A vector with K elements of marginal probabilities with the
  #               the i-th entry refering to P(X_i = 1).
  #   ...: arguments that can be passed to the Ipfp function such as tol, iter,
  #        print, compute.cov.
  #     
  # Returns: A list generated by the Ipfp function whose elements are
  #   joint.proba: Resulting multivariate joint-probabilities.
  #   stp.crit: Final value of the stopping criterion.
  #   evol.stp.crit: Evolution of the stopping criterion over the iterations.
  #   conv: boolean Indicating whether the algorithm converged to a solution.
  #   check.margins: A list returning, for each margin, the absolute maximum 
  #                  deviation between the desired and generated margin.
  #   label: The names of the variables.
  
  # checking if marginal probabilities are provided
  if (is.null(marg.probs) == TRUE) {
    stop("Marginal probabilities need to be provided!\n")
  }
  
  # checking if odds ratio or correlations are provided
  if (is.null(odds) == TRUE & is.null(corr) == TRUE) {
    stop("Correlation or odds ratios must be provided in a K by K matrix!")
  }
  
  # generation of the pairwise probabilities
  if (is.null(odds) == TRUE) {
    pair.proba <- Corr2PairProbs(corr, marg.probs)    
  } else {  
    pair.proba <- Odds2PairProbs(odds, marg.probs)    
  }
  
  # retrieving and storing the label of the variables
  lab = rownames(pair.proba)
  
  # generating the Ipfp seed
  K <- length(marg.probs)
  K2 <- 2^K
  seed <- array(rep(1, K2), dim = rep(2, K))
  
  # initialize the Ipfp constraints
  target.data <- list()
  target.list <- list()
  
  # fixing the marginal probabilities (Ipfp constraints)
  for (k in 1:K) {
    target.data[[k]] <- c(marg.probs[k], 1 - marg.probs[k])
    target.list[[k]] <- k
  }
  
  # fixing the pairwise probabilities (Ipfp constraints)
  index <- K  
  # ... loop over the rows
  for (i in 1:(K - 1)) {      
    # ... loop over the columns
    for (j in (i + 1):K) {      
      index <- index + 1
      target.data[[index]] <- array(c(pair.proba[i, j], NA, NA, NA), 
                                    dim = c(2, 2))
      target.list[[index]] <- c(i, j)      
    }    
  }
  
  # obtain the joint-probabilities
  joint.proba <- Ipfp(seed = seed, target.list = target.list, 
                      target.data = target.data, na.target = TRUE, ...)
  
  # gathering the results
  results <- list("joint.proba" = joint.proba$p.hat, 
                  "stp.crit" = joint.proba$stp.crit, "conv" = joint.proba$conv,
                  "check.margins" = joint.proba$check.margins,
                  "label" = lab)
  
  # returning the generated joint-distribution
  return(results)
  
}

RMultBinary <- function(n = 1, mult.bin.dist, target.values = NULL) {
  # Simulate multivariate Bernoulli distribution.
  # 
  # This function generates a sample from a multinomial distribution of K 
  # dependent binary (Bernoulli) variables (X_1, X_2, ..., X_K) defined by an
  # array (of 2^K cells) of joint-probabilities.
  #
  # Author: T. Suesse
  #  
  # Args:
  #   n: Desired sample size. Default = 1.
  #   mult.bin.dist: A list detailing the multivariate distribution containing
  #                  - joint.proba:
  #                    an array detailing the joint-probabilities of the K 
  #                    binary variables. The array has K dimensions of size 2, 
  #                    referring to the 2 possible outcome of the considered
  #                    variable). Hence, the total number of elements is 2^K;
  #                  - var.labels (optionnal)
  #                    the names of the K variables.
  #                  This list can be generated via the ObtainMultBinaryDist
  #                  function.
  #   target.values: A list describing the possibles outcomes of each binary
  #                  variable. By default = {0, 1}.
  #     
  # Returns: A list whose elements are
  #   binary.sequences: The generated K x n random sequence.
  #   possible.binary.sequences: The possible binary sequences, i.e. the domain.
  #   chosen.random.index: The index of the random draws in the domain.
  
  if (is.null(mult.bin.dist) == TRUE) {
    stop("A list detailing the multivariate binary distribution to be simulated
         must be provide! The list must contain at least an array detailing the
         joint-probabilities (see ObtainMultBinaryDist)\n")
  }
  
  if (is.null(mult.bin.dist$joint.proba) == TRUE) {
    stop("the mult.bin.dist list must at least contain an element joint.proba
         (an array detailing the joint-probabilities). See 
         ObtainMultBinaryDist to generate an appropriate structure.")   
  }
  
  # retrieving the variables names
  if (is.null(mult.bin.dist$label) == FALSE) {
    var.labels = mult.bin.dist$label    
  } else {
    var.labels <- NULL
  }
  
  array.prob = mult.bin.dist$joint.proba
  dims <- dim(array.prob)
  K <- length(dims)
  
  # initialization of the possible outcome for each attribute
  if (is.null(target.values) == TRUE) {
    target.values <- vector("list", K)    
    target.values[] <- list(c(1,0))
  }
  
  # gives matrix of size 2^K times n 
  y <- rmultinom(n = n, size = 1, prob = c(array.prob))
  
  # gives observations 1:2^K
  y <- apply(y, 2, which.max)
  
  # transform this number to binary sequence
  domain <- expand.grid(target.values)
  binary.sequences <- as.matrix(domain[y, ])
  rownames(binary.sequences) <- NULL  
  colnames(binary.sequences) <- var.labels
  
  # returning the n random sequences of size K
  return(list("binary.sequences" = binary.sequences, "chosen.random.index" = y, 
              "possible.binary.sequences" = domain))
  
}
