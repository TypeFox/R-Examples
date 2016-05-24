# File mipfp/R/models.R
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
# This file provides a function to estimate a contingency table using several
# model-based approaches.
# ------------------------------------------------------------------------------

ObtainModelEstimates <- function(seed, target.list, target.data, method = "ml", 
                                 tol.margins = 1e-10, 
                                 replace.zeros = 1e-10, ...) {
  # Estimates N-way tables using max likelihood, min chi2 or least squares.
  # 
  # This function provides several alternative estimating methods to the IPFP 
  # when estimating multiway table subject to known constrains/totals: Maximum 
  # likelihood method (ML), minimum chi-squared (CHI2) and weighted least 
  # squares (LSQ).
  # 
  # Author: T. Suesse - modified by J. Barthelemy
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must 
  #         be non-negative.
  #   target.list: A list of the target margins provided in target.data. Each 
  #                component of the list is an array whose cells indicates which
  #                dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each 
  #                component of the list is an array storing a margin. The list 
  #                order must follow the one defined in target.list. Note that 
  #                the cells of the arrays must be non-negative.
  #   method: by default "ml" (Maximum likelihood), other options "chi2" (
  #           minimum chi-squared) and "lsq" (least squares).
  #   replace.zeros: constant that is added to zero cell counts, as the 
  #                  procedures require strictly positive cell counts.
  #   tol.margins: The tolerance for margins consistency.
  #   ...: Additional parameters that can be passed to control the optimisation
  #        process.
  #     
  # Returns: A mipfp object consisting of a list whose elements are
  #   call: A call object in which all the specified arguments are given by
  #         their full names.
  #   method: The selected method for estimation.
  #   pi.hat: Array of the estimated table probabilities.
  #   xi.hat: Array of the estimated table frequencies.  
  #   error.margins: for each list element of target.data, error.margins shows 
  #                  the maximum absolute value of A * pi.hat - margins.vector.
  #                  The elements should approximate zero, otherwise margins are
  #                  not met.
  #   solnp.res: For optimisation it uses the R package Rsolnp and solnp is the 
  #              corresponding object returned by Rsolnp.
  
  # checking if a seed is provided
  if (is.null(seed) == TRUE) {
    stop('Error: no seed specified!')
  }
  
  # checking if target are provided
  if (is.null(target.data) == TRUE | is.null(target.list) == TRUE) {
    stop('Error: target.data and/or target.list not specified!')
  }
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min)) < 0.0 | min(seed) < 0.0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }
  
  # checking if NA in target cells
  if (is.na(min(sapply(target.data, min))) == TRUE)  {
    stop('Error: NA values present in the margins')
  }
  
  # checking if NA in seed
  if (is.na(min(seed)) == TRUE) {
    stop('Error: NA values present in the seed!')
  }
  
  # checking the strict positivy of replace.zeros
  if (replace.zeros <= 0.0) {
    stop('Error: replace.zeros must be strictly positive!')
  }
  
  # checking the margins consistency if no missing values in the targets
  error.margins <- TRUE
  if (length(target.data) > 1) {
    for (m in 2:length(target.data)) {      
      if (abs(sum(target.data[[m-1]]) - sum(target.data[[m]])) > tol.margins) {   
        error.margins <- FALSE
        warning('Target not consistents - shifting to probabilities!
                  Check input data!\n')
        break  
      }      
    }
  }

  # if margins are not consistent, shifting from frequencies to probabilities
  if (error.margins == FALSE) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }  
  
  # create vector version of the seed compatible with the function 
  # 'MarginalMatrix' from cmm that requires that the last index moves fastest
  seed.vector <- Array2Vector(seed)
    
  # dimensions of the problem
  n.sets.margins <- length(target.list)    
  K <- dim(seed)
  K.length <- length(K)
  K.prod <- prod(K)
  n <- sum(seed.vector)
  target.n = sum(target.data[[1]])
  
  # scaling input to probabilities and removing 0 cells
  seed.vector[seed.vector == 0] <- replace.zeros * n
  seed.prop.vector <- seed.vector / n  
        
  # generation of the constraints matrix A such that A' * pi.hat = target.data
  marg.list <- ComputeA(dim.arr = dim(seed), target.data = target.data, 
                        target.list = target.list)
  A <- t(marg.list$marginal.matrix)
  margins.vector <- marg.list$margins
  dim.A <- dim(A)
  
  # defining some functions used by the optimisation
  # ... log-likelihood function
  fun.loglik <- function(p) {
    return(- sum(seed.vector * log(p)))
  }
  
  # ... chi-square function
  fun.chisq <- function(p) {
    return(sum((p - seed.prop.vector)^2 / p))
  }
  
  # ... least-square function
  fun.lsq <- function(p) {
    return(sum((p - seed.prop.vector)^2 / seed.prop.vector))
  }
  
  # ... equality constraints
  eqfun1 <- function(p) {    
    return(A %*% p - margins.vector)
  }
  
  # switching to the desired user method
  # ... converting method to a integer index
  method.num <- switch(method, ml = 1, chi2 = 2, lsq = 3, 4)
  if (method.num == 4) {
    warning("'method' must be 'ml', 'chi2' or 'lsq', switching to 'ml'!")
    method.num <- 1
    method <- "ml"
  }
  # ... index determine the method: 1 = ML, 2 = CHI2 and 3 = LSQ
  switch(method.num,
    fun <- fun.loglik,
    fun <- fun.chisq,
    fun <- fun.lsq
  )
  
  # calling of "solnp" of package "Rsolnp"
  rsolnp <- try(Rsolnp::solnp(pars=seed.prop.vector, fun = fun,
                              LB = rep(0 + replace.zeros^1.5, K.prod), 
                              UB = rep(1 - replace.zeros^1.5, K.prod), 
                              eqfun = eqfun1, eqB = c(rep(0, dim(A)[1])),                              
                              control = list(trace = 0, ...)))
  
  # assessing solnp's convergence, returning an error if no convergence
  conv <- TRUE
  if (inherits(rsolnp, "try-error") | rsolnp$convergence > 0 ){
    conv <- FALSE
    warning("No reliable solutions found by solnp!
            Check the delta and tol parameters.
            replace.zeros might be too low!\n")
  }
    
  # saving solution: probabilities and frequencies
  pi.hat <- rsolnp$pars
  pi.hat.array <- Vector2Array(pi.hat, dim.out = K)  
  xi.hat.array <- pi.hat.array * target.n
  
  # copying the names and dimnames
  names(pi.hat.array) <- names(seed)
  names(xi.hat.array) <- names(seed)
  dimnames(pi.hat.array) <- dimnames(seed)
  dimnames(xi.hat.array) <- dimnames(seed)
  
  # computing final max difference between generated and target margins
  error.margins <- vector(mode = "numeric", length = n.sets.margins)  
  for (j in 1:n.sets.margins) {
    error.margins[j] = max(abs(target.data[[j]] 
                          - apply(xi.hat.array, target.list[[j]], sum)))
    if (is.null(names(dimnames(seed))) == FALSE) {
      names(error.margins)[j] <- paste(names(dimnames(seed))[target.list[[j]]],
                                       collapse = ".")  
    }
  }      
  
  # gathering the results
  results <- list("x.hat" = xi.hat.array, "p.hat" = pi.hat.array,
                  "error.margins" = error.margins, "solnp.res" = rsolnp,
                  "conv" = conv)  

  # adding the method applied
  results$method <- method
  
  # adding the calling expression
  results$call <- match.call()
  
  # updating the class attribute  
  class(results) <- c("list", "mipfp")
  
  # returning results
  return(results)

}
