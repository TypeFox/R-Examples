# File mipfp/R/ipfpMultiDim.R
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
# This file provides the functions Ipfp and Ipfp.Covar that respectively
# implements the iterative proportional fitting procedure and a Delta method
# that computes the covariance matrix of the estimate produced by Ipfp.
# ------------------------------------------------------------------------------

Ipfp <- function(seed, target.list, target.data, print = FALSE, iter = 1000, 
                 tol = 1e-10, tol.margins = 1e-10, na.target = FALSE) {
  # Update an array using the iterative proportional fitting procedure.
  #
  # Author: J. Barthelemy
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must
  #         be non-negative.
  #   target.list: A list of the target margins provided in target.data. Each
  #                component of the list is an array whose cells indicates
  #                which dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each
  #                component of the list is an array storing a margin.
  #                The list order must follow the one defined in target.list. 
  #                Note that the cells of the arrays must be non-negative, but
  #                can contains NA values.
  #   print: Verbose parameter: if TRUE prints the current iteration number
  #          and the value of the stopping criterion.
  #   iter: The maximum number of iteration allowed; must be greater than 0.
  #   tol: If the maximum absolute difference between two iteration is lower
  #        than the value specified by tol, then ipfp has reached convergence
  #        (stopping criterion); must be greater than 0.
  #   tol.margins: The tolerance for margins consistency.
  #   na.target: If set to TRUE, allows the targets to have NA cells. In that
  #              case the margins consistency is not checked.
  #
  # Returns: A mipfp object consisting of a list whose elements are
  #   call: A call object in which all the specified arguments are given by
  #         their full names.
  #   method: The selected method for estimation.
  #   stp.crit: The final value of the stopping criterion.
  #   evol.stp.crit: Evolution of the stopping criterion over the iterations.
  #   conv: A boolean indicating whether the algorithm converged to a solution.
  #   error.margins: A list returning, for each margin, the absolute maximum 
  #                  deviation between the target and generated margin.
  
  # checking if a seed is provided
  if (is.null(seed) == TRUE) {
    stop('Error: no seed specified!')
  }
  
  # checking if target are provided
  if (is.null(target.data) == TRUE | is.null(target.data) == TRUE) {
    stop('Error: target.data and/or target.data not specified!')
  }
  
  # checking if NA in target cells if na.target is set to FALSE
  if (is.na(min(sapply(target.data, min))) == TRUE & na.target == FALSE)  {
    stop('Error: NA values present in the margins - use na.target = TRUE!')
  }
  
  # checking if NA in seed
  if (is.na(min(seed)) == TRUE) {
    stop('Error: NA values present in the seed!')
  }
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min), na.rm = na.target) < 0 | min(seed) < 0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }  
  
  # checking the strict positiviy of tol and iter
  if (iter < 1 | tol <= 0.0) {
    stop('Error: tol and iter must be strictly positive!')
  }
  
  # checking if NA allowed and requesting the covariance matrices
  if (na.target == TRUE) { #& compute.cov == TRUE) {
    warning('Missing values allowed in the target margins.
             Computation of the covariance matrices set to FALSE!')
    compute.cov <- FALSE
  }
  
  # checking the margins consistency if no missing values in the targets
  error.margins <- TRUE  
  if (na.target == FALSE) {
    if (length(target.data) > 1) {
      for (m in 2:length(target.data)) {      
        if (abs(sum(target.data[[m-1]]) - sum(target.data[[m]])) > 
            tol.margins) {
          error.margins <- FALSE
          warning('Target not consistents - shifting to probabilities!
                  Check input data!\n')
          break
        }      
      }
    }
  } else {
    if (print == TRUE) {
      cat('NOTE: Missing values present in target cells. ')
      cat('Margins consistency not checked!\n')  
    }        
  }
  
  # if margins are not consistent, shifting from frequencies to probabilities
  if (error.margins == FALSE) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }
  
  if (print == TRUE & error.margins == TRUE & na.target == FALSE) {
    cat('Margins consistency checked!\n')
  }
    
  # initial value is the seed
  result <- seed  
  converged <- FALSE
  tmp.evol.stp.crit <- vector(mode="numeric", length = iter)
  
  # ipfp iterations
  for (i in 1:iter) {
    
    if (print == TRUE) {
      cat('... ITER', i, '\n')
    } 
    
    # saving previous iteration result (for testing convergence)
    result.temp <- result
            
    # loop over the constraints
    for (j in 1:length(target.list)) {
      # ... extracting current margins      
      temp.sum <- apply(result, target.list[[j]], sum)
      # ... computation of the update factor, taking care of 0 and NA cells   
      update.factor <- ifelse(target.data[[j]] == 0 | temp.sum == 0, 0,
                              target.data[[j]] / temp.sum)
      if (na.target == TRUE) {
        update.factor[is.na(update.factor)] <- 1;
      }
      # ... apply the update factor
      result <- sweep(result, target.list[[j]], update.factor, FUN = "*")
    }
    
    # stopping criterion
    stp.crit <- max(abs(result - result.temp))
    tmp.evol.stp.crit[i] <- stp.crit
    if (stp.crit < tol) {
      converged <- TRUE
      if (print == TRUE) {
        cat('       stoping criterion:', stp.crit, '\n')
        cat('Convergence reached after', i, 'iterations!\n')
      } 
      break
    }
    
    if (print == TRUE) {
      cat('       stoping criterion:', stp.crit, '\n')
    }
    
  }
    
  # checking the convergence
  if (converged == FALSE) {
    warning('IPFP did not converged after ', iter, ' iteration(s)! 
            This migh be due to 0 cells in the seed, maximum number 
            of iteration too low or tolerance too small\n')
  }        
  
  # computing final max difference between generated and target margins
  diff.margins <- vector(mode = "numeric", length = length(target.list))
  if (na.target == FALSE) {
    for (j in 1:length(target.list)) {
      diff.margins[j] = max(abs(target.data[[j]] 
                                - apply(result, target.list[[j]], sum)))
      if (is.null(names(dimnames(seed))) == FALSE) {
        names(diff.margins)[j] <- paste(names(dimnames(seed))[target.list[[j]]],
                                        collapse = ".")
      }
    }
  }
  
  # storing the evolution of the stopping criterion
  evol.stp.crit <- tmp.evol.stp.crit[1:i]
  
  # computing the proportions
  result.prop <- result / sum(result)
  
  # gathering the results in a list
  results.list <- list("x.hat" = result, "p.hat" = result.prop, 
                       "conv" = converged, "error.margins" = diff.margins, 
                       "evol.stp.crit" = evol.stp.crit)
  
  # adding the method applied
  results.list$method <- "ipfp"
  
  # adding the calling expression
  results.list$call <- match.call()
  
  # updating the class attribute
  class(results.list) <- c("list", "mipfp")      
  
  # returning the result
  return(results.list)
  
}
