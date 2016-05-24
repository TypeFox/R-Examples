################################################################################
# PARAMETERS
################################################################################
# G:          Integer: Number of generations to simulate
#
# H0:         Matrix:  The haplotypes of the initial population.
#                      The number of loci is the number of columns of H0.
#
# N0:         Vector:  The i'th element is the count of the 
#                      haplotype H0[i, ] in the initial population.
#                      sum(N0) is the size of initial population
#
# mutmodel:   mutmodel
#
# alpha:      Vector:  Vector of length G of growth factors (1 correspond to expected 
#                      constant population size)
#             Scalar:  If a scalar, the same value is used for all generations
#
# SNP         Bool:    Make alleles modulus 2 to immitate SNPs
#
# save_generations:    NULL:    No intermediate population will be saved.
#             Vector:  A vector of the generation numbers to save
#
# progress:   Bool:    Show progress of simulations
# trace:      Bool:    Show trace of simulations (more verbose than progress)

############################################
# RETURN VALUE
############################################
# pars:               The parameters used to create the result
# saved.populations:  intermediate populations if chosen by setting 
#                     corresponding values of save_generations (see above)
# population:         the haplotypes and counts of the end population
# pop.sizes:          the population size for each generation
# expected.pop.sizes: the expected population sizes

fwsim <- function(
  G, 
  H0, 
  N0,  
  mutmodel, 
  alpha = 1.0,
  SNP = FALSE,
  save_generations = NULL, 
  progress = TRUE, 
  trace = FALSE, ...) {

  if (is.null(G) || length(G) != 1L || !is.integer(G) || G <= 0L) {
    stop("G must be an integer >= 1L (note the postfix L)")
  }

  if (is.null(H0)) {
    stop("H0 must be a vector or matrix of integers")
  }

  if (!is.matrix(H0)) {
    if (!is.integer(H0)) {
      stop("H0 must be a vector or matrix of integers")
    }
    
    H0_names <- NULL
    
    if (!is.null(names(H0))) {
      H0_names <- names(H0)
    }
    
    H0 <- matrix(H0, nrow = 1)
    
    if (!is.null(H0_names)) {
      colnames(H0) <- H0_names
    }
  } 
  
  if (ncol(H0) == 0L || nrow(H0) == 0L || !is.integer(H0)) {
    stop("H0 must be a vector or matrix of integers")
  }
    
  r <- ncol(H0)  
  
  if (is.null(N0) || !is.integer(N0) || length(N0) != nrow(H0) || any(N0 <= 0L)) {
    stop("N0 must be a integer vector with length corresponding to the number of rows of H0 and have elements >= 1L (note the postfix L)")
  }
  
  #if (is.null(mutmodel) || !is(mutmodel, "mutmodel")) {
  #  stop("mutmodel must be a mutmodel created with init_mutmodel")
  #}
  
  if (is.null(mutmodel)) {
    stop("mutmodel must be a mutmodel created with init_mutmodel or numeric vector of mutation probabilities")
  }
  
  if (!is(mutmodel, "mutmodel")) {
    if (length(mutmodel) <= 0L || !is.numeric(mutmodel) || any(mutmodel < 0) || any(mutmodel > 1)) {
      stop("If mutmodel is a numeric vector, it must have elements between 0 and 1 (mutation probabilities)")
    }
    
    mutpars <-  matrix(c(mutmodel / 2, mutmodel / 2), ncol = length(mutmodel), byrow = TRUE)
    
    if (!is.null(names(mutmodel))) {
      colnames(mutpars) <- names(mutmodel)
    } else if (!is.null(colnames(H0))) {
      colnames(mutpars) <- colnames(H0)
    }
    
    mutmodel <- init_mutmodel(modeltype = 1L, mutpars = mutpars)
  }
    
  if (is.null(ncol(mutmodel$mutpars)) || ncol(mutmodel$mutpars) != r) {
    stop("mutmodel and H0 each specifies different number of loci")
  }

  if (!is.numeric(alpha) || (length(alpha) != 1L && length(alpha) != G)) {
    stop("The growth rate, alpha, must be a numeric vector of size 1L or G")
  }
  if (length(alpha) == 1L) alpha <- rep(alpha, G)
  
  if (is.null(SNP) || !is.logical(SNP)) {
    stop("SNP must be logical/boolean")
  }
  
  if (!is.null(save_generations) && length(save_generations) > 0L) {
    if (!is.integer(save_generations)) stop("save_generations must be a vector of integers.")
    
    save_generations <- sort(unique(save_generations))
    
    if (any(save_generations <= 0L | save_generations >= G)) stop("If not null, save_generations must be numbers between 0L and G, both excluded.")
    
    new.gs <- rep(0L, G - 1L)
    new.gs[save_generations] <- 1L
    save_generations <- c(new.gs, 0L)
  } else {
    save_generations <- rep(0L, G) # For easier handling in C++
  }
  
  if (is.null(progress) || !is.logical(progress) || length(progress) != 1L) {
    stop("progress must one logical/boolean")
  }
  
  if (is.null(trace) || !is.logical(trace) || length(trace) != 1L) {
    stop("trace must one logical/boolean")
  }
  #progress <- as.integer(progress)

	res <- Cpp_fwpopsim(G, H0, N0, alpha, mutmodel, SNP, save_generations, progress, trace)
	
	#################
	expected.pop.sizes <- numeric(G + 1L)
  expected.pop.sizes[1L] <- sum(N0)

  for (generation in 2L:(G+1L)) {
    expected.pop.sizes[generation] <- alpha[generation - 1L]*expected.pop.sizes[generation - 1L]
  }
  
  expected.pop.sizes <- expected.pop.sizes[-1L]
  
	res$expected_pop_sizes <- expected.pop.sizes
	#################
	
	colnames(res$population) <- c(colnames(mutmodel$mutpars), "N")
	colnames(res$pars$H0) <- colnames(mutmodel$mutpars)
	
	for (g in seq_along(res$saved_populations)) {
	  if (is.matrix(res$saved_populations[[g]])) {
	    colnames(res$saved_populations[[g]]) <- c(colnames(mutmodel$mutpars), "N")
	  }
	}
	
	res$population <- as.data.frame(res$population)
	res$size_model <- "stochastic"
	
	class(res) <- c("fwsim", class(res))

  return(res)
}


fwsim_fixed <- function(
  G, 
  H0, 
  N0,  
  mutmodel, 
  SNP = FALSE,
  save_generations = NULL, 
  progress = TRUE, 
  trace = FALSE, ...) {

  if (is.null(G) || length(G) != 1L || !is.integer(G) || G <= 0L) {
    stop("G must be an integer >= 1L (note the postfix L)")
  }

  if (is.null(H0)) {
    stop("H0 must be a vector or matrix of integers")
  }

  if (!is.matrix(H0)) {
    if (!is.integer(H0)) {
      stop("H0 must be a vector or matrix of integers")
    }
    
    H0_names <- NULL
    
    if (!is.null(names(H0))) {
      H0_names <- names(H0)
    }
    
    H0 <- matrix(H0, nrow = 1)
    
    if (!is.null(H0_names)) {
      colnames(H0) <- H0_names
    }
  } 
  
  if (ncol(H0) == 0L || nrow(H0) == 0L || !is.integer(H0)) {
    stop("H0 must be a vector or matrix of integers")
  }
    
  r <- ncol(H0)  
  
  if (is.null(N0) || !is.integer(N0) || length(N0) != nrow(H0) || any(N0 <= 0L)) {
    stop("N0 must be a integer vector with length corresponding to the number of rows of H0 and have elements >= 1L (note the postfix L)")
  }
  
  #if (is.null(mutmodel) || !is(mutmodel, "mutmodel")) {
  #  stop("mutmodel must be a mutmodel created with init_mutmodel")
  #}
  
  if (is.null(mutmodel)) {
    stop("mutmodel must be a mutmodel created with init_mutmodel or numeric vector of mutation probabilities")
  }
  
  if (!is(mutmodel, "mutmodel")) {
    if (length(mutmodel) <= 0L || !is.numeric(mutmodel) || any(mutmodel < 0) || any(mutmodel > 1)) {
      stop("If mutmodel is a numeric vector, it must have elements between 0 and 1 (mutation probabilities)")
    }
    
    mutpars <-  matrix(c(mutmodel / 2, mutmodel / 2), ncol = length(mutmodel), byrow = TRUE)
    
    if (!is.null(names(mutmodel))) {
      colnames(mutpars) <- names(mutmodel)
    } else if (!is.null(colnames(H0))) {
      colnames(mutpars) <- colnames(H0)
    }
    
    mutmodel <- init_mutmodel(modeltype = 1L, mutpars = mutpars)
  }
    
  if (is.null(ncol(mutmodel$mutpars)) || ncol(mutmodel$mutpars) != r) {
    stop("mutmodel and H0 each specifies different number of loci")
  }

  if (is.null(SNP) || !is.logical(SNP)) {
    stop("SNP must be logical/boolean")
  }
  
  if (!is.null(save_generations) && length(save_generations) > 0L) {
    if (!is.integer(save_generations)) stop("save_generations must be a vector of integers.")
    
    save_generations <- sort(unique(save_generations))
    
    if (any(save_generations <= 0L | save_generations >= G)) stop("If not null, save_generations must be numbers between 0L and G, both excluded.")
    
    new.gs <- rep(0L, G - 1L)
    new.gs[save_generations] <- 1L
    save_generations <- c(new.gs, 0L)
  } else {
    save_generations <- rep(0L, G) # For easier handling in C++
  }
  
  if (is.null(progress) || !is.logical(progress) || length(progress) != 1L) {
    stop("progress must one logical/boolean")
  }
  
  if (is.null(trace) || !is.logical(trace) || length(trace) != 1L) {
    stop("trace must one logical/boolean")
  }
  #progress <- as.integer(progress)

	res <- Cpp_fwpopsim_fixed(G, H0, N0, mutmodel, SNP, save_generations, progress, trace)
	
	#################
	expected.pop.sizes <- rep(sum(N0), G)  
	res$expected_pop_sizes <- expected.pop.sizes
	#################
	
	colnames(res$population) <- colnames(mutmodel$mutpars)
	colnames(res$pars$H0) <- colnames(mutmodel$mutpars)
	
	for (g in seq_along(res$saved_populations)) {
	  if (is.matrix(res$saved_populations[[g]])) {
	    colnames(res$saved_populations[[g]]) <- c(colnames(mutmodel$mutpars), "N")
	  }
	}
	
	res$population <- as.data.frame(res$population)	
	res$size_model <- "fixed"
	
	class(res) <- c("fwsim", class(res))

  return(res)
}



