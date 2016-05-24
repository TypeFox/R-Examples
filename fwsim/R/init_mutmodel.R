################################################################################
# PARAMETERS
################################################################################
# modeltype:  Integer: 1 = SMM (traditional)
#                      2 = LMM (Arne Jochens)
#                      3 = EMM (Svante)
#
# mutpars:    Matrix:  Mutation parameters for each locus
#
#                      Mutmodel 1 (SMM): 2 parameters per locus
#                        P(i -> i-1) = mu_d
#                        P(i -> i+1) = mu_u
#                        P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
#                                    = 1 - mu_d - mu_u
#
#                        mutpars[1, locus]: mu_d
#                        mutpars[2, locus]: mu_u
#
#                      Mutmodel 2 (LMM): 6 parameters per locus
#                        P(i -> i-1) = gamma_d / (1 + exp(alpha_d*(beta_d - i)))
#                        P(i -> i+1) = gamma_u / (1 + exp(alpha_u*(beta_u - i)))
#                        P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
#
#                        mutpars[1, locus]: gamma_d
#                        mutpars[2, locus]: alpha_d
#                        mutpars[3, locus]: beta_d
#                        mutpars[4, locus]: gamma_u
#                        mutpars[5, locus]: alpha_u
#                        mutpars[6, locus]: beta_u
#
#                      Mutmodel 3 (EMM): 4 parameters per locus
#                        P(i -> i-1) = 1/((1 + exp(a + b*i))*(1 + exp(alpha + beta*i)))
#                        P(i -> i+1) = exp(alpha + beta*i)/((1 + exp(a + b*i))*(1 + exp(alpha + beta*i)))
#                        P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
#                                    = exp(a + b*i)/(1 + exp(a + b*i))
#
#                        mutpars[1, locus]: a
#                        mutpars[2, locus]: b
#                        mutpars[3, locus]: alpha
#                        mutpars[4, locus]: beta
#
#             Vector:  If a vector, the same values are used for all loci.

############################################
# RETURN VALUE
############################################
# mutmodel:           the mutation model

init_mutmodel <- function(
  modeltype = 1L, 
  mutpars = NULL, ...) {

  if (is.null(modeltype) || length(modeltype) != 1L || !is.integer(modeltype) || !(modeltype %in% c(1L:3L))) {
    stop("modeltype must be 1L, 2L, or 3L")
  }
   
  if (is.null(mutpars)) {
    stop("mutpars must be a matrix or a vector")
  }

  if (is.null(nrow(mutpars)) | is.null(ncol(mutpars))) {
    stop("Expected a matrix with mutation parameters")
  }  
  
  if (is.null(colnames(mutpars))) {
    colnames(mutpars) <- paste("Locus", 1L:ncol(mutpars), sep = "")
  }
  
  if (modeltype == 1L) {
    rownames(mutpars) <- c("mu_d", "mu_u")
  } else if (modeltype == 2L) {
    rownames(mutpars) <- c("gamma_d", "alpha_d", "beta_d", "gamma_u", "alpha_u", "beta_u")
  } else if (modeltype == 3L) {
    rownames(mutpars) <- c("a", "b", "alpha", "beta")
  } else {
    stop("Unexpected model type")
  }
  
  ans <- list(modeltype = modeltype, mutpars = mutpars)
  
  class(ans) <- c("mutmodel", class(ans))
  
  return(ans)  
}

