#' Generate a list of hyperparameters
#' 
#' Generates a list of hyperparameters for use in \code{hcmm_impute}. Specifying
#' only \code{hcmmdat} or \code{q} AND \code{cx} will generate default values
#' (see citation). 
#'
#' @param hcmmdat An \code{hcmm_data} object
#' @param q The number of continuous variables
#' @param cx A length p vector (where p is the number of categorical variables). cx[j]
#' is the number of distinct values taken by X[,j]
#' @param alpha_a,alpha_b Gamma prior on top-level concentration parameter,
#' where E(alpha) = alpha_a/alpha_b
#' @param beta_x_a,beta_x_b Gamma prior on X model concentration parameter
#' @param beta_y_a,beta_y_b Gamma prior on Y model concentration parameter
#' @param tau_a,tau_b Gamma prior on coefficient precision parameters
#' @param v,w Degree of freedom parameters in the hierarchical inverse-Wishart/Wishart prior
#' @param Sigma0 Centering matrix in the hierarchical inverse-Wishart/Wishart prior
#' @param gamma Parameter of the symmetric Dirichlet priors in the product multinomial kernel. 
#' (Should be a length p vector.) 
#' @param sigma2_0beta Variance of the prior on B0
#'
#' @return A list of hyperparameters
#' @export
#'
hcmm_hyperpar = function( hcmmdat = NULL,
                          q = ncol(hcmmdat$Y),
                          cx = hcmmdat$cx,
                          alpha_a = 0.5,
                          alpha_b = 0.5,
                          beta_x_a = 0.5,
                          beta_x_b = 0.5,
                          beta_y_a = 0.5,
                          beta_y_b = 0.5,
                          tau_a = 0.5,
                          tau_b = 0.5,
                          v = q+1,
                          w = q+2,
                          Sigma0 = diag(1, q)/v,
                          gamma = 1/cx,
                          sigma2_0beta = 10) {
  
  
    return(list(alpha_a = alpha_a,
                alpha_b = alpha_b,
                beta_x_a = beta_x_a,
                beta_x_b = beta_x_b,
                beta_y_a = beta_y_a,
                beta_y_b = beta_y_b,
                tau_a = tau_a,
                tau_b = tau_b,
                v = v,
                w = w,
                sigma2_0beta = sigma2_0beta,
                Sigma0 = Sigma0,
                gamma = gamma))
}