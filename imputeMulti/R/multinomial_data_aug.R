

#' @title Data Augmentation algorithm for multinomial data
#' @description Implement the Data Augmentation algorithm for multvariate multinomial data given
#' observed counts of complete and missing data (\eqn{Y_obs} and \eqn{Y_mis}). Allows for specification
#' of a Dirichlet conjugate prior. 
#' @param x_y A \code{data.frame} of observed counts for complete observations.
#' @param z_Os_y A \code{data.frame} of observed marginal-counts for incomplete observations.
#' @param enum_comp A \code{data.frame} specifying a vector of all possible observed patterns.
#' @param conj_prior A string specifying the conjugate prior. One of 
#' \code{c("none", "data.dep", "flat.prior", "non.informative")}.
#' @param alpha The vector of counts \eqn{\alpha} for a \eqn{Dir(\alpha)} prior. Must be specified if 
#' \code{conj_prior} is either \code{c("data.dep", "flat.prior")}. If \code{flat.prior}, specify 
#' as a scalar. If \code{data.dep}, specify as a vector with key matching \code{enum_comp}.
#' @param burnin A scalar specifying the number of iterations to use as a burnin. Defaults 
#' to \code{500}.
#' @param post_draws An integer specifying the number of draws from the posterior distribution.
#'  Defaults to \code{1000}.
#' @param verbose Logical. If \code{TRUE}, provide verbose output on each iteration.
#' @return An object of class \code{\link{mod_imputeMulti-class}}.
#' @seealso \code{\link{multinomial_em}}, \code{\link{multinomial_impute}}
#' 
#' @examples \dontrun{
#'  data(tract2221)
#'  x_y <- multinomial_stats(tract2221[,1:4], output= "x_y")
#'  z_Os_y <- multinomial_stats(tract2221[,1:4], output= "z_Os_y")
#'  x_possible <- multinomial_stats(tract2221[,1:4], output= "possible.obs")
#'
#'  imputeDA_mle <- multinomial_data_aug(x_y, z_Os_y, x_possible, n_obs= nrow(tract_2221),
#'                      conj_prior= "none", verbose= TRUE)
#' }
#' 
#' @export
multinomial_data_aug <- function(x_y, z_Os_y, enum_comp, 
                                 conj_prior= c("none", "data.dep", "flat.prior", "non.informative"), 
                                 alpha= NULL, burnin= 500, post_draws= 1000,
                                 verbose= FALSE) {
  # check some errors
  conj_prior <- match.arg(conj_prior, several.ok= FALSE)
  if (conj_prior %in% c("data.dep", "flat.prior") & is.null(alpha) ) {
    stop("Please supply argument alpha as prior.")
  }
  
  mc <- match.call()
  z_p <- ncol(z_Os_y)
  count_p <- ncol(enum_comp)
  
  # 01. Merge in prior if supplied; calculate if requested
  #----------------------------------------------
  enum_comp <- check_prior(conj_prior= conj_prior, alpha= alpha, verbose= verbose,
                           outer= FALSE, enum_comp= enum_comp)
  
  # 02. I and P Steps
  #----------------------------------------------
  iter <- 0
  while (iter < burnin) {
    # I Step
    log_lik <- log_lik0 <- 0
    enum_comp$counts <- 0
    
    for (s in 1:nrow(z_Os_y)) { 
      # random allocation of observed marginal counts to complete pattern y
      # (x_y| z_Os_y, theta) = \sum_s (Xsy|Zsy, gamma)
      # (Xsy|Zy_theta) ~ M(Zsy, gamma)
      comp_ind <- marg_complete_compare(z_Os_y[s, -z_p], enum_comp[, 1:count_p], 
                                    marg_to_complete= TRUE) # pattern match to complete
      
      b_Os_y <- sum(enum_comp$theta_y[unlist(comp_ind)])
      
      E_Xsy_Zy_theta <- as.vector(stats::rmultinom(1, size= z_Os_y$counts[s], 
                                  prob= enum_comp$theta_y[unlist(comp_ind)] / b_Os_y)) # normalized probability
      # expected count += random draw based on marginally-observed
      enum_comp$counts[unlist(comp_ind)] <- enum_comp$counts[unlist(comp_ind)] + E_Xsy_Zy_theta
      
      # update log-lik
      if (b_Os_y > 0) {
        log_lik <- log_lik + z_Os_y$counts[s] * log(b_Os_y)
      }
    }
    # expected count += observed counts
    enum_comp$counts[as.integer(rownames(x_y))] <- enum_comp$counts[as.integer(rownames(x_y))] + x_y$counts
    # update log-lik
    log_lik <- log_lik + sum(ifelse(enum_comp$theta_y[as.integer(rownames(x_y))] == 0, 0,
                                    x_y$counts * log(enum_comp$theta_y[as.integer(rownames(x_y))])))
    
    # P Step
    if (conj_prior == "none") {
      # in case of random zeros: use non-informative prior
      enum_comp$theta_y1 <- as.vector(gtools::rdirichlet(n=1, alpha= enum_comp$counts + 1))
    } else {
      enum_comp$theta_y1 <- as.vector(gtools::rdirichlet(n=1, alpha= enum_comp$counts + enum_comp$alpha))
    }
    
    ### update iteration; print likelihood if verbose
    iter <- iter + 1
    if (verbose) {
      cat("Iteration", iter, ": log-likelihood =", sprintf("%.10f", log_lik), "... \n")
    }
    
    enum_comp$theta_y <- enum_comp$theta_y1
    log_lik0 <- log_lik
  }
  
  # 03. if iter >= max_iter, exit
  # MLE for theta_y is taken to be the mean of n= post_draws draws from the 
  # posterior distribution
  #----------------------------------------------
  # update log-lik for prior
  if (conj_prior != "none") {
    log_lik <- log_lik + sum(ifelse(enum_comp$alpha == 0 | enum_comp$theta_y == 0, 0,
                                    enum_comp$alpha * log(enum_comp$theta_y)))
  }
  
  if (conj_prior == "none") {
    # in case of random zeros: use non-informative prior
    theta_post <- gtools::rdirichlet(n= post_draws, alpha= enum_comp$counts + 1) 
    enum_comp$theta_y <- colMeans(theta_post)
  } else {
    theta_post <- gtools::rdirichlet(n= post_draws, alpha= enum_comp$counts + enum_comp$alpha)
    enum_comp$theta_y <- colMeans(theta_post)
  }
  
  enum_comp$theta_y1 <- NULL
  enum_comp$counts <- NULL
  
  mod <- methods::new("mod_imputeMulti",
             method= "DA",
             mle_call= mc,
             mle_iter= iter,
             mle_log_lik= log_lik,
             mle_cp= conj_prior,
             mle_x_y= enum_comp)
  
  # mod <- list(method= "DA",
  #             mle_call= mc,
  #             mle_iter= iter,
  #             mle_log_lik= log_lik,
  #             mle_cp= conj_prior,
  #             mle_x_y= enum_comp)
  # class(mod) <- "mod_imputeMulti"
  
  return(mod)
}