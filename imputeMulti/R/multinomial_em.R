

#' @title EM algorithm for multinomial data
#' @description Implement the EM algorithm for multvariate multinomial data given
#' observed counts of complete and missing data (\eqn{Y_obs} and \eqn{Y_mis}). Allows for specification
#' of a Dirichlet conjugate prior.
#' @param x_y A \code{data.frame} of observed counts for complete observations.
#' @param z_Os_y A \code{data.frame} of observed marginal-counts for incomplete observations.
#' @param enum_comp A \code{data.frame} specifying a vector of all possible observed patterns.
#' @param n_obs An integer specifying the number of observations in the original data.
#' @param conj_prior A string specifying the conjugate prior. One of
#' \code{c("none", "data.dep", "flat.prior", "non.informative")}.
#' @param alpha The vector of counts \eqn{\alpha} for a \eqn{Dir(\alpha)} prior. Must be specified if
#' \code{conj_prior} is either \code{c("data.dep", "flat.prior")}. If \code{flat.prior}, specify
#' as a scalar. If \code{data.dep}, specify as a vector with key matching \code{enum_comp}.
#' @param tol A scalar specifying the convergence criteria. Defaults to \code{5e-7}
#' @param max_iter An integer specifying the maximum number of allowable iterations. Defaults
#' to \code{10000}.
#' @param verbose Logical. If \code{TRUE}, provide verbose output on each iteration.
#' @return An object of class \code{\link{mod_imputeMulti-class}}.
#' @seealso \code{\link{multinomial_data_aug}}, \code{\link{multinomial_impute}}
#' 
#' @examples \dontrun{
#'  data(tract2221)
#'  x_y <- multinomial_stats(tract2221[,1:4], output= "x_y")
#'  z_Os_y <- multinomial_stats(tract2221[,1:4], output= "z_Os_y")
#'  x_possible <- multinomial_stats(tract2221[,1:4], output= "possible.obs")
#'
#'  imputeEM_mle <- multinomial_em(x_y, z_Os_y, x_possible, n_obs= nrow(tract_2221),
#'                      conj_prior= "none", verbose= TRUE)
#' }
#' 
#' @export
multinomial_em <- function(x_y, z_Os_y, enum_comp, n_obs,
                           conj_prior= c("none", "data.dep", "flat.prior", "non.informative"),
                           alpha= NULL, tol= 5e-7, max_iter= 10000,
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

  # 02. E and M Steps
  #----------------------------------------------
  iter <- 0
  while (iter < max_iter) {
    # E Step
    log_lik <- log_lik0 <- 0
    enum_comp$counts <- 0

    for (s in 1:nrow(z_Os_y)) {
      # allocate observed marginal counts proportionally to complete patterns
      # E(x_y| z_Os_y, theta) = \sum_s [E_Xsy_Zy_theta]
      # E_Xsy_Zy_theta = (z_Os_y * theta_y) / b_Os_y
      comp_ind <- marg_complete_compare(z_Os_y[s, -z_p], enum_comp[, 1:count_p],
                                    marg_to_complete= TRUE) # pattern match to complete

      b_Os_y <- sum(enum_comp$theta_y[unlist(comp_ind)])
      E_Xsy_Zy_theta <- z_Os_y$counts[s] * enum_comp$theta_y[unlist(comp_ind)] / b_Os_y # normalize

      # expected count += proportional marginally-observed
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

    # M Step
    if (conj_prior == "none") {
      enum_comp$theta_y1 <- enum_comp$counts / n_obs
    } else {
      D <- nrow(enum_comp)
      alpha_0 <- sum(enum_comp$alpha)
      enum_comp$theta_y1 <- (enum_comp$counts + enum_comp$alpha - 1) / (n_obs + alpha_0 - D)
    }

    # update iteration; print likelihood if verbose
    iter <- iter + 1
    if (verbose) {
      cat("Iteration", iter, ": log-likelihood =", sprintf("%.10f", log_lik), "\n Convergence Criteria =",
                sprintf("%.10f", supDist(enum_comp$theta_y, enum_comp$theta_y1)), "... \n")
    }

  # 03. check convergence to exit and return
  #----------------------------------------------
    if (supDist(enum_comp$theta_y, enum_comp$theta_y1) < tol |
        abs(log_lik - log_lik0) < tol * 100) {
      # update log-lik for prior
      if (conj_prior != "none") {
        log_lik <- log_lik + sum(ifelse(enum_comp$alpha == 0 | enum_comp$theta_y == 0, 0,
                                        enum_comp$alpha * log(enum_comp$theta_y)))
      }
      enum_comp$theta_y1 <- NULL
      enum_comp$counts <- NULL

      mod <- methods::new("mod_imputeMulti",
                 method= "EM",
                 mle_call= mc,
                 mle_iter= iter,
                 mle_log_lik= log_lik,
                 mle_cp= conj_prior,
                 mle_x_y= enum_comp)
      
      # mod <- list(method= "EM", 
      #             mle_call= mc,
      #             mle_iter= iter,
      #             mle_log_lik= log_lik,
      #             mle_cp= conj_prior,
      #             mle_x_y= enum_comp)
      # 
      # class(mod) <- "mod_imputeMulti"

      return(mod)
    } else {
    enum_comp$theta_y <- enum_comp$theta_y1
    log_lik0 <- log_lik
    }
  }
  # 04. if iter >= max_iter, exit
  #----------------------------------------------
  # update log-lik for prior
  if (conj_prior != "none") {
    log_lik <- log_lik + sum(ifelse(enum_comp$alpha == 0 | enum_comp$theta_y == 0, 0,
                                    enum_comp$alpha * log(enum_comp$theta_y)))
  }
  enum_comp$theta_y1 <- NULL
  enum_comp$counts <- NULL

  mod <- methods::new("mod_imputeMulti",
             method= "EM",
             mle_call= mc,
             mle_iter= iter,
             mle_log_lik= log_lik,
             mle_cp= conj_prior,
             mle_x_y= enum_comp)
  
  # mod <- list(method= "EM", 
  #             mle_call= mc,
  #             mle_iter= iter,
  #             mle_log_lik= log_lik,
  #             mle_cp= conj_prior,
  #             mle_x_y= enum_comp)
  # class(mod) <- "mod_imputeMulti"

  return(mod)
}
