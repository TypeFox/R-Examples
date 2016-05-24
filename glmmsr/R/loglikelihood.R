#' Find the log-likelihood function
#' @param modfr a model frame, the output of \code{find_modfr_glmm}
#' @inheritParams glmm
#' @export
find_lfun_glmm <- function(modfr, method, control = NULL)
{
  con <- find_control_with_defaults(control, method)
  find_lfun_glmm_internal(modfr, method, con)
}


find_lfun_glmm_internal <- function(modfr, method, control)
{
  devfun_lme4 <- lme4::mkGlmerDevfun(fr = modfr$fr, X = modfr$X,
                                     reTrms = modfr$reTrms, family = modfr$family,
                                     control = lme4_control())
  nAGQ_lme4 <- ifelse(method == "AGQ", control$nAGQ, 1)
  devfun_lme4 <- lme4::updateGlmerDevfun(devfun_lme4, modfr$reTrms, nAGQ = nAGQ_lme4)
  lfun_lme4 <- function(x) { -devfun_lme4(x) / 2 }
  switch(method,
         Laplace = lfun_lme4,
         AGQ = lfun_lme4,
         SR = find_lfun_SR(modfr, devfun_lme4,
                           nSL = control$nSL,
                           nAGQ = control$nAGQ),
         IS = find_lfun_IS(modfr, devfun_lme4, nIS = control$nIS),
         stop(cat("method", method, "not available"))
  )
}

find_lfun_SR <- function(modfr, devfun_lme4, nSL, nAGQ) {
  check_modfr_SR(modfr)

  n_fixed <- ncol(modfr$X)
  n_random <- length(modfr$reTrms$theta)

  factorization_terms <- find_factorization_terms(modfr)
  calibration_pars <- calibration_parameters()
  calibration_pars$n_quadrature_points <- nAGQ
  calibration_pars$n_sparse_levels <- nSL
  calibration_pars$family <- modfr$family$family
  calibration_pars$link <- modfr$family$link
  beliefs <- cluster_graph(factorization_terms)

  complexity_limit <- 5e5
  if(beliefs$width^(2 * nSL) > complexity_limit) {
    stop(paste("The sequential reduction approximation with", nSL,
               "sparse levels \n is too difficult to compute in this case. \n",
               "Consider reducing nSL, or using a different approximation method."),
         call. = FALSE)
  }

  lfun <- function(pars) {
    if(length(pars) != (n_fixed + n_random)) {
      stop("cannot compute loglikelihood for parameter of length ",
           length(pars), " != ", n_fixed + n_random,
           call. = FALSE)
    }
    calibration_pars$theta <- pars[1:n_random]
    calibration_pars$beta <- pars[-(1:n_random)]
    normal_approx <- compute_normal_approx(pars, devfun_lme4)
    l_SR <- beliefs$compute_log_normalizing_constant(normal_approx$mean,
                                                     normal_approx$precision,
                                                     calibration_pars)
    if(l_SR > normal_approx$l_laplace / 2)
      stop("Suspect sequential reduction approximation is unstable at that parameter value",
           call. = FALSE)
    l_SR
  }
  lfun
}

compute_normal_approx <- function(pars, devfun_lme4)
{
  l_laplace <- -devfun_lme4(pars)/2
  PR <- get("pp", environment(devfun_lme4))
  L_tot <- Matrix::expand(PR$L())
  L <- L_tot$L
  P <- L_tot$P
  precision <- Matrix::t(P)%*%Matrix::tcrossprod(L)%*%P
  mean <- PR$delu
  return(list(l_laplace = l_laplace, mean = mean, precision = precision))
}
