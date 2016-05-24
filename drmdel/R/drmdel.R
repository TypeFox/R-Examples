# ##############################
# This software is written by Song Cai and published under GPLv3.
#
# Version 1.3.1, December 31, 2014.
# ##############################

negLDL <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate negative log dual empirical likelihood for a
# given value of parameter.
#
# logDualLWrapper prototype:
# void logDualLWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_val /*output*/)

  logDL <- .C("logDualLWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_val=double(1))

  return(-logDL$ldl_val)

}

# User specified basis function version
negLDLUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate negative log dual empirical likelihood for a
# given value of parameter.
#
# basis_func must be a function.
#
# logDualLUfWrapper prototype:
# SEXP logDualLUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
#     SEXP h_func, SEXP env, SEXP x)


  val <- .Call("logDualLUfWrapper", as.double(n_total), as.double(n_samples),
        as.double(m), as.double(d), as.double(par),
        basis_func, new.env(), as.double(x))

  return(-val)

}


negLDLGr <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate the gradient of the negative log dual empirical
# likelihood for a given value of parameter.
#
# logDualLDLGrWrapper prototype:
# void logDualLGrWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_gradient /*output*/)


  logDLGr <- .C("logDualLGrWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_gr=double(m*(d+1)))

  return(-logDLGr$ldl_gr)

}

# User specified basis function version
negLDLGrUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate the gradient of the negative log dual empirical
# likelihood for a given value of parameter.
#
# basis_func must be a function.
#
# logDualLGrUfWrapper prototype:
# SEXP logDualLGrUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
#     SEXP h_func, SEXP env, SEXP x /*input*/)

  ldl_gr <- .Call("logDualLGrUfWrapper", as.double(n_total),
                  as.double(n_samples), as.double(m), as.double(d),
                  as.double(par), basis_func, new.env(), as.double(x))

  return(-ldl_gr)

}

gen_par_pos <- function(m, d){

  par_pos_alpha <- seq(from=1, by=(d+1), length.out=m) 
  par_pos_beta <- 1:(m*d) + rep(1:m, each=d)

  return(list(alpha=par_pos_alpha, beta=par_pos_beta))

}

#g_null_jac_full <- function(g_null_jac, par_null, par_alpha,
                            #par_dim, par_dim_null, par_pos) {
## par_dim is the dimension of DRM parameter $theta$, which
##   equals m*(d+1).
## par_dim_null is the dimension of the null mapping, i.e.
##   dim(gamma).
##
## Jacobian matrix for g_null_full (null mapping from (alpha,
## gammma) to (theta_1, theta_2, ..., theta_m))

  #par_dim_null_full <- m + par_dim_null
  #g_null_full_jac <- matrix(rep(0, par_dim*par_dim_null_full),
                                #nrow=par_dim,
                                #ncol=par_dim_null_full)

  ##par_gamma <- par_null_full[-(1:m)]  # get gamma
  #g_null_full_jac[par_pos$beta, (m+1):par_dim_null_full] <- g_null_jac(par_null_full[-(1:m)])
  #g_null_full_jac[par_pos$alpha, 1:m] <- diag(m) 

  #return(g_null_full_jac)
#}

#g_null_full <- function(g_null, par_null, par_alpha,
                        #par_dim, par_pos) {

  #par <- numeric(par_dim)

  #par_beta <- g_null(par_null)

  #par[par_pos$alpha] <- par_alpha
  #par[par_pos$beta] <- par_beta

  #return(par)

#}

negLDL_null <- function(par_null_full, g_null, g_null_jac=NULL,
                        par_pos, par_dim, par_dim_null=NULL,
                        x, n_total, n_samples, m, model, d) {
# g_null is a null mapping from gamma to beta, but
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#   alpha: the first m elements of par_null_full
#   gamma: the rest (from (m+1)th to the last) elements of
#     par_null_full
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par[par_pos$beta] <- g_null(par_null_full[-(1:m)])  # get beta

  return(negLDL(par=par, x=x, n_total=n_total, n_samples=n_samples, m=m,
         model=model, d=d))

}

# User specified basis function version
negLDLUf_null <- function(par_null_full, g_null, g_null_jac=NULL,
                          par_pos, par_dim, par_dim_null=NULL,
                          x, n_total, n_samples, m, basis_func, d) {
# basis_func must be a function.
#
# g_null is a null mapping from gamma to beta, but
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#   alpha: the first m elements of par_null_full
#   gamma: the rest (from (m+1)th to the last) elements of
#     par_null_full
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par[par_pos$beta] <- g_null(par_null_full[-(1:m)])  # get beta

  return(negLDLUf(par=par, x=x, n_total=n_total, n_samples=n_samples, m=m,
         basis_func=basis_func, d=d))

}

negLDLGr_null <- function(par_null_full, g_null, g_null_jac,
                          par_pos, par_dim, par_dim_null,
                          x, n_total, n_samples, m, model, d) {
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#
# g_null is a null mapping from gamma to beta.
# g_null_jac is the jacobian matrix of g_null, i.e. a matrix
#   of dimension m*d by dim(gamma).
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par_gamma <- par_null_full[-(1:m)]  # get gamma
  par[par_pos$beta] <- g_null(par_gamma)

  par_full_gr <- negLDLGr(par=par, x=x, n_total=n_total, n_samples=n_samples,
                          m=m, model=model, d=d)

  # Jacobian matrix for g_null_full (null mapping from (alpha,
  # gammma) to (theta_1, theta_2, ..., theta_m))
  par_dim_null_full <- m + par_dim_null
  g_null_full_jac <- matrix(rep(0, par_dim*par_dim_null_full),
                                nrow=par_dim,
                                ncol=par_dim_null_full)

  if (par_dim_null > 0) {
    g_null_full_jac[par_pos$beta, (m+1):par_dim_null_full] <- g_null_jac(par_gamma)
  }
  g_null_full_jac[par_pos$alpha, 1:m] <- diag(m) 

  return(as.numeric(t(g_null_full_jac) %*% par_full_gr))

}

# User specified basis function version
negLDLGrUf_null <- function(par_null_full, g_null, g_null_jac, 
                            par_pos, par_dim, par_dim_null,
                            x, n_total, n_samples, m, basis_func, d) {
# basis_func must be a function.
#
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#
# g_null is a null mapping from gamma to beta.
# g_null_jac is the jacobian matrix of g_null, i.e. a matrix
#   of dimension m*d by dim(gamma).
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par_gamma <- par_null_full[-(1:m)]  # get gamma
  par[par_pos$beta] <- g_null(par_gamma)

  par_full_gr <- negLDLGrUf(par=par, x=x, n_total=n_total, n_samples=n_samples,
                            m=m, basis_func=basis_func, d=d)

  # Jacobian matrix for g_null_full (null mapping from (alpha,
  # gammma) to (theta_1, theta_2, ..., theta_m))
  par_dim_null_full <- m + par_dim_null
  g_null_full_jac <- matrix(rep(0, par_dim*par_dim_null_full),
                                nrow=par_dim,
                                ncol=par_dim_null_full)

  if (par_dim_null > 0) {
    g_null_full_jac[par_pos$beta, (m+1):par_dim_null_full] <- g_null_jac(par_gamma)
  }
  g_null_full_jac[par_pos$alpha, 1:m] <- diag(m) 

  return(as.numeric(t(g_null_full_jac) %*% par_full_gr))

}

negLDLHessian <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate the hessian of the negative log dual empirical
# likelihood for a given value of parameter.
#
# logDualLDLHessianWrapper prototype:
# void logDualLHessianWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_hessian /*output*/)


  logDLHessian <- .C("logDualLHessianWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_hessian=double(m*(d+1)*m*(d+1)))

  return( matrix(-(logDLHessian$ldl_hessian), m*(d+1), m*(d+1),
                 byrow=TRUE) )

}

# User specified basis function version
negLDLHessianUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate the hessian of the negative log dual empirical
# likelihood for a given value of parameter.
#
# basis_func must be a function.
#
# logDualLDLHessianUfWrapper prototype:
# SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)

  ldl_hessian <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                       as.double(n_samples), as.double(m), as.double(d),
                       as.double(par), basis_func, new.env(), as.double(x))

  return(-ldl_hessian)

}

Wmat <- function(n_total, n_samples, m, d) {
# Calculate the W matrix needed for calculating a consistently estimated
# asymptotic variance of the MELE of parameters in DRM.
#
# WmatWrapper prototype:
# void WmatWrapper(double * restrict n_total, /*inputs*/
#   double * restrict n_samples, /*inputs*/
#   double * restrict m, double * restrict d,
#   double * restrict W /*output*/)


  Wmat <- .C("WmatWrapper", as.double(n_total),
             as.double(n_samples), as.double(m), as.double(d),
             W=double(m*(d+1)*m*(d+1)))

  return( matrix(Wmat$W, m*(d+1), m*(d+1), byrow=TRUE) )

}

meleVarHat1 <- function(mele, x, n_total, n_samples, m, model, d, W) {
# Calculate an consistently estimated asymptotic variance of the MELE of
# parameters in DRM.
#
# W -- the matrix obtained using the Wmat function.
#
# logDualLDLHessianWrapper prototype:
# void logDualLHessianWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_hessian /*output*/)

  hessianDim <- m*(d+1)

  logDLHessianMdele <- .C("logDualLHessianWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(mele), as.double(model), as.double(x),
              ldl_hessian=double(hessianDim*hessianDim))

  U <- matrix(-(logDLHessianMdele$ldl_hessian)/n_total, hessianDim, hessianDim,
              byrow=TRUE)

  return( solve(U) - W )

}

# User specified basis function version
meleVarHatUf1 <- function(mele, x, n_total, n_samples, m, basis_func, d, W) {
# Calculate an consistently estimated asymptotic variance of the MELE of
# parameters in DRM.
#
# W -- the matrix obtained using the Wmat function.
#
# basis_func must be a function.
#
# logDualLDLHessianUfWrapper prototype:
# SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)

  ldl_hessian <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                       as.double(n_samples), as.double(m), as.double(d),
                       as.double(mele), basis_func, new.env(), as.double(x))

  U <- -ldl_hessian/n_total

  return( solve(U) - W )

}

meleVarHat <- function(mele, x, n_total, n_samples, m, model, d) {
# Calculate an consistently estimated asymptotic variance of the MELE of
# parameters in DRM.
#
# WmatWrapper prototype:
# void WmatWrapper(double * restrict n_total, /*inputs*/
#   double * restrict n_samples, /*inputs*/
#   double * restrict m, double * restrict d,
#   double * restrict W /*output*/)
#
# logDualLDLHessianWrapper prototype:
# void logDualLHessianWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_hessian /*output*/)

  # Calculating W matrix
  Wmat <- .C("WmatWrapper", as.double(n_total),
             as.double(n_samples), as.double(m), as.double(d),
             W=double(m*(d+1)*m*(d+1)))

  W <- matrix(Wmat$W, m*(d+1), m*(d+1), byrow=TRUE)

  # Calculate the Hessian of logDL
  hessianDim <- m*(d+1)

  logDLHessianMdele <- .C("logDualLHessianWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(mele), as.double(model), as.double(x),
              ldl_hessian=double(hessianDim*hessianDim))

  # Observed information matrix
  U <- matrix(-(logDLHessianMdele$ldl_hessian)/n_total, hessianDim, hessianDim,
              byrow=TRUE)

  # Estimated asymptotic variance of the MELE
  return( solve(U) - W )

}

# User specified basis function version
meleVarHatUf <- function(mele, x, n_total, n_samples, m, basis_func, d) {
# Calculate an consistently estimated asymptotic variance of the MELE of
# parameters in DRM.
#
# WmatWrapper prototype:
# void WmatWrapper(double * restrict n_total, /*inputs*/
#   double * restrict n_samples, /*inputs*/
#   double * restrict m, double * restrict d,
#   double * restrict W /*output*/)
#
# basis_func must be a function.
#
# logDualLDLHessianUfWrapper prototype:
# SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)

  # Calculating W matrix
  Wmat <- .C("WmatWrapper", as.double(n_total),
             as.double(n_samples), as.double(m), as.double(d),
             W=double(m*(d+1)*m*(d+1)))

  W <- matrix(Wmat$W, m*(d+1), m*(d+1), byrow=TRUE)

  # Calculate the Hessian of logDL
  ldl_hessian <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                       as.double(n_samples), as.double(m), as.double(d),
                       as.double(mele), basis_func, new.env(), as.double(x))

  # Observed information matrix
  U <- -ldl_hessian/n_total
  #U <- matrix(-ldl_hessian/n_total, m*(d+1), m*(d+1), byrow=TRUE)

  # Estimated asymptotic variance of the MELE
  return( solve(U) - W )

}

pEst <- function(x, n_total, n_samples, m, model, d, mele) {
# Extract estimated probabilities given MELE
#
# probEstWrapper prototype:
# void probEstWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict normalize,
#     double * restrict pEst /*output*/)

  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  pEstimates <- .C("probEstWrapper", as.double(n_total),
                      as.double(n_samples), as.double(m), as.double(d),
                      as.double(mele), as.double(model), as.double(x),
                      1.0,
                      p_est=double((m+1)*n_total))

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))

  return(data.frame(k=population_indicator, x=x, p_est=pEstimates$p_est))

}

# User specified basis function version
pEstUf <- function(x, n_total, n_samples, m, basis_func, d, mele) {
# Extract estimated probabilities given MELE
#
# basis_func must be a function.
#
# probEstUfWrapper prototype:
# SEXP probEstUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x, SEXP normalize /*input*/)

  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  p_est <- .Call("probEstUfWrapper", as.double(n_total),
                    as.double(n_samples), as.double(m), as.double(d),
                    as.double(mele), basis_func, new.env(), as.double(x),
                    1.0)

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))

  return(data.frame(k=population_indicator, x=x, p_est=p_est))

}

pBlEst <- function(x, n_total, n_samples, m, model, d, mele) {
# Extract estimated probabilities given MELE
#
# probBlEstWrapper prototype:
# void probBlEstWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict normalize,
#     double * restrict pBlEst /*output*/)


  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  pBlEstimates <- .C("probBlEstWrapper", as.double(n_total),
                      as.double(n_samples), as.double(m), as.double(d),
                      as.double(mele), as.double(model), as.double(x),
                      1.0,
                      p_bl_est=double(n_total))

  return(data.frame(x=x, p_est=pBlEstimates$p_bl_est))

}

# User specified basis function version
pBlEstUf <- function(x, n_total, n_samples, m, basis_func, d, mele) {
# Extract estimated probabilities given MELE
#
# basis_func must be a function.
#
# probBlEstUfWrapper prototype:
# SEXP probBlEstUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x, SEXP normalize /*input*/)

  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  p_bl_est <- .Call("probBlEstUfWrapper", as.double(n_total),
                       as.double(n_samples), as.double(m), as.double(d),
                       as.double(mele), basis_func, new.env(), as.double(x),
                       1.0)

  return(data.frame(x=x, p_est=p_bl_est))

}

cdfEst <- function(x, n_samples, p_est) {
# Estimate distribution functions F_k's, k = 0, 1, ..., m.

  # order of x
  x_order <- order(x)
  x_sort <- x[x_order]

  # extract p_est for k^{th} sample
  m <- length(n_samples) - 1
  n_total <- sum(n_samples)
  cum_p <- NULL
  for (i in 1:(m+1)) {
    pe_tmp <- p_est[p_est$k==(i-1), 3]
    cum_p <- c(cum_p, cumsum(pe_tmp[x_order]))
  }

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))
  x_data <- rep(x=x_sort, m+1)

  return(data.frame(k=population_indicator, x=x_data, cdf_est=cum_p))

}

quantEst <- function(k, p, cdf_est, n_samples,
                     interpolation=TRUE, adj=FALSE,
                     adj_val=NULL) {
# Estimate the quantile of the k[i]^th, k[i] = 0, 1, ..., m, population at
#   probability p[i] without giving covariance estimates.

  # Arguments handling and checking
  if ((length(k) != length(p)) && (length(k) != 1 ) && (length(p) != 1)) {
    stop("The lengths of vector 'k' and 'p' must be the same, or one of the vector must have length one!")
  }
  nK <- max(length(k), length(p))
  if (length(k)==1) k <- rep(k, nK)
  if (length(p)==1) p <- rep(p, nK)

  if (!is.logical(interpolation)) {
    stop("The argument 'interpolation' must be a logical variable (either TRUE or FALSE)!")
  }

  if (!is.logical(adj)) {
    stop("The argument 'adj' must be a logical variable (either TRUE or FALSE)!")
  }

  if (adj==TRUE) {

    if (is.null(adj_val)) {

      adj_val <- -1/(2*n_samples[(k+1)])

    } else {

      if (!is.numeric(adj_val)) {
        stop("The argument 'adj_val' must either be NULL or a numerical value vector!")
      } else if (length(adj_val) != 1 && length(adj_val) != nK) {
        stop("The length of the numerical argument 'adj_val' must either be 1 or the same as the length of the argument 'k' or, when the length of 'k' is 1, the same as the length of the argument 'p')!")
      }

      if (length(adj_val) == 1) adj_val <- rep(adj_val, nK)

    }

  }

  # Extract sorted data and useful information about DRM 
  x_sort <- cdf_est[cdf_est$k==0,]$x

  qe <- numeric(nK)
  for (i in 1:nK) {
    cdf_est_tmp <- cdf_est[cdf_est$k==k[i], 3]
    if (adj==TRUE) {
      cdf_est_tmp <- cdf_est_tmp + adj_val[i]
    }
    if (p[i] >= tail(cdf_est_tmp, 1)) {
      pos_tmp <- length(x_sort) - 1
    } else {
      pos_tmp <- length(which(cdf_est_tmp < p[i]))
    }

    if (interpolation==TRUE) {

      # linear interpolation
      qe[i] <- x_sort[pos_tmp] + (x_sort[pos_tmp + 1] - x_sort[pos_tmp]) *
        (p[i] - cdf_est_tmp[pos_tmp]) /
          (cdf_est_tmp[pos_tmp+1] - cdf_est_tmp[pos_tmp])

    } else {

      # without linear interpolation
      qe[i] <- x_sort[pos_tmp + 1]

    }

  }

  return(qe)

}

bwEst <- function(k, n_samples, p_est, cdf_est, interpolation=TRUE) {
# Estimate bandwith for density estimation for distribution functions F_k, k = 0, 1, ..., m.

  nK <- length(k)

  # Estimated means
  sd_est <- numeric(nK)
  x <- p_est[p_est$k==0,2]
  for (i in 1:nK) {
    pe_tmp <- p_est[p_est$k==k[i],3]
    mu_tmp <- sum(x * pe_tmp)
    sd_est[i] <- sqrt( sum(x^2 * pe_tmp) - mu_tmp^2 )
  }

  q_est <- quantEst(k=c(k, k), 
                    p=c(rep(0.25, nK), rep(0.75, nK)),
                    cdf_est, n_samples,
                    interpolation=interpolation,
                    adj=FALSE, adj_val=NULL)
  iqr <- q_est[(nK+1):(2*nK)] - q_est[1:nK]

  bw <- numeric(nK)
  const <- 1.06 / ( (sum(n_samples))^(1/5) )
  for (i in 1:nK) {
    bw[i] <- const * min(sd_est[i], iqr[i]/1.34) 
  }

  bwnames <- NULL 
  for (i in 1:nK) {
    bwnames <- c(bwnames, paste("F", k[i], sep=""))
  }
  names(bw) <- bwnames
  rm(bwnames)

  return(bw)
}

summaryDRMFEst <- function(n_samples, p_est, cdf_est, interpolation=TRUE) {
# Estimated mean, variance, IQR, of distribution functions F_k's, k = 0, 1, ..., m.

  m <- length(n_samples) - 1
  n_total <- sum(n_samples)

  # Estimated means
  mu_est <- numeric(m+1)
  var_est <- numeric(m+1)
  x <- p_est[p_est$k==0,2]
  for (i in 1:(m+1)) {
    pe_tmp <- p_est[p_est$k==(i-1),3]
    mu_est[i] <- sum(x * pe_tmp)
    var_est[i] <- sum(x^2 * pe_tmp) - mu_est[i]^2
  }
  sd_est <- sqrt(var_est)

  k <- 0:m
  q1 <- quantEst(k, 0.25, cdf_est, n_samples,
                 interpolation=interpolation, adj=FALSE,
                 adj_val=NULL)
  q3 <- quantEst(k, 0.75, cdf_est, n_samples,
                 interpolation=interpolation, adj=FALSE,
                 adj_val=NULL)
  iqr <- q3 - q1


  result <- data.frame(mean=mu_est, var=var_est, sd=sd_est,
                       Q1=q1, Q3=q3, IQR=iqr)

  rnames <- c("F0")
  for (i in 1:m) {
    rnames <- c(rnames, paste("F", i, sep=""))
  }
  row.names(result) = rnames 

  return(result)
}

displayPar <- function(par, m) {

  d <- length(par)/m - 1

  par_out <- matrix(par, d+1, m) 
  par_out <- t(par_out)
  par_out <- data.frame(par_out)
  rnames <- NULL
  for (i in 1:m) {
    #rnames <- c(rnames, paste("theta[", i, "]", sep=""))
    rnames <- c(rnames, paste("F", i, sep=""))
  }
  cnames <- paste("alpha[]")
  for (i in 1:d) {
    cnames <- c(cnames, paste("beta[,", i, "]", sep=""))
  }
  row.names(par_out) <- rnames
  colnames(par_out) <- cnames

  return(par_out)

}

drmdel <- function(x, n_samples, basis_func, g_null=NULL,
                   g_null_jac=NULL, par_dim_null=NULL,
                   ...)
{
  # setting useful constants
  n_total <- sum(n_samples)
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  rho <- n_samples/n_total
  m <- length(n_samples) - 1

  if (is.function(basis_func)) {
    d <- length(basis_func(x[1]))
  } else if (basis_func %in% 1:4) {
    d <- 1
  } else if (basis_func %in% 5:6) {
    d <- 2
  } else if (basis_func %in% 7:10) {
    d <- 3
  } else if (basis_func == 11){
    d <- 4
  } else {
    stop("Parameter 'basis_func' must be a function of a single variable or an integer between 1 and 11!")
  }

  drm_info <- list(m=m, d=d, basis_func=basis_func,
                   n_samples=n_samples, n_total=n_total, rho=rho)

  par_dim <- m*(d+1)

  par_pos <- gen_par_pos(m, d)

  # extracting ... arguments for function optim()
  dot_args <- list(...)
  # setting default method as "BFGS" for function optim()
  if (is.null(dot_args$method)) dot_args <- c(dot_args, list(method="BFGS"))
  # setting default maximum iteration number to 1000 for function optim()
  if (is.null(dot_args$control)) {
    #dot_args <- c(dot_args, list(control=list(maxit=10000,
                                              #reltol=.Machine$double.eps)))
    dot_args <- c(dot_args, list(control=list(maxit=10000)))
  } else {

    if (is.null(dot_args$control$maxit)) {
      dot_args$control <- c(dot_args$control, list(maxit=10000))
    }

    #if (is.null(dot_args$control$reltol)) {
      #dot_args$control <- c(dot_args$control,
                            #list(reltol=.Machine$double.eps))
    #}

  }
  # safeguard one-dimensional case
  if (par_dim==1) {
    dot_args$method="Brent"
    if (is.null(dot_args$lower)) dot_args <- c(dot_args, list(lower=-100))
    if (is.null(dot_args$upper)) dot_args <- c(dot_args, list(upper=100))
  }

  # calculating MELE
  # About initial values: it is better to always set initial values to zeros
  # for all parameters to ensure that, at the initial value, the dual
  # log empirical likelihood is finite!
  par_init <- rep(0, par_dim)  # 

  if (is.function(basis_func)) {

    drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDLUf,
                                     gr=negLDLGrUf),
                                dot_args,
                                list(x=x, n_total=n_total,
                                     n_samples=n_samples, m=m,
                                     basis_func=basis_func,
                                     d=d)))
    mele <- drm_opt$par
    mnames <- NULL
    for (i in 1:m) {
      mnames <- c(mnames, paste("alpha[", i, "]", sep=""))
      for (j in 1:d) {
        mnames <- c(mnames, paste("beta[", i, ",", j, "]", sep=""))
      }
    }
    names(mele) <- mnames

    negldl <- drm_opt$value

    # ##### Estimate information matrix #####
    logDLHessianMdele <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                               as.double(n_samples), as.double(m),
                               as.double(d), as.double(mele), basis_func,
                               new.env(), as.double(x))

    # Observed information matrix
    info_mat <- matrix(-logDLHessianMdele/n_total, par_dim, par_dim,
                       byrow=TRUE)

    rm(logDLHessianMdele)
    # ##########

    # estimate dF_{k}(x_kj)
    p_est <- pEstUf(x=x, n_total=n_total, n_samples=n_samples, m=m,
                    basis_func=basis_func, d=d, mele=mele)

    # estimate F_{k}(x_kj)
    cdf_est <- cdfEst(x=x, n_samples=n_samples, p_est=p_est)

    if (!is.null(g_null)) {

      if (is.null(par_dim_null)) {
        stop("One must provide 'par_dim_null', the dimension of the null parameter!")
      }

      par_init_null_full <- rep(0, (m+par_dim_null))
      par_dim_null_full <- par_dim_null + m

      dot_args_null <- dot_args
      # safeguard one-dimensional case
      if (par_dim_null_full==1) {
        dot_args_null$method="Brent"
        if (is.null(dot_args_null$lower)) {
          dot_args_null <- c(dot_args_null, list(lower=-100))
        }
        if (is.null(dot_args_null$upper)) {
          dot_args_null <- c(dot_args_null, list(upper=100))
        }
      }

      if (!is.null(g_null_jac) || par_dim_null==0) {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDLUf_null,
                                              gr=negLDLGrUf_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              g_null_jac=g_null_jac,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              par_dim_null=par_dim_null,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              basis_func=basis_func, d=d)))
      } else {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDLUf_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              basis_func=basis_func, d=d)))
      }

      negldl_null <- drm_opt_null$value
      mele_null <- list(alpha=drm_opt_null$par[1:m],
                        gamma=drm_opt_null$par[-(1:m)])
      delr <- -2*(negldl - negldl_null)  # DEL ratio statistics
      df <- par_dim - par_dim_null_full
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  mele_null=mele_null, negldl_null=negldl_null,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    } else {

      delr <- -2*negldl  # DEL ratio statistics
      df <- m*d
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    }

  } else {

    drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDL, gr=negLDLGr),
                                dot_args,
                                list(x=x, n_total=n_total,
                                     n_samples=n_samples, m=m,
                                     model=basis_func,
                                     d=d)))
    mele <- drm_opt$par
    mnames <- NULL
    for (i in 1:m) {
      mnames <- c(mnames, paste("alpha[", i, "]", sep=""))
      for (j in 1:d) {
        mnames <- c(mnames, paste("beta[", i, ",", j, "]", sep=""))
      }
    }
    names(mele) <- mnames

    negldl <- drm_opt$value

    # ##### Estimate information matrix #####
    logDLHessianMdele <- .C("logDualLHessianWrapper", as.double(n_total),
                as.double(n_samples), as.double(m), as.double(d),
                as.double(mele), as.double(basis_func), as.double(x),
                ldl_hessian=double(par_dim*par_dim))

    # Observed information matrix
    info_mat <- matrix(-(logDLHessianMdele$ldl_hessian)/n_total, par_dim,
                       par_dim, byrow=TRUE)

    rm(logDLHessianMdele)
    # ##########

    # estimate dF_{k}(x_kj)
    p_est <- pEst(x=x, n_total=n_total, n_samples=n_samples, m=m,
                  model=basis_func, d=d, mele=mele)

    # estimate F_{k}(x_kj)
    cdf_est <- cdfEst(x=x, n_samples=n_samples, p_est=p_est)

    if (!is.null(g_null)) {

      if (is.null(par_dim_null)) {
        stop("One must provide 'par_dim_null', the dimension of the null parameter!")
      }

      par_init_null_full <- rep(0, (m+par_dim_null))
      par_dim_null_full <- par_dim_null + m

      dot_args_null <- dot_args
      # safeguard one-dimensional case
      if (par_dim_null_full==1) {
        dot_args_null$method="Brent"
        if (is.null(dot_args_null$lower)) {
          dot_args_null <- c(dot_args_null, list(lower=-100))
        }
        if (is.null(dot_args_null$upper)) {
          dot_args_null <- c(dot_args_null, list(upper=100))
        }
      }

      if (!is.null(g_null_jac) || par_dim_null==0) {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDL_null,
                                              gr=negLDLGr_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              g_null_jac=g_null_jac,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              par_dim_null=par_dim_null,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              model=basis_func, d=d)))
      } else {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDL_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              model=basis_func, d=d)))
      }

      negldl_null <- drm_opt_null$value
      mele_null <- list(alpha=drm_opt_null$par[1:m],
                        gamma=drm_opt_null$par[-(1:m)])
      delr <- -2*(negldl - negldl_null)  # DEL ratio statistics
      df <- par_dim - par_dim_null_full
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  mele_null=mele_null, negldl_null=negldl_null,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    } else {

      delr <- -2*negldl  # DEL ratio statistics
      df <- m*d
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    }

  }

}

summaryDRM <- function(drmfit)
{
  space3 <- "   "
  cat("Basic information about the DRM:\n")
  cat(space3, "Number of samples (m+1):", drmfit$drm_info$m+1, "\n")
  if (is.function(drmfit$drm_info$basis_func)) {
    cat(space3, "Basis function:\n")
    print(drmfit$drm_info$basis_func)
  } else {
    cat(space3, "Basis function:", drmfit$drm_info$basis_func, "\n")
  }
  cat(space3, "Dimension of the basis function (d):", drmfit$drm_info$d, "\n")
  cat(space3, "Sample sizes:", drmfit$drm_inf$n_samples, "\n")
  cat(space3, "Sample proportions (rho):",
      format(drmfit$drm_inf$rho, digits=3), "\n")

  cat("\n")
  mele_display <- displayPar(drmfit$mele, drmfit$drm_info$m)
  cat("Maximum empirical likelihood estimator (MELE):\n")
  print(format(mele_display, digits=3))

  cat("\n")
  if (is.null(drmfit$mele_null)) {
    cat("Default hypothesis testing problem:\n")
    cat(space3, "H_0: All distribution functions, F_k's, are equal.\n")
    cat(space3, "H_1: One of the distribution functions is different from the others.\n")
  }
  cat("Dual empirical likelihood ratio statistc (DELR):",
      format(drmfit$delr, nsmall=3), "\n")
  cat("Degree of freedom:", drmfit$df, "\n")
  cat("p-value:", format(drmfit$p_val, digits=3), "\n")

  cat("\n")
  summaryStatF <- summaryDRMFEst(drmfit$drm_info$n_samples, drmfit$p_est,
                                 drmfit$cdf_est, interpolation=TRUE)
  cat("Summary statistics of the estimated F_k's (mean, var -- variance, sd -- standard deviation, Q1 -- first quartile, Q3 -- third quartile, IQR -- inter-quartile range):\n")
  print(format(summaryStatF, digits=3))

}

meleCov <- function(drmfit)
{

  Wmat <- .C("WmatWrapper", as.double(drmfit$drm_info$n_total),
             as.double(drmfit$drm_info$n_samples),
             as.double(drmfit$drm_info$m), as.double(drmfit$drm_info$d),
             W=double((drmfit$drm_info$m*(drmfit$drm_info$d+1))^2))

  W <- matrix(Wmat$W, drmfit$drm_info$m*(drmfit$drm_info$d+1),
              drmfit$drm_info$m*(drmfit$drm_info$d+1), byrow=TRUE)

  return(solve(drmfit$info_mat) - W)

}

densityDRM <- function(k, drmfit, interpolation=TRUE, ...)
{
  # Arguments handling and checking
  if (!is.numeric(k) || length(k) != 1) {
    stop("The argument 'k' must be an integer in the set of {0, 1, ..., m}!")
  }

  dot_args <- list(...)
  if (!is.null(dot_args$x)) {
    warning("'x' is not allowed to be specified manually. It is automatically set to the observed data points.")
    message("'x' now is by replaced by the observed data points ...")
  }
  dot_args$x <- drmfit$p_est[drmfit$p_est$k==k,2]

  if (!is.null(dot_args$weights)) {
    stop("'weights' is not allowed to be specified manually. It is automatically set to be the estimated dF_{k}(x_{kj})'s under the DRM.")
    message("'weights' now is replaced by the estimated dF_{k}(x_{kj})'s under the DRM...")
  }
  dot_args$weights <- drmfit$p_est[drmfit$p_est$k==k,3]

  if (is.null(dot_args$bw)) {
    dot_args$bw <- bwEst(k, drmfit$drm_info$n_samples,
                         drmfit$p_est, drmfit$cdf_est,
                         interpolation=interpolation)
  }

  return(do.call(density, dot_args))

}

cdfDRM <- function(k, x=NULL, drmfit, interpolation=TRUE)
{
  # Arguments handling and checking
  nK <- length(k)

  if (!is.logical(interpolation)) {
    stop("The argument 'interpolation' must be a logical variable (either TRUE or FALSE)!")
  }

  # CDF estimation
  cdf_est <- list(nK)
  if (is.null(x)) {

    for (i in 1:nK) {
      cdf_est[[i]] <- drmfit$cdf_est[drmfit$cdf_est$k==k[i],2:3]
      row.names(cdf_est[[i]]) <- NULL
    }

  } else if (is.list(x) || is.numeric(x)) {

    if (is.list(x)) {
      if (nK != length(x)) {
        stop("The lengths of vector 'k' and list 'p' must be the same.")
      }
    }

    if (is.numeric(x)) {
      x_tmp <- x
      x <- list(nK)
      for (i in 1:nK) {
        x[[i]] <- x_tmp
      }
      rm(x_tmp)
    }

    x_sort <- drmfit$cdf_est[drmfit$cdf_est$k==0,]$x

    if (interpolation==TRUE) { # linear interpolation

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]

        cdf_est[[i]] <- numeric(length(x[[i]]))
        for (j in 1:length(x[[i]])) {
          if (x[[i]][j] <= x_sort[1]) {
            cdf_est[[i]][j] <- cdf_est_tmp[1]
          } else if (x[[i]][j] >= tail(x_sort,1)) {
            cdf_est[[i]][j] <- tail(cdf_est_tmp,1)
          } else {
            pos_tmp <- length(which(x_sort <= x[[i]][j]))
            cdf_est[[i]][j] <- cdf_est_tmp[pos_tmp] +
              (cdf_est_tmp[pos_tmp + 1] - cdf_est_tmp[pos_tmp]) *
              (x[[i]][j] - x_sort[pos_tmp]) /
                (x_sort[pos_tmp + 1] - x_sort[pos_tmp])
          }
        }

        cdf_est[[i]] <- data.frame(x=x[[i]], cdf_est=cdf_est[[i]])

      }

    } else { # without linear interpolation

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]

        cdf_est[[i]] <- numeric(length(x[[i]]))
        for (j in 1:length(x[[i]])) {
          if (x[[i]][j] <= x_sort[1]) {
            cdf_est[[i]][j] <- cdf_est_tmp[1]
          } else if (x[[i]][j] >= tail(x_sort,1)) {
            cdf_est[[i]][j] <- tail(cdf_est_tmp,1)
          } else {
            pos_tmp <- length(which(x_sort <= x[[i]][j]))
            cdf_est[[i]][j] <- cdf_est_tmp[pos_tmp]
          }
        }

        cdf_est[[i]] <- data.frame(x=x[[i]], cdf_est=cdf_est[[i]])

      }

    }

  } else {

    stop("The argument 'x' can be a list of the same length as the argument 'k', each component of which is a numerical vector, or 'x' can be a single numerical vector, or NULL (default). Other types are not allowed!")

  }

  # assign names to cdf_est components
  cdfnames <- NULL 
  for (i in 1:nK) {
    cdfnames <- c(cdfnames, paste("F", k[i], sep=""))
  }
  names(cdf_est) <- cdfnames
  rm(cdfnames)

  return(cdf_est)

}

quantileDRM <- function(k, p, drmfit, cov=TRUE, interpolation=TRUE,
                        adj=FALSE, adj_val=NULL, bw=NULL, show_bw=FALSE)
{
  # Arguments handling and checking
  if ((length(k) != length(p)) && (length(k) != 1 ) && (length(p) != 1)) {
    stop("The lengths of vector 'k' and 'p' must be the same, or one of the vector must have length one!")
  }
  nK <- max(length(k), length(p))
  if (length(k)==1) k <- rep(k, nK)
  if (length(p)==1) p <- rep(p, nK)

  if (!is.logical(cov)) {
    stop("The argument 'cov' must be a logical variable (either TRUE or FALSE)!")
  }

  if (!is.logical(interpolation)) {
    stop("The argument 'interpolation' must be a logical variable (either TRUE or FALSE)!")
  }

  if (!is.logical(adj)) {
    stop("The argument 'adj' must be a logical variable (either TRUE or FALSE)!")
  }

  # Extract basic DRM information --- sample sizes
  n_samples <- drmfit$drm_info$n_samples
  if (adj==TRUE) {

    if (is.null(adj_val)) {

      adj_val <- -1/(2*n_samples[(k+1)])

    } else {

      if (!is.numeric(adj_val)) {
        stop("The argument 'adj_val' must either be NULL or a numerical value vector!")
      } else if (length(adj_val) != 1 && length(adj_val) != nK) {
        stop("The length of the numerical argument 'adj_val' must either be 1 or the same as the length of the argument 'k' or, when the length of 'k' is 1, the same as the length of the argument 'p')!")
      }

      if (length(adj_val) == 1) adj_val <- rep(adj_val, nK)

    }

  }

  if (!is.logical(show_bw)) {
    stop("The argument 'show_bw' must be a logical variable (either TRUE or FALSE)!")
  }

  # Extract basic DRM information
  m <- drmfit$drm_info$m
  d <- drmfit$drm_info$d
  mele <- drmfit$mele
  n_total <- drmfit$drm_info$n_total
  basis_func <- drmfit$drm_info$basis_func
  rho <- drmfit$drm_info$rho[k+1]

  # Extract sorted data from a drm-fit object
  x_sort <- drmfit$cdf_est[drmfit$cdf_est$k==0,]$x

  # Estimate quantiles and the corresponding covariance matrix
  if (cov==FALSE) {

    qe <- numeric(nK)
    if (interpolation==TRUE) {

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]
        if (adj==TRUE) {
          cdf_est_tmp <- cdf_est_tmp + adj_val[i]
        }

        if (p[i] >= tail(cdf_est_tmp, 1)) {
          pos_tmp <- length(x_sort) - 1
        } else {
          pos_tmp <- length(which(cdf_est_tmp < p[i]))
        }

        # linear interpolation
        qe[i] <- x_sort[pos_tmp] + (x_sort[pos_tmp + 1] - x_sort[pos_tmp]) *
          (p[i] - cdf_est_tmp[pos_tmp]) /
            (cdf_est_tmp[pos_tmp+1] - cdf_est_tmp[pos_tmp])

      }

    } else {

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]
        if (adj==TRUE) {
          cdf_est_tmp <- cdf_est_tmp + adj_val[i]
        }

        if (p[i] >= tail(cdf_est_tmp, 1)) {
          pos_tmp <- length(x_sort) - 1
        } else {
          pos_tmp <- length(which(cdf_est_tmp < p[i]))
        }

        # without linear interpolation
        qe[i] <- x_sort[pos_tmp + 1]

      }

    }

    return(list(est=qe))

  } else {

    # More arguments handling and checking
    if (is.null(bw)) {
      bw <- bwEst(k, drmfit$drm_info$n_samples,
                  drmfit$p_est, drmfit$cdf_est,
                  interpolation=interpolation)
    } else {
      if (is.numeric(bw) && length(bw) == 1) {
        bw <- rep(bw, nK)
      } else if (!is.numeric(bw) || length(bw) != nK) {
        stop("The argument 'bw' must be a numerical variable, and its length must be either one or the same as the length of the longer one of the arguments 'k' and 'p'!")
      }

      bwnames <- NULL 
      for (i in 1:nK) {
        bwnames <- c(bwnames, paste("F", k[i], sep=""))
      }
      names(bw) <- bwnames
      rm(bwnames)
    }

    # Estimate quantiles
    qe <- numeric(nK)  # quantile estimator
    length_vec <- numeric(nK)  # length of data that are <= the quantile estimates
    density_est <- numeric(nK)  ## density estimators at estimated quantiles

    if (interpolation==TRUE) {

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]
        if (adj==TRUE) {
          cdf_est_tmp <- cdf_est_tmp + adj_val[i]
        }

        if (p[i] >= tail(cdf_est_tmp, 1)) {
          pos_tmp <- length(x_sort) - 1
        } else {
          pos_tmp <- length(which(cdf_est_tmp < p[i]))
        }

        # linear interpolation
        qe[i] <- x_sort[pos_tmp] + (x_sort[pos_tmp + 1] - x_sort[pos_tmp]) *
          (p[i] - cdf_est_tmp[pos_tmp]) /
            (cdf_est_tmp[pos_tmp+1] - cdf_est_tmp[pos_tmp])

        if (p[i] >= cdf_est_tmp[pos_tmp+1]) {
          length_vec[i] <- pos_tmp+1
        } else {
          length_vec[i] <- pos_tmp
        }

        # density estimation
        density_est[i] <-
          density(x=drmfit$p_est[drmfit$p_est$k==k[i], 2],
                  bw=bw[i], adjust=1, kernel="gaussian",
                  weights=drmfit$p_est[drmfit$p_est$k==k[i], 3],
                  n=1, from=qe[i], to=qe[i])$y

      }

    } else {

      for (i in 1:nK) {

        cdf_est_tmp <- drmfit$cdf_est[drmfit$cdf_est$k==k[i], 3]
        if (adj==TRUE) {
          cdf_est_tmp <- cdf_est_tmp + adj_val[i]
        }

        if (p[i] >= tail(cdf_est_tmp, 1)) {
          pos_tmp <- length(x_sort) - 1
        } else {
          pos_tmp <- length(which(cdf_est_tmp < p[i]))
        }

        # without linear interpolation
        length_vec[i] <- pos_tmp + 1
        #qe[i] <- x_sort[pos_tmp + 1]
        qe[i] <- x_sort[length_vec[i]]

        # density estimation
        density_est[i] <-
          density(x=drmfit$p_est[drmfit$p_est$k==k[i], 2],
                  bw=bw[i], adjust=1, kernel="gaussian",
                  weights=drmfit$p_est[drmfit$p_est$k==k[i], 3],
                  n=1, from=qe[i], to=qe[i])$y

      }

    }

    # Estimate the covariance matrix of the quantile estimators
    qe_cov <- matrix(rep(0, nK*nK), nK, nK)
    BHat <- matrix(rep(0, nK*m*(d+1)), nK, m*(d+1))

    if (is.function(basis_func)) {

      for (i in 1:nK) {

        B_tmp <- .Call("BEstUfWrapper", as.double(k[i]),
                       as.double(length_vec[i]), as.double(n_samples),
                       as.double(m), as.double(d), as.double(mele),
                       basis_func, new.env(), as.double(x_sort))

        BHat[i,] <- B_tmp/n_total

        rm(B_tmp)

      }

      for (i in 1:nK) {

        aHat_tmp <- .Call("aEstUfWrapper", as.double(k[i]), as.double(k[i]),
                          as.double(length_vec[i]), as.double(n_samples),
                          as.double(m), as.double(d), as.double(mele),
                          basis_func, new.env(), as.double(x_sort))

        qe_cov[i,i] <- 1/rho[i]*(p[i] - p[i]^2) -
          1/(rho[i]^2) * 
            (aHat_tmp/n_total - 
             crossprod(BHat[i,], solve(drmfit$info_mat, BHat[i,])))
        qe_cov[i,i] <- qe_cov[i,i]/(density_est[i]^2)

        rm(aHat_tmp)

        j <- i+1 
        while (j <= nK) {

          aHat_tmp <- .Call("aEstUfWrapper", as.double(k[i]), as.double(k[j]),
                            as.double(min(length_vec[i], length_vec[j])),
                            as.double(n_samples), as.double(m), as.double(d),
                            as.double(mele), basis_func, new.env(),
                            as.double(x_sort))

          if (k[i]==k[j]) {
            qe_cov[i,j] <- 1/rho[i]*(p[i] - p[i]*p[j]) -
              1/(rho[i]*rho[j]) * 
                (aHat_tmp/n_total
                 - crossprod(BHat[i,], solve(drmfit$info_mat, BHat[j,])))
          } else {
            qe_cov[i,j] <- -1/(rho[i]*rho[j]) * 
              (aHat_tmp/n_total
               - crossprod(BHat[i,], solve(drmfit$info_mat, BHat[j,])))
          }

          qe_cov[i,j] <- qe_cov[i,j]/(density_est[i]*density_est[j])
          qe_cov[j,i] <- qe_cov[i,j]

          rm(aHat_tmp)

          j <- j+1

        }

      }

    } else {

      for (i in 1:nK) {

        B_tmp <- .C("BEstWrapper", as.double(k[i]), as.double(length_vec[i]),
                    as.double(n_samples), as.double(m), as.double(d),
                    as.double(mele), as.double(basis_func),
                    as.double(x_sort), B=double(m*(d+1)))

        BHat[i,] <- (B_tmp$B)/n_total

        rm(B_tmp)

      }

      for (i in 1:nK) {

        aHat_tmp <- .C("aEstWrapper", as.double(k[i]), as.double(k[i]),
                       as.double(length_vec[i]), as.double(n_samples),
                       as.double(m), as.double(d), as.double(mele),
                       as.double(basis_func), as.double(x_sort), a=double(1))

        qe_cov[i,i] <- 1/rho[i]*(p[i] - p[i]^2) -
          1/(rho[i]^2) * 
            (aHat_tmp$a/n_total - 
             crossprod(BHat[i,], solve(drmfit$info_mat, BHat[i,])))
        qe_cov[i,i] <- qe_cov[i,i]/(density_est[i]^2)

        rm(aHat_tmp)

        j <- i+1
        while (j <= nK) {

          aHat_tmp <- .C("aEstWrapper", as.double(k[i]), as.double(k[j]),
                         as.double(min(length_vec[i], length_vec[j])),
                         as.double(n_samples), as.double(m), as.double(d),
                         as.double(mele), as.double(basis_func),
                         as.double(x_sort), a=double(1))

          if (k[i]==k[j]) {
            qe_cov[i,j] <- 1/rho[i]*(p[i] - p[i]*p[j]) -
              1/(rho[i]*rho[j]) * 
                (aHat_tmp$a/n_total
                 - crossprod(BHat[i,], solve(drmfit$info_mat, BHat[j,])))
          } else {
            qe_cov[i,j] <- -1/(rho[i]*rho[j]) * 
              (aHat_tmp$a/n_total
               - crossprod(BHat[i,], solve(drmfit$info_mat, BHat[j,])))
          }

          qe_cov[i,j] <- qe_cov[i,j]/(density_est[i]*density_est[j])
          qe_cov[j,i] <- qe_cov[i,j]

          rm(aHat_tmp)

          j <- j+1

        }

      }

    }

    if (show_bw==TRUE) {

      return(list(est=qe, cov=qe_cov, bw=bw))

    } else {

      return(list(est=qe, cov=qe_cov))

    }

  }

}

quantileCompWald <- function(quantileDRMObject, n_total,
                             pairwise=TRUE, p_adj_method="none",
                             A=NULL, b=NULL) {

  # Arguments handling and checking
  if (is.list(quantileDRMObject)) {

    if(is.null(quantileDRMObject$est)) {
      stop("The argument quantileDRMObject is not a proper output from quantileDRM() function!")
    }

    if(is.null(quantileDRMObject$cov)) {
      stop("To perform Wald-type test, the covariance matrix of the quantile estimators must be estimated. Specify 'cov=TRUE' in quantileDRM() function!")

    }

  } else {

    stop("The argument quantileDRMObject is not a proper output from quantileDRM() function!")

  }

  nK <- length(quantileDRMObject$est)

  if (!is.logical(pairwise)) {
    stop("The argument 'pairwise' must be a logical variable (either TRUE or FALSE)!")
  }

  if ((pairwise==TRUE) && (nK <= 1)) {
    warning("Only one quantile estimate available, no pairwise comparison to be done! The argument 'pairwise' is set to be FALSE!")
    pairwise <- FALSE
  }

  if (pairwise==FALSE && is.null(A) && is.null(b)) {
    stop("Nonthing to be done! Specify 'pairwise=TRUE' for pairwise comparison or 'A' and 'b' for linear hypothesis about the quantiles.")
  }

  if ( (is.null(A) && !is.null(b)) || (!is.null(A) && is.null(b)) ) {
    stop("To test a linear hypothesis, both 'A' - the left-hand side matrix in the linear hypothesis, and 'b' - the right-hand side vector in the linear hypothesis, must be specified.")
  }

  if (!is.null(A)) {

    if (!is.matrix(A)) {

      stop("The argument 'A' must be a non-singular matrix whose number of columns equals length(quantileDRMObject$est)!")

    } else {

      # Checking the dimension of A
      if (ncol(A) != nK) {
        stop("The number of columns of the matrix 'A' must equal length(quantileDRMObject$est)!")
      }
      
      if (nrow(A) > nK) {
        stop("The number of rows of the matrix 'A' must be no bigger than length(quantileDRMObject$est)!")
      }

      # Checking the singularity of A
      if (is_exact_singular(A)) {
        stop("The argument 'A' is exactly singular. It must be a non-singular matrix!")
      } else {
        rcond_A <- rcond(A)
        if (rcond_A < 1.1*.Machine$double.eps) {
          stop(paste("The argument 'A' is computationally singular: reciprocal condition number = ", rcond_A, ". It must be a computationally non-singular matrix!", sep=""))
        }
      }

    }

    if (is.matrix(b)) {
      if (ncol(b) != 1) {
        stop("The argument 'b' must be a vector or a matrix of a single column whose length equals length(quantileDRMObject$est)!")
      } else {
        b <- as.vector(b)
      }
    }

    if (!is.vector(b)) {
      stop("The argument 'b' must be a vector or a matrix of a single column whose length equals length(quantileDRMObject$est)!")

    } else {

      if (length(b) != nrow(A)) {
        stop("The length of 'b' must equal the number of rows of the matrix 'A'!")
      }

    }

  }

  result <- list()

  if ((pairwise==TRUE) && (nK > 1)) {

    p_val_pair <- numeric(0)

    for (i in 1:(nK-1)) {

      for (j in (i+1):nK) {

        # Wald statistic
        wald_stat_tmp <- 
          n_total * (quantileDRMObject$est[i] - quantileDRMObject$est[j])^2 /
          ( rbind(c(-1,1))%*% quantileDRMObject$cov[c(i, j), c(i, j)] %*%
           cbind(c(-1,1)) ) 

        # p-value
        p_val_pair <- c( p_val_pair, (1 - pchisq(wald_stat_tmp, df=1)) )

        rm(wald_stat_tmp)

      }

    }

    # p-value adjustment for multiple comparisons
    if (p_adj_method != "none") {
      p_val_pair <- p.adjust(p_val_pair, method=p_adj_method)
    }
    
    # create output p-value matrix (lower triangular)
    p_val_pair_mat <- matrix(rep(NA, nK*nK), nK, nK)
    p_val_pair_mat[lower.tri(p_val_pair_mat)] <- p_val_pair

    p_val_names <- character(0)
    for (i in 1:nK) {
      p_val_names <- c(p_val_names, paste("q", i, sep=""))
    }
    rownames(p_val_pair_mat) <- p_val_names
    colnames(p_val_pair_mat) <- p_val_names

    result <- c(result, list(p_val_pair=p_val_pair_mat))

  }

  if (!is.null(A)) {

    # Define test-statistic
    vec_tmp <- A %*% quantileDRMObject$est - cbind(b)
    wald_stat_tmp <- as.numeric( n_total *
      crossprod( vec_tmp, solve( (A %*% quantileDRMObject$cov %*% t(A)),
                                vec_tmp ) ) )

    # p-value
    df_tmp <- nrow(A)
    p_val_linear_test <- 1 - pchisq(wald_stat_tmp, df=df_tmp)

    rm(vec_tmp, df_tmp, wald_stat_tmp)

    result <- c(result, list(p_val=p_val_linear_test))

  }

  return(result)

}

is_exact_singular <- function(x) {

  # Arguments handling and checking
  if (is.vector(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("The argument of the function is_exact_singluar() must be a numerical vector or matrix!")

  m <- nrow(x)
  n <- ncol(x)
  storage.mode(x) = "double"

  # LAPACK subroutine DGETRF (double precision) computes an LU factorization
  # of a general M-by-N matrix A using partial pivoting with row interchanges.
  # It also provides test for exact singularity.
  #
  # subroutine dgetrf(integer M, 
  #                   integer N,
  #                   double precision, dimension( lda, * ) A,
  #                   integer LDA,
  #                   integer, dimension( * ) IPIV,
  #                   integer INFO 
  #                  )
  #
  # Input: M, N, A, LDA
  # Output: A, IPIV, INFO

  lu_decomp <- .C("dgetrfCWrapper", as.double(m), as.double(n),
                  as.double(x), as.double(m), double(min(m, n)),
                  info=double(1))

  return(lu_decomp$info > 0)

}

trueParNormal <- function(m, mu, sigma) {
# m -- number of samples - 1
# mu -- mean vector of length m+1
# sigma -- standard deviation (sd) vector of length m+1

  par <- numeric(m*3)
  for (i in 1:m) {
    par[1 + 3*(i-1)] <- log(sigma[1]/sigma[i+1]) + (mu[1])^2/2/(sigma[1])^2 -
      (mu[i+1])^2/2/(sigma[i+1])^2  # alpha

    par[2 + 3*(i-1)] <- mu[i+1]/(sigma[i+1])^2 - mu[1]/(sigma[1])^2 # beta 1

    par[3 + 3*(i-1)] <- 1/2/(sigma[1])^2 - 1/2/(sigma[i+1])^2 # beta 2
  }

  return(par)

}

#trueParGamma <- function(m, shape, scale) {
trueParGamma <- function(m, shape, rate) {
# m -- number of samples - 1
# shape -- shape vector of length m+1
# scale -- scale vector of length m+1
# rate -- rate vector of length m+1
# gamma density:
  #f(x) = ( 1/(scale^shape * Gamma(shape)) ) * x^(shape) * exp(-x/scale)
  #or
  #f(x) = ( rate^shape / Gamma(shape) ) * x^(shape) * exp(-rate * x)

  par <- numeric(m*3)

  # scale formulation
  #for (i in 1:m) {
    #par[1 + 3*(i-1)] <- (shape[1] * log(scale[1]) + log(gamma(shape[1]))) -
    #(shape[i+1] * log(scale[i+1]) + log(gamma(shape[i+1])))  # alpha

    #par[2 + 3*(i-1)] <- 1/scale[1] - 1/scale[i+1] # beta 1

    #par[3 + 3*(i-1)] <- shape[i+1] - shape[1] # beta 2
  #}

  # rate formulation
  for (i in 1:m) {
    par[1 + 3*(i-1)] <- (shape[i+1] * log(rate[i+1]) - log(gamma(shape[i+1]))) -
      (shape[1] * log(rate[1]) - log(gamma(shape[1]))) # alpha

    par[2 + 3*(i-1)] <- rate[1] - rate[i+1] # beta 1

    par[3 + 3*(i-1)] <- shape[i+1] - shape[1] # beta 2
  }

  return(par)

}
