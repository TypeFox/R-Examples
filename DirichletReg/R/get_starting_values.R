get_starting_values <- function(Y, X.mats, Z.mat, repar, base, weights){

  ops <- options(warn = -1L)
  on.exit(options(ops))

  if(!repar){###################################################### COMMON MODEL

    ### collinearity check begin
    exclude_par <- lapply(X.mats, function(list_el){
      lin_coef <- lm.fit(x = list_el, y = runif(nrow(list_el)))[["coefficients"]]
      if(any(na_pos <- is.na(lin_coef))){
        temp_names <- names(lin_coef)
        lin_coef <- seq_along(lin_coef)
        names(lin_coef) <- temp_names
        return(lin_coef[na_pos])
      } else {
        return(NULL)
      }
    })
    ### collinearity check end

    beta.LL <- function(x, y, X, w){
      b <- matrix(x, ncol = 2L)
      if(ncol(X) > 1L){
        LL <- w * dbeta(y, exp(X%*%b[,1]), exp(X%*%b[,2]), log=TRUE)
      } else {
        LL <- w * dbeta(y, unlist(exp(X*x[1])), unlist(exp(X*x[2])), log=TRUE)
      }
      return(LL)
    }

    beta.LL.deriv <- function(x, y, X, w){
      b <- matrix(x, ncol=2L)
      grad <- matrix(0.0, nrow=nrow(X), ncol=prod(dim(b)))
      element <- 1L
      if(ncol(X) > 1L){
        for(cc in seq_len(ncol(b)))for(rr in seq_len(nrow(b))){
          grad[,element] <- w * X[,rr]*(psigamma(exp(X%*%b[,1L]+X%*%b[,2L]))-psigamma(exp(X%*%b[,cc]))+log(y))
          element <- element + 1L
        }
      } else {
        for(cc in seq_len(ncol(b))){
          grad[,element] <- w * X*(psigamma(exp(X%*%x[1L]+X%*%x[2L]))-psigamma(exp(X%*%x[cc]))+log(y))
          element <- element + 1L
        }
      }
      return(grad)
    }

    unidim_fit <- lapply(seq_len(ncol(Y)), function(i){
      if(is.null(exclude_par[[i]])){
        correctX <- X.mats[[i]]
      } else {
        correctX <- X.mats[[i]][ , -exclude_par[[i]], drop = FALSE]
      }
      #suppressWarnings(
        maxBFGS(beta.LL, beta.LL.deriv,
          start        = rep(0, 2*ncol(correctX)),
          tol          = 1e-05,
          finalHessian = FALSE,
          X            = correctX,
          y            = Y[,i],
          w            = weights)$estimate[seq_len(ncol(correctX))]
      #)
    })

    for(cmp in seq_len(ncol(Y))){
      if(is.null(exclude_par[[cmp]])){
        break
      } else {
        for(NA_vars in seq_along(exclude_par[[cmp]])) unidim_fit[[cmp]] <- append(unidim_fit[[cmp]], NA, exclude_par[[cmp]][NA_vars] - 1L)
      }
    }

    unidim_fit <- unlist(unidim_fit)

  } else {#################################################### ALTERNATIVE MODEL

    Y_logr <- log(Y[,-base,drop=FALSE]/(Y[,base,drop=TRUE]))
    unidim_fit <- as.numeric(lm(Y_logr ~ X.mats[[1L]] - 1, weights = weights)[["coefficients"]])

    epsilon <- matrix(1.0, nrow = nrow(Y), ncol = ncol(Y))

    start_var <- 1L
    n_mean_par <- ncol(X.mats[[1L]])

    for(i in seq_len(ncol(Y))){
      if(i == base){
        next
      } else {
        epsilon[, i] <- as.numeric(exp(X.mats[[1L]] %*% unidim_fit[seq.int(start_var, start_var + n_mean_par - 1L)]))
        start_var <- start_var + n_mean_par
      }
    }

    MU <- epsilon / rowSums(epsilon)

    log_phi <- optimize(function(x){ sum(weights*ddirichlet(Y, MU * exp(x), log = TRUE)) }, c(-20, 20), maximum = TRUE)[["maximum"]]
    gammas <- as.numeric(lm(I(rep(log_phi, nrow(Y))) ~ Z.mat - 1, weights = weights)[["coefficients"]])

    unidim_fit <- c(unidim_fit, gammas)

  }

  return(as.numeric(unlist(unidim_fit)))

}
