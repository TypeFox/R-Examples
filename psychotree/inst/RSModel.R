### S4 StatModel model
RSModel <- function (reltol = 1e-10, deriv = c("sum", "diff"),
                     hessian = TRUE, maxit = 100L) {
  new("StatModel",
      capabilities = new("StatModelCapabilities"),
      name = "RSM",
      dpp = ModelEnvFormula,
      fit = function (object, weights = NULL, ...){
        y <- object@get("response")
        z <- RSModel.fit(y = y, weights = weights, reltol = reltol,
                         deriv = deriv, hessian = hessian, maxit = maxit)
        z$ModelEnv <- object
        z$addargs <- list(...)
        z
    }
  )
}

## methods needed for mob()
reweight.RSModel <- function (object, weights, ...) {
  fit <- RSModel(reltol = object$reltol)@fit
  do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

estfun.RSModel <- function (x, ...) {

  ## get relevant informations
  dat <- x$data                    # completely cleaned (downcoded, null cats treatment, weights) data.
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  n <- nrow(dat)
  m <- ncol(dat)
  o <- max(dat, na.rm = TRUE)
  npar <- x$df
  ptot <- rowSums(dat, na.rm = TRUE) + 1 # +1 because gamma of score 0 is in row 1.
  ctot <- t(apply(dat, 1, tabulate, nbins = o))
  
  ## helper variables for first derivative calculation
  mv <- 1:m
  ov <- 1:o
  beta_index <- rep.int(mv,  rep.int(o, m))
  tau_index <- rep.int(ov, m)

  ## calculate gradient (matrix with 1:(m-1) cols for beta_j, j = 2, m, and o - 1 cols for tau^*_k, k = 2, ..o)
  if (!x$na) {
    
    ## select zero and first derivatives of gamma functions
    ## transform derivatives with item-category parameters to derivatives of RSM-parameters
    gamma0 <- x$esf[[1]][ptot]
    gamma1_pcm <- x$esf[[2]]
  
    ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
    gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m + o)
    for (j in mv) gamma1[, j] <- gamma1_pcm[, beta_index == j, drop = FALSE] %*% ov
    for (k in ov) gamma1[, k + m] <- rowSums(gamma1_pcm[, tau_index == k, drop = FALSE])
    gamma1 <- apply(gamma1, 2, "[", ptot)

    ## finally: the gradient
    agrad <- matrix(0, nrow = n, ncol = m + o)
    agrad[, mv] <- weights * (- dat + (gamma1 / gamma0)[, mv, drop = FALSE])
    agrad[, m + ov] <- weights * (- ctot + (gamma1 / gamma0)[, m + ov, drop = FALSE])

  } else {

    ## return value
    agrad <- matrix(0, nrow = n, ncol = m + o)

    ## observed NA patterns
    na_patterns <- factor(apply(is.na(dat), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop through na patterns, select derivatives and calculate gradient
    for(i in seq_len(nlevels(na_patterns))) {

      ## parse NA patterns
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      mv_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      mv_i <- if(length(mv_i) < 1) mv else mv[-mv_i]
      m_i <- length(mv_i)
      mv_i_1 <- 1:m_i

      ## calculate/fetch necessary stuff for gradient
      weights_i <- weights[na_i]
      ptot_i <- ptot[na_i]
      gamma0 <- x$esf[[i]][[1]][ptot_i]
      gamma1_pcm <- x$esf[[i]][[2]]
      beta_index_i <- rep.int(mv_i, rep.int(o, m_i))
      tau_index_i <- rep.int(ov, m_i)

      ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
      gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m_i + o)
      for (j in mv_i_1) gamma1[, j] <- gamma1_pcm[, beta_index_i == mv_i[j], drop = FALSE] %*% ov
      for (k in ov) gamma1[, k + m_i] <- rowSums(gamma1_pcm[, tau_index_i == k, drop = FALSE])
      gamma1 <- apply(gamma1, 2, "[", ptot_i)
      if (!is.matrix(gamma1)) gamma1 <- matrix(gamma1, nrow = 1)

      ## finally: the gradient for NA group i
      agrad[na_i, mv_i] <- weights_i * (- dat[na_i, mv_i, drop = FALSE] + (gamma1 / gamma0)[, mv_i_1, drop = FALSE])
      agrad[na_i, m_i + ov] <- weights_i * (- ctot[na_i, , drop = FALSE] + (gamma1 / gamma0)[, m_i + ov, drop = FALSE])
    }

  }

  ## collect and return matrix of initial size with gradients plugged in.
  grad <- matrix(0, ncol = npar, nrow = length(weights_org))
  grad[weights_org > 0, ] <- agrad[, -c(1, m + 1)]

  grad
}
