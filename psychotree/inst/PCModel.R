### S4 StatModel model
PCModel <- function (nullcats = c("keep", "downcode", "ignore"), reltol = 1e-10,
                     deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L) {
  new("StatModel",
      capabilities = new("StatModelCapabilities"),
      name = "PCM",
      dpp = ModelEnvFormula,
      fit = function (object, weights = NULL, ...){
        y <- object@get("response")
        z <- PCModel.fit(y = y, weights = weights, nullcats = nullcats,
                         reltol = reltol, deriv = deriv, hessian = hessian, maxit = maxit)
        z$ModelEnv <- object
        z$addargs <- list(...)
        z
    }
  )
}

## methods needed for mob()
reweight.PCModel <- function (object, weights, ...) {
  fit <- PCModel(reltol = object$reltol)@fit
  do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

estfun.PCModel <- function (x, ...) {

  ## get relevant informations
  dat <- x$data                    # completely cleaned (downcoded, null cats treatment, weights) data.
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  n <- nrow(dat)
  m <- ncol(dat)
  oj_vec <- x$categories
  oj <- sapply(x$categories, length)
  npar_all <- sum(oj)
  npar_ident <- x$df
  ptot <- rowSums(dat, na.rm = TRUE) + 1 # +1 because gamma of score 0 is in row 1.

  ## helper variables
  parindex <- unlist(oj_vec)
  itemindex <- rep.int(1:m, oj)

  ## calculate gradient
  if (!x$na) {

    ## select gamma zero and first derivatives with ptot
    gamma0 <- x$esf[[1]][ptot]
    gamma1 <- apply(x$esf[[2]], 2, "[", ptot)

    ## construct data matrix ('selection' matrix, 0/1, cols = parameters)
    if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
      est_par <- !unlist(x$nullcats)
      gamma1 <- gamma1[, est_par, drop = FALSE]
    }
    xmat <- matrix(FALSE, nrow = n, ncol = npar_all)
    for (i in 1:n) xmat[i, ] <- dat[i, itemindex] == parindex 

    ## calculate gradient
    agrad <- weights * (- xmat + (gamma1 / gamma0))

  } else {

    ## return value & helper variables
    agrad <- matrix(0, nrow = n, ncol = npar_all)
    mv <- 1:m

    ## observed NA patterns 
    na_patterns <- factor(apply(is.na(dat), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop through na patterns, select derivatives and calculate gradient
    for(i in seq_len(nlevels(na_patterns))) {

      ## parse NA patterns and setup necessary stuff for gradient calculation
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      n_na_i <- length(na_i)
      mv_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      mv_i <- if(length(mv_i) < 1) mv else mv[-mv_i]
      oj_i <- oj[mv_i]
      oj_vec_i <- oj_vec[mv_i]
      weights_i <- weights[na_i]
      ptot_i <- ptot[na_i]
      dat_i <- dat[na_i, , drop = FALSE]

      ## select gamma zero and first derivatives with ptot_i
      gamma0_i <- x$esf[[i]][[1]][ptot_i]
      gamma1_i <- apply(x$esf[[i]][[2]], 2, "[", ptot_i)
      if (!is.matrix(gamma1_i)) gamma1_i <- matrix(gamma1_i, nrow = 1)

      ## construct data matrix ('selection' matrix, 0/1, cols = parameters) for NA group i
      parindex_i <- unlist(oj_vec_i)
      itemindex_i <- rep.int(mv_i, oj_i)

      if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
        est_par_i<- !unlist(x$nullcat[mv_i])
        gamma1_i <- gamma1_i[, est_par_i, drop = FALSE]
      }

      xmat_i <- matrix(FALSE, nrow = n_na_i, ncol = sum(oj_i))
      for (i in 1:n_na_i) xmat_i[i, ] <- dat_i[i, itemindex_i] == parindex_i

      ## finally: the gradient for NA group i
      agrad[na_i, itemindex %in% mv_i] <- weights_i * (- xmat_i + (gamma1_i/ gamma0_i))
    }

  }  

  ## collect and return matrix of initial size with gradients plugged in.
  grad <- matrix(0, ncol = npar_ident, nrow = length(weights_org))
  grad[weights_org > 0, ] <- agrad[, -1, drop = FALSE]

  grad
}
