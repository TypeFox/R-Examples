`pcusp.old` <-
function (y, alpha, beta)
Vectorize(function(x) integrate(dcusp, -Inf, x, alpha = alpha,
    beta = beta)$value)(y)

`pcusp` <- # new pcusp the utilizes the c-code underlying cusp.nc because the older pcusp uses integrate from within R
  function (y, alpha, beta, subdivisions = 100, rel.tol = .Machine$double.eps^0.25,
            abs.tol = rel.tol, stop.on.error = TRUE, aux = NULL, keep.order = TRUE)
  {
    x = y
    limit <- as.integer(subdivisions)
    if (limit < 1 || (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps,
                                                    5e-29)))
      stop("invalid parameter values")

    # remember the order of alpha and beta
    idx = order(alpha, beta)

    # First compute the normalizing constants
    lower <- -Inf
    upper <- Inf
    inf <- 2 # integrate over (-infinity, infinity)
    bound <- 0
    wk <- .External("cuspnc", as.double(alpha), as.double(beta),
                    as.double(bound), as.integer(inf), as.double(abs.tol),
                    as.double(rel.tol), limit = limit)

    # Second compute the unnormalized probabilities
    lower <- -Inf
    upper <- Inf
    inf <- -1 # integrate over (-infinity, bound)
    pvals = value = abs.error = subdivs = ierr = matrix(NA, length(x), length(idx))
    for (i in 1:length(x)) {
      bound <- x[i]
      wk2 <- .External("cuspnc", as.double(alpha), as.double(beta),
                       as.double(bound), as.integer(inf), as.double(abs.tol),
                       as.double(rel.tol), limit = limit)
      pvals[i, ] = wk2$value / wk$value
      value[i, ] = wk2$value
      abs.error[i, ] = wk2$abs.error
      subdivs[i, ] = wk2$subdivisions
      ierr[i, ] = wk2$ierr
    }
    if (!keep.order) drop(pvals[,seq_along(idx)[idx]]) else drop(pvals) # somehow this is inverted... ???
  }
