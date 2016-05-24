gvcmcatfitridge <-
function (
x,
y,
weights,
family,
control,
offset = rep(0, nobs),
...
)

{

# definitions zu design + names
    start = NULL
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y) else names(y)
    conv <- FALSE
    nobs <- nrow(x)
    nlik <- if (control$scaled.lik) nobs else 1
    nvars <- ncol(x)
    EMPTY <- nvars == 0

# check gewichte/offset
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    if (length(offset)==1)
        offset <- rep.int(offset, nobs)

# family defs
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)

    coefold <- if (nvars >= nobs || control$start.ml==TRUE) control$oml else rep(1, nvars) # rep(1, nvars) fur brd
    if (any(is.na(coefold))) coefold[which(is.na(coefold))] <- 0

    mustart <- etastart <- NULL
    if (is.null(start)) {  # mustart  # kein start value
        suppressWarnings(eval(family$initialize))
        etastart <- family$linkfun(mustart)
        start <- coefold # nur fur 1. iteration penalty!
    }     else {  # start value
        start <- if (length(start) != nvars)   # 
                    {stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                      nvars, paste(deparse(xnames), collapse = ", ")),
                      domain = NA)}
                 else { # 
                    isnogood <- is.na(start)
                    if (any(isnogood)) { start[isnogood] <- 0 }
                    start }
        eval(family$initialize)
        coefold <- start
        etastart <- offset + as.vector(if (NCOL(x) == 1) {x * start} else {x %*% start})
        mustart <- linkinv(etastart) # mukeep
    }

# was wenn keine Kovariablen?
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }

# was wenn Kovariablen?!
    if (!EMPTY) {

    # startwert eta => offset!!!
      eta <- etastart

    # start wert mu aus startwert eta, check mu
      mu <- mustart # linkinv(eta)
      if (!(validmu(mu) && valideta(eta)))
          stop("cannot find valid starting values: please specify some",
              call. = FALSE)

    # startwert dev
      devold <- sum(dev.resids(y, mu, weights))

    # startwert converged, boundary
      boundary <- conv <- FALSE

    # Matrix A
        acoefs <- diag(2, nrow=ncol(x), ncol=ncol(x))
        lambda <- .2
        A <- function(beta, acoefs=acoefs, lambda=lambda) {
          nbs <- t(acoefs)%*%beta
          fd <- as.integer(drop(nbs)!=0)
          appro <- fd * lambda * rep(2,length(nbs))
          crossprod(t(acoefs)*sqrt(appro))
          }

    # P-IRLS
      for (i in 1L:control$maxi) {

            good <- weights > 0  # index fuer beos, die beruecksichtigt werden

            # initials fuer fitting alg.
            varmu <- variance(mu)[good] # v.i <- variance(mu.i)/weights
            if (any(is.na(varmu)))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")

            mu.eta.val <- mu.eta(eta) # d.i <- mu.eta(eta.i)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")

            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ", i)
                break
            }

            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good]) # w.wurzel.i <- as.vector(d.i/sqrt(v.i))

            A.lambda <- nlik * A(start, acoefs=acoefs, lambda=lambda)
            x.star <- w*x #  in glm: as.vector(w*x)
            y.schlange <- as.vector(w*z)  # w.wurzel.i * (eta.i + (y - mu.i)/d.i)

            # fitting algorithm (mit beta.i = start)
                p.imat.new <- crossprod(x.star) + A.lambda
                chol.pimat.new <- chol(p.imat.new)
                inv.pimat.new <- chol2inv(chol.pimat.new)
                start.new <- control$g * drop(inv.pimat.new %*% t(x.star) %*%
                    y.schlange) + (1 - control$g) * start # start als matrix?!
            # ende fitting alg.

            # definitions output
            start <- start.new   # in glm: start[fit$pivot] <- fit$coefficients

            # check ergebnis fitting alg: ceofs (+ rank in glm)
            if (any(!is.finite(start))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", i), domain = NA)
                break
            }

            # berechne werte mu, eta, dev fuer alg in neuer iteration
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            boundary <- FALSE

            # check dev
            if (!is.finite(dev)) {
                # 1. stop, falls coefold null
#                if (is.null(coefold))
                if (any(coefold==NaN) || any(coefold==Inf) || any(is.na(coefold)))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                # 2. probiere stepsize zu verkleinern
                warning("step size truncated due to divergence",
                  call. = FALSE)
                ii <- 1
                    while (!is.finite(dev)) {
                      if (ii > control$maxi)
                        stop("inner loop 1; cannot correct step size",
                          call. = FALSE)
                      ii <- ii + 1
                      start <- (start + coefold)/2
                      eta <- drop(x %*% start)
                      mu <- linkinv(eta <- eta + offset)
                      dev <- sum(dev.resids(y, mu, weights))
                    }
                # 3. Meldung, dass Lsg am Rand
                boundary <- TRUE
            }

            # check eta, mu
            if (!(valideta(eta) && validmu(mu))) {
                # 1.
#                if (is.null(coefold))
                if (any(coefold==NaN) || any(coefold==Inf) || any(is.na(coefold)))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                # 2.
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxi)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                # 3.
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
            }

            # converged?!
            # if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            if (sum(abs(start - coefold))/sum(abs(coefold)) <= control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            } else {  # wenn nein weiter mit neuen coef werten
                devold <- dev
                coef <- coefold <- start
            }

      } # ende PIRLS
        # "output" PIRLS: start, mu, eta, dev, w, x.star, inv.pimat.new
        # recalculation of start, mu, eta, dev with rounded values of start
        
        start <- coef <- round(start, control$accuracy)
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))

      # check result fitting loop
      if (!conv)
          warning("Convergence warning for lambda = ", lambda[1], "\n", call. = FALSE)
      if (boundary)
          warning("Algorithm stopped at boundary value",
              call. = FALSE)
      #eps <- 10 * .Machine$double.eps
      if (family$family == "binomial") {
          if (any(mu > 1 - control$epsilon) || any(mu < control$epsilon))
              warning("Fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
      }
      if (family$family == "poisson") {
          if (any(mu < control$epsilon))
              warning("Fitted rates numerically 0 occurred",
                call. = FALSE)
      }

     
      names(coef) <- xnames # rownames

} # Ende Fall mit Kovariablen

# output
list(coefficients = coef)

}

