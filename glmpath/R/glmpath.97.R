step.length <- function(corrector, lambda2, min.lambda, max.arclength,
                        frac.arclength, add.newvars, backshoot, h0 = NULL,
                        eps = .Machine$double.eps)
  {
    active <- corrector$active
    force.active <- corrector$force.active
    lambda <- corrector$lambda - min.lambda
    Lambda2 <- c(0, rep(lambda2, length(active) - 1))
    C <- corrector$corr
    b <- corrector$b[active]
    xw <- corrector$xw
    xwx <- t(xw) %*% xw[, active]
    xwx.in <- solve(xwx[active, ] + diag(Lambda2))
    db <- drop(xwx.in %*% c(rep(0, length(force.active)),
                            sign(C[active[-force.active]])))
    newa <- NULL
    if (!backshoot) {
      Cmax <- max(abs(C))
      inactive <- seq(C)[-active]
      ninact <- length(inactive)
      a <- drop(xwx[inactive, ] %*% db)
      gam <- c((Cmax - C[inactive]) / (1 - a), (Cmax + C[inactive]) / (1 + a))
      ha <- min(gam[gam > eps], lambda)
      hd <- -b[-force.active] / db[-force.active]
      h <- min(hd[hd > eps], ha)
      if (h == ha && h < lambda) {
        ii <- which(gam > eps)
        nii <- length(ii)
        oo <- order(gam[ii])
        oo <- oo[1:min(nii, add.newvars)]
        oo <- ii[oo]
        oo <- unique(ifelse(oo <= ninact, oo, oo - ninact))
        newa <- inactive[oo]
      }
      h <- min(h * frac.arclength, max.arclength / sum(abs(db)))
    } else {
      hd <- b[-force.active] / db[-force.active]
      ii <- hd > eps & hd < -h0
      if (any(ii)) h <- -max(hd[ii])
      else h = 0
    }
    list(h = -h, db = db, newa = newa, arclength = h * sum(abs(db)))
  }

predictor1 <- function(b, step)
  {
    b - step$db * step$h
  }

corrector1 <- function(x, y, family, weight, offset, active, tmpa,
                       force.active, lambda, lambda2, b0, bshoot.threshold,
                       relax.lambda, trace, no.iter = FALSE,
                       eps = .Machine$double.eps)
  {
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    dev.resids <- family$dev.resids
    p <- length(tmpa)
    if (!no.iter) {
      b2 <- c(pmax(b0[tmpa], 0), -pmin(b0[tmpa], 0))
      xa <- x[, tmpa, drop = FALSE]
      penalty <- rep(1, p)
      penalty[force.active] <- 0
      distr <- switch(family$family, gaussian = 0, binomial = 1, poisson = 2)
      z <- c(as.vector(xa), y, weight, offset, penalty)
      mz <- c(length(y), distr, lambda, lambda2)
      sol <- .C('solve_glmpath',
                as.integer(2 * p),
                as.double(b2),
                as.double(rep(0, 2 * p)),
                as.double(rep(1e300, 2 * p)),
                as.integer(0),
                as.double(z),
                as.double(mz))
      if (sol[[5]] != 0) {
        cat('Convergence warning\n')
      }
      b0[tmpa] <- sol[[2]][1:p] - sol[[2]][(p + 1):(2 * p)]
    }
    Lambda2 <- c(0, rep(lambda2, length(b0) - 1))
    eta <- drop(x %*% b0) + offset
    mu <- linkinv(eta)
    mu.eta.eta <- mu.eta(eta)
    w <- (weight * mu.eta.eta^2 / variance(mu))^0.5
    z <- (y - mu) / mu.eta.eta
    xw <- x * w
    wz <- w * z
    corr <- drop(t(wz) %*% xw) - Lambda2 * b0
    if (p > length(force.active)) {
      i <- which(abs(corr) >= lambda * (1 - relax.lambda))
      newa <- i[!i%in%tmpa]
      newactive <- i[!i%in%active]
      i <- which(abs(b0[active[-force.active]]) < eps)
      inactive <- active[-force.active][i]
      active <- active[!active%in%inactive]
      active <- c(active, newactive)
      b0[-active] <- 0
    } else {
      lambda <- max(abs(corr[-active]))
      c1 <- which(abs(corr) == lambda)
      newa <- newactive <- c1
      inactive <- NULL
      active <- c(active, c1)
    }
    df <- length(active) - length(newactive)
    dev <- sum(dev.resids(y, mu, weight))
    backshoot <- ifelse(any(abs(b0[newactive]) > bshoot.threshold),
                        TRUE, FALSE)
    list(b = b0, active = active, force.active = force.active,
         newactive = newactive, newa = newa, inactive = inactive, corr = corr,
         lambda = lambda, xw = xw, df = df, dev = dev, backshoot = backshoot)
  }

glmpath <- function(x, y, data, nopenalty.subset = NULL, family = binomial,
                    weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5,
                    max.steps = 10 * min(n, m), max.norm = 100 * m,
                    min.lambda = (if (m >= n) 1e-6 else 0), max.vars = Inf,
                    max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
                    bshoot.threshold = 0.1, relax.lambda = 1e-8,
                    standardize = TRUE, eps = .Machine$double.eps,
                    trace = FALSE)
  {
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = 'function', envir = parent.frame())
    if (is.function(family))
        family <- family()
    else if (family$family == 'gaussian') family <- gaussian()
    else if (family$family == 'binomial') family <- binomial()
    else if (family$family == 'poisson') family <- poisson()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (!missing(data)) {
      x <- data$x
      y <- data$y
    }
    n <- nrow(x)
    m <- ncol(x)
    no.iter <- FALSE
    if (family$family == 'gaussian' && family$link == 'identity')
      no.iter <- TRUE
    if (family$family == 'binomial') {
      if (any(y == -1)) {
        y[y == -1] <- 0
        cat('y=-1 converted to y=0.\n')
      }
    }
    if (length(weight) != n) {
      stop('Length of the weight vector != sample size')
    }
    if (any(weight < 0)) {
      stop('Negative weights are not allowed.')
    }
    if (length(offset) != n) {
      stop('Length of the offset vector != sample size')
    }
    if (frac.arclength > 1 || frac.arclength <= 0) {
      frac.arclength <- 1
      cat('frac.arclength should be in (0,1]. frac.arclength is set to 1.\n')
    }
    if (max.arclength < Inf && frac.arclength < 1) {
      frac.arclength <- 1
      cat(paste('frac.arclength<1 can be used only if max.arclength=Inf.',
                'frac.arclength is set to 1.\n'))
    }
    n.repeat <- n.repeat1 <- ceiling(1 / frac.arclength)
    one <- rep(1, n)
    meanx <- drop(one %*% x) / n
    x <- scale(x, meanx, FALSE)
    if (standardize) {
      sdx <- sqrt(drop(one %*% (x^2)) / (n - 1))
      ignore <- sdx < eps
      if (any(ignore)) sdx[ignore] <- 1
    } else {
      sdx <- rep(1, m)
    }
    x <- scale(x, FALSE, sdx)
    x <- cbind(1, x)
    if (is.null(dimnames(x)[[2]])) {
      xnames <- c('Intercept', paste('x', seq(m), sep=''))
    } else {
      xnames <- c('Intercept', dimnames(x)[[2]][-1])
    }
    if (length(nopenalty.subset) >= m) {
      stop('cannot force in all the variables')
    } else {
      nopenalty.subset <- c(1, nopenalty.subset + 1)
    }
    force.active <- c(1:length(nopenalty.subset))
    lam.vec <- step.len <- rep(0, max.steps)
    bmat.pred <- bmat.corr <- cmat <- matrix(0, max.steps, m + 1)
    new.df <- df <- dev <- rep(0, max.steps)
    new.A <- rep(FALSE, max.steps)
    actions <- vector('list', length = max.steps)
    backshoot <- FALSE
    b <- rep(0, (m + 1))
    corrector <- corrector1(x, y, family, weight, offset, nopenalty.subset,
                            nopenalty.subset, force.active, 0, lambda2, b,
                            bshoot.threshold, relax.lambda, trace)
    if (trace) cat('The model begins with', xnames[nopenalty.subset])
    k <- 1
    b <- bmat.pred[k, ] <- bmat.corr[k, ] <- corrector$b
    cmat[k, ] <- corrector$corr
    lam.vec[k] <- lambda <- corrector$lambda
    new.df[k] <- df[k] <- corrector$df
    dev[k] <- corrector$dev
    new.A[k] <- TRUE
    active <- corrector$active
    actions[[k]] <- active[-1] - 1
    names(actions[[k]]) <- xnames[active[-1]]
    if (trace) {
      cat(paste('\nLambda=', lambda, 'lets the first factor in.\nStep', k, ':'))
      cat(paste('\t', xnames[active[-force.active]], 'added'))
    }
    if (max.steps <= 1) {
      stop('Increase max.steps.')
    }
    if (max.norm <= sum(abs(b))) {
      stop('Increase max.norm.')
    }
    if (lambda <= min.lambda) {
      stop('Decrease min.lambda')
    }
    if (max.vars <= 1) {
      stop('Increase max.vars.')
    }
    while(TRUE) {
      if (!backshoot) {
        arclength <- max.arclength
        if (new.A[k]) {
          frac.arclength1 <- frac.arclength
        } else {
          frac.arclength1 <- 1
          if (n.repeat1 > 1 && max.arclength == Inf) {
            arclength <- step$arclength
          }
        }
        k <- k + 1
        if (trace) cat(paste('\nStep', k, ':'))
        step <- step.length(corrector, lambda2, min.lambda, arclength,
                            frac.arclength1, add.newvars, backshoot)
        b[active] <- predictor1(b[active], step)
        bmat.pred[k, ] <- b
        step.len[k-1] <- h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        tmpa <- c(active, step$newa)
      } else {
        if (trace) cat(paste('\nStep', k, ':'))
        step <- step.length(corrector, lambda2, min.lambda, Inf, 1,
                            add.newvars, backshoot, h)
        step.len[k-1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
      }
      corrector <- corrector1(x, y, family, weight, offset, active, tmpa,
                              force.active, lambda, lambda2, b,
                              bshoot.threshold, relax.lambda, trace, no.iter)
      newa <- corrector$newa
      while(length(newa) > 0) {
        if (trace) cat(paste('\nRepeating step', k, ':'))
        tmpa <- c(tmpa, newa)
        corrector <- corrector1(x, y, family, weight, offset, active, tmpa,
                                force.active, lambda, lambda2, b,
                                bshoot.threshold, relax.lambda, trace, no.iter)
        newa <- corrector$newa
      }
      newaction <- c(corrector$newactive, -corrector$inactive)
      if (length(newaction) > 0 && length(corrector$active) <= n) {
        if (corrector$backshoot && !backshoot) {
          if (trace) cat('\nOvershooting occurred: increasing lambda again')
          backshoot <- TRUE
          n.repeat1 <- 1
        } else {
          active <- corrector$active
          b <- corrector$b
          actions[[k]] <- sign(newaction) * (abs(newaction) - 1)
          names(actions[[k]]) <- xnames[abs(newaction)]
          new.df[k] <- corrector$df
          new.A[k] <- TRUE
          if (trace) {
            na <- newaction[newaction > 0]
            ina <- -newaction[newaction < 0]
            if (length(na) > 0) cat(paste('\t', xnames[na], 'added'))
            if (length(ina) > 0) cat(paste('\t', xnames[ina], 'dropped'))
          }
          backshoot <- FALSE
          n.repeat1 <- n.repeat
        }
      } else if (length(corrector$active) <= n) {
        active <- corrector$active
        b <- corrector$b
        backshoot <- FALSE
        n.repeat1 <- max(n.repeat1 - 1, 1)
      }
      if (!backshoot) {
        bmat.corr[k, ] <- b
        cmat[k, ] <- corrector$corr
        df[k] <- corrector$df
        dev[k] <- corrector$dev
        if (lambda <= min.lambda || k == max.steps ||
            length(corrector$active) > min(n, max.vars) ||
            sum(b[-nopenalty.subset]) >= max.norm) {
          if (length(corrector$active) > min(n, max.vars)) {
            k <- k - 1
          }
          if (trace) {
            if (lambda <= min.lambda) {
              cat('\nLambda=', min.lambda, '\n')
            } else if (k == max.steps) {
              cat(paste('\nMaximum steps (', max.steps, ') taken.\n'))
            } else if (length(corrector$active) > min(n, max.vars)) {
              cat('\nNumber of active variables has reached its maximum.\n')
            } else {
              cat('\n|beta| >=', max.norm, '\n')
            }
          }
          break
        }
      }
    }
    bmat.pred <- bmat.pred[1:k, ]
    bmat.corr <- bmat.corr[1:k, ]
    cmat <- cmat[1:k, ]
    bmat.pred[, -1] <- scale(bmat.pred[, -1], FALSE, sdx)
    bmat.corr[, -1] <- scale(bmat.corr[,-1], FALSE, sdx)
    bmat.pred[, 1] <- (bmat.pred[, 1] -
                       bmat.pred[, -1, drop = FALSE] %*% meanx)
    bmat.corr[, 1] <- (bmat.corr[, 1] -
                       bmat.corr[, -1, drop = FALSE] %*% meanx)
    dimnames(bmat.pred) <- dimnames(bmat.corr) <- dimnames(cmat) <-
      list(seq(k), xnames)
    df <- df[1:k]
    dev <- dev[1:k]
    if (family$family == 'gaussian')
      Dev <- sum(weight) * (log(2 * pi * dev / sum(weight)) + 1) + 2
    else Dev <- dev
    aic <- Dev + 2 * df
    bic <- Dev + log(n) * df
    if (length(nopenalty.subset) > 1) {
      names(nopenalty.subset) <- xnames[nopenalty.subset]
      nopenalty.subset <- nopenalty.subset[-1] - 1
    } else {
      nopenalty.subset <- NULL
    }
    object <- list(call = call, lambda = lam.vec[1:k], lambda2 = lambda2,
                   step.length = abs(step.len[1:(k - 1)]), corr = cmat,
                   new.df = new.df[1:k], df = df, deviance = dev, aic = aic,
                   bic = bic, b.predictor = bmat.pred, b.corrector = bmat.corr,
                   new.A = new.A[1:k], actions = actions[1:k], meanx = meanx,
                   sdx = sdx, xnames = xnames, family = family,
                   weight = weight, offset = offset,
                   nopenalty.subset = nopenalty.subset,
                   standardize = standardize)
    class(object) <- 'glmpath'
    object
  }

plot.glmpath <- function(x, xvar = c('norm', 'lambda', 'step'),
                         type = c('coefficients', 'aic', 'bic'),
                         plot.all.steps = FALSE, xlimit = NULL,
                         predictor = FALSE, omit.zero = TRUE, breaks = TRUE,
                         mar = NULL, eps = .Machine$double.eps, main = NULL,
                         ...)
  {
    object <- x
    ii <- object$new.A
    if (plot.all.steps) {
      ii[!ii] <- TRUE
    } else {
      ii[length(ii)] <- TRUE
    }
    lam <- object$lambda[ii]
    xvar <- match.arg(xvar)
    type <- match.arg(type)
    coef.pred <- scale(object$b.predictor[ii, -1], FALSE, 1 / object$sdx)
    coef.corr <- scale(object$b.corrector[ii, -1], FALSE, 1 / object$sdx)
    xnames <- object$xnames[-1]
    if (omit.zero) {
      c1 <- drop(rep(1, nrow(coef.corr)) %*% abs(coef.corr))
      nonzero <- c1 > eps
      xnames <- xnames[nonzero]
      coef.pred <- coef.pred[, nonzero, drop = FALSE]
      coef.corr <- coef.corr[, nonzero, drop = FALSE]
    }
    m <- ncol(coef.pred)
    k <- nrow(coef.pred)
    s <- switch(xvar,
      norm = if (is.null(object$nopenalty.subset)) apply(abs(coef.corr), 1, sum)
       else (apply(abs(coef.corr[, -object$nopenalty.subset, drop = FALSE]),
                   1, sum)),
      lambda = lam,
      step = seq(k))
    if (xvar != 'lambda') {
      if (is.null(xlimit)) xlimit <- max(s)
      else if (xlimit <= min(s)) stop('Increase xlimit.')
      xi <- s <= xlimit
    } else {
      if (is.null(xlimit)) xlimit <- min(s)
      else if (xlimit >= max(s)) stop('Decrease xlimit.')
      xi <- s >= xlimit
    }
    k <- max(which(xi))
    xname <- switch(xvar, norm = '|beta|', lambda = 'lambda', step = 'step')
    if (!is.null(mar)) par(mar=mar)
    if (type == 'aic') {
      aic <- object$aic[ii][xi]
      plot(s[xi], aic, xlab = xname, ylab = 'AIC', type = 'b', pch = 16,
           cex = 0.3, ...)
      if (is.null(main)) title('AIC', line = 2.5)
      else title(main, line = 2.5)
    } else if (type == 'bic') {
      bic <- object$bic[ii][xi]
      plot(s[xi], bic, xlab = xname, ylab = 'BIC', type = 'b', pch = 16,
           cex = 0.3, ...)
      if (is.null(main)) title('BIC', line = 2.5)
      else title(main, line = 2.5)
    } else {
      ylab <- ifelse(object$standardize, 'Standardized coefficients',
                     'Coefficients')
      matplot(s[xi], coef.corr[xi, ], xlab = xname, type = 'b', pch = '*',
              ylab = ylab, lty = 1, ...)
      if (is.null(main)) title('Coefficient path', line = 2.5)
      else title(main, line = 2.5)
      abline(h = 0, lty = 3)
      axis(4, at = coef.corr[k, ], labels = xnames, cex = 0.8, adj = 0, las = 1)
      if (predictor) {
        for (i in 1:m) {
          segments(s[xi][-k], coef.corr[xi, ][-k, i], s[xi][-1],
                   coef.pred[xi, ][-1, i], lty = 2, col = i)
        }
      }
    }
    if (breaks) {
      new <- object$new.A[ii] & xi
      axis(3, at = s[new], labels = object$new.df[ii][new], cex = 0.8)
      abline(v = s[new])
    }
  }

predict.glmpath <- function(object, newx, newy, s,
                            type = c('link', 'response', 'loglik',
                              'coefficients'),
                            mode = c('step', 'norm.fraction', 'norm',
                              'lambda.fraction', 'lambda'), weight = NULL,
                            offset = NULL, eps = .Machine$double.eps, ...)
  {
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(newx) && type != 'coefficients') {
      warning('no newx argument; type switched to coefficients')
      type <- 'coefficients'
    }
    if (missing(newy) && type == 'loglik') {
      warning('no newy argument; type switched to coefficients')
      type <- 'coefficients'
    }
    if (type != 'coefficients') {
      if (is.vector(newx)) newx <- matrix(newx, nrow = 1)
    }
    if (any(object$offset != 0) && type != 'coefficients') {
      if (is.null(offset)) {
        stop('offset missing')
      }
      if (length(offset) != nrow(newx)) {
        stop('Length of the offset vector != sample size in newx')
      }
    }
    b <- object$b.corrector
    std.b <- scale(b[, -1], FALSE, 1 / object$sdx)
    if (!is.null(object$nopenalty.subset))
      std.b <- std.b[, -object$nopenalty.subset, drop = FALSE]
    k <- nrow(b)
    steps <- seq(k)
    if (missing(s)) {
      s <- steps[object$new.A]
      if (mode != 'step') {
        warning('no s argument; mode switched to step')
        mode <- 'step'
      }
    }
    sb <- switch(mode, step = {
      if (any(s < 1) || any(s > k))
        stop('Argument s out of range')
      steps
    }, norm.fraction = {
      if (any(s > 1) || any(s < 0))
        stop('Argument s out of range')
      bnorm <- apply(abs(std.b), 1, sum)
      bnorm / bnorm[k]
    }, norm = {
      bnorm <- apply(abs(std.b), 1, sum)
      if (any(s > bnorm[k]) || any(s < bnorm[1]))
        stop('Argument s out of range')
      bnorm
    }, lambda.fraction = {
      if (any(s > 1) || any(s < 0))
        step('Argument s out of range')
      lam <- object$lambda
      lam[lam < eps] <- eps
      lam <- log(lam)
      (lam - min(lam)) / (max(lam) - min(lam))
    }, lambda = {
      lam <- object$lambda
      if (any(s > lam[1]) || any(s < lam[k]))
        stop('Argument s out of range')
      lam
    })
    sfrac <- (s - sb[1]) / (sb[k] - sb[1])
    sb <- (sb - sb[1]) / (sb[k] - sb[1])
    usb <- unique(sb)
    useq <- match(usb, sb)
    sb <- sb[useq]
    b <- b[useq, ]
    coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- (((sb[right] - sfrac) * b[left, , drop = FALSE] +
              (sfrac - sb[left]) * b[right, , drop = FALSE]) /
             (sb[right] - sb[left]))
    newb[left == right, ] <- b[left[left == right], ]
    if (type != 'coefficients') {
      if (missing(offset)) offset <- rep(0, nrow(newx))
      fit <- cbind(1, newx) %*% t(newb) + offset
      if (type != 'link') {
        fit <- object$family$linkinv(fit)
        if (type == 'loglik') {
          if (is.null(weight)) weight <- rep(1, length(newy))
          loglik <- matrix(0, length(newy), length(s))
          for (k in 1:length(s)) {
            loglik[, k] <- (-object$family$dev.resids(newy, fit[, k], weight) /
                            2)
          }
          fit <- loglik
        }
      }
      dimnames(fit) <- list(seq(nrow(newx)), s)
    } else {
      fit <- newb
      dimnames(fit) <- list(s, object$xnames)
    }
    attr(fit, 's') <- s
    attr(fit, 'fraction') <- sfrac
    attr(fit, 'mode') <- mode
    fit
}

cv.glmpath <- function(x, y, data, family = binomial, weight = rep(1, n),
                       offset = rep(0, n), nfold = 10,
                       fraction = seq(0, 1, length = 100),
                       type = c('loglik', 'response'),
                       mode = c('norm', 'lambda'), plot.it = TRUE, se = TRUE,
                       ...)
  {
    type <- match.arg(type)
    mode <- match.arg(mode)
    if (is.character(family))
        family <- get(family, mode = 'function', envir = parent.frame())
    if (is.function(family))
        family <- family()
    else if (family$family == 'gaussian') family <- gaussian()
    else if (family$family == 'binomial') family <- binomial()
    else if (family$family == 'poisson') family <- poisson()
    if (!missing(data)) {
      x <- data$x
      y <- data$y
    }
    if (family$family == 'binomial') {
      uy <- unique(y)
      if (all(uy == c(1, -1)) | all(uy == c(-1, 1))) {
        y[y == -1] <- 0
        cat('y=-1 converted to y=0 \n')
      }
    }
    n <- length(y)
    folds <- split(sample(seq(n)), rep(1:nfold, length = n))
    errormat <- matrix(0, length(fraction), nfold)
    for (i in seq(nfold)) {
        omit <- folds[[i]]
        fit <- glmpath(x[-omit, ], y[-omit], family = family,
                       weight = weight[-omit], offset = offset[-omit], ...)
        pred <- switch(mode, norm = {
          predict(fit, x[omit, , drop = FALSE], y[omit], offset[omit],
                  s = fraction, type = type, mode = 'norm.fraction',
                  weight = weight[omit])
        }, lambda = {
          predict(fit, x[omit, , drop = FALSE], y[omit], offset[omit],
                  s = fraction, type = type, mode = 'lambda.fraction',
                  weight = weight[omit])
        })
        if (length(omit) == 1) pred <- matrix(pred, nrow = 1)
        if (type == 'loglik') errormat[, i] <- -apply(pred, 2, mean)
        else {
          if (family$family == 'binomial') pred <- ifelse(pred > 0.5, 1, 0)
          errormat[, i] <- apply((y[omit] - pred)^2, 2, mean)
        }
        cat('CV Fold', i, '\n')
    }
    cv.error <- apply(errormat, 1, mean)
    cv.se <- sqrt(apply(errormat, 1, var) / nfold)
    object <- list(fraction = fraction, cv.error = cv.error, cv.se = cv.se,
                   folds = folds)
    if (plot.it) {
      lab <- switch(type, loglik='Cross-validated minus log-likelihood',
                    response = 'Cross-validation errors')
      plot(fraction, cv.error, type='l',
           ylim = range(c(cv.error - cv.se, cv.error + cv.se)),
           xlab = switch(mode, norm='Norm fraction',
             lambda='log(lambda) fraction'), ylab = lab, main = lab)
      if (se) segments(fraction, cv.error - cv.se, fraction, cv.error + cv.se)
    }
    invisible(object)
}

print.glmpath <- function(x, ...)
  {
    cat('Call:\n')
    dput(x$call)
    actions <- x$actions
    xn <- x$xnames[-1]
    k <- length(actions)
    for (i in 1:k) {
      if (length(actions[[i]]) > 0) {
        cat('Step', i, ':')
        for (ii in actions[[i]]) {
          if (ii > 0) cat(' ', xn[ii])
          else cat(' -', xn[-ii])
        }
        cat('\n')
      }
    }
  }

summary.glmpath <- function(object, ...)
  {
    cat('Call:\n')
    dput(object$call)
    ii <- object$new.A
    ii[length(ii)] <- TRUE
    M <- data.frame(Df = object$df[ii], Deviance = object$deviance[ii],
                    AIC = object$aic[ii], BIC = object$bic[ii])
    dimnames(M)[[1]] <- paste('Step', which(ii), sep = ' ')
    M
  }

bootstrap.path <- function(x, y, data, B, index = NULL,
                           path = c('glmpath','coxpath'),
                           method = c('aic','bic'), trace = FALSE, ...)
  {
    path <- match.arg(path)
    method <- match.arg(method)
    if (!missing(data)) x <- data$x
    n <- nrow(x)
    p <- ncol(x)
    if (!is.null(index)) B <- nrow(index)
    else index <- matrix(sample(c(1:n), n * B, replace=T), nrow = B)
    beta <- matrix(0, B, p)
    if (path == 'glmpath') {
      if (!missing(data)) y <- data$y
      fit <- glmpath(x, y, ...)
      s <- switch(method, aic = which.min(fit$aic), bic = which.min(fit$bic))
      beta0 <- fit$b.corrector[s, -1] * fit$sdx
      for (b in 1:B) {
        bx <- x[index[b, ], ]
        by <- y[index[b, ]]
        fit <- glmpath(bx, by, ...)
        s <- switch(method, aic = which.min(fit$aic), bic = which.min(fit$bic))
        beta[b, ] <- fit$b.corrector[s, -1] * fit$sdx
        if (trace) cat(b)
      }
    } else {
      time <- data$time
      status <- data$status
      fit <- coxpath(data, trace = FALSE, ...)
      s <- switch(method, aic = which.min(fit$aic), bic = which.min(fit$bic))
      beta0 <- fit$b.corrector[s, ] * fit$sdx
      for (b in 1:B) {
        bx <- x[index[b, ], ]
        btime <- time[index[b, ]]
        bstatus <- status[index[b, ]]
        fit <- coxpath(list(x = bx, time = btime, status = bstatus), ...)
        s <- switch(method, aic = which.min(fit$aic), bic = which.min(fit$bic))
        beta[b, ] <- fit$b.corrector[s, ] * fit$sdx
        if (trace) cat(b)
      }
    }
    dimnames(beta) <- list(seq(B), dimnames(x)[[2]])
    attr(beta, 'coefficients') <- beta0
    class(beta) <- 'bootpath'
    beta
  }

plot.bootpath <- function(x, type = c('histogram', 'pairplot'), mfrow = NULL,
                          mar = NULL, ...)
  {
    type <- match.arg(type)
    beta0 <- attr(x, 'coefficients')
    if (!is.null(mar)) par(mar = mar)
    if (type == 'histogram') {
      p <- ncol(x)
      if (!is.null(mfrow)) par(mfrow = mfrow)
      else par(mfrow = c(2, ceiling(p / 2)))
      for (j in 1:p) {
        hist(x[, j], main = dimnames(x)[[2]][j], xlab = 'Bootstrap coefficient',
             ylab = 'Frequency', freq=T, ...)
        if (sum(x[, j] == 0) > 0) segments(0, 0, 0, sum(x[, j] == 0), lwd = 3)
        abline(v = beta0[j], col = 2)
      }
    } else {
      bpanel <- function(x, y) {
        abline(v = 0, h = 0, lwd = 2, col = 3)
        points(x, y)
        points(x[1], y[1], pch = 16, col = 2)
      }
      pairs(rbind(beta0, x), panel = bpanel, ...)
    }
  }
