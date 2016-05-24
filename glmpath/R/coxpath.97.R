step.length.cox <- function(corrector, lambda2, x, d, rslist, wlist,
                            min.lambda, max.arclength, frac.arclength,
                            add.newvars, backshoot, approx.Gram, h0 = NULL,
                            eps = .Machine$double.eps)
  {
    active <- corrector$active
    force.active <- corrector$force.active
    lambda <- corrector$lambda - min.lambda
    fk <- length(force.active)
    k <- length(active) - fk
    Lambda2 <- rep(lambda2, length(active))
    C <- corrector$corr
    b <- corrector$b[active]
    eta <- corrector$eta
    wsum <- corrector$wsum
    n <- length(d)
    p1 <- length(active)
    p2 <- ncol(x)
    A <- rep(0, n)
    if (!approx.Gram) AA <- matrix(0, n, n)
    rset <- w <- NULL
    for (i in 1:sum(d == 1)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1, length(rset0)), wlist[[i]]) * eta[rset]
      w1 <- w / wsum[i]
      A[rset] <- A[rset] + w1 - w1^2
      if (!approx.Gram) {
        k <- length(rset)
        AA[1:k, 1:k] <- AA[1:k, 1:k] - outer(w1[rset], w1[rset])
      }
    }
    if (approx.Gram) dhdb <- t(x[, active, drop = FALSE] * A) %*% x
    else {
      diag(AA) <- A
      dhdb <- t(x[, active, drop = FALSE]) %*% AA %*% x
    }
    if (is.null(force.active)) {
      if (length(active) == 1) {
        db <- sign(C[active]) / (dhdb[, active] + lambda2)
      } else {
        db <- (solve(dhdb[, active, drop = FALSE] + diag(Lambda2)) %*%
               sign(C[active]))
      }
    } else {
      db <- (solve(dhdb[, active, drop = FALSE] + diag(Lambda2)) %*%
             c(rep(0, length(force.active)), sign(C[active[-force.active]])))
    }
    newa <- NULL
    if (!backshoot) {
      Cmax <- max(abs(C))
      inactive <- seq(C)[-active]
      ninact <- length(inactive)
      a <- drop(t(db) %*% dhdb[, inactive, drop = FALSE])
      gam <- c((Cmax - C[inactive]) / (1 - a), (Cmax + C[inactive]) / (1 + a))
      ha <- min(gam[gam > eps], lambda)
      if (is.null(force.active)) hd <- -b / db
      else hd <- -b[-force.active] / db[-force.active]
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
      if (is.null(force.active)) hd <- b / db
      else hd <- b[-force.active]/db[-force.active]
      ii <- hd > eps & hd < -h0
      if (any(ii)) h <- -max(hd[ii])
      else h = 0
    }
    list(h = -h, db = db, newa = newa, arclength = h * sum(abs(db)))
  }

predictor.cox <- function(b, step)
  {
    b - step$db * step$h
  }

corrector.cox <- function(x, d, rslist, wlist, rept, method, active, tmpa,
                          force.active, lambda, lambda2, b0,
                          bshoot.threshold, relax.lambda, trace,
                          eps = .Machine$double.eps)
  {
    nobs <- nrow(x)
    p <- length(tmpa)
    fa <- length(force.active)
    if (p > 0) {
      b2 <- c(pmax(b0[tmpa], 0), -pmin(b0[tmpa], 0))
      xa <- x[, tmpa, drop = FALSE]
      penalty <- rep(1, p)
      penalty[force.active] <- 0
      z <- c(as.vector(xa), d, rept, penalty, rep(0, nobs), rep(0, nobs))
      mz <- c(nobs, method, lambda, lambda2, 0)
      sol <- .C('solve_coxpath',
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
      i <- (p + 2) * nobs + p
      eta <- sol[[6]][i + c(1:nobs)]
      i <- i + nobs
      wsum <- sol[[6]][i + c(1:nobs)][d == 1]
      lp <- sol[[7]][5]
    }
    else {
      eta <- rep(1, nobs)
      wsum <- rep(0, sum(d))
    }
    rset <- NULL
    a <- d == 1
    for (i in 1:sum(a)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1, length(rset0)), wlist[[i]]) * eta[rset]
      if (wsum[i] == 0) wsum[i] <- sum(w)
      a[rset] <- a[rset] - w / wsum[i]
    }
    if (p == 0) lp <- -sum(log(wsum))
    corr <- apply(x * a, 2, sum) - lambda2 * b0
    if (p > fa) {
      i <- which(abs(corr) >= lambda * (1 - relax.lambda))
      newa <- i[!i%in%tmpa]
      newactive <- i[!i%in%active]
      if (fa == 0) {
        i <- which(abs(b0[active]) < eps)
        inactive <- active[i]
      } else {
        i <- which(abs(b0[active[-force.active]]) < eps)
        inactive <- active[-force.active][i]
      }
      active <- active[!active%in%inactive]
      active <- c(active, newactive)
      b0[-active] <- 0
    } else {
      if (fa == 0) lambda <- max(abs(corr))
      else lambda <- max(abs(corr[-active]))
      c1 <- which(abs(corr) == lambda)
      newa <- newactive <- c1
      inactive <- NULL
      active <- c(active, c1)
    }
    df <- length(active) - length(newactive)
    backshoot <- ifelse(any(abs(b0[newactive]) > bshoot.threshold),
                        TRUE, FALSE)
    list(eta = eta, wsum = wsum, b = b0, lp = lp, active = active,
         force.active = force.active, newactive = newactive, newa = newa,
         inactive = inactive, corr = corr, lambda = lambda, df = df,
         backshoot = backshoot)
  }

coxpath <- function(data, nopenalty.subset = NULL,
                    method = c('breslow', 'efron'), lambda2 = 1e-5,
                    max.steps = 10 * min(n, m), max.norm = 100 * m,
                    min.lambda = (if (m >= n) 1e-3 else 0), max.vars = Inf,
                    max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
                    bshoot.threshold = 0.1, relax.lambda = 1e-7,
                    approx.Gram = FALSE, standardize = TRUE,
                    eps = .Machine$double.eps, trace = FALSE)
  {
    call <- match.call()
    method <- match.arg(method)
    mthd <- switch(method, breslow = 1, efron = 2)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- nrow(x)
    m <- ncol(x)
    o <- order(status)
    oo <- o[order(time[o], decreasing = TRUE)]
    x <- x[oo, ]
    time <- time[oo]
    status <- status[oo]
    complete <- which(status == 1)
    nnc <- length(complete)
    rept <- rep(0, n)
    for (i in complete) rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    rslist <- wlist <- vector('list', length = nnc)
    for (i in 1:nnc) {
      if (i == 1) {
        ii <- time >= time[complete[i]]
        rslist[[i]] <- which(ii)
      } else if (rept[complete[i]] >= rept[complete[i] - 1]) {
        ii <- (time >= time[complete[i]]) & (time < time[complete[i - 1]])
        rslist[[i]] <- which(ii)
      }
      wlist[[i]] <- rep(1, sum(ii))
      if (mthd == 2) {
        if (rept[complete[i]] > 0) {
          tie <- time[ii] == time[complete[i]] & status[ii] == 1
          di <- max(rept[ii][tie])
          wlist[[i]][tie] <- wlist[[i]][tie] - (di - rept[complete[i]]) / di
        }
      }
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
    if (is.null(dimnames(x)[[2]])) xnames <- paste('x', seq(m), sep='')
    else xnames <- dimnames(x)[[2]]
    lam.vec <- step.len <- rep(0, max.steps)
    bmat.pred <- bmat.corr <- cmat <- matrix(0, nrow = max.steps, ncol = m)
    lp <- df <- new.df <- rep(0, max.steps)
    new.A <- rep(FALSE, max.steps)
    actions <- vector('list', length = max.steps)
    backshoot <- FALSE
    b <- rep(0, m)
    if (!is.null(nopenalty.subset)) {
      force.active <- c(1:length(nopenalty.subset))
      names(nopenalty.subset) <- xnames[nopenalty.subset]
    } else {
      force.active <- NULL
    }
    corrector <- corrector.cox(x, status, rslist, wlist, rept, mthd,
                               nopenalty.subset, nopenalty.subset,
                               force.active, 0, lambda2, b,
                               bshoot.threshold, relax.lambda, trace)
    if (trace && !is.null(force.active))
      cat('The model begins with',xnames[nopenalty.subset], '\n')
    k <- 1
    b <- bmat.pred[k, ] <- bmat.corr[k, ] <- corrector$b
    cmat[k, ] <- corrector$corr
    lam.vec[k] <- lambda <- corrector$lambda
    new.df[k] <- df[k] <- corrector$df
    lp[k] <- corrector$lp
    new.A[k] <- TRUE
    actions[[k]] <- active <- corrector$active
    names(actions[[k]]) <- xnames[active]
    if (trace) {
      cat(paste('Lambda=', lambda, 'lets the first factor in.\nStep', k, ':'))
      if (is.null(force.active))
        cat(paste('\t', xnames[active],'added'))
      else cat(paste('\t', xnames[active[-force.active]], 'added'))
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
        step <- step.length.cox(corrector, lambda2, x, status, rslist, wlist,
                                min.lambda, arclength, frac.arclength1,
                                add.newvars, backshoot, approx.Gram)
        b[active] <- predictor.cox(b[active], step)
        bmat.pred[k, ] <- b
        step.len[k - 1] <- h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        tmpa <- c(active, step$newa)
        if (is.null(force.active)) a <- abs(b[tmpa])
        else a <- abs(b[tmpa[-force.active]])
      } else {
        if (trace) cat(paste('\nStep', k, ':'))
        step <- step.length.cox(corrector, lambda2, x, status, rslist, wlist,
                                min.lambda, Inf, 1, add.newvars, backshoot,
                                approx.Gram, h)
        step.len[k - 1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        if (is.null(force.active)) a <- abs(b[tmpa])
        else a <- abs(b[tmpa[-force.active]])
      }
      corrector <- corrector.cox(x, status, rslist, wlist, rept, mthd,
                                 active, tmpa, force.active, lambda, lambda2,
                                 b, bshoot.threshold, relax.lambda, trace)
      newa <- corrector$newa
      while(length(newa) > 0) {
        if (trace) cat(paste('\nRepeating step', k, ':'))
        tmpa <- c(tmpa, newa)
        if (is.null(force.active)) a <- abs(b[tmpa])
        else a <- abs(b[tmpa[-force.active]])
        corrector <- corrector.cox(x, status, rslist, wlist, rept, mthd,
                                   active, tmpa, force.active, lambda, lambda2,
                                   b, bshoot.threshold, relax.lambda, trace)
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
          actions[[k]] <- newaction
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
      }
      else if (length(corrector$active) <= n) {
        active <- corrector$active
        b <- corrector$b
        backshoot <- FALSE
        n.repeat1 <- max(n.repeat1 - 1, 1)
      }
      if (!backshoot) {
        bmat.corr[k, ] <- b
        cmat[k, ] <- corrector$corr
        lp[k] <- corrector$lp
        df[k] <- corrector$df
        if (lambda <= min.lambda || k == max.steps ||
            length(corrector$active) > min(n, max.vars) ||
            sum(corrector$a) >= max.norm) {
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
    cmat <- cmat[1:k, ,drop=F]
    bmat.pred <- scale(bmat.pred, FALSE, sdx)
    bmat.corr <- scale(bmat.corr, FALSE, sdx)
    dimnames(bmat.pred) <- dimnames(bmat.corr) <- dimnames(cmat) <-
      list(seq(k), xnames)
    df <- df[1:k]
    lp <- lp[1:k]
    aic <- -2 * lp + 2 * df
    bic <- -2 * lp + log(n) * df
    object <- list(call = call, lambda = lam.vec[1:k], lambda2 = lambda2,
                   step.length = abs(step.len[1:(k-1)]), corr = cmat,
                   new.df = new.df[1:k], df = df, loglik = lp, aic = aic,
                   bic = bic, b.predictor = bmat.pred, b.corrector = bmat.corr,
                   new.A = new.A[1:k], actions = actions[1:k], meanx = meanx,
                   sdx = sdx, xnames = xnames, method = method,
                   nopenalty.subset = nopenalty.subset,
                   standardize = standardize)
    class(object) <- 'coxpath'
    object
  }

plot.coxpath <- function(x, xvar = c('norm', 'lambda', 'step'),
                         type = c('coefficients', 'aic', 'bic'),
                         plot.all.steps = FALSE, xlimit = NULL,
                         predictor = FALSE, omit.zero = TRUE, breaks = TRUE,
                         mar = NULL, main = NULL, eps = .Machine$double.eps,
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
    coef.pred <- scale(object$b.predictor[ii, ], FALSE, 1 / object$sdx)
    coef.corr <- scale(object$b.corrector[ii, ], FALSE, 1 / object$sdx)
    xnames <- object$xnames
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
    if (!is.null(mar)) par(mar = mar)
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
          segments(s[xi][-k], coef.corr[xi, ][-k,i], s[xi][-1],
                   coef.pred[xi, ][-1,i], lty = 2, col = i)
        }
      }
    }
    if (breaks) {
      new <- object$new.A[ii] & xi
      axis(3, at = s[new], labels = object$new.df[ii][new], cex = 0.8)
      abline(v = s[new])
    }
  }

predict.coxpath <- function(object, data, s,
                            type = c('coefficients', 'loglik', 'lp', 'risk',
                              'coxph'),
                            mode = c('step', 'norm.fraction', 'norm',
                              'lambda.fraction', 'lambda'),
                            eps = .Machine$double.eps, ...)
  {
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(data) && type != 'coefficients') {
        warning('No data argument; type switched to coefficients')
        type <- 'coefficients'
    }
    if (!missing(s)) {
      if (length(s) > 1 && type == 'coxph') {
        warning('Length(s) > 1. Only the first element is used.')
        s <- s[1]
      }
    }
    b <- object$b.corrector
    std.b <- scale(b, FALSE, 1 / object$sdx)
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
    b <- b[useq, ,drop=F]
    coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- (((sb[right] - sfrac) * b[left, , drop = FALSE] +
              (sfrac - sb[left]) * b[right, , drop = FALSE]) /
             (sb[right] - sb[left]))
    newb[left == right, ] <- b[left[left == right], ]
    coef <- newb
    if (type == 'coefficients') {
      fit <- coef
      dimnames(fit) <- list(s, object$xnames)
    } else if (type == 'loglik') {
      fit <- logplik(data$x, data$time, data$status, t(coef), object$method)
      names(fit) <- s
    } else if (type == 'lp' || type == 'risk') {
      b0 <- coef %*% object$meanx
      fit <- scale(data$x %*% t(coef), b0, FALSE)
      if (type == 'risk') fit <- exp(fit)
      dimnames(fit) <- list(seq(nrow(data$x)), s)
    } else {
      coef <- drop(coef)
      active <- abs(coef) > eps
      coef <- coef[active]
      x <- data$x[, active, drop = FALSE]
      time <- data$time
      status <- data$status
      fit <- coxph(Surv(time, status) ~ x, method = object$method)
      junk <- logplik(x, time, status, coef, object$method, TRUE)
      w <- junk$w
      dmat <- junk$dmat
      oo <- junk$oo
      a <- sum(active)
      info <- matrix(0, a, a)
      for (i in 1:sum(status == 1)) {
        ind <- dmat[, i] > 0
        xr <- x[oo[ind], , drop = FALSE]
        wr <- w[ind, i]
        v1 <- xr * wr
        v2 <- apply(v1, 2, sum)
        info <- info + t(xr) %*% v1 - outer(v2, v2)
      }
      fit$coefficients <- coef
      fit$var <- solve(info)
      fit$loglik <- c(fit$loglik[1], junk$loglik)
      fit$iter <- fit$residuals <- NULL
      fit$linear.predictors <- junk$eta - sum(coef * object$meanx[active])
      fit$method <- object$method
      fit$assign <- seq(a)
      fit$wald.test <- sum(coef*(info %*% coef))
    }
    attr(fit, 's') <- s
    attr(fit, 'fraction') <- sfrac
    attr(fit, 'mode') <- mode
    return(fit)
  }

logplik <- function(x, time, status, b, method = c('breslow', 'efron'),
                    return.all = FALSE)
  {
    method <- match.arg(method)
    n <- length(time)
    o <- order(status, decreasing=T)
    oo <- o[order(time[o])]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    complete <- which(status == 1)
    nnc <- length(complete)
    if (nnc == 0) {
      stop('No complete observation. Failed to compute partial likelihood.')
    }
    dmat <- matrix(0, n, nnc)
    for (i in 1:nnc) {
      dmat[time >= time[complete[i]], i] <- 1
      if (method == 'efron') {
        if (rept[complete[i]] > 0) {
          tie <- time == time[complete[i]] & status == 1
          di <- max(rept[tie])
          dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]]) / di
        }
      }
    }
    eta <- x %*% b
    eeta <- exp(eta)
    k <- ncol(eta)
    loglik <- rep(0, k)
    for (i in 1:k) {
      w <- dmat * eeta[oo, i]
      wsum <- apply(w, 2, sum)
      loglik[i] <- sum(eta[oo, i][status == 1]) - sum(log(wsum))
    }
    if (return.all) {
      return(list(loglik = loglik, w = scale(w, F, wsum), eta = eta,
                  dmat = dmat, oo = oo))
    } else {
      return(loglik)
    }
  }

cv.coxpath <- function(data, method = c('breslow', 'efron'), nfold = 5,
                       fraction = seq(0, 1, length = 100),
                       mode = c('norm', 'lambda'), plot.it = TRUE, se = TRUE,
                       ...)
  {
    method <- match.arg(method)
    mode <- match.arg(mode)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- length(time)
    folds <- split(sample(seq(n)), rep(1:nfold, length = n))
    errormat <- matrix(0, length(fraction), nfold)
    for (i in seq(nfold)) {
      omit <- folds[[i]]
      trdata <- list(x = x[-omit, ], time = time[-omit],
                     status = status[-omit])
      tsdata <- list(x = x[omit, ], time = time[omit], status = status[omit])
      fit <- coxpath(trdata, method = method, ...)
      pred <- switch(mode, norm = {
        predict(fit, tsdata, fraction, 'loglik', 'norm.fraction')
      }, lambda = {
        predict(fit, tsdata, fraction, 'loglik', 'lambda.fraction')
      })
      if (length(omit) == 1) pred <- matrix(pred, nrow = 1)
      errormat[, i] <- -pred / length(omit)
      cat('CV Fold', i, '\n')
    }
    cv.error <- apply(errormat, 1, mean)
    cv.se <- sqrt(apply(errormat, 1, var) / nfold)
    object <- list(fraction = fraction, cv.error = cv.error, cv.se = cv.se,
                   folds = folds)
    if (plot.it) {
      plot(fraction, cv.error, type = 'l',
           ylim = range(c(cv.error - cv.se, cv.error + cv.se)),
           xlab = switch(mode, norm = 'Norm fraction',
             lambda = 'log(lambda) fraction'),
           ylab = 'Minus log-partial-likelihood',
           main = 'Cross-validated minus log-partial-likelihood')
      if (se) segments(fraction, cv.error - cv.se, fraction, cv.error + cv.se)
    }
    invisible(object)
}

print.coxpath <- function(x, ...)
  {
    cat('Call:\n')
    dput(x$call)
    actions <- x$actions
    xn <- x$xnames
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

summary.coxpath <- function(object, ...)
  {
    cat('Call:\n')
    dput(object$call)
    ii <- object$new.A
    ii[length(ii)] <- TRUE
    M <- data.frame(Df = object$df[ii], Log.p.lik = object$loglik[ii],
                    AIC = object$aic[ii], BIC = object$bic[ii])
    dimnames(M)[[1]] <- paste('Step', which(ii), sep=' ')
    M
  }
