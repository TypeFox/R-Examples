`drm` <-
  function (formula, family = binomial, data = sys.parent(), weights, 
            offset, subset = NULL, na.action, start = NULL, link = "cum", 
            dep = "I", Ncond = TRUE, Lclass = 2, dropout = FALSE, drop.x = NULL, 
            save.profiles = TRUE, pmatrix = NULL, print.level = 2, iterlim = 200, 
            ...) 
{
  if (length(dep) > 1) {
    assoc <- dep[2:length(dep)]
    dep <- dep[[1]]
  }
  else {
    dep <- unlist(dep)
    assoc <- NULL
  }
  match.arg(dep, c("I", "N", "L", "B", "D", "NL", "NB", "ND", 
                   "M", "NM", "LM", "NLM", "M2", "NM2"))
  call <- match.call()
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  m <- match.call(expand = FALSE)
  m$family <- m$dep <- m$start <- m$link <- m$dropout <- m$pmatrix <- NULL
  m$print.level <- m$Ncond <- m$Lclass <- m$iterlim <- m$save.profiles <- m$... <- NULL
  if (missing(na.action)) 
    m$na.action <- "na.include"
  Terms <- if (missing(data)) 
    terms(formula, c("cluster", "Time"))
  else terms(formula, c("cluster", "Time"), data = data)
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  cluster <- attr(Terms, "specials")$cluster
  Time <- attr(Terms, "specials")$Time
  if (!length(cluster)) 
    stop("Cluster not specified. Use e.g. y ~ x + cluster(id)")
  if (length(cluster) > 1) 
    stop("Only one cluster allowed. Use e.g. y ~ x + cluster(id)")
  if (!length(Time)) 
    if (length(grep("M", dep))) 
      stop("Time not specified for the M-structure:\n use e.g. y ~ x + cluster(id) + Time(time)")
  if (length(Time) > 1) 
    stop("Only one Time-dimension allowed.\n Use e.g. y ~ x + cluster(id) + Time(time)")
  tempc <- untangle.specials(Terms, "cluster", 1:10)
  tempt <- untangle.specials(Terms, "Time", 1:10)
  cl <- strata(m[, tempc$vars], shortlabel = TRUE)
  ord <- attr(Terms, "order")[c(tempc$terms, tempt$terms)]
  if (any(ord > 1)) 
    stop("cluster() or Time() can not be used in an interaction")
  if (length(grep("M", dep)) > 0) {
    ord <- c(cluster, Time)
    note <- paste("Note:  ordering of the responses within each cluster is by Time")
  }
  else {
    if (any(is.na(m))) {
      ord <- c(cluster, attr(Terms, "response"))
      note <- paste("Note:  ordering of the responses within each cluster is by response value (with NA last)")
    }
    else {
      ord <- cluster
      note <- NULL
    }
  }
  mord <- eval(parse(text = paste("order(", paste("m[,", ord, 
                       "]", collapse = ","), ")")))
  m <- m[mord, ]
  X <- model.matrix(Terms, m)
  dropx <- match(c(tempc$vars, tempt$vars), dimnames(X)[[2]])
  marg.names <- X[, dropx, drop = FALSE]
  X <- X[, -(dropx), drop = FALSE]
  if (any(is.na(X))) 
    warning("Missing covariate values detected")
  cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
  dimnames(X)[[1]] <- sort(as.numeric(as.factor(cl)))
  nrep <- table(table(cluster))
  if (length(nrep) > 1) 
    stop(paste("Can handle only balanced data frame.", "For unbalanced data frames, fill NA's to obtain equal cluster size", 
               sep = "\n"))
  nrep <- as.numeric(names(nrep))
  wt <- model.extract(m, weights)
  if (length(wt) > 0) {
    nwt <- table(unlist(tapply(wt, cluster, table)))
    if (length(nwt > 1)) 
      stop("Only equal weights per cluster allowed")
    wt <- wt[!duplicated(cluster)]
  }
  xint <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  n <- nrow(X)
  pc <- ncol(X) - 1
  if (!length(wt)) 
    wt <- 1
  offset <- model.extract(m, offset)
  if (length(offset) == 1) 
    offset <- rep(offset, n)
  if (length(offset)==0)
    offset <- rep(0,n)
  y <- model.extract(m, "response")
  if (any(is.na(y))) 
    cat(paste("Note: Missing response values; Assuming", 
              if (dropout) 
              "MNAR"
              else ("MAR"), "\n"))
  nclass <- length(levels(factor(y)))
  if (xint == 0 && nclass > 2) 
    stop("Can't fit without intercept(s)")
  mis <- apply(cbind(y, X), 1, function(i) any(is.na(i)))
  x <- X[!mis, drop = FALSE, ]
  y1 <- ifelse(as.numeric(as.factor(y)) > 1, 1, 0)[!mis]
  if (is.null(start)) {
    parm1 <- glm.fit(x, y1, family = family)$coefficients
    if (nclass > 2) {
      theta <- matrix(0, nrow = length(parm1), ncol = nclass - 
                      2)
      for (i in 2:(nclass - 1)) {
        y1 <- ifelse(as.numeric(as.factor(y)) > i, 1, 
                     0)[!mis]
        theta[, i - 1] <- glm.fit(x, y1, family = family)$coefficients
      }
      parms <- c(as.vector(-parm1[1]), (-theta[1, ]), parm1[-1])
      names(parms)[1:(nclass - 1)] <- paste("(Intercept)", 
                                            1:(nclass - 1), sep = "")
      if (link == "bcl") {
        parms[1:(nclass - 1)] <- log(table(y)[nclass]) - 
          log(table(y)[-nclass])
        bclreg <- as.vector(-theta[-1, ])
        names(bclreg) <- rep(names(parm1)[-1], nclass - 
                             2)
        parms <- c(-parms, bclreg)
        names(parms)[-(1:(nclass - 1))] <- paste(names(parms)[-(1:(nclass - 
                                                                   1))], rep(1:(nclass - 1), rep(nrow(theta) - 
                                                                                                 1, nclass - 1)), sep = "")
      }
      if (link == "acl") {
        parms[1:(nclass - 1)] <- log(table(y)[-nclass]) - 
          log(table(y)[-1])
      }
    }
    else (parms <- parm1)
  }
  else {
    parms <- rep(0, ncol(x) + (nclass - 2))
    names(parms)[1:(nclass - 1)] <- paste("(Intercept)", 
                                          if (nclass > 2) 
                                          1:(nclass - 1)
                                          else (""), sep = "")
    names(parms)[nclass:(ncol(x) + nclass - 2)] <- dimnames(x)[[2]][-1]
  }
  npar <- length(parms)
  if (length(grep("N", dep)) > 0) 
    parms <- c(parms, nu = 1)
  if (length(grep("M", dep)) > 0) {
    if (length(grep("M2", dep)) > 0) {
      if (nclass > 2) 
        stop("M2 currently only allowed for binary response")
      if (nrep == 2) 
        stop("M2 only for 3 or more repeated measures: Use M instead")
      parms <- c(parms, tau12 = 1, tau13 = 1, tau123 = 1)
    }
    else {
      mm <- expand.grid(c(1:nclass), c(1:nclass))
      tn <- apply(mm - 1, 1, paste, collapse = "")
      p2 <- rep(1, length(tn[-grep("0", tn)]))
      names(p2) <- paste("tau", if (nclass > 2) 
                         tn[-grep("0", tn)], sep = "")
      parms <- c(parms, p2)
    }
  }
  if (length(grep("L", dep)) > 0) {
    if (Lclass > 2) {
      if (nclass > 2) 
        stop("Lclass > 2 for structure L currently implemented only for binary responses")
      if (length(grep("LM", dep)) > 0) 
        stop("Lclass > 2 for structure LM currently not implemented")
      else {
        p2 <- c(rep(0.01/(Lclass - 2), Lclass - 2), 0.99, 
                c(1:(Lclass - 1))/Lclass)
        names(p2) <- c(paste("nu", 1:(Lclass - 1), sep = ""), 
                       paste("kappa", 0:(Lclass - 2), sep = ""))
      }
    }
    else {
      p2 <- c(nu1 = 0.99, rep(0.5, nclass - 1))
      if (nclass == 2) 
        names(p2) <- c("nu1", "kappa")
      else (names(p2) <- c("nu1", paste("kappa", 1:(nclass - 
                                                    1), sep = "")))
    }
    parms <- c(parms, p2)
  }
  if (length(grep("B", dep)) > 0) {
    p2 <- rep(c(1, 0), (nclass - 1))
    names(p2) <- paste(c("xi", "xi"), rep((nclass - 1):0, 
                                          nclass - 1), sep = "")
    parms <- c(parms, p2)
  }
  if (length(grep("D", dep)) > 0) {
    p2 <- c(1e-04, rep(5, (nclass - 1)))
    names(p2) <- c(paste("xi", 0:(nclass - 1), sep = ""))
    parms <- c(parms, p2)
  }
  else (parms <- parms)
  Y <- matrix(as.numeric(factor(y)), ncol = nrep, byrow = TRUE)
  N <- apply(Y - 1, 1, sum, na.rm = TRUE)
  if (dropout) {
    yd <- sort(na.omit(unique(y)))
    if (!is.na(match("(drop.x)", names(m)))) {
      dropc <- matrix(m[, "(drop.x)"], ncol = nrep, byrow = TRUE)
      dropx <- expand.grid(list(yd, yd, sort(na.omit(unique(dropc)))))
      droplab <- dropx
      droplab[, 1:2] <- sapply(1:2, function(i) as.numeric(factor(droplab[, 
                                                                          i])))
      droplab <- apply(droplab, 1, paste, collapse = "")
      names(dropx) <- c(paste(names(m)[1], c("cur", "prev"), 
                              sep = "."), paste(deparse(substitute(drop.x)), 
                                "prev", sep = "."))
      misn <- lapply(seq(nrow(Y)), function(i, Y, dropc) c(Y[i, 
                                                             which(is.na(Y[i, ]) == TRUE)[1] - 1], dropc[i, 
                                                                          which(is.na(Y[i, ]) == TRUE)[1] - 1]), Y = Y, 
                     dropc = dropc)
    }
    else {
      dropc <- rbind(rep("", nrep))
      dropx <- expand.grid(list(yd, yd))
      droplab <- apply(sapply(dropx, function(i) as.numeric(factor(i))), 
                       1, paste, collapse = "")
      names(dropx) <- c(paste(names(m)[1], c("cur", "prev"), 
                              sep = "."))
      misn <- lapply(seq(nrow(Y)), function(i, Y, dropc) c(Y[i, 
                                                             which(is.na(Y[i, ]) == TRUE)[1] - 1]), Y = Y, 
                     dropc = dropc)
    }
    dropx <- eval(parse(text = paste("model.matrix(~", paste(names(dropx), 
                          collapse = "+"), ",data=dropx)", collapse = "")))
    dimnames(dropx)[[2]][1] <- paste(dimnames(dropx)[[2]][1], 
                                     "d", sep = "")
    dropstart <- length(Y[, 1][is.na(apply(Y, 1, sum))])/(length(Y[, 
                                                                   1]) * (nrep - 1))
    if (dropstart == 0) 
      stop("Data complete; set dropout=FALSE")
    dropstart <- c(log(dropstart/(1 - dropstart)), rep(0, 
                                                       ncol(dropx) - 1))
    names(dropstart) <- dimnames(dropx)[[2]]
    parms <- c(parms, dropstart)
    cat("Calculating kronecker products...\n")
    w <- lapply(1:nrow(Y), function(i, Y, nclass, dep) {
      kroneckerd.drm(Y[i, ], nclass = nclass, dep = dep)
    }, nclass = nclass, Y = Y, dep = dep)
    if (length(grep("M", dep)) > 0) {
      dims <- dim(kronecker.drm(Y[1, ], nclass = nclass, 
                                dep = dep))
      Nf <- sapply(w, function(i) rep(0, dim(i)[3]))
    }
    else (Nf <- sapply(w, function(i) rep(0, ncol(i))))
    for (i in 1:length(Nf)) if (N[i] == 0) 
      Nf[[i]][1] <- 1
  }
  else {
    cat("Calculating kronecker products...\n")
    w <- apply(Y, 1, function(i, Y, nclass, dep) {
      kronecker.drm(i, nclass = nclass, dep = dep)
    }, nclass = nclass, dep = dep)
    if (length(grep("M", dep)) > 0) {
      dims <- dim(kronecker.drm(Y[1, ], nclass = nclass, 
                                dep = dep))
      w <- array(w, dim = c(dims, nrow(Y)))
    }
    Nf <- apply(Y - 1, 1, sum, na.rm = TRUE)
    drop.cov <- dropx <- droplab <- NULL
  }
  assoc <- getass.drm(dep = dep, assoc = assoc, parms = parms, 
                      nrep = nrep, npar = npar, nclass = nclass, data = data, 
                      subset = subset, cluster = cl)
  parms <- assoc$parms
  equal <- assoc$equal
  ass <- assoc$tfunction
  asscov <- assoc$ass
  constant <- assoc$constant
  if (length(grep("M2", dep)) > 0) {
    if (nclass > 2) 
      stop("structure M2 currently implemented only for binary response")
  }
  if (!is.null(assoc$tfunction)) {
    tlen <- cumsum(sapply(assoc$tfunction, function(i) length(formals(i))))
    tlen <- cbind(c(0, tlen[-length(tlen)]) + 1, tlen)
    tlen <- lapply(1:nrow(tlen), function(i, tlen) tlen[i, 
                                                        1]:tlen[i, 2], tlen = tlen)
  }
  if (!is.null(start)) {
    if (length(start) != length(parms)) 
      stop("Wrong number of initial parameter values")
    parms <- structure(start, names = names(parms))
  }
  cat("Initial parameter values\n")
  print(round(parms, 3))
  cat("\n")
  if (any(duplicated(names(parms)))) 
    stop(paste("Coefficient-names conflict with default association parameter names.\n", 
               "Provide different variable names for the regression coefficients"))
  f <- rep(list(1:nclass), if (length(grep("M", dep)) > 0) {
    if (length(grep("M2", dep)) > 0) 3 else 2
  } else nrep)
  mm <- expand.grid(f)
  tm <- sapply(1:(nclass - 1), function(i, mm) {
    apply(mm, 1, function(k, i) {
      length(k[k == i + 1])
    }, i = i)
  }, mm = mm)
  if (length(asscov$assmodel) == 0) 
    asscov$alab <- rep(1, n/nrep)
  Xlab <- apply(matrix(apply(cbind(offset, X, rep(asscov$alab, rep(nrep, n/nrep))), 1, paste, collapse = ""),
                       ncol = nrep, byrow = TRUE), 1, paste, collapse = "")
  lab <- match(dimnames(X)[[1]], grep("FALSE", paste(duplicated(Xlab))))
  X.org <- X
  offset.org <- offset
  X <- X[!is.na(lab), drop = FALSE, ]
  offset <- offset[!is.na(lab)]
  dimnames(X) <- NULL
  lab <- as.numeric(factor(Xlab, levels = unique(Xlab)))
  if (length(asscov$assmodel) > 0) {
    if (length(grep("tau", unlist(sapply(asscov$assmodel, 
                                         dimnames)))) > 0) {
      if (length(ass) > 0) 
        stop("Covariates for the M-association parameters currently not implemented")
    }
  }
  if (dropout) {
    if (length(asscov$assmodel) > 0) 
      stop("Both dropout=TRUE and covariates for the association currently not implemented")
    .loglik <- get("loglikd.drm")
  }
  else (.loglik <- get("logliks.drm"))
  tim <- system.time(llinit <- .loglik(parms, y = Y, X = X, 
                                       dep = dep, asscov = asscov, tlen = tlen, equal = equal, 
                                       npar = npar, nclass = nclass, nrep = nrep, tm = tm, w = w, 
                                       lab = lab, inv = family$linkinv, ass = ass, link = link, 
                                       constant = constant, Ncond = Ncond, Lclass = Lclass, 
                                       offset = offset, wt = wt, Nf = Nf, misn = misn, drop.cov = dropc, 
                                       dropx = dropx, droplab = droplab))
  if (llinit == .Machine$double.xmax) 
    stop("Poor initial regression parameter values: set new with `start'-argument")
  cat(paste("Initial log likelihood", -round(llinit, 2), ";", 
            "Function evaluation time:", round(tim[3], 2), "sec.\n"))
  cat("Starting numerical optimisation of log likelihood...\n\n")
  est <- nlm(f = .loglik, p = parms, hessian = TRUE, y = Y, 
             X = X, dep = dep, asscov = asscov, npar = npar, nclass = nclass, 
             nrep = nrep, tm = tm, w = w, print.level = print.level, 
             lab = lab, equal = equal, inv = family$linkinv, ass = ass, 
             link = link, constant = constant, Ncond = Ncond, Lclass = Lclass, 
             tlen = tlen, offset = offset, wt = wt, Nf = Nf, misn = misn, 
             dropx = dropx, droplab = droplab, drop.cov = dropc, iterlim = iterlim, 
             ...)
  estimate <- est$estimate
  names(estimate) <- names(parms)
  if (qr(est$hessian)$rank < ncol(est$hessian)) {
    warning("Apparently a singular Hessian matrix. Dumping out estimate and Hessian")
    list(call = call, coefficients = estimate, hessian = est$hessian)
  }
  else {
    if (is.null(pmatrix)) {
      v <- expand.grid(rep(list(1:nclass), nrep))
      cat(paste("Creating matrix for all possible profiles (", 
                nrow(v), "x", nrow(v), ")...\n", sep = ""))
      w <- apply(v, 1, function(i, nclass, dep) {
        kronecker.drm(i, nclass = nclass, dep = dep)
      }, nclass = nclass, dep = dep)
      if (length(grep("M", dep)) > 0) {
        w <- array(w, dim = c(dims, nrow(v)))
      }
    }
    else (w <- get(pmatrix))
    cat(paste("Estimating all possible profiles ...\n", sep = ""))
    prof <- getpath.drm(p = estimate, v = v, X = X, dep = dep, 
                        tlen = tlen, npar = npar, nclass = nclass, equal = equal, 
                        nrep = nrep, tm = tm, w = w, lab = lab, inv = family$linkinv, 
                        ass = ass, link = link, constant = constant, Ncond = Ncond, 
                        Lclass = Lclass, offset = offset, wt = wt, save.profiles = save.profiles, 
                        asscov = asscov)
    cat("\nDone.\n")
    if (save.profiles) {
      if (!is.null(pmatrix)) 
        v <- expand.grid(rep(list(1:nclass), nrep))
      dimnames(prof) <- list(apply(v, 1, paste, collapse = ""), 
                             unique(marg.names[, 1]))
      if (!is.null(note)) 
        attributes(prof)$note <- note
    }
    if (nclass == 2) 
      mu <- matrix(family$linkinv(offset.org + X.org %*% 
                                  estimate[1:npar]), ncol = 1)
    else {
      if (link == "bcl") {
        if (npar == (nclass - 1)) 
          e.eta <- exp(sapply(estimate[1:(nclass - 1)], 
                              function(i, offs) i + offs, offs = offset.org))
        else (e.eta <- exp(offset.org + sapply(0:(nclass - 
                                                  2), function(i, x, p, nclass) {
                                                    x %*% p[c(1 + i, (ncol(x) - 1) * i + (nclass:(ncol(x) + 
                                                                                                  nclass - 2)))]
                                                  }, x = X.org, p = estimate, nclass = nclass)))
        mu <- e.eta/(1 + c(e.eta %*% rep(1, nclass - 
                                         1)))
      }
      else {
        eta <- offset.org
        if (npar > (nclass - 1)) 
          eta <- offset.org + (X.org[, -1, drop = FALSE] %*% 
                               estimate[nclass:npar])
        if (link == "cum") {
          mu <- matrix(family$linkinv(sapply(estimate[1:(nclass - 
                                                         1)], function(i, eta) i - eta, eta = eta)), 
                       ncol = (nclass - 1))
        }
        if (link == "acl") {
          eta <- sapply(estimate[1:(nclass - 1)], function(i, 
                                                           eta) i + eta, eta = eta)
          tri <- lower.tri(matrix(1, ncol = (nclass - 
                                             1), nrow = (nclass - 1)), TRUE)
          mu <- exp(eta %*% tri)/c(1 + exp(eta %*% tri) %*% 
                                   rep(1, (nclass - 1)))
        }
      }
    }
    dimnames(mu) <- list(apply(marg.names, 1, paste, collapse = ";"), 
                         names(estimate)[1:(nclass - 1)])
    mu <- mu[match(seq(nrow(mu)), mord), ]
    if (length(grep("L", dep)) > 0) {
      if (nclass == 2) 
        mucond <- mu
      else {
        if (link != "cum") 
          mucond <- cbind(mu[, -1], 1 - apply(mu, 1, 
                                              sum))
        else (mucond <- cbind(mu[, 2:ncol(mu)], 1) - 
              mu)
      }
      assp <- estimate[(npar + 1):length(estimate)]
      if (!is.null(equal)) {
        if (length(constant[[1]] > 0)) {
          assc <- c(equal, constant[[1]])
          assc[constant[[1]]] <- constant[[2]]
          assc[-constant[[1]]] <- assp[equal]
          assp <- assc
        }
        else (assp <- assp[equal])
      }
      if (is.na(match("nu1", names(estimate))) | Lclass > 
          2) {
        cat("Note: conditional probabilities not calculated; See ?getass.drm for an example\n")
        mucond <- NULL
      }
      else {
        L.p <- assp[match("nu1", names(estimate)[(npar + 
                                                  1):length(estimate)]):(match("nu1", names(estimate)[(npar + 
                                                                                                       1):length(estimate)]) + (nclass - 1))]
        mucond <- (t(t(mucond)/(L.p[1] + (1 - L.p[1]) * 
                                L.p[-1])))
        mucond2 <- (t(t(mucond) * L.p[-1]))
        if (nclass > 2 && link == "cum") {
          mucond <- t(apply(cbind(1 - apply(mucond, 1, 
                                            sum), mucond), 1, cumsum))[, -nclass, drop = FALSE]
          mucond2 <- t(apply(cbind(1 - apply(mucond2, 
                                             1, sum), mucond2), 1, cumsum))[, -nclass, 
                                                                            drop = FALSE]
        }
        dimnames(mucond) <- dimnames(mucond2) <- dimnames(mu)
        mucond <- list("L=1" = mucond, "L=0" = mucond2)
      }
    }
    else (mucond <- NULL)
    covmat <- solve(est$hessian)
    structure(list(call = call, coefficients = estimate, 
                   cov.scaled = covmat, deviance = 2 * est$minimum, 
                   aic = if (dropout) NA else 2 * (est$minimum + length(estimate)), 
                   df.residual = length(y) - length(estimate), df.null = length(y) - 
                   (nclass - 1), code = est$code, niter = est$iterations, 
                   fitted.marginals = mu, fitted.profiles = prof, fitted.conditionals = mucond, 
                   terms = Terms), class = "drm")
  }
}
