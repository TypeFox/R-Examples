aftgee <- function(formula, data, subset, id = NULL, contrasts = NULL,
                   weights = NULL, margin = NULL,
                   corstr="independence",
                   binit = "srrgehan", B = 100,
                   control = aftgee.control()
                   ) {
  scall <- match.call()
  mnames <- c("", "formula", "data", "weights", "margin", "subset", "na.action", "id")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  ##  if (is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  y <- model.extract(m, "response")
  if (ncol(y) > 2) 
      stop("aftgee only supports Surv object with right censoring", call. = FALSE)
  N <- NROW(y)
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts) 
  weights <- model.extract(m, weights)
  if (is.null(weights))
      weights <- rep(1, N)
  id <- model.extract(m, id)
  if (is.null(id))
      id <- 1:nrow(y)
  margin <- model.extract(m, margin)
  if (is.null(margin))
      margin <- rep(1, N)
  out <- NULL
  if (sum(y[,2]) == nrow(y) & corstr %in% c("indep", "independence")) {
      warning("Response is uncensored and correlation structure is independence,
ordinary least squares is used", call. = FALSE)
      out <- lm(log(y[,1]) ~ x - 1)
      out$coef.init <- out$coef.res <- out$coefficients
      out$coefficients <- cbind(out$coefficients, out$coefficients)
      out$var.res <- vcov(out)
  }
  else {
      out <- aftgee.fit(y = y, x = x, id = id, corstr = corstr, B = B, binit = binit,
                        weights = weights, margin = margin, control = control)
  } 
  out$call <- scall
  out$y <- y[,1]
  out$x <- x
  rownames(out$coefficients) <- names(out$coef.res) <- names(out$coef.init) <- colnames(out$x)
  out$intercept <- FALSE
  if (sum(x[,1]) == nrow(x)) 
      out$intercept <- TRUE
  colnames(out$coefficients) <- c("binit", "AFTGEE")
  class(out) <- "aftgee"
  out
}

aftgee.fit <- function(y, x, id, corstr="independence",
                       weights = rep(1, nrow(x)),
                       margin = rep(1, nrow(x)),
                       B = 100, binit = "lm",
                       control = aftgee.control()) {
  x <- as.matrix(x)
  n <- length(unique(id))
  rm <- NULL
  rmName <- NULL
  ## include.intercept <- 0
  ## if ("(Intercept)" %in% colnames(x)) {
  ##     x <- x[, -which(colnames(x) == "(Intercept)")]
  ##     include.intercept <- 1
  ## }
  p <- ncol(x)
  firstBeta <- firstSd <- firstSdMat <- firstconvergence <- NULL
  clsize <- unlist(lapply(split(id, id), length))
  N <- sum(clsize)
  if (is.numeric(binit)) {
      if (length(binit) != p) 
          stop("binit value length does not match with numbers of covariates", call. = FALSE)
      firstBeta <- binit      
  }
  if (!(is.numeric(binit))) {
      if (!(binit %in% c("lm", "srrgehan"))) 
          stop("Invalid binit value method", call. = FALSE)
  }
  Y <- log(y[,1])
  Y <- ifelse(Y == -Inf, 0, Y)
  delta <- y[,2]
  if (!(is.numeric(binit))) {
      if ("(Intercept)" %in% colnames(x)) {
          xtemp <- as.matrix(x[,-(colnames(x) != "(Intercept)")])
      } else {
          xtemp <- x
      }
      if (binit == "lm") {
          linfit <- summary(lm(Y ~ xtemp - 1))
          first <- list(beta = linfit$coef[,1], sd = linfit$coef[,2])
          firstBeta <- first$beta
          firstSd <- first$sd
          firstconvergence <- first$convergence
      }

      if (binit == "srrgehan") {
          first <- aftsrr.fit(Y = Y, delta = delta, X = xtemp,
                              id = id, weights = weights, B = 0, control = control)
          firstBeta <- first$beta
          firstSdMat <- NA
          firstconvergence <- first$convergence
      }
     if ("(Intercept)" %in% colnames(x)) {
         firstBeta <- c(mean(Y - xtemp %*% firstBeta), firstBeta)
      }
  }
  binitValue <- list(beta = firstBeta, sd = firstSd, sdMat = firstSdMat)
  result <- aftgee.est(Y, x, delta, binitValue$beta, id, corstr,
                       rep(1, length(Y)), margin, weights, control)
  ## variance estimation
  sample <- zout <- NULL
  if (B > 0) {
    sample <- matrix(0, nrow = B, ncol = length(result$beta))
    for (i in 1:B){
      Z <- as.vector(rep(rexp(n,1), time = clsize))
      zout <- cbind(zout, Z)
      ## sample[i,] <- aftgee.est(Y, x, delta, result$iniBeta, id, corstr,
      ##                          Z, margin, weights, control)$beta

      sample[i,] <- aftgee.est(y = Y, x = x, delta = delta, beta = result$iniBeta, id, corstr,
                               Z, margin, weights, control)$beta
  }
    ## if (include.intercept == 1) {
    ##     result$beta <- c(eRes(e = Y - x %*% result$beta, delta = delta)[[3]], result$beta)
    ##     binitValue$beta <- c(eRes(e = Y - x %*% binitValue$beta, delta = delta)[[3]],
    ##                          binitValue$beta)
    ##     sample <- cbind(apply(sample, 1, function(beta) eRes(e = Y - x %*% beta, delta)[[3]]),
    ##                     sample)
    ## }
    vhat <- var(sample)
  }
  if (B == 0) {
    vhat <- NULL
  }
  ini.beta <- c(binitValue$beta)
  ini.sd <- c(binitValue$sd)
  ini.sdMat <- c(binitValue$sdMat)
  fit <- list(coefficients = cbind(ini.beta, result$beta),
              coef.res = result$beta,
              var.res = vhat,
              varMargin = result$gamma,
              alpha = result$alpha,
              coef.init = ini.beta,
              sd.init = ini.sd,
              var.init.mat = ini.sdMat,
              binit = binit,
              conv = result$convergence,
              ini.conv = firstconvergence,
              bhat = sample,
              zout = zout, 
              conv.step = result$convStep)
  class(fit) <- "aftgee.fit"
  fit
}

aftgee.control <- function(maxiter = 50,
                           reltol = 0.001,
                           trace = FALSE) {
  list(maxiter = maxiter,
       reltol = reltol,
       trace = trace)
}

aftgee.est <- function(y, x, delta, beta, id, corstr = "independence", Z = rep(1, length(y)),
                       margin = rep(1, length(id)), weights = rep(1, length(y)),
                       control = aftgee.control()) {
    xmat <- as.matrix(x) 
    nobs <- length(y)
    for (i in 1:control$maxiter) {
        betaprev <- beta
        eres <- NULL
        eres2 <- NULL
        if (sum(margin == 1) == nobs) {
            e <- y - xmat %*% beta
            eres <- eRes(e, delta = delta, z = Z * weights)
            yhat <- delta * y + (1 - delta) * (eres[[1]] + xmat %*% beta)
            yhatZ <- sqrt(Z) * yhat
            xmatZ <- sqrt(Z) * xmat
            geefit <- geese.fit(xmatZ, yhatZ, id, corstr = corstr, weights =  weights)
        }
        if (sum(margin == 1) != nobs) {
            e <- y - xmat %*% beta
            ## if (res == FALSE) {
            ##     eres <- eRes(e, delta, Z * weights)
            ##     yhat <- delta * y + (1 - delta) * (eres[[1]] + xmat %*% beta)
            ##     yhatZ <- sqrt(Z * weights) * yhat
            ##     xmatZ <- sqrt(Z * weights) * xmat
            ##     geefit <- geese.fit(xmatZ, yhatZ, id, zsca = model.matrix(~factor(margin) -1 ), corstr = corstr)
            ## }
            ## else { # (res == TRUE) {
            er1 <- NULL
            er2 <- NULL
            for (m in unique(margin)) {
                temp <- eRes(e[margin == m], delta[margin == m], Z[margin == m])
                temp[[2]] <- ifelse(delta[margin == m] == 1, e[margin == m]^2, temp[[2]])
                eres2[m] <- mean(temp[[2]], na.rm = TRUE)
                dum <- cumsum(ifelse(margin == m, 1, 0))
                er1temp <- temp[[1]][ifelse(margin == m, dum, NA)]
                er1 <- rbind(er1, er1temp)
            }
            er1 <- as.vector(er1)
            er1 <- er1[!is.na(er1)]
            yhat <- delta * y + (1 - delta) * (er1 + xmat %*% beta)
            yhatZ <- sqrt(Z * weights) * yhat
            xmatZ <- sqrt(Z * weights) * xmat
            er2 <- as.matrix(eres2[margin])
            ## geefit <- geese.fit(xmat, yhat, id, zsca = er2, scale.fix = TRUE, corstr = corstr, weights = Z * weights)
            geefit <- geese.fit(xmatZ, yhatZ, id, zsca = er2, scale.fix = TRUE, corstr = corstr)
        }
        beta <- geefit$beta
        if (control$trace) {
            cat("\n beta:\n")
            ## cat(beta)
            print(as.numeric(beta))
        }
        convStep = i
        if (max(abs(beta - betaprev)/abs(beta)) <= control$reltol) break
    } ## end i for 1:maxiter
    iniBeta <- geefit$beta
    if ("(Intercept)" %in% colnames(x)) {
    ##    beta <- c(eRes(e = y - as.matrix(x[,-1]) %*% geefit$beta[-1], delta = delta)[[3]],
      ##            geefit$beta[-1])
beta <- c(mean(yhat - as.matrix(x[,-1]) %*% geefit$beta[-1]), geefit$beta[-1])
    } else {
        beta <- geefit$beta
    }
    alpha <- geefit$alpha
    gamma <- eres2
    convergence <- ifelse(i == control$maxiter, 1, 0)
    out <- list(beta = beta, alpha = alpha, gamma = gamma, iniBeta = iniBeta,
                convergence = convergence, convStep = convStep)
    return(out)
}

eRes <- function(e, delta, z = rep(1, length(e)))
{
  nobs <- length(e)
  ord <- order(e)
  ei <- e[ord]
  deltai <- delta[ord]
  zi <- z[ord]
  dummy <- 1:nobs
  repeats <- table(ei)
  Shat <- survfit(Surv(ei, deltai) ~ 1, weights = zi)$surv
  Shat <- rep(Shat, repeats)
  edif <- c(diff(ei), 0)  ## diff(ei) gives 1 less terms
  ehat <- rev(cumsum(rev(edif * Shat)))
  inpt <- mean(ehat)
  ehat2 <- rev(cumsum(rev(ei * edif * Shat)))
  ehat <- ehat/Shat + ei    ## +ei because there was a diff() in edif
  ehat2 <- 2 * ehat2/Shat + ei^2
  ehat[is.na(ehat)] = ei[is.na(ehat)]
  ehat2[is.na(ehat2)] = ei[is.na(ehat2)]^2
  ehat2[which(ehat2 < 0)] <- NaN
  eres <- ehat
  eres2 <- ehat2
  eres[dummy[ord]] <- ehat  ## puting it back to the original order
  eres2[dummy[ord]] <- ehat2
  return(list(eres, eres2, inpt))
}
