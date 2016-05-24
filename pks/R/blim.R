## Fitting the basic local independence model (BLIM) by MDML
blim <- function(K, N.R, method = c("MD", "ML", "MDML"), R = as.binmat(N.R),
                 P.K = rep(1/nstates, nstates),
                 beta = rep(0.1, nitems),
                  eta = rep(0.1, nitems),
                 betafix = rep(NA, nitems), etafix = rep(NA, nitems),
                 betaequal = NULL, etaequal = NULL,
                 errtype = c("both", "error", "guessing"),
                 errequal = FALSE, randinit = FALSE, incradius = 0,
                 tol = 1e-07, maxiter = 10000, zeropad = 12) {

  K       <- as.matrix(K)
  N.R     <- setNames(as.integer(N.R), names(N.R))  # convert to named int
  N       <- sum(N.R)
  nitems  <- ncol(K)
  npat    <- nrow(R)
  nstates <- nrow(K)

  ## Uniformly random initial values
  if (randinit) {
      beta <- runif(nitems)                       # constraint: beta + eta < 1
       eta <- runif(nitems)
      beta <- ifelse(beta + eta < 1, beta, 1 - beta)
       eta <- ifelse(beta + eta < 1,  eta, 1 -  eta)
         x <- c(0, sort(runif(nstates - 1)), 1)
       P.K <- x[-1] - x[-length(x)]               # constraint: sum(P.K) == 1
  }

  ## Equality restrictions
  betaeq <- etaeq <- diag(nitems)
  if (!is.null(betaequal)) for (i in betaequal) betaeq[i, i] <- 1
  if (!is.null( etaequal)) for (i in  etaequal)  etaeq[i, i] <- 1

  errtype <- match.arg(errtype)                        # overrides arguments
  if (errtype == "error") {
    etafix <- rep(0, nitems)
    warning("errtype is deprecated, use etafix = rep(0, nitems) instead")
  }
  if (errtype == "guessing") {
    betafix <- rep(0, nitems)
    warning("errtype is deprecated, use betafix = rep(0, nitems) instead")
  }
  if (errequal) {
    betaeq <- etaeq <- matrix(1, nitems, nitems)
    warning("errequal is deprecated, use betaequal or etaequal instead")
  }
  beta[!is.na(betafix)] <- betafix[!is.na(betafix)]    # overrides arguments
   eta[!is.na( etafix)] <-  etafix[!is.na( etafix)]

  names(P.K) <- if(is.null(rownames(K))) as.pattern(K) else rownames(K)
  names(beta) <- names(eta) <-
    if (is.null(colnames(K))) {
      make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
                  sep = "")
    } else colnames(K)
  dimnames(betaeq) <- dimnames(etaeq) <- list(names(eta), names(eta))

  ## Assigning state K given response R
  d.RK <- if (length(which(c(betafix, etafix) == 0)) > 0) {
    apply(K, 1, function(k) {
      RwoK <- t(R) & !k
      idx <- which(RwoK, arr.ind=TRUE)
      RwoK[idx[idx[, "row"] %in% which(etafix == 0), ]] <- NA
    
      KwoR <- k & !t(R)
      idx <- which(KwoR, arr.ind=TRUE)
      KwoR[idx[idx[, "row"] %in% which(betafix == 0), ]] <- NA
      colSums(RwoK) + colSums(KwoR)
    })
  } else
    apply(K, 1, function(k) colSums(xor(t(R), k)))
  d.min <- apply(d.RK, 1, min, na.rm = TRUE)             # minimum discrepancy
  i.RK  <- (d.RK <= (d.min + incradius)) & !is.na(d.RK)

  ## Minimum discrepancy distribution
  disc.tab <- xtabs(N.R ~ d.min)
  disc     <- as.numeric(names(disc.tab)) %*% disc.tab / N

  iter    <- 0
  maxdiff <- 2 * tol
  em      <- switch(method <- match.arg(method), MD = 0, ML = 1, MDML = 1)
  md      <- switch(method, MD = 1, ML = 0, MDML = 1)
  beta.num <- beta.denom <- eta.num <-  eta.denom <- beta
  while ((maxdiff > tol) && (iter < maxiter) &&
         ((md*(1 - em) != 1) || (iter == 0))) {
    pi.old   <- P.K
    beta.old <- beta
    eta.old  <- eta

    P.R.K <- apply(K, 1, function(k) apply(
           beta^((1 - t(R))*k) * (1 - beta)^(t(R)*k) *
            eta^(t(R)*(1 - k)) * (1 - eta)^((1 - t(R))*(1 - k)),
           2, prod))
    P.R    <- as.numeric(P.R.K %*% P.K)
    P.K.R  <- P.R.K * outer(1/P.R, P.K)         # prediction of P(K|R)
    mat.RK <- i.RK^md * P.K.R^em
    m.RK   <- (mat.RK / rowSums(mat.RK)) * N.R  # m.RK = E(M.RK) = P(K|R)*N(R)

    ## Distribution of knowledge states
    P.K <- colSums(m.RK) / N

    ## Careless error and guessing parameters
    for (j in seq_len(nitems)) {
      beta.num[j]   <- sum(m.RK[which(R[, j] == 0), which(K[, j] == 1)])
      beta.denom[j] <- sum(m.RK[, which(K[, j] == 1)])
       eta.num[j]   <- sum(m.RK[which(R[, j] == 1), which(K[, j] == 0)])
       eta.denom[j] <- sum(m.RK[, which(K[, j] == 0)])
    }
    beta <- drop(betaeq %*% beta.num / betaeq %*% beta.denom)
     eta <- drop( etaeq %*%  eta.num /  etaeq %*%  eta.denom)
    beta[is.na(beta)] <- 0
     eta[is.na( eta)] <- 0
    beta[!is.na(betafix)] <- betafix[!is.na(betafix)]  # reset fixed parameters
     eta[!is.na( etafix)] <-  etafix[!is.na( etafix)]

    maxdiff <- max(abs(c(P.K, beta, eta) - c(pi.old, beta.old, eta.old)))
    iter <- iter + 1
  }
  if(iter >= maxiter) warning("iteration maximum has been exceeded")

  ## Mean number of errors
  P.Kq <- numeric(nitems)
  for(j in seq_len(nitems))
    P.Kq[j] <- sum(P.K[which(K[,j] == 1)])
  nerror <- c("careless error" = sum(beta * P.Kq),
                 "lucky guess" = sum( eta * (1 - P.Kq)))

  ## If there are missing response patterns, create complete R and N.R
  if(npat < 2^nitems && nitems <= zeropad) {
    N.Rincomp <- N.R
    R   <- expand.grid(rep(list(0:1), nitems), KEEP.OUT.ATTRS=FALSE)
    N.R <- setNames(integer(nrow(R)), as.pattern(R)) # named int filled w/zeros
    R   <- as.binmat(N.R)                            # named int again
    N.R[names(N.Rincomp)] <- N.Rincomp
  }

  ## Recompute predictions and likelihood
  P.R.K <- apply(K, 1, function(k) apply(
         beta^((1 - t(R))*k) * (1 - beta)^(t(R)*k) *
          eta^(t(R)*(1 - k)) * (1 - eta)^((1 - t(R))*(1 - k)),
         2, prod))
  P.R <- as.numeric(P.R.K %*% P.K)
  if (sum(P.R) < 1) P.R <- P.R/sum(P.R)      # if no zero padding: normalize
  loglik <- sum(log(P.R) * N.R, na.rm=TRUE)

  ## Number of parameters
  npar <- nstates - 1 + qr(betaeq)$rank - sum(!is.na(betafix)) +
                        qr( etaeq)$rank - sum(!is.na( etafix))

  ## Goodness of fit, df = number of patterns or persons
  fitted <- setNames(N*P.R, names(N.R))
  G2     <- 2*sum(N.R*log(N.R/fitted), na.rm=TRUE)
# df     <- min(2^nitems - 1, N) - npar        # number of patterns or persons
  df     <- min(if(nitems <= zeropad) 2^nitems - 1 else npat, N) - npar
  gof    <- c(G2=G2, df=df, pval = 1 - pchisq(G2, df))

  z <- list(discrepancy=c(disc), P.K=P.K, beta=beta, eta=eta,
    disc.tab=disc.tab, K=K, N.R=N.R, nitems=nitems, nstates=nstates,
    npatterns=npat, ntotal=N, nerror=nerror, npar=npar,
    method=method, iter=iter, loglik=loglik, fitted.values=fitted,
    goodness.of.fit=gof)
  class(z) <- "blim"
  z
}


print.blim <- function(x, P.Kshow = FALSE, errshow = TRUE,
  digits=max(3, getOption("digits") - 2), ...){
  cat("\nBasic local independence models (BLIMs)\n")
  cat("\nNumber of knowledge states:", x$nstates)
  cat("\nNumber of response patterns:", x$npatterns)
  cat("\nNumber of respondents:", x$ntotal)

  method <- switch(x$method,
            MD = "Minimum discrepancy",
            ML = "Maximum likelihood",
          MDML = "Minimum discrepancy maximum likelihood")
  cat("\n\nMethod:", method)
  cat("\nNumber of iterations:", x$iter)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval, digits=digits), "\n", sep="")

  cat("\nMinimum discrepancy distribution (mean = ",
    round(x$discrepancy, digits=digits), ")\n", sep="")
  disc.tab <- x$disc.tab
  names(dimnames(disc.tab)) <- NULL
  print(disc.tab)
  cat("\nMean number of errors (total = ",
    round(sum(x$nerror), digits=digits), ")\n", sep="")
  print(x$nerror)
  if(P.Kshow){
    cat("\nDistribution of knowledge states\n")
    printCoefmat(cbind("P(K)"=x$P.K), digits=digits, cs.ind=1, tst.ind=NULL,
      zap.ind=1)
  }
  if(errshow){
    cat("\nError and guessing parameters\n")
    printCoefmat(cbind(beta=x$beta, eta=x$eta), digits=digits, cs.ind=1:2,
      tst.ind=NULL, zap.ind=1:2)
  }
  cat("\n")
  invisible(x)
}


## Log-likelihood for blim objects
logLik.blim <- function(object, ...){
  if(length(list(...)))
    warning("extra arguments discarded")
  p <- object$npar
  val <- object$loglik
  attr(val, "df") <- p
  attr(val, "nobs") <- object$npatterns
  class(val) <- "logLik"
  val
}


## Number of observations
nobs.blim <- function(object, ...) object$npatterns


## Residuals for BLIMs
residuals.blim <- function(object, type=c("deviance", "pearson"), ...){

  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  type <- match.arg(type)
  wts <- object$ntotal
  y <- object$N.R / wts
  mu <- object$fitted.values/wts
  res <- switch(type,
    deviance = if(object$goodness['df'] > 0){
        d.res <- sqrt(pmax(dev.resids(y, mu, wts), 0))
        ifelse(y > mu, d.res, -d.res)  # sign
      }
      else rep.int(0, length(mu)),
    pearson = (y - mu) * sqrt(wts)/sqrt(mu)
  )
  if(!is.null(object$na.action)) res <- naresid(object$na.action, res)
  res
}


## Diagnostic plot for BLIMs
plot.blim <- function(x,
  xlab="Predicted response probabilities", ylab="Deviance residuals", ...){

  xres <- resid(x)
  mu   <- x$fitted.values/x$ntotal
  plot(mu, xres, xlab = xlab, ylab = ylab, type="n", ...)
  abline(h = 0, lty = 2)
  panel.smooth(mu, xres)
}


## Simulate responses from BLIM
simulate.blim <- function(object, nsim = 1, seed = NULL, ...){
     P.K <- object$P.K
    beta <- object$beta
     eta <- object$eta
      tK <- t(as.matrix(object$K))
       N <- object$ntotal
  nitems <- nrow(tK)

  state.id <- sample(seq_along(P.K), N, replace=TRUE, prob=P.K)  # draw states

  P.1.K <- tK*(1 - beta) + (1 - tK)*eta               # P(resp = 1 | K)
  R     <- matrix(0, N, nitems)                       # response matrix
  for(i in seq_len(N))
    R[i,] <- rbinom(nitems, 1, P.1.K[, state.id[i]])  # draw a response

  as.pattern(R, freq = TRUE)
}


## Convert binary matrix to vector of response patterns
as.pattern <- function(R, freq = FALSE, as.letters = FALSE, as.set = FALSE){
  if(freq){
    N.R <- table(apply(R, 1, paste, collapse=""))
    setNames(as.integer(N.R), names(N.R))          # convert to named int
  }else
    if(as.letters | as.set){
      nitems <- ncol(R)
      item.names <- 
       make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
                     sep="")
      lett <- apply(R, 1, function(r) paste(item.names[which(r == 1)],
                    collapse=""))
      lett[lett == ""] <- "0"

      if(as.set){
        # Separate elements in lett by "_", remove leading "_",
        # then strsplit along "_" (trailing "_" are ignored by strsplit)
        setfam <- as.set(lapply(strsplit(
          gsub("^_(.+)", "\\1", gsub("([0-9]*)", "\\1_", unname(lett))),
          "_"), as.set))
        if (set_contains_element(setfam, set("0")))
          setfam[[set("0")]] <- set()  # proper empty set
        setfam  # return family of sets, class set
      }else
        lett    # return letters, class character
    }else
      unname(apply(R, 1, paste, collapse=""))
}


## Convert vector of response patterns to named binary matrix
as.binmat <- function(N.R, uniq = TRUE, col.names = NULL){
  if (is.set(N.R)) {
    states <- lapply(N.R, as.character)
    items <- sort(unique(unlist(states)))
    R <- matrix(0, length(N.R), length(items),
                dimnames=list(NULL,
                              if(is.null(col.names)) items else col.names))
    for (i in seq_len(nrow(R))) R[i, states[[i]]] <- 1
  } else {
    pat <- if(is.null(names(N.R))) N.R else names(N.R)
    R   <- if(uniq) strsplit(pat, "") else strsplit(rep(pat, N.R), "")
    R   <- do.call(rbind, R)

    colnames(R) <- 
      if(is.null(col.names)){
        nitems <- ncol(R)
        make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
          sep="")
      }else
        col.names
  }
  storage.mode(R) <- "integer"
  R
}


anova.blim <- function(object, ..., test = c("Chisq", "none")){
  ## Adapted form MASS::anova.polr and stats::anova.glmlist

  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
      stop('anova is not implemented for a single "blim" object')
  mlist <- list(object, ...)
  nmodels <- length(mlist)
  names(mlist) <- sapply(match.call()[-1],
      function(s) paste(deparse(s), collapse="\n"))[seq_len(nmodels)]

  if (any(!sapply(mlist, inherits, "blim")))
      stop('not all objects are of class "blim"')

  ns <- sapply(mlist, function(x) length(x$fitted))
  if (any(ns != ns[1]))
      stop("models were not all fitted to the same size of dataset")

  dfs <- sapply(mlist, function(x) x$goodness.of.fit["df"])
  lls <- sapply(mlist, function(x) x$goodness.of.fit["G2"])
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, pchisq(x2[-1], df[-1], lower.tail = FALSE))

  out <- data.frame(Resid.df = dfs, Deviance = lls, Df = df, Chisq = x2,
                    Prob = pr)
  dimnames(out) <- list(1:nmodels, c("Resid. Df", "Resid. Dev", "Df",
                                     "Deviance", "Pr(>Chi)"))
  if (test == "none") out <- out[, -ncol(out)]

  structure(out,
            heading = c("Analysis of Deviance Table\n",
                        paste0("Model ", format(1L:nmodels), ": ",
                               names(mlist), collapse = "\n")),
            class = c("anova", "data.frame"))
}


deviance.blim <- function(object, ...) object$goodness.of.fit["G2"]


coef.blim <- function(object, ...){
  c(setNames(object$beta, paste("beta", names(object$beta), sep=".")),
    setNames(object$eta,  paste( "eta", names(object$eta),  sep=".")),
    object$P.K)
}

