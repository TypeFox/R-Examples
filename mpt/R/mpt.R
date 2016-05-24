# Apr/12/2016 anova.mpt() now works with stats::print.anova()
#
# Mar/18/2016 better fix using eval(..., as.list(spec$par)); this renders more
#             expressions fittable via EM including base functions
#
# Aug/27/2015 BUG FIX: mpt(..., method = "EM") could fail when a symbol in the
#             model was also used as an object name in the work space; fixed
#             by using evalq() instead of eval(); problem with evalq() is that
#             it hides base functions like sqrt()
#
# Mar/11/2015 add multinomial constant to logLik
#
# Sep/10/2014 new infrastructure, mptspec(), mpt(..., method = "BFGS")
#
# Dec/15/2013 simplify extraction of EM constants (a, b, c)
#
# Jan/24/2013 BUG FIX: typo in vcov.mpt(), hence wrong standard errors
#             (reported by Rainer Alexandrowicz and Bartosz Gula)


## Fit MPT model via maximum likelihood (BFGS or EM)
mpt <- function(spec, data, start = NULL, method = c("BFGS", "EM"),
         treeid = "treeid", freqvar = "freq",
         optimargs =
           if(method == "BFGS") list(control =
             list(reltol = .Machine$double.eps^(1/1.2), maxit = 1000))
           else list())
{
  stopifnot(class(spec) == "mptspec")

  ## Either, 'data' is a dataframe
  if(is.data.frame(data)) {
    y <- data[, freqvar]
    tid <- if(length(treeid) == length(y)) factor(treeid)
           else if(length(treeid) == 1 && treeid %in% names(data))
             factor(data[, treeid])
           else if(length(names(spec$prob)) == length(y))  # read from spec
             factor(gsub("^(.+)\\..*", "\\1", names(spec$prob)))
           else rep(1, length(y))
    data <- matrix(y, nrow=1L)

  ## Or a matrix/vector of frequencies
  } else {
    ## sanity checking and reordering of data
    if(is.null(dim(data))) data <- matrix(data, nrow=1L,
                                          dimnames=list(NULL, names(data)))
    if(!is.null(dnam <- colnames(data)) & !is.null(snam <- names(spec$prob))){
      if(!all(snam == dnam)) warning("variable names do not match")
      # if(!all(snam %in% dnam)) {
      #   warning("variable names do not match")
      # } else {
      #   data <- data[, snam, drop = FALSE]
      # }
    }
    tid <- if(length(treeid) == NCOL(data)) factor(treeid)
           else if(length(names(spec$prob)) == NCOL(data))
             factor(gsub("^(.+)\\..*", "\\1", names(spec$prob)))
           else if(!is.null(colnames(data)))
             factor(gsub("^(.+)\\..*", "\\1", colnames(data))) # before 1st dot
           else rep(1, NCOL(data))
  }
  if(NCOL(data) != length(spec$prob))
    stop("number of response categories and model equations do not match")

  ## for fitting only sums are needed
  y <- colSums(data)
  
  method <- match.arg(method)
  ## determine number of parameters and starting values
  if(is.null(start)) {
    start <- spec$par[is.na(spec$par)]  # FIX ME: is.na still necessary?
    start[] <- if (method == "EM") 0.5 else 0  # completely ad hoc
  } else {
    ## do sanity checking of starting values/names/etc.
    if(is.null(names(start))) names(start) <- names(spec$par[is.na(spec$par)])
    if (method == "BFGS") start <- qlogis(start)  # logit transform
  }

  if (method == "BFGS") {

    ## set up log-likelihood and gradient
    nll <- function(par) -sum(y * log(spec$par2prob(plogis(par))))
    grad <- function(par) {
      yp <- drop(y/spec$par2prob(plogis(par)))
      dp <- spec$par2deriv(plogis(par))$deriv
      -drop(dp %*% yp) * dlogis(par)  # FIX ME: dlogis(par) optional? Ask Z.
    }
  
    optArgs <- list(par=start, fn=nll, gr=grad, method="BFGS")
    optArgs <- c(optArgs, as.list(optimargs))
    opt <- do.call(optim, optArgs)
    # opt <- optim(start, nll, gr = grad, method = "BFGS",
    #   control = list(reltol = .Machine$double.eps^(1/1.2), maxit = 1000))
    coef <- opt$par
    loglik <- -opt$value
    pcat <- spec$par2prob(plogis(coef))
    aa <- bb <- cc <- NULL
 
  # } else if (method == "BFGS") {
  #   opt <- optim(start, nll, gr = grad, method = "BFGS",
  #     control = list(reltol = .Machine$double.eps^(1/1.2), maxit = 1000))

  #   coef <- plogis(opt$par)
  #   vc <- solve(H(coef))
  #   loglik <- -sum(log(spec$par2prob(coef)) * y)

  } else {  # EM

    ## Get constants for EM algorithm
    terms <- sapply(lapply(spec$prob, as.character), strsplit, "\\+")  # "+"
    terms <- lapply(terms, function(x) gsub("[[:space:]]", "", x))

    aa <- bb <- array(NA, c(max(sapply(terms, length)),  # max paths to categ
                            length(terms),               # n categories
                            length(start)))              # n pars
    cc <- matrix(1, dim(aa)[1], dim(aa)[2])

    for(j in 1:dim(aa)[2]){
      for(i in 1:sapply(terms, length)[j]){
        pterms <- strsplit(terms[[j]][i], "\\*")[[1]]
        cc[i, j] <- prod(sapply(parse(text=pterms), eval, as.list(spec$par)),
                         na.rm=TRUE)

        for(s in seq_along(start)){
          tname <- names(start)[s]

          aa[i, j, s] <- sum(grepl(paste0("^", tname, "$"), pterms))
          powix <- grepl(paste0("^", tname, "\\^[0-9]+"), pterms)
          aa[i, j, s] <- sum(aa[i, j, s],
            as.numeric(gsub(paste0("^", tname, "\\^([0-9]+)"), "\\1",
                            pterms)[powix]))

          ## Brackets () are optional
          bb[i, j, s] <- sum(grepl(paste0("^\\(?1-", tname, "\\)?$"), pterms))
          powix <- grepl(paste0("^\\(1-", tname, "\\)\\^[0-9]+"), pterms)
          bb[i, j, s] <- sum(bb[i, j, s],
            as.numeric(gsub(paste0("^\\(1-", tname, "\\)\\^([0-9]+)"), "\\1",
                            pterms)[powix]))
        }
      }
    }
    dimnames(aa)[[3]] <- dimnames(bb)[[3]] <- as.list(names(start))

    ## Call mptEM
    optArgs <- list(theta=start, data=y, a=aa, b=bb, c=cc)
    optArgs <- c(optArgs, as.list(optimargs))
    opt <- do.call(mptEM, optArgs)
    # opt <- mptEM(start, y, aa, bb, cc, ...)
    coef <- opt$theta
    loglik <- opt$loglik
    pcat <- opt$pcat
  }

  snam <- if(!is.null(names(spec$prob))) names(spec$prob)
          else if(!is.null(colnames(data))) colnames(data)
          else paste(tid, unlist(lapply(rle(as.character(tid))$lengths,
                                        seq_len)), sep=".")
  ncat    <- table(tid)
  nobs    <- sum(ncat - 1)
  ntrees  <- length(ncat)
# n       <- setNames(tapply(y, tid, sum)[as.character(tid)], snam)
  nbytree <- tapply(y, tid, sum)
  n       <- setNames(nbytree[as.character(tid)], snam)
  fitted  <- n*pcat
  G2      <- 2*sum(y*log(y/fitted), na.rm=TRUE)
  df      <- nobs - length(coef)
  gof     <- c(G2=G2, df=df, pval = pchisq(G2, df, lower.tail=FALSE))

  rval <- list(
    coefficients = coef,
    loglik = sum(lfactorial(nbytree)) - sum(lfactorial(y)) + loglik,
    nobs = nobs,         # nrow(data),
    # df = length(start),
    fitted = fitted,
    goodness.of.fit = gof,
    ntrees = ntrees,
    n = n,
    y = setNames(y, snam),
    pcat = setNames(pcat, snam),
    treeid = tid,
    a = aa, b = bb, c = cc,
    spec = spec,
    method = method,
    optim = opt
  )
  class(rval) <- "mpt"
  return(rval)
}


## EM algorithm
mptEM <- function(theta, data, a, b, c, maxit = 1000, tolerance = 1e-8,
                  stepsize = 1, verbose = FALSE){
  nbranch <- dim(a)[1]
  pbranch <- matrix(NA_real_, nbranch, length(data))
  loglik0 <- -Inf
  theta1  <- theta

  iter <- 1
  while(iter < maxit){
    if(verbose) print(c(iter, loglik0))

    ## E step
    for(i in seq_len(nbranch))
      for(j in seq_along(data))
        pbranch[i, j] <- c[i,j] * prod(theta^a[i,j,] * (1 - theta)^b[i,j,])
    
    pcat    <- colSums(pbranch, na.rm=TRUE)
    loglik1 <- sum(data*log(pcat))
    if(loglik1 - loglik0 < tolerance) break  # stop if converged
    loglik0 <- loglik1
    m       <- t(data*t(pbranch)/pcat)
    
    ## M step
    for(s in seq_along(theta))
      theta1[s] <-
        sum(a[,,s]*m, na.rm=TRUE)/sum((a[,,s] + b[,,s])*m, na.rm=TRUE)
    theta   <- theta - stepsize*(theta - theta1)
    iter    <- iter + 1
  }
  if(iter >= maxit) warning("iteration maximum has been exceeded")
  out <- list(theta=theta, loglik=loglik1, pcat=pcat, pbranch=pbranch,
              iter=iter)
  out
}


coef.mpt <- function(object, logit = FALSE, ...){
  coef <- object$coefficients
  if (logit) {
    nm <- paste0("logit(", names(coef), ")")
    if (object$method == "EM")
      setNames(qlogis(coef), nm)
    else
      setNames(coef, nm)
  } else {
    if (object$method == "EM")
      coef
    else
      plogis(coef)
  }
}


vcov.mpt <- function(object, logit = FALSE, what = c("vcov", "fisher"),
                     ...){
  what <- match.arg(what)
  coef <- coef(object, logit=logit)

  # Negative Hessian (information) on probability scale.
  # Which one is correct? Should be very similar. Stick with I.obs.
  # %*% is slightly faster than sum( * )
  H <- function(par, y = object$y, spec = object$spec,
                type = c("observed", "estimated", "expected")){
    pp  <- spec$par2prob(par)
    yp  <- drop(y/pp)
    dp  <- spec$par2deriv(par)$deriv
    d2p <- spec$par2deriv(par)$deriv2
    npar <- length(par)
    H <- matrix(NA, npar, npar)
    for (i in seq_len(npar))
      for (j in i:npar)
        switch(EXPR = match.arg(type),
          observed =
          H[i, j] <- yp %*% (dp[i, ]*dp[j, ]/pp - d2p[i, j, ]),
        # H[i, j] <- sum(yp * (dp[i, ]*dp[j, ]/pp - d2p[i, j, ]))

          estimated =
          H[i, j] <- yp %*% (dp[i, ]*dp[j, ]/pp),
        # H[i, j] <- sum(y*dp[i, ]*dp[j, ]/pp^2)
        # H[i, j] <- sum(yp*dp[i, ]*dp[j, ]/pp)
        # H[i, j] <- sum(yp/pp * dp[i, ]*dp[j, ])

          # Only correct for single-tree models; for joint MPT models,
          # calculate per-tree info and add them.
          expected =
          H[i, j] <- sum(y)*sum(dp[i, ]*dp[j, ]/pp)
        )
    H[lower.tri(H)] <- t(H)[lower.tri(H)]
    dimnames(H) <- list(names(par), names(par))
    H
  }

  if (logit) {  # delta method
    # Note that plogis(x)*(1 - plogis(x)) == dlogis(x)
    # Automatic dimnames for tcrossprod(...)
    # Conjecture: the smaller the matrix the more reliable its inverse
    # dhessian <- diag(1/(plogis(coef) * (1 - plogis(coef))), length(coef))
    # dhessian %*% solve(H(plogis(coef))) %*% dhessian
    # solve(tcrossprod(dlogis(coef)) * H(plogis(coef)))  # possibly singular
    if (what == "vcov") 1/tcrossprod(dlogis(coef)) * solve(H(plogis(coef)))
    else                  tcrossprod(dlogis(coef)) *       H(plogis(coef))
  } else {
    if (what == "vcov") solve(H(coef))
    else                      H(coef)
  }
}


## Based on stats::confint.default
confint.mpt <- function(object, parm, level = 0.95, logit = TRUE, ...)
{
  cf <- coef(object, logit=logit)
  pnames <- names(cf)
  if (missing(parm)) 
      parm <- pnames
  else if (is.numeric(parm)) 
      parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100*a, trim=TRUE, scientific=FALSE, digits=3), "%")
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object, logit=logit)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}


# logLik.mpt <- function(object, ...)
#                    structure(object$loglik, df = object$df, class = "logLik")


print.mpt <- function(x, digits = max(3, getOption("digits") - 3),
                      logit=FALSE, ...){
  cat("\nMultinomial processing tree (MPT) models\n\n")
  cat("Parameter estimates:\n")
  print.default(format(coef(x, logit=logit), digits=digits), print.gap=2,
                quote = FALSE)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval, digits=digits), "\n", sep="")
  cat("\n")
  invisible(x)
}


anova.mpt <- function (object, ..., test = c("Chisq", "none")){
  ## Adapted form MASS::anova.polr and stats::anova.glmlist

  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
      stop('anova is not implemented for a single "mpt" object')
  mlist <- list(object, ...)
  nmodels <- length(mlist)
  names(mlist) <- c(deparse(substitute(object)),
                    as.character(substitute(...[]))[2:nmodels])

  if (any(!sapply(mlist, inherits, "mpt")))
      stop('not all objects are of class "mpt"')

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


## Log-likelihood for mpt objects
logLik.mpt <- function(object, ...){
  if(length(list(...)))
      warning("extra arguments discarded")
  p <- length(object$coefficients)
  val <- object$loglik
  attr(val, "df") <- p
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}


## Number of observations
nobs.mpt <- function(object, ...) object$nobs


## Residuals for mpt models
residuals.mpt <- function(object, type=c("deviance", "pearson"), ...){

  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  type <- match.arg(type)
  wts <- object$n
  y <- object$y / wts
  mu <- object$pcat
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


## Diagnostic plot for mpt models
plot.mpt <- function(x, showNames = TRUE,
                     xlab="Predicted response probabilities",
                     ylab="Deviance residuals", ...){
  xres <- resid(x)
  mu   <- x$pcat
  plot(mu, xres, xlab = xlab, ylab = ylab, type="n", ...)
  abline(h = 0, lty = 2)
  if(showNames){
    text(mu, xres, names(xres), cex=0.8)
    panel.smooth(mu, xres, cex=0)
  }else{
    panel.smooth(mu, xres)
  }
}


# ## Covariance matrix for MPT model parameters
# vcov.mpt <- function(object, ..., what = c("vcov", "fisher")){
#   a       <- object$a
#   b       <- object$b
#   y       <- object$y
#   pcat    <- object$pcat
#   pbranch <- object$pbranch
#   theta   <- coef(object)
# 
#   ## as(Theta), bs(Theta)
#   as.t <- bs.t <- numeric(length(theta))
#   for(s in seq_along(theta)){
#     for(j in seq_along(pcat)){
#       as.t[s] <- as.t[s] + y[j]*sum(a[,j,s]*pbranch[,j]/pcat[j], na.rm=TRUE)
#       bs.t[s] <- bs.t[s] + y[j]*sum(b[,j,s]*pbranch[,j]/pcat[j], na.rm=TRUE)
#     }
#   }
#   
#   ## d as(Theta)/d t, d bs(Theta)/d t
#   das.t <- dbs.t <- matrix(0, length(theta), length(theta))
#   for(s in seq_along(theta)){
#     for(r in seq_along(theta)){
#       for(j in seq_along(pcat)){
#         das.t[s, r] <- das.t[s, r] + y[j] * (
#         sum(a[,j,s] * pbranch[,j] *
#           sum((a[,j,r]/theta[r] - b[,j,r]/(1 - theta[r])) * pbranch[,j],
#             na.rm = TRUE) /
#           pcat[j]^2, na.rm = TRUE) -
#         sum(a[,j,s] *
#           (a[,j,r]/theta[r] - b[,j,r]/(1 - theta[r])) * pbranch[,j] / pcat[j],
#           na.rm = TRUE)
#         )
#   
#         dbs.t[s, r] <- dbs.t[s, r] + y[j] * (
#         sum(b[,j,s] * pbranch[,j] *
#           sum((a[,j,r]/theta[r] - b[,j,r]/(1 - theta[r])) * pbranch[,j],
#             na.rm = TRUE) /
#           pcat[j]^2, na.rm = TRUE) -
#         sum(b[,j,s] *
#           (a[,j,r]/theta[r] - b[,j,r]/(1 - theta[r])) * pbranch[,j] / pcat[j],
#           na.rm = TRUE)
#         )
#       }
#     }
#   }
#   
#   ## I(Theta)
#   info.t <- das.t/theta - dbs.t/(1 - theta) +
#             diag(as.t/theta^2 + bs.t/(1 - theta)^2)
#   dimnames(info.t) <- list(names(theta), names(theta))
#   what <- match.arg(what)
#   if (what == "vcov") solve(info.t) else info.t
# }


summary.mpt <- function(object, ...){
  x <- object
  coef <- coef(x, logit=TRUE)
  pcoef <- coef(x, logit=FALSE)

  ## Catch vcov error, so there are at least some NA's in the summary
  s.err <- tryCatch(sqrt(diag(vcov(x, logit=TRUE))),
                    error = function(e) rep(NA, length(coef)))

  tvalue <- coef / s.err
  pvalue <- 2 * pnorm(-abs(tvalue))
  dn <- c("Estimate", "Logit Estim.", "Std. Error")
  coef.table <- cbind(pcoef, coef, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(pcoef), c(dn, "z value", "Pr(>|z|)"))

  aic <- AIC(x)
  ans <- list(ntrees=x$ntrees, coefficients=coef.table, aic=aic,
              gof=x$goodness.of.fit, X2=sum(resid(x, "pearson")^2))
  class(ans) <- "summary.mpt"
  return(ans)
}


print.summary.mpt <- function(x, digits = max(3, getOption("digits") - 3),
                              cs.ind = 2:3, ...){
  cat("\nCoefficients:\n")
  printCoefmat(x$coef, digits=digits, cs.ind=cs.ind, ...)
  # cat("\nGoodness of fit:\n")
  cat("\nLikelihood ratio G2:", format(x$gof[1], digits=digits), "on",
    x$gof[2], "df,", "p-value:", format(x$gof[3], digits=digits), "\n")
  cat("Pearson X2: ", format(x$X2, digits=digits), ",    ",
      "AIC: ", format(x$aic, digits=max(4, digits + 1)), sep="")
  cat("\n")
  cat("Number of trees:", x$ntrees, "\n")
  invisible(x)
}


## Simulate responses from mpt model
simulate.mpt <- function(object, nsim, seed, pool = TRUE, ...){

  if(pool){
    tid  <- object$treeid
    freq <- unlist( lapply(unique(tid),
      function(i) rmultinom(1, object$n[tid == i], object$pcat[tid == i])) )
    names(freq) <- tid
  }else{
    stop("individual response simulation not yet implemented")
  }
  setNames(freq, names(object$fitted))
}


deviance.mpt <- function(object, ...) object$goodness.of.fit["G2"]


predict.mpt <- function(object, newdata = NULL, type = c("freq", "prob"),
                        ...){
  type <- match.arg(type)
  if(type == "prob") object$pcat
  else
    if(is.null(newdata)) fitted(object)
    else {
      stopifnot(length(newdata) == length(object$pcat))
      tid <- object$treeid
      object$pcat * setNames(tapply(newdata, tid, sum)[as.character(tid)],
                             names(object$pcat))
    }
}

