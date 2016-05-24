# 2004/MAY/04: removing cov2cor (it's there anyway)
#              changing 'print.coefmat' to 'printCoefmat'
#              'print.matrix' doesn't work anymore (for strans)
#
# 2004/JUL/20: trying to make it fit for CRAN:
#   removing beta2eba.R (now in linear2btl.R)
#   making summary.eba generic: add '...', rename x to object
#   changing all 'T' to 'TRUE'
#
# 2005/APR/26: added new functions:
#   pcX            paired-comparison design matrix
#   wald.test      testing linear hypotheses (Cp = 0) in EBA models
#   group.test     groupwise testing in EBA models
#   residuals.eba  deviance and Pearson residuals
#   plot.eba       diagnostic plot
#   BUG FIX: df in the imbalance test corrected
#
# 2008/MAR/28: added new functions:
#   anova.eba      anova function for eba models
#   vcov.eba       vcov function for eba models
#
# 2009/MAY/24: moved functions strans and print.strans to extra file
#
# 2010/JAN/07: new function simulate.eba
#
# 2010/MAY/10: replace fdHess() by nlme::fdHess()
#
# 2011/FEB/28: new zap.ind, tst.ind in printCoefmat in print.summary.eba and
#              in print.group.test
#
# 2011/MAR/01:
#   removed dependencies on $estimate and $se components in objects of class
#     eba
#   used seq_along and seq_len where possible
#   uscale: utility scale extractor function
#   cov.u gains norm argument for normalized scale values


OptiPt <- function(M, A = 1:I, s = rep(1/J, J), constrained = TRUE){
  # parameter estimation for BTL/Pretree/EBA models
  # M: paired-comparison matrix
  # A: model specification list(c(1,6), c(2,6), c(3,7),...)
  # s: starting vector (optional)
  # constrained: constrain parameters to be positive
  # author: Florian Wickelmaier (wickelmaier@web.de)
  #
  # Reference: Wickelmaier, F. & Schmid, C. (2004). A Matlab function
  #   to estimate choice-model parameters from paired-comparison data.
  #   Behavior Research Methods, Instruments, and Computers, 36, 29--40.

  I <- ncol(M)         # number of alternatives/stimuli
  J <- max(unlist(A))  # number of eba parameters

  idx1 <- idx0 <- matrix(0, I*(I - 1)/2, J)  # index matrices
  rdx  <- 1
  for(i in seq_len(I - 1)){           # for(i in 1:(I-1)){
    for(j in (i + 1):I){
      idx1[rdx, setdiff(A[[i]], A[[j]])] <- 1
      idx0[rdx, setdiff(A[[j]], A[[i]])] <- 1
      rdx <- rdx + 1
    }
  }

  y1 <- t(M)[lower.tri(t(M))]  # response vectors
  y0 <- M[lower.tri(M)]
  n  <- y1 + y0
  names(y1) <- names(y0) <- names(n) <- NULL
  logL.sat  <- sum(dbinom(y1, n, y1/n, log=TRUE))  # logLik of the sat. model

  if(constrained){  # minimization
    out <- nlm(L.constrained, s, y1=y1, m=n, i1=idx1, i0=idx0)  # constrained
  }else{
    out <- nlm(L, s, y1=y1, m=n, i1=idx1, i0=idx0)            # unconstrained
  }

  p <- out$estimate  # optimized parameters
  names(p) <- 1:J
  hes  <- nlme::fdHess(p, L, y1, n, idx1, idx0)$Hessian  # numerical Hessian
  cova <- solve(rbind(cbind(hes, 1), c(rep(1, J), 0)))[1:J, 1:J]
  dimnames(cova) <- list(names(p), names(p))
  logL.eba <- -out$min                     # likelihood of the specified model

  fitted <- matrix(0, I, I)                        # fitted PCM
  fitted[lower.tri(fitted)] <- n/(1 + idx0%*%p / idx1%*%p)
  mu <- as.numeric(1/(1 + idx0%*%p / idx1%*%p))    # predicted probabilities
  fitted <- t(fitted)
  fitted[lower.tri(fitted)] <- n/(1 + idx1%*%p / idx0%*%p)
  dimnames(fitted) <- dimnames(M)

  G2   <- 2 * (logL.sat - logL.eba)  # G2 goodness-of-fit statistic
  df   <- I*(I - 1)/2 - (J - 1)
  pval <- 1 - pchisq(G2, df)
  gof  <- c(G2, df, pval)
  names(gof) <- c("-2logL", "df", "pval")
  X2   <- sum((M - fitted)^2 / fitted, na.rm=TRUE)

  u <- numeric()  # scale values
  for(i in seq_len(I)) u <- c(u, sum(p[A[[i]]]))
  names(u) <- colnames(M)

  z <- list(coefficients=p, estimate=p, fitted=fitted,
            logL.eba=logL.eba, logL.sat=logL.sat, goodness.of.fit=gof,
            u.scale=u, hessian=-hes, cov.p=cova, chi.alt=X2, A=A, y1=y1,
            y0=y0, n=n, mu=mu, nobs=length(n))
  class(z) <- "eba"
  z
}


eba <- OptiPt  # wrapper for OptiPt


summary.eba <- function(object, ...){
  x      <- object
  I      <- length(x$A)
  J      <- length(coef(x))      # length(x$estimate)
  y      <- c(x$y1, x$y0)
  coef   <- coef(x)              # x$estimate
  s.err  <- sqrt(diag(vcov(x)))  # x$se
  tvalue <- coef / s.err
  pvalue <- 2 * pnorm(-abs(tvalue))
  dn     <- c("Estimate", "Std. Error")
  coef.table <- cbind(coef, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef), c(dn, "z value", "Pr(>|z|)"))

  tests <- rbind(
   # mean poisson model vs. saturated poisson model (on y)
    c(df1 <- 1,
      df2 <- I*(I - 1),
      l.1 <- sum(dpois(y, mean(y), log=TRUE)),
      l.2 <- sum(dpois(y, y, log=TRUE)),
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # EBA model vs. saturated binomial model
    c(df1 <- J - 1,
      df2 <- I*(I - 1)/2,
      l.1 <- x$logL.eba,
      l.2 <- x$logL.sat,
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # Null model vs. EBA model
    c(df1 <- 0,
      df2 <- J - 1,
      l.1 <- sum(dbinom(x$y1, x$n, 1/2, log=TRUE)),
      l.2 <- x$logL.eba,
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # mean poisson model vs. saturated poisson model (on n)
    c(df1 <- 1,
      df2 <- I*(I - 1)/2,
      l.1 <- sum(dpois(x$n, mean(x$n), log=TRUE)),
      l.2 <- sum(dpois(x$n, x$n, log=TRUE)),
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1))
  )
  rownames(tests) <- c("Overall", "EBA", "Effect", "Imbalance")
  colnames(tests) <- c("Df1","Df2","logLik1","logLik2","Deviance","Pr(>Chi)")

  aic <- -2*x$logL.eba + 2*(length(coef)-1)
  ans <- list(coefficients=coef.table, aic=aic, logL.eba=x$logL.eba,
    logL.sat=x$logL.sat, tests=tests, chi.alt=x$chi.alt)
  class(ans) <- "summary.eba"
  return(ans)
}


print.eba <- function(x, digits=max(3, getOption("digits")-3),
  na.print="", ...){
  cat("\nElimination by aspects (EBA) models\n\n")
  cat("Parameter estimates:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
      quote = FALSE)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (-2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval,digits=digits), "\n", sep="")
  cat("\n")
  invisible(x)
}


print.summary.eba <- function(x, digits=max(3, getOption("digits") - 3),
  na.print="", signif.stars=getOption("show.signif.stars"), ...){
  cat("\nParameter estimates:\n")
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
  cat("\nModel tests:\n")
  printCoefmat(x$tests, digits = digits, signif.stars = signif.stars,
    zap.ind = 1:2, tst.ind = 3:5, ...)
  cat("\nAIC: ", format(x$aic, digits=max(4, digits + 1)), "\n")
  cat("Pearson X2:", format(x$chi.alt, digits=digits))
  cat("\n")
  invisible(x)
}


L <- function(p, y1, m, i1, i0)
  -sum(dbinom(y1, m, 1/(1 + i0 %*% p/i1 %*% p), log=TRUE))


L.constrained <- function(p, y1, m, i1, i0)  # constrain search space
  ifelse(all(p > 0), -sum(dbinom(y1, m, 1/(1+i0%*%p/i1%*%p), log=TRUE)), 1e20)


uscale <- function(object, norm = "sum", log = FALSE){
  ## Extract utility scale from eba object
  
  uscale <- object$u.scale

  if(is.null(norm))
    uscale <- uscale
  else if(norm == "sum")
    uscale <- uscale/sum(uscale)
  else if(as.numeric(norm) %in% seq_along(uscale))
    uscale <- uscale/uscale[as.numeric(norm)]
  else
    stop(sprintf(
         "normalization has to be 'sum' or a number form 1 to %i or 'NULL'",
         length(uscale)))

  if(log) log(uscale) else uscale
}


cov.u <- function(object, norm = "sum"){
  ## Covariance matrix of the utility scale

  x     <- object
  A     <- x$A
  cov.p <- x$cov.p
  cov.u <- matrix(0, length(A), length(A))
  for(i in seq_along(A)){                  # for(i in 1:length(A)){
    for(j in seq_along(A)){                #   for(j in 1:length(A)){
      cell <- 0
      for(k in seq_along(A[[i]]))          #     for(k in 1:length(A[[i]]))
        for(l in seq_along(A[[j]]))        #       for(l in 1:length(A[[j]]))
          cell <- cell + cov.p[A[[i]][k], A[[j]][l]]
      cov.u[i, j] <- cell
    }
  }
  colnames(cov.u) <- rownames(cov.u) <- names(x$u.scale)

  ## Normalization
  if(is.null(norm))
    cov.u
  else if(norm == "sum")
    cov.u/sum(x$u.scale)^2
  else if(as.numeric(norm) %in% seq_len(nrow(cov.u)))
    cov.u/x$u.scale[as.numeric(norm)]^2
  else
    stop(sprintf(
         "normalization has to be 'sum' or a number form 1 to %i or 'NULL'",
         length(uscale)))
}


pcX <- function(nstimuli, omitRef=TRUE){
  ## Paired comparison design matrix

  X <- matrix(0, choose(nstimuli, 2), nstimuli)
  count <- 1
  for(i in seq_len(nstimuli - 1)){     # for(i in 1:(nstimuli - 1)){
    for(j in (i + 1):nstimuli){
      X[count, i] <- 1
      X[count, j] <- -1
      count <- count + 1
    }
  }
  if(omitRef) X[,-1]
  else X
}


group.test <- function(groups, A=1:I, s=rep(1/J,J), constrained=TRUE){
  # groups: 3d array of group matrices (one matrix per group)
  # BUG FIX: combinatorial constant is added to the pooled models!

  pool <- apply(groups, 1:2, sum)  # pooled data matrix
  I <- ncol(pool)                  # number of stimuli
  J <- max(unlist(A))              # number of eba parameters
  ngroups <- dim(groups)[3]        # number of groups

  eba.p <- OptiPt(pool, A, s, constrained)  # EBA for pooled data
  ebas <- NULL  # list of eba models (one per group)
  for(i in seq_len(ngroups))
    ebas[[i]] <- OptiPt(groups[,,i], A, s, constrained)

  C1 <- sum(log(choose(eba.p$n, eba.p$y1)))
  C2 <- 0
  for(i in seq_len(ngroups))
    C2 <- C2 + sum(log(choose(ebas[[i]]$n, ebas[[i]]$y1)))
  C <- C2 - C1  # combinatorial constant

  logL.eba.group <- 0
  y <- n <- NULL
  for(i in seq_len(ngroups)){          # for(i in 1:ngroups){
    logL.eba.group <- logL.eba.group + ebas[[i]]$logL.eba
    y <- c(y, c(ebas[[i]]$y1, ebas[[i]]$y0))
    n <- c(n, ebas[[i]]$n)
  }
  logL.sat.group <- 0
  for(i in seq_len(ngroups))
    logL.sat.group <- logL.sat.group + ebas[[i]]$logL.sat

  tests <- rbind(
   # mean poisson model vs. saturated poisson group model (on y)
    c(df1 <- 1,
      df2 <- ngroups * I * (I - 1),
      l.1 <- sum(dpois(y, mean(y), log=TRUE)),
      l.2 <- sum(dpois(y, y, log=TRUE)),
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # EBA group model vs. saturated binomial group model
    c(df1 <- ngroups * (J - 1),
      df2 <- ngroups * I * (I - 1)/2,
      l.1 <- logL.eba.group,
      l.2 <- logL.sat.group,
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # EBA pool model vs. EBA group model
    c(df1 <- J - 1,
      df2 <- ngroups * (J - 1),
      l.1 <- eba.p$logL.eba + C,  # add constant to obtain correct logLik
      l.2 <- logL.eba.group,
      dev <- 2 * (l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # Null model vs. EBA pool model
    c(df1 <- 0,
      df2 <- J - 1,
      l.1 <- sum(dbinom(eba.p$y1, eba.p$n, 1/2, log=TRUE)) + C,  # add C
      l.2 <- eba.p$logL.eba + C,                                 # add C
      dev <- 2 * (l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # mean poisson model vs. saturated poisson group model (on n)
    c(df1 <- 1,
      df2 <- ngroups * (I * (I - 1)/2),
      l.1 <- sum(dpois(n, mean(n), log=TRUE)),
      l.2 <- sum(dpois(n, n, log=TRUE)),
      dev <- 2 * (l.2 - l.1), 1 - pchisq(dev, df2 - df1))
  )
  rownames(tests) <- c("Overall", "EBA.g", "Group", "Effect", "Imbalance")
  colnames(tests) <- c("Df1","Df2","logLik1","logLik2","Deviance","Pr(>Chi)")

  z <- list(tests=tests)
  class(z) <- "group.test"
  z
}


print.group.test <- function(x, digits=max(3,getOption("digits")-3),
  na.print="", signif.stars=getOption("show.signif.stars"), ...){
  cat("\nTesting for group effects in EBA models:\n")
  cat("\n")
  printCoefmat(x$tests, digits = digits, signif.stars = signif.stars,
    zap.ind = 1:2, tst.ind = 3:5, ...)
  cat("\n")
  invisible(x)
}


wald.test <- function(object, C, u.scale=TRUE){
  # Wald test of linear hypothesis Cp = 0 for EBA models
  # u.scale=TRUE: test on the u.scale values
  # u.scale=FALSE: test on the EBA parameters

  if(u.scale){
    p   <- uscale(object, norm=NULL)   # object$u.scale
    COV <-  cov.u(object, norm=NULL)
  }else{
    p   <- coef(object)                # object$estimate
    COV <- vcov(object)                # object$cov.p
  }
  if(!is.matrix(C)) stop("C is not a matrix")
  if(dim(C)[2] != length(p))
    stop("column number of C and length of p do not agree")
  if(u.scale) colnames(C) <- names(object$u.scale)

  W <- t(C%*%p) %*% solve( C%*%COV%*%t(C) ) %*% (C%*%p)
  z <- list(W=W, df=qr(C)$rank, pval=1 - pchisq(W, qr(C)$rank), C=C)
  class(z) <- "wald.test"
  z
}


print.wald.test <- function(x, digits=max(3,getOption("digits")-4), ...){
  cat("\nWald Test: Cp = 0\n\n")
  cat("C:\n")
  print(x$C)
  cat("\nW = ", x$W, ", df = ", x$df, ", p-value = ", x$pval, "\n", sep='')
  cat("\n")
  invisible(x)
}


residuals.eba <- function(object, type=c("deviance", "pearson"), ...){

  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
              (1-y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))

  type <- match.arg(type)
  wts <- object$n
  y <- object$y1 / wts
  mu <- object$mu
  res <- switch(type,
    deviance = if(object$goodness['df'] > 0){
        d.res <- sqrt(pmax(dev.resids(y, mu, wts), 0))
        ifelse(y > mu, d.res, -d.res)  # sign
      }
      else rep.int(0, length(mu)),
    pearson = (y - mu) * sqrt(wts)/sqrt(mu * (1 - mu))
  )
  if(!is.null(object$na.action)) res <- naresid(object$na.action, res)
  res
}


plot.eba <- function(x, xlab="Predicted choice probabilities",
  ylab="Deviance residuals", ...){
  plot(x$mu, resid(x), xlab = xlab, ylab = ylab, ...)
  abline(h = 0, lty = 2)
  panel.smooth(x$mu, resid(x))
}


anova.eba <- function(object, ..., test=c("Chisq", "none")){
  # Adapted form MASS::anova.polr and stats::anova.glmlist
  # Also works for objects of class eba.order

  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
      stop('anova is not implemented for a single "eba" object')
  mlist <- list(object, ...)
  nmodels <- length(mlist)
  names(mlist) <- sapply(match.call()[-1],
      function(s) paste(deparse(s), collapse="\n"))[seq_len(nmodels)]

  ## All models must be either eba or eba.order
  if (any(!sapply(mlist, inherits, "eba")) &&
      any(!sapply(mlist, inherits, "eba.order")))
      stop('not all objects are of class "eba" or of class "eba.order"')

  ns <- sapply(mlist, function(x) length(x$mu))
  if (any(ns != ns[1]))
      stop("models were not all fitted to the same size of dataset")

  dfs <- sapply(mlist, function(x) x$goodness.of.fit["df"])
  lls <- sapply(mlist, function(x) x$goodness.of.fit["-2logL"])
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


vcov.eba <- function(object, ...) object$cov.p


## Log-likelihood for eba objects
logLik.eba <- function(object, ...){
    if(length(list(...)))
        warning("extra arguments discarded")
    p <- length(object$estimate) - 1
    val <- object$logL.eba
    attr(val, "df") <- p
    attr(val, "nobs") <- object$nobs
    class(val) <- "logLik"
    val
}


## Number of observations
nobs.eba <- function(object, ...) object$nobs


simulate.eba <- function(object, nsim, seed, pool = TRUE, ...){
  ## Simulate responses from eba model

  if(pool){
    n  <- length(object$mu)
    y1 <- rbinom(n, size=object$n, prob=object$mu)
    y0 <- object$n - y1

    mat <- matrix(0, nrow(object$fitted), ncol(object$fitted))
    mat[lower.tri(mat)] <- y1
    mat <- t(mat)
    mat[lower.tri(mat)] <- y0
    dimnames(mat) <- dimnames(object$fitted)
  }else{
    stop("individual response simulation not yet implemented")
  }
  mat
}


deviance.eba <- function(object, ...) object$goodness.of.fit["-2logL"]

