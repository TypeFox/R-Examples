betabin <-
  function(data, start = c(.5,.5),
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"),
           vcov = TRUE, corrected = TRUE, gradTol = 1e-4, ...)
{
  m <- match.call(expand.dots = FALSE)
  call <- match.call()
  m$method <- NULL
  m[[1]] <- as.name("list")
  m <- eval.parent(m)
  doFit <- TRUE # A little trick:
  if("doFit" %in% names(m$...)) doFit <- (m$...)$doFit
  if(is.data.frame(m$data)) m$data <- as.matrix(m$data)
  if(!is.matrix(m$data))
    stop("'data' is not a matrix or data.frame")
  if(NCOL(m$data) != 2 || NROW(m$data) < 3)
    stop("'data' should have 2 columns and > 3 rows")
  method <- match.arg(method)
  pGuess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  name <- c("mu", "gamma")
  if(any(start < 1e-3) || any(start > 1- 1e-3))
    stop("start has to be in the open interval (0, 1)")
  bbRho <- bbEnvir(parent.frame(), X=as.matrix(m$data),
                   corrected = corrected, pGuess = pGuess,
                   start = start)
  if(!doFit) return(bbRho)
  ## optimize log-likelihood:
  fit <- optim(getParBB(bbRho), fn = function(par) setParBB(bbRho, par),
               method = "L-BFGS-B", hessian = FALSE,
               lower = bbRho$lbounds, upper = bbRho$ubounds,
               control = list(parscale = c(.01, .01)))
  ## test for adequate convergence:
  if(all(fit$par > bbRho$lbounds) && all(fit$par < bbRho$ubounds)) {
    grad <- grad(function(par) setParBB(bbRho, par), fit$par)
    if(max(abs(grad)) > gradTol)
      warning(sprintf("Optimizer terminated with max|gradient|: %e",
                      max(abs(grad))), call. = FALSE) }
  else
    warning("Parameters at boundary occurred", call. = FALSE)
  ## collect results:
  coef <- fit$par
  names(coef) <- name
  res <- list(coefficients = coef,
              convergence = fit$convergence,
              message = fit$message, counts = fit$counts,
              call = match.call(), data = m$data,
              method = method, corrected = corrected)
  ## make likelihood ratio tests:
  res$logLik <- -fit$value + bbRho$Factor
  res$logLikNull <-
        sum(dbinom(bbRho$x, bbRho$n, prob = pGuess, log = TRUE))
  res$logLikMu <-
    with(bbRho, sum(dbinom(x, n, prob = sum(x)/sum(n),
                           log = TRUE)))
  ## compute variance-covariance matrix of the parameters:
  if(vcov) {
    res$vcov <- matrix(NA, 2, 2)
    names(res$vcov) <- list(name, name)
    if(all(fit$par > bbRho$lbounds) && all(fit$par < bbRho$ubounds))
      res$vcov <-
        solve(hessian(function(par) setParBB(bbRho, par), coef))
  }
  class(res) <- c("betabin")
  return(res)
}

bbEnvir <- function(parent, X, corrected, pGuess, start)
### Constructs an environment, rho for the computation of
### beta-binomial and chance-corrected beta-binomial models.
### arg: parent: parent environment
###      X: data matrix with at least 2 col and 4 rows
###      corrected: boolean, fit the chance-corrected model?
###      pGuess: the guessing probability; used for the
###              chance-corrected model
###      start: starting values; two element vector.
{
    rho <- new.env(parent = parent)
    rho$corrected <- corrected
    rho$X <- X
    rho$n <- X[,2]
    rho$x <- X[,1]
    rho$N <- nrow(X)
    rho$start <- rho$par <- start
    rho$lbounds <- c(1e-6, 1e-6)
    if(corrected) rho$lbounds <- c(1e-6, 1e-3)
    rho$ubounds <- c(1 - 1e-6, 1 - 1e-6)
    rho$nllAux <- function(x, a, b)
        log(sum(choose(x[1], 0:x[1]) * (1 - pGuess)^(0:x[1]) *
                pGuess^(-(0:x[1])) *
                beta(a + 0:x[1], b + x[2] - x[1])))
    ## The term choose(x[1], 0:x[1]) does not have to be computed
    ## every time.
    rho$Factor <- sum(lchoose(rho$n, rho$x))
    if(corrected)
        rho$Factor <- rho$Factor +
            sum((rho$n-rho$x) * log(1-pGuess)) + sum(rho$x * log(pGuess))
    rho
}

getParBB <- function(rho) rho$par
setParBB <- function(rho, par)
### Set parameters of the model and return the negative log-likelihood
### at those parameters.
{
  if(!missing(par))
    rho$par <- par
  with(rho, {
    a <- par[1] * (1 - par[2]) / par[2]
    b <- (1 - par[2]) * (1 - par[1]) / par[2] ## }
    if(corrected)
      N * lbeta(a, b) - sum(apply(X, 1, nllAux, a=a, b=b))
    else
      N * lbeta(a, b) - sum(lbeta(a + x, b - x + n))
  })
}

print.betabin <-
  function(x, digits = max(3, getOption("digits") - 3), ... )
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.betabin <-
  function(object, level = 0.95, ...)
{
  if(is.null(object$vcov))
    stop("summary needs vcov in object")
  ## Construct output table:
  table <- array(NA, dim = c(5, 4))
  rownames(table) <- c("mu", "gamma", "pc", "pd", "d-prime")
  colnames(table) <- c("Estimate", "Std. Error", "Lower", "Upper")
  table[1:2,1] <- coef(object)
  mu <- table[1,1]
  table[1:2,2] <- se <- sqrt(diag(object$vcov))
  object$level <- level
  a <- (1 - level)/2
  ci <- coef(object) + se %o% qnorm(c(a, 1 - a))
  ci[ci < 0] <- 0
  ci[ci > 1] <- 1
  table[1:2, 3:4] <- ci
  muCI <- table[1, 3:4]
  ## Transform estimate to pc, pd and d-prime scales:
  if(!is.na(mu)) {
    if(object$corrected) obj <- rescale(pd = mu, method = object$method)
    else obj <- rescale(pc = mu, method = object$method)
    table[3:5,1] <- unlist(coef(obj))
  }
  ## Transform std.err:
  if(all(!is.na(c(mu, se[1])))) {
    if(object$corrected)
      obj <- rescale(pd = mu, std.err = se[1], method = object$method)
    else
      obj <- rescale(pc = mu, std.err = se[1], method = object$method)
    table[3:5,2] <- unlist(obj$std.err)
  }
  ## Transform Wald CI on mu-scale to pc, pd and d-prime scales:
  if(all(!is.na(muCI))) {
    if(object$corrected)
      intervals <- rescale(pd = muCI, method = object$method)
    else
      intervals <- rescale(pc = muCI, method = object$method)
    table[3:5,3] <- unlist(coef(intervals)[1,])
    table[3:5,4] <- unlist(coef(intervals)[2,])
  }
  object$LR.OD <- 2 * with(object, logLik - logLikMu)
  object$p.value.OD <-
    pchisq(object$LR.OD, df = 1, lower.tail = FALSE)
  object$LR.null <- 2 * with(object, logLik - logLikNull)
  object$p.value.null <-
    pchisq(object$LR.null, df = 2, lower.tail = FALSE)
  object$coefficients <- table
  class(object) <- "summary.betabin"
  object
}

print.summary.betabin <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(x$corrected)
    cat(paste("\nChance-corrected beta-binomial model for the ", x$method,
          " protocol\nwith ", round(100 * x$level, 3),
          " percent confidence intervals\n\n", sep = ""))
  else
    cat(paste("\nBeta-binomial model for the ", x$method,
          " protocol\nwith ", round(100 * x$level, 3),
          " percent confidence intervals\n\n", sep = ""))
  printCoefmat(x$coefficients, tst.ind=integer(0), cs.ind=1:4,
               digits=digits, P.values=FALSE, has.Pvalue=FALSE, ...)
  cat("\nlog-likelihood: ", round(x$logLik, digits), "\n")
  cat("LR-test of over-dispersion, G^2:", round(x$LR.OD, digits),
      "df:", 1, "p-value:", format.pval(x$p.value.OD, digits=digits),
      "\n")
  cat("LR-test of association, G^2:", round(x$LR.null, digits), "df:",
      2, "p-value:", format.pval(x$p.value.null, digits=digits), "\n")
  invisible(x)
}

vcov.betabin <- function(object, ...) {
  object$vcov
}

logLik.betabin <- function(object, ...) {
  val <- object$logLik
  names(val) <- NULL
  attr(val, "nobs") <- sum(object$data[,2])
  attr(val, "df") <- 2
  class(val) <- "logLik"
  val
}

## Potential functions:
## profile.betabin <- function() {
##
## }
##
## plot.profile.betabin <- function() {
##
## }
##
## confint.betabin <- function() {
##
## }

