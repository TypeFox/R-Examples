betareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    link.phi = NULL, type = c("ML", "BC", "BR"),
                    control = betareg.control(...),
                    model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(!(min(Y) > 0 & max(Y) < 1)) stop("invalid dependent variable, all observations must be in (0, 1)")

  ## convenience variables
  n <- length(Y)

  ## type of estimator
  type <- match.arg(type)

  ## links
  if(is.character(link)) link <- match.arg(link)
  if(is.null(link.phi)) link.phi <- if(simple_formula) "identity" else "log"
  if(is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## in mean part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## in precision part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for mean)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(mean = offsetX, precision = offsetZ)

  ## call the actual workhorse: betareg.fit()
  rval <- betareg.fit(X, Y, Z, weights, offset, link, link.phi, type, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, precision = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, precision = Z)

  class(rval) <- "betareg"
  return(rval)
}

betareg.control <- function(phi = TRUE,
  method = "BFGS", maxit = 5000, hessian = FALSE, trace = FALSE, start = NULL,
  fsmaxit = 200, fstol = 1e-8, ...)
{
  rval <- list(phi = phi, method = method, maxit = maxit, hessian = hessian, trace = trace, start = start,
    fsmaxit = fsmaxit, fstol = fstol)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}

betareg.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
  link = "logit", link.phi = "log", type = "ML", control = betareg.control())
{
  ## response and regressor matrix
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(mean = offset, precision = rep.int(0, n))
  if(is.null(z)) {
    m <- 1L
    z <- matrix(1, ncol = m, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
    phi_const <- TRUE
  } else {
    m <- NCOL(z)
    if(m < 1L) stop("dispersion regression needs to have at least one parameter")
    phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[, 1L]), rep.int(1, n)))
  }

  ## link processing
  if(is.character(link)) {
    linkstr <- link
    if(linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      ## add dmu.deta potentially needed for BC/BR
      linkobj$dmu.deta <- make.dmu.deta(linkstr)
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) -log(-log(mu)),
        linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
        mu.eta = function(eta) {
          eta <- pmin(eta, 700)
          pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
        },
        dmu.deta = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
        valideta = function(eta) TRUE,
        name = "loglog"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link
    linkstr <- link$name
    if(type != "ML" && is.null(linkobj$dmu.deta)) warning("link needs to provide dmu.deta component for BC/BR")
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta
  if(is.character(link.phi)) {
    phi_linkstr <- link.phi
    phi_linkobj <- make.link(phi_linkstr)
    phi_linkobj$dmu.deta <- make.dmu.deta(phi_linkstr)
  } else {
    phi_linkobj <- link.phi
    phi_linkstr <- link.phi$name
    if(type != "ML" && is.null(phi_linkobj$dmu.deta)) warning("link.phi needs to provide dmu.deta component for BC/BR")
  }
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.eta <- phi_linkobj$mu.eta
  phi_dmu.deta <- phi_linkobj$dmu.deta
  ## y* transformation
  ystar <- qlogis(y)

  ## control parameters
  ocontrol <- control
  phi_full <- control$phi
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  control$phi <- control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL

  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, linkfun(y), weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    sigma2 <- sum(weights * res^2)/((sum(weights) - k) * (dlink)^2)
    phi_y <- weights * yhat * (1 - yhat)/(sum(weights) * sigma2) - 1/n
    phi <- rep.int(0, ncol(z))
    phi[1L] <- suppressWarnings(phi_linkfun(sum(phi_y)))
    ## i.e., start out from the fixed dispersion model as described
    ## in Ferrari & Cribari-Neto (2004) (and differing from Simas et al. 2009)
    ## An alternative would be
    ##   phi <- lm.wfit(z, phi_linkfun(phi_y), weights)$coefficients
    ## but that only works in general if all(phi_y > 0) which is not necessarily
    ## the case.
    ##
    ## Additionally, sum(phi_y) might not even be > 0 which should be caught.
    if(!isTRUE(phi_linkinv(phi[1L]) > 0)) {
      warning("no valid starting value for precision parameter found, using 1 instead")
      phi[1L] <- 1
    }
    start <- list(mean = beta, precision = phi)
  }
  if(is.list(start)) start <- do.call("c", start)

  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par, deriv = 0L) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = m) + k]
    eta <- as.vector(x %*% beta + offset[[1L]])
    phi_eta <- as.vector(z %*% gamma + offset[[2L]])
    mu <- linkinv(eta)
    phi <- phi_linkinv(phi_eta)
    mustar <- if(deriv >= 1L) digamma(mu * phi) - digamma((1 - mu) * phi) else NULL
    psi1 <- if(deriv >= 2L) trigamma(mu * phi) else NULL
    psi2 <- if(deriv >= 2L) trigamma((1 - mu) * phi) else NULL
    list(
      beta = beta,
      gamma = gamma,
      eta = eta,
      phi_eta = phi_eta,
      mu = mu,
      phi = phi,
      mustar = mustar,
      psi1 = psi1,
      psi2 = psi2
    )
  }

  ## objective function
  loglikfun <- function(par, fit = NULL) {
    ## extract fitted parameters
    if(is.null(fit)) fit <- fitfun(par)
    alpha <- fit$mu * fit$phi
    beta <- (1 - fit$mu) * fit$phi
    
    ## compute log-likelihood
    if(any(!is.finite(fit$phi)) | any(alpha > 1e300) | any(beta > 1e300)) NaN else { ## catch extreme cases without warning
      ll <- suppressWarnings(dbeta(y, alpha, beta, log = TRUE))
      if(any(!is.finite(ll))) NaN else sum(weights * ll) ## again: catch extreme cases without warning
    }
  }

  ## gradient (by default) or gradient contributions (sum = FALSE)
  gradfun <- function(par, sum = TRUE, fit = NULL) {
    ## extract fitted means/precisions
    if(is.null(fit)) fit <- fitfun(par, deriv = 1L)
    mu <- fit$mu
    phi <- fit$phi
    eta <- fit$eta
    phi_eta <- fit$phi_eta
    mustar <- fit$mustar

    ## compute gradient contributions
    rval <- cbind(
      phi * (ystar - mustar) * mu.eta(eta) * weights * x,
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
        phi_mu.eta(phi_eta) * weights * z
    )
    if(sum) colSums(rval) else rval
  }

  ## analytical Hessian (expected information) or covariance matrix (inverse of Hessian)
  hessfun <- function(par, inverse = FALSE, fit = NULL) {
    ## extract fitted means/precisions
    if(is.null(fit)) fit <- fitfun(par, deriv = 2L)
    mu <- fit$mu
    phi <- fit$phi
    eta <- fit$eta
    phi_eta <- fit$phi_eta
    mustar <- fit$mustar
    psi1 <- fit$psi1
    psi2 <- fit$psi2

    ## auxiliary transformations
    a <- psi1 + psi2
    b <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
    ## compute elements of W
    wbb <- phi^2 * a * mu.eta(eta)^2
    wpp <- b * phi_mu.eta(phi_eta)^2
    wbp <- phi * (mu * a - psi2) * mu.eta(eta) * phi_mu.eta(phi_eta)
    ## compute elements of K
    kbb <- if(k > 0L) crossprod(sqrt(weights) * sqrt(wbb) * x) else crossprod(x)
    kpp <- if(m > 0L) crossprod(sqrt(weights) * sqrt(wpp) * z) else crossprod(z)
    kbp <- if(k > 0L & m > 0L) crossprod(weights * wbp * x, z) else crossprod(x, z)

    ## put together K (= expected information)
    K <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
    if(!inverse) K else chol2inv(chol(K))
    ## previously computed K^(-1) via partitioned matrices, but this appears to be
    ## slower - even for moderately sized problems
    ##   kbb1 <- if(k > 0L) chol2inv(qr.R(qr(sqrt(weights) * sqrt(wbb) * x))) else kbb
    ##   kpp1 <- if(m > 0L) solve(kpp - t(kbp) %*% kbb1 %*% kbp) else kpp
    ##   vcov <- cbind(rbind(kbb1 + kbb1 %*% kbp %*% kpp1 %*% t(kbp) %*% kbb1,
    ##     -kpp1 %*% t(kbp) %*% kbb1), rbind(-kbb1 %*% kbp %*% kpp1, kpp1))
  }

  ## compute biases and adjustment for bias correction/reduction
  biasfun <- function(par, fit = NULL, vcov = NULL) {
    if (is.null(fit)) fit <- fitfun(par, deriv = 2L)
    InfoInv <- if(is.null(vcov)) try(hessfun(par, inverse = TRUE), silent = TRUE) else vcov
    mu <- fit$mu
    phi <- fit$phi
    eta <- fit$eta
    phi_eta <- fit$phi_eta
    D1 <- mu.eta(eta)
    D2 <- phi_mu.eta(phi_eta)
    D1dash <- dmu.deta(eta)
    D2dash <- phi_dmu.deta(phi_eta)
    Psi2 <- fit$psi2
    dPsi1 <-  psigamma(mu * phi, 2)       ## potentially move to fitfun() when we add support for
    dPsi2 <-  psigamma((1 - mu) * phi, 2) ## observed information (as opposed to expected)
    kappa2 <- fit$psi1 + Psi2
    kappa3 <- dPsi1 - dPsi2
    Psi3 <- psigamma(phi, 1)
    dPsi3 <- psigamma(phi, 2)
    ## PQsum produces the the adustments to the score functions and is suggested for iteration
    PQsum <- function(t) {
      if (t <= k)  {
        Xt <- x[,t]
        bb <- if (k > 0L)
          crossprod(x, weights * phi^2 * D1 * (phi * D1^2 * kappa3 + D1dash * kappa2) * Xt * x)
        else
          crossprod(x)
        bg <- if ((k > 0L) & (m > 0L))
          crossprod(x, weights * phi * D1^2 * D2 * (mu * phi * kappa3 + phi * dPsi2 + kappa2) * Xt * z)
        else
          crossprod(x, z)
        gg <- if (m > 0L)
          crossprod(z, weights * phi * D1 * D2^2 * (mu^2 * kappa3 - dPsi2 + 2 * mu * dPsi2) * Xt * z) +
            crossprod(z, weights * phi * D1 * D2dash * (mu * kappa2 - Psi2) * Xt * z)
        else
          crossprod(z)
      } else {
        Zt <- z[, t - k]
        bb <- if (k > 0L)
          crossprod(x, weights * phi * D2 * (phi * D1^2 * mu * kappa3 + phi * D1^2 * dPsi2 + D1dash * mu * kappa2 - D1dash * Psi2) * Zt * x)
        else
          crossprod(x)
        bg <- if ((k > 0L) & (m > 0L))
          crossprod(x, weights * D1 * D2^2 * (phi * mu^2 * kappa3 + phi * (2 * mu - 1) * dPsi2 + mu * kappa2 - Psi2) * Zt * z)
        else
          crossprod(x, z)
        gg <- if (m > 0L)
          crossprod(z, weights * D2^3 * (mu^3 * kappa3 + (3 * mu^2 - 3 * mu + 1) * dPsi2 - dPsi3) * Zt * z) +
            crossprod(z, weights * D2dash * D2 * (mu^2 * kappa2 + (1 - 2 * mu) * Psi2 - Psi3) * Zt * z)
        else
          crossprod(z)
      }
      pq <- rbind(cbind(bb, bg), cbind(t(bg), gg))
      sum(diag(InfoInv %*% pq))/2
    }
    if (inherits(InfoInv, "try-error")) {
      bias <- adjustment <- rep.int(NA_real_, k + m)
    }
    else {
      adjustment <- sapply(1:(k + m), PQsum)
      bias <- - InfoInv %*% adjustment
    }
    list(bias = bias, adjustment = adjustment)
  }


  ## optimize likelihood
  opt <- optim(par = start, fn = loglikfun, gr = gradfun,
    method = method, hessian = hessian, control = control)
  par <- opt$par

  ## conduct further (quasi) Fisher scoring to move ML derivatives
  ## even further to zero or conduct bias reduction
  ## (suppressed if fsmaxit = 0 or if only numerical optim result desired)
  if(type == "BR" & fsmaxit <= 0) warning("BR cannot be performed with fsmaxit <= 0")
  step <- .Machine$integer.max
  iter <- 0
  if(fsmaxit > 0 & !(hessian & type == "ML"))
  {
    for (iter in 1:fsmaxit) {
      stepPrev <- step
      stepFactor <- 0
      testhalf <- TRUE
      while (testhalf & stepFactor < 11) {
        fit <- fitfun(par, deriv = 2L)
        scores <- gradfun(par, fit = fit)
	InfoInv <- try(hessfun(par, fit = fit, inverse = TRUE))
	if(failedInv <- inherits(InfoInv, "try-error")) {
          warning("failed to invert the information matrix: iteration stopped prematurely")
          break
        }
        bias <- if(type == "BR") biasfun(par, fit = fit, vcov = InfoInv)$bias else 0
        par <- par + 2^(-stepFactor) * (step <- InfoInv %*% scores - bias)
        stepFactor <- stepFactor + 1
        testhalf <- drop(crossprod(stepPrev) < crossprod(step))
      }
      if (failedInv | (all(abs(step) < fstol))) {
        break
      }
    }
  }

  ## check whether both optim() and manual iteration converged IK:
  ## modified the condition a bit... optim might fail to converge but
  ## if additional iteration are requested Fisher scoring might get
  ## there
  if((fsmaxit == 0 & opt$convergence > 0) | iter >= fsmaxit) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }

  ## conduct single bias correction (if BC selected) else do not
  ## estimate the first order biases
  if(type == "BC") {
    bias <- as.vector(biasfun(par)$bias)
    par <- par - bias
  }
  else {
    bias <- rep.int(NA_real_, k + m)
  }

  ## extract fitted values/parameters
  fit <- fitfun(par, deriv = 2L)
  beta <- fit$beta
  gamma <- fit$gamma
  eta <- fit$eta
  mu <- fit$mu
  phi <- fit$phi

  ## log-likelihood/gradients/covariance matrix at optimized parameters
  ll <- loglikfun(par, fit = fit)
  ## No need to evaluate ef below.
  ## ef <- gradfun(par, fit = fit, sum = FALSE)
  vcov <- if (hessian & (type == "ML")) solve(-as.matrix(opt$hessian)) else hessfun(fit = fit, inverse = TRUE)

  ## R-squared
  pseudor2 <- if(var(eta) * var(ystar) <= 0) NA else cor(eta, linkfun(y))^2

  ## names
  names(beta) <- colnames(x)
  names(gamma) <- if(phi_const & phi_linkstr == "identity") "(phi)" else colnames(z)
  rownames(vcov) <- colnames(vcov) <- names(bias) <- c(colnames(x),
    if(phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"))

  ## set up return value
  rval <- list(
    coefficients = list(mean = beta, precision = gamma),
    residuals = y - mu,
    fitted.values = structure(mu, .Names = names(y)),
    type = type,
    optim = opt,
    method = method,
    control = ocontrol,
    scoring = iter,
    start = start,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(mean = if(identical(offset[[1L]], rep.int(0, n))) NULL else offset[[1L]],
      precision = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    df.null = nobs - 2,
    df.residual = nobs - k - m,
    phi = phi_full,
    loglik = ll,
    vcov = vcov,
    bias = bias,
    pseudo.r.squared = pseudor2,
    link = list(mean = linkobj, precision = phi_linkobj),
    converged = converged
  )
  return(rval)
}

print.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients$mean)) {
      cat(paste("Coefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in mean model)\n\n")
    if(x$phi) {
      if(length(x$coefficients$precision)) {
        cat(paste("Phi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
        print.default(format(x$coefficients$precision, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      } else cat("No coefficients (in precision model)\n\n")
    }
  }

  invisible(x)
}

summary.betareg <- function(object, phi = NULL, type = "sweighted2", ...)
{
  ## treat phi as full model parameter?
  if(!is.null(phi)) object$phi <- phi

  ## residuals
  type <- match.arg(type, c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type

  ## extend coefficient table
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$precision)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], precision = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$precision) <- names(object$coefficients$precision)
  object$coefficients <- cf

  ## number of iterations
  object$iterations <- c("optim" = as.vector(tail(na.omit(object$optim$count), 1)), "scoring" = as.vector(object$scoring))

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.betareg"
  object
}

print.summary.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    types <- c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
    Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals",
      "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$mean)) {
      cat(paste("\nCoefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in mean model)\n")

    if(x$phi) {
      if(NROW(x$coefficients$precision)) {
        cat(paste("\nPhi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
        printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
      } else cat("\nNo coefficients (in precision model)\n")
    }

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat("\nType of estimator:", x$type, switch(x$type,
      "ML" = "(maximum likelihood)",
      "BC" = "(bias-corrected)",
      "BR" = "(bias-reduced)"))
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    if(x$iterations[2L] > 0) {
      cat(paste("\nNumber of iterations:", x$iterations[1L],
        sprintf("(%s) +", x$method), x$iterations[2L], "(Fisher scoring) \n"))
    } else {
      cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
    }
  }

  invisible(x)
}

predict.betareg <- function(object, newdata = NULL,
  type = c("response", "link", "precision", "variance", "quantile"),
  na.action = na.pass, at = 0.5, ...)
{
  type <- match.arg(type)

  if(type == "quantile") {
    qfun <- function(at, mu, phi) {
      rval <- sapply(at, function(p) qbeta(p, mu * phi, (1 - mu) * phi))
      if(length(at) > 1L) {
        if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
	  dimnames = list(unique(names(rval)), NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      } else {
        rval <- drop(rval)
      }
      rval   
    }
  }

  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values
      },
      "link" = {
        object$link$mean$linkfun(object$fitted.values)
      },
      "precision" = {
        gamma <- object$coefficients$precision
        z <- if(is.null(object$x)) model.matrix(object, model = "precision") else object$x$precision
        offset <- if(is.null(object$offset$precision)) rep.int(0, NROW(z)) else object$offset$precision
        object$link$precision$linkinv(drop(z %*% gamma + offset))
      },
      "variance" = {
        gamma <- object$coefficients$precision
        z <- if(is.null(object$x)) model.matrix(object, model = "precision") else object$x$precision
        offset <- if(is.null(object$offset$precision)) rep.int(0, NROW(z)) else object$offset$precision
        phi <- object$link$precision$linkinv(drop(z %*% gamma + offset))
        object$fitted.values * (1 - object$fitted.values) / (1 + phi)
      },
      "quantile" = {
        mu <- object$fitted.values
        gamma <- object$coefficients$precision
        z <- if(is.null(object$x)) model.matrix(object, model = "precision") else object$x$precision
        offset <- if(is.null(object$offset$precision)) rep.int(0, NROW(z)) else object$offset$precision
        phi <- object$link$precision$linkinv(drop(z %*% gamma + offset))
        qfun(at, mu, phi)
      }
    )
    return(rval)

  } else {

    tnam <- switch(type,
      "response" = "mean",
      "link" = "mean",
      "precision" = "precision",
      "variance" = "full",
      "quantile" = "full")

    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    newdata <- newdata[rownames(mf), , drop = FALSE]
    offset <- list(mean = rep.int(0, nrow(mf)), precision = rep.int(0, nrow(mf)))

    if(type %in% c("response", "link", "variance", "quantile")) {
      X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
      if(!is.null(object$call$offset)) offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
      if(!is.null(off.num <- attr(object$terms$mean, "offset"))) {
        for(j in off.num) offset[[1L]] <- offset[[1L]] + eval(attr(object$terms$mean, "variables")[[j + 1L]], newdata)
      }
    }
    if(type %in% c("precision", "variance", "quantile")) {
      Z <- model.matrix(object$terms$precision, mf, contrasts = object$contrasts$precision)
      if(!is.null(off.num <- attr(object$terms$precision, "offset"))) {
        for(j in off.num) offset[[2L]] <- offset[[2L]] + eval(attr(object$terms$precision, "variables")[[j + 1L]], newdata)
      }
    }

    rval <- switch(type,
      "response" = {
        object$link$mean$linkinv(drop(X %*% object$coefficients$mean + offset[[1L]]))
      },
      "link" = {
        drop(X %*% object$coefficients$mean + offset[[1L]])
      },
      "precision" = {
        object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + offset[[2L]]))
      },
      "variance" = {
        mu <- object$link$mean$linkinv(drop(X %*% object$coefficients$mean + offset[[1L]]))
        phi <- object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + offset[[2L]]))
        mu * (1 - mu) / (1 + phi)
      },
      "quantile" = {
        mu <- object$link$mean$linkinv(drop(X %*% object$coefficients$mean + offset[[1L]]))
        phi <- object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + offset[[2L]]))
        qfun(at, mu, phi)
      }
    )
    return(rval)

  }
}

coef.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...) {
  cf <- object$coefficients

  model <- if(is.null(phi)) {
    if(missing(model)) ifelse(object$phi, "full", "mean") else match.arg(model)
  } else {
    if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
    ifelse(phi, "full", "mean")
  }

  switch(model,
    "mean" = {
      cf$mean
    },
    "precision" = {
      cf$precision
    },
    "full" = {
      nam1 <- names(cf$mean)
      nam2 <- names(cf$precision)
      cf <- c(cf$mean, cf$precision)
      names(cf) <- c(nam1, if(identical(nam2, "(phi)")) "(phi)" else paste("(phi)", nam2, sep = "_"))
      cf
    }
  )
}

vcov.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...) {
  vc <- object$vcov
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$precision)

  model <- if(is.null(phi)) {
    if(missing(model)) ifelse(object$phi, "full", "mean") else match.arg(model)
  } else {
    if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
    ifelse(phi, "full", "mean")
  }

  switch(model,
    "mean" = {
      vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
    },
    "precision" = {
      vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$precision)
      vc
    },
    "full" = {
      vc
    }
  )
}

bread.betareg <- function(x, phi = NULL, ...) {
  vcov(x, phi = phi) * x$nobs
}

estfun.betareg <- function(x, phi = NULL, ...) {
  ## extract response y and regressors X and Z
  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  xmat <- if(is.null(x$x)) model.matrix(x, model = "mean") else x$x$mean
  zmat <- if(is.null(x$x)) model.matrix(x, model = "precision") else x$x$precision
  offset <- x$offset
  if(is.null(offset[[1L]])) offset[[1L]] <- rep.int(0, NROW(xmat))
  if(is.null(offset[[2L]])) offset[[2L]] <- rep.int(0, NROW(zmat))
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  phi_full <- if(is.null(phi)) x$phi else phi

  ## extract coefficients
  beta <- x$coefficients$mean
  gamma <- x$coefficients$precision

  ## compute y*
  ystar <- qlogis(y)

  ## compute mu*
  eta <- as.vector(xmat %*% beta + offset[[1L]])
  phi_eta <- as.vector(zmat %*% gamma + offset[[2L]])
  mu <- x$link$mean$linkinv(eta)
  phi <- x$link$precision$linkinv(phi_eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ## compute scores of beta
  rval <- phi * (ystar - mustar) * as.vector(x$link$mean$mu.eta(eta)) * wts * xmat

  ## combine with scores of phi
  if(phi_full) {
    rval <- cbind(rval,
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
      as.vector(x$link$precision$mu.eta(phi_eta)) * wts * zmat)
    colnames(rval) <- names(coef(x, phi = phi_full))
  }
  attr(rval, "assign") <- NULL

  ##
  if(x$type == "BR") {
    nobs <- nrow(xmat)
    k <- ncol(xmat)
    m <- ncol(zmat)
    InfoInv <- x$vcov
    D1 <- x$link$mean$mu.eta(eta)
    D2 <- x$link$precision$mu.eta(phi_eta)
    D1dash <- x$link$mean$dmu.deta(eta)
    D2dash <- x$link$precision$dmu.deta(phi_eta)
    Psi2 <- psigamma((1 - mu) * phi, 1)
    dPsi1 <-  psigamma(mu * phi, 2)
    dPsi2 <-  psigamma((1 - mu) * phi, 2)
    kappa2 <- psigamma(mu * phi, 1) + Psi2
    kappa3 <- dPsi1 - dPsi2
    Psi3 <- psigamma(phi, 1)
    dPsi3 <- psigamma(phi, 2)
    PQ <- function(t) {
      prodfun <- function(mat1, mat2) {
        sapply(seq_len(nobs), function(i) tcrossprod(mat1[i,], mat2[i,]), simplify = FALSE)
      }
      if (t <= k)  {
        Xt <- xmat[,t]
        bb <- if (k > 0L) {
          bbComp <- wts * phi^2 * D1 * (phi * D1^2 * kappa3 + D1dash * kappa2) * Xt * xmat
          prodfun(xmat, bbComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, k, k))
        bg <- if ((k > 0L) & (m > 0L)) {
          bgComp <- wts * phi * D1^2 * D2 * (mu * phi * kappa3 + phi * dPsi2 + kappa2) * Xt * zmat
          prodfun(xmat, bgComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, k, m))
        gg <- if (m > 0L) {
          ggComp <- wts * phi * D1 * D2^2 * (mu^2 * kappa3 - dPsi2 + 2 * mu * dPsi2) * Xt * zmat +
            wts * phi * D1 * D2dash * (mu * kappa2 - Psi2) * Xt * zmat
          prodfun(zmat, ggComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, m, m))
      } else {
        Zt <- zmat[, t - k]
        bb <- if (k > 0L) {
          bbComp <- wts * phi * D2 * (phi * D1^2 * mu * kappa3 + phi * D1^2 * dPsi2 + D1dash * mu * kappa2 - D1dash * Psi2) * Zt * xmat
          prodfun(xmat, bbComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, k, k))
        bg <- if ((k > 0L) & (m > 0L)) {
          bgComp <- wts * D1 * D2^2 * (phi * mu^2 * kappa3 + phi * (2 * mu - 1) * dPsi2 + mu * kappa2 - Psi2) * Zt * zmat
          prodfun(xmat, bgComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, k, m))
        gg <- if (m > 0L) {
          ggComp <- wts * D2^3 * (mu^3 * kappa3 + (3 * mu^2 - 3 * mu + 1) * dPsi2 - dPsi3) * Zt * zmat +
            wts * D2dash * D2 * (mu^2 * kappa2 + (1 - 2 * mu) * Psi2 - Psi3) * Zt * zmat
          prodfun(zmat, ggComp)
        }
        else
          sapply(1:nobs, function(x) matrix(0, m, m))
      }
      sapply(seq_len(nobs), function(i)
             sum(diag(InfoInv %*% rbind(cbind(bb[[i]], bg[[i]]), cbind(t(bg[[i]]), gg[[i]]))))/2,
             simplify = TRUE)
    }
    if (inherits(InfoInv, "try-error")) {
      adjustment <- rep.int(NA_real_, k + m)
    }
    else
      adjustment <- sapply(1:(k + m), PQ)
    rval <- rval + adjustment
  }
  return(rval)
}

coeftest.betareg <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)

logLik.betareg <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}

terms.betareg <- function(x, model = c("mean", "precision"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.betareg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.Formula(formula$formula)) formula$call$formula <- formula$formula <-
    formula(formula$formula, collapse = TRUE)
  formula$terms <- formula$terms$full
  NextMethod()
}

model.matrix.betareg <- function(object, model = c("mean", "precision"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

residuals.betareg <- function(object,
  type = c("sweighted2", "deviance", "pearson", "response", "weighted", "sweighted"), ...)
{
  ## raw response residuals and desired type
  res <- object$residuals
  type <- match.arg(type)
  if(type == "response") return(res)

  ## extract fitted information
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  mu <- fitted(object)
  wts <- weights(object)
  if(is.null(wts)) wts <- 1
  phi <- predict(object, type = "precision")

  res <- switch(type,

    "pearson" = {
      sqrt(wts) * res / sqrt(mu * (1 - mu) / (1 + phi))
    },

    "deviance" = {
      ll <- function(mu, phi)
        (lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
        (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y))
      sqrt(wts) * sign(res) * sqrt(2 * abs(ll(y, phi) - ll(mu, phi)))
    },

    "weighted" = {
      ystar <- qlogis(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(phi * v)
    },

    "sweighted" = {
      ystar <- qlogis(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(v)
    },

    "sweighted2" = {
      ystar <- qlogis(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(v * (1 - hatvalues(object)))
    })

  return(res)
}

cooks.distance.betareg <- function(model, ...)
{
    h <- hatvalues(model)
    k <- length(model$coefficients$mean)
    res <- residuals(model, type = "pearson")
    h * (res^2)/(k * (1 - h)^2)
}

update.betareg <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
