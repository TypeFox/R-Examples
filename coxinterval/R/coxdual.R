coxdual <- function(formula, data = parent.frame(), subset, init = NULL,
                    formula.coxph = NULL, init.coxph = FALSE, control, ...)
{
  ## extract model frame and perform input validation
  cl <- match.call(expand.dots = FALSE)
  ## set parameters controlling model fit
  control <- if (missing(control)) coxdual.control(...)
             else do.call(coxdual.control, control)
  datargs <- c("formula", "data", "subset")
  mf <- cl[c(1, match(datargs, names(cl), nomatch = 0))]
  mf[[1]] <- as.name("model.frame")
  specials <- c("trans", "cluster", "strata", "tt")
  ftrm <- if (missing(data)) terms(formula, specials)
          else terms(formula, data = data, specials)
  if (with(attr(ftrm, "specials"), length(c(strata, tt))))
    stop("The 'strata' and 'tt' terms not supported.")
  ## store column indices of terms in model frame
  irsp <- attr(ftrm, "response")
  ityp <- attr(ftrm, "specials")$type
  itrn <- attr(ftrm, "specials")$trans
  icls <- attr(ftrm, "specials")$cluster
  if (length(itrn) != 1) stop("Model requires exactly one 'trans' term.")
  if (length(icls) != 1) stop("Model requires exactly one 'cluster' term.")
  icov <- (1:(length(attr(ftrm, "variables")) - 1))[-c(irsp, itrn, icls)]
  mf$formula <- ftrm
  mf$na.action <- as.name("na.coxdual")
  suppressWarnings(mf <- eval(mf, parent.frame()))
  mt <- attr(mf, "terms")
  mm <- model.matrix(mt, mf)
  ## find covariates model matrix
  if (length(icov)) {
    asgn <- frame.assign(mf, mt, mm)
    jcov <- subset.data.frame(asgn, frame %in% icov)$matrix
  }
  else jcov <- 1
  if (!inherits(mf[, irsp], "Surv")) stop("Response is not a 'Surv' object.")
  if (attr(mf[, irsp], "type") != "counting")
    stop("Response is not a 'counting'-type 'Surv' object.")
  ## in case NA action didn't apply trans attributes
  if (!is.null(attr(mf[, itrn], "states"))) {
    attr(mf, "states") <- attr(mf[, itrn], "states")
    attr(mf, "types") <- attr(mf[, itrn], "types")
  }
  states <- attr(mf, "states")
  if (length(states) != 3 | length(attr(mf, "types")) != 3)
    stop("Invalid state transitions in the model 'trans' term.")
  ## sort data
  ord <- order(mf[, icls], mf[, itrn][, 1], mf[, itrn][, 2])
  mf <- mf[ord, ]
  mm <- mm[ord, ]
  d <- coxdual.data(mf[, icls], mf[, irsp][, 1], mf[, irsp][, 2],
                    mf[, itrn][, 1], mf[, itrn][, 2], mf[, irsp][, 3],
                    mm[, jcov], states, control$sieve, control$eps)
  n <- nrow(d$z)
  ncov <- length(jcov)
  if (control$sieve) {
    part <- mapply(function(s, r) s[r >= control$risk.min],
                   d$supp, d$risk, SIMPLIFY = FALSE)
    npart <- mapply(function(s, k) max(1, round(length(s)/k)), part,
                    with(control, sieve.const * n^sieve.rate), SIMPLIFY = FALSE)
    npart <- mapply(function(s, k) seq(k, length(s), k),
                    part, npart, SIMPLIFY = FALSE)
    part <- mapply(function(s, k) c(floor(min(d$u)), s[k[-c(1, length(k))]],
                                    ceiling(max(d$v))), part, npart,
                   SIMPLIFY = FALSE)
  }
  else part <- lapply(d$supp, function(x) c(min(d$u), x))
  npart <- sapply(part, length)
  ## type-specific partition
  tvec <- do.call("c", part)
  names(tvec) <- NULL
  ## common partition
  svec <- sort(unique(tvec))
  ## strictly dual right censoring?
  if (with(d, all(left[contrib == 1] == right[contrib == 1]))) censor <- "right"
  else censor <- "interval"
  ## fit right-censored data alternatives with survival's coxph
  if (censor == "right" & is.null(formula.coxph) & init.coxph)
    formula.coxph <- list(cl$formula)
  else if (class(formula.coxph) == "formula")
    formula.coxph <- list(formula.coxph)
  fit.coxph <- list()
  if (!is.null(formula.coxph)) {
    for (i in 1:length(formula.coxph)) {
      fit.coxph[[i]] <- cl
      fit.coxph[[i]]$formula <- update.formula(fit.coxph[[i]]$formula,
                                               formula.coxph[[i]])
      temp <- gsub("^trans\\(", "strata\\(",
                   attr(terms(fit.coxph[[i]]$formula), "term.labels"))
      ## dispense with extraneous terms for null model
      if (!length(icov)) temp <- c(temp[grep("^strata\\(", temp)], "1")
      fit.coxph[[i]]$formula <-
        update.formula(fit.coxph[[i]]$formula,
                       as.formula(paste("~", paste(temp, collapse = " + "))))
      temp <- list(formula = fit.coxph[[i]]$formula, data = data,
                   na.action = "na.omit")
      if (!missing(subset)) temp <- c(temp, subset)
      invisible(capture.output(fit.coxph[[i]] <- try(do.call("coxph", temp))))
      if (inherits(fit.coxph[[i]], "try-error"))
        fit.coxph[[i]] <- list(call = temp, error = fit.coxph[[i]][1])
      else {
        temp <- if (is.null(fit.coxph[[i]]$na.action)) 1:nrow(mf)
                else -fit.coxph[[i]]$na.action
        fit.coxph[[i]] <-
          list(call = fit.coxph[[i]]$call,
               n = length(unique(mf[temp, icls])),
               m = fit.coxph[[i]]$n,
               p = length(fit.coxph[[i]]$coefficients),
               coef = fit.coxph[[i]]$coefficients,
               var = fit.coxph[[i]]$var,
               basehaz = basehaz(fit.coxph[[i]], centered = FALSE),
               iter = fit.coxph[[i]]$iter,
               loglik = fit.coxph[[i]]$loglik,
               na.action = fit.coxph[[i]]$na.action)
        if (censor == "right" & i == 1) {
          if (length(icov))
            rownames(fit.coxph[[i]]$var) <-
              colnames(fit.coxph[[i]]$var) <- colnames(mm)[jcov]
          names(fit.coxph[[i]]$basehaz) <- c("hazard", "time", "trans")
        }
        fit.coxph[[i]]$call$data <- cl$data
        class(fit.coxph[[i]]) <- "coxinterval"
      }
    }
  }
  else fit.coxph <- NULL
  ## initial values
  if (is.null(init)) init <- list()
  init.coxph <- init.coxph & !is.null(fit.coxph[[1]]$basehaz)
  if (init.coxph) {
    init$coef <- fit.coxph[[1]]$coef
    init$basehaz <- fit.coxph[[1]]$basehaz
  }
  if (is.null(init$coef) | !length(icov)) init$coef <- rep(0, ncov)
  else {
    if (length(init$coef) == 1) init$coef <- rep(init$coef, ncov)
    else if (length(init$coef) != ncov)
      stop("Initial value needs ", ncov, " regression coefficients.")
  }
  if (is.null(init$basehaz)) {
    init$basehaz <- tvec / max(d$v)
    basehaz <- NULL
  }
  else {
    if (ncol(init$basehaz) != 3)
      stop("Invalid initial baseline cumulative intensity data frame.")
    if (!is.data.frame(init$basehaz)) init$basehaz <- data.frame(init$basehaz)
    if (!all(names(init$basehaz) == c("hazard", "time", "trans")))
      init$basehaz <- init$basehaz[c("hazard", "time", "trans")]
    else
      names(init$basehaz) <- c("hazard", "time", "trans")
    init$basehaz <- init$basehaz[with(init$basehaz, order(trans, time)), ]
    if (length(unique(init$basehaz$trans)) != length(states))
      stop("Invalid initial baseline cumulative intensity transition types.")
    basehaz <- init$basehaz
  }
  if (!is.null(basehaz)) {
    basehaz <- step2jump(basehaz, stratum = "trans")
    if (any(basehaz$hazard < 0))
      stop("Initial baseline cumulative intensity must be nondecreasing.")
    basehaz <- jump2step(basehaz, stratum = "trans")
    init$basehaz <-
      do.call("c", mapply(linapprox,
                          with(basehaz, split(data.frame(time, hazard), trans)),
                          part, SIMPLIFY = FALSE))
  }
  init$basehaz <- data.frame(hazard = init$basehaz, time = tvec,
                             trans = rep(c(1, 2, 12), times = npart))
  basehaz <- if (control$sieve) lin2const(init$basehaz, stratum = "trans")
             else step2jump(init$basehaz, stratum = "trans")
  list(init = init, basehaz = basehaz)
  fit <- .C("coxdual",
            coef = as.double(init$coef),
            basehaz = as.double(basehaz$hazard),
            as.integer(ncov),
            as.integer(npart),
            as.double(tvec),
            as.double(svec),
            as.integer(length(svec)),
            as.double(as.matrix(d$z)),
            as.integer(n),
            as.double(d$left),
            as.double(d$right),
            as.double(d$u),
            as.double(d$v),
            as.integer(d$contrib),
            as.integer(d$absorb),
            var = as.double(rep(0, ncov^2)),
            loglik = as.double(rep(0, control$iter.max + 1)),
            as.double(control$eps),
            as.integer(control$iter.max),
            as.double(control$coef.typ),
            as.double(control$coef.max),
            zerocoef = as.integer(length(icov) == 0),
            sieve = as.integer(control$sieve),
            iter = as.integer(0),
            maxnorm = as.double(0),
            gradnorm = as.double(0),
            cputime = as.double(0),
            flag = as.integer(0),
            NAOK = TRUE)
  if (fit$flag == 1 & length(icov) > 0)
    stop("Parameter estimation failed; coefficient Hessian not invertible.")
  if (with(fit, any(is.na(coef), is.nan(coef), is.na(basehaz), is.nan(basehaz))))
    stop("Parameter estimation failed.")
  if (fit$flag == 2)
    stop("Variance estimation failed; profile information not invertible.")
  if (with(fit, any(is.na(diag(var)), is.nan(diag(var)))))
    stop("Variance estimation failed.")
  if (with(fit, iter == control$iter.max & maxnorm > control$eps))
    warning("Maximum iterations reached before convergence.")
  names(fit$coef) <- names(init$coef) <- colnames(mm)[jcov]
  var <- matrix(fit$var, ncov)
  rownames(var) <- colnames(var) <- colnames(mm)[jcov]
  init$init.coxph <- init.coxph
  init$basehaz$trans <- as.factor(init$basehaz[, 3])
  levels(init$basehaz[, 3]) <- attr(mf, "types")
  basehaz <- data.frame(hazard = fit$basehaz, time = init$basehaz$time,
                        trans = init$basehaz$trans)
  basehaz <- if (control$sieve) const2lin(basehaz, stratum = "trans")
             else jump2step(basehaz, stratum = "trans")
  rownames(basehaz) <- rownames(init$basehaz) <- NULL
  censor.rate <-
    with(d, cbind("Status and survival" = sum(contrib != 0 & absorb),
                  "Status only" = sum(contrib != 0 & !absorb),
                  "Survival only" = sum(contrib == 0 & absorb),
                  "Neither" = sum(contrib == 0 & !absorb))) / n
  rownames(censor.rate) <- "Observation rate"
  if (length(fit.coxph) == 1) fit.coxph <- fit.coxph[[1]]
  fit <- list(call = cl, censor = censor, n = n, m = nrow(mf),
              p = ncov * (length(icov) > 0), coef = fit$coef, var = var,
              basehaz = basehaz, init = init,
              loglik = n * with(fit, loglik[1:(iter + 1)]),
              iter = fit$iter, maxnorm = fit$maxnorm, gradnorm = fit$gradnorm,
              cputime = fit$cputime, coxph = fit.coxph,
              na.action = attr(mf, "na.action"), censor.rate = censor.rate,
              control = control, nullmodel = fit$zerocoef)
  if (control$data) fit$data <- d
  if (length(fit$coxph) == 1) fit$coxph <- fit$coxph[[1]]
  class(fit) <- c("coxdual", "coxinterval")
  fit
}
