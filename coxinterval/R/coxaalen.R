coxaalen <- function(formula, data = parent.frame(), subset, init = NULL,
                     formula.timereg = NULL, init.timereg = FALSE, control, ...)
{
  ## extract model frame and perform input validation
  cl <- match.call(expand.dots = FALSE)
  if (!is.loaded("coxaalen", "coxinterval"))
    stop("Required CPLEX-dependent libraries are unavailable.")
  ## set parameters controlling model fit
  control <- if (missing(control)) coxaalen.control(...)
             else do.call(coxaalen.control, control)
  datargs <- c("formula", "data", "subset")
  mf <- cl[c(1, match(datargs, names(cl), nomatch = 0))]
  mf[[1]] <- as.name("model.frame")
  specials <- c("prop", "strata", "cluster", "tt")
  ftrm <- if (missing(data)) terms(formula, specials)
          else terms(formula, data = data, specials)
  if (with(attr(ftrm, "specials"), length(c(strata, cluster, tt))))
    stop("The 'strata', 'cluster' and 'tt' terms not supported")
  ## store column indices of terms in model frame
  irsp <- attr(ftrm, "response")
  iprp <- attr(ftrm, "specials")$prop
  iadd <- (1:(length(attr(ftrm, "variables")) - 1))[-c(irsp, iprp)]
  if (!length(iprp)) stop("Model requires 'prop' terms")
  mf$formula <- ftrm
  mf$na.action <- as.name("na.omit")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mm <- model.matrix(mt, mf)
  ## find proportional and additive terms in model matrix
  asgn <- frame.assign(mf, mt, mm)
  jprp <- subset.data.frame(asgn, frame %in% iprp)$matrix
  ## additive term always includes intercept
  jadd <- if (length(iadd)) c(1, subset.data.frame(asgn, frame %in% iadd)$matrix)
          else 1
  if (!inherits(mf[, irsp], "Surv")) stop("Response is not a 'Surv' object")
  ## Surv converts interval2 to interval
  if (attr(mf[, irsp], "type") != "interval")
    stop("Response is not an 'interval2'-type 'Surv' object")
  ## fit right-censored data alternatives with timereg's cox.aalen
  fit.timereg <- list()
  if (requireNamespace("timereg", quietly = TRUE) & !is.null(formula.timereg)) {
    if (class(formula.timereg) == "formula")
      formula.timereg <- list(formula.timereg)
    for (i in 1:length(formula.timereg)) {
      if (!missing(subset)) warning("Model alternatives not based on subset")
      keep <- 1:nrow(data)
      if (!is.null(omit <- attr(mf, "na.action"))) keep <- keep[-omit]
      fit.timereg[[i]] <- cl
      fit.timereg[[i]]$formula <-
        update.formula(fit.timereg[[i]]$formula, formula.timereg[[i]])
      ## cox.aalen arguments, constructed to avoid NA-related errors
      temp <- list(formula = fit.timereg[[i]]$formula, data = data[keep, ],
                   robust = 0, silent = 1, max.timepoint.sim = nrow(data))
      invisible(capture.output(fit.timereg[[i]] <-
        try(do.call(cox.aalen, temp))))
      if (inherits(fit.timereg[[i]], "try-error"))
        fit.timereg[[i]] <- list(call = temp, error = fit.timereg[1])
      else {
        temp <- rownames(fit.timereg[[i]]$gamma)
        fit.timereg[[i]] <- list(call = fit.timereg[[i]]$call,
                                 n = length(keep),
                                 p = nrow(fit.timereg[[i]]$gamma),
                                 coef = as.vector(fit.timereg[[i]]$gamma),
                                 var = fit.timereg[[i]]$var.gamma,
                                 basehaz = as.data.frame(fit.timereg[[i]]$cum))
        names(fit.timereg[[i]]$coef) <- temp
      }
      fit.timereg[[i]]$call$data <- cl$data
      class(fit.timereg[[i]]) <- "coxinterval"
    }
  }
  else fit.timereg <- NULL
  n <- nrow(mf)
  nadd <- length(jadd)
  nprp <- length(jprp)
  tmax <- max(mf[, irsp][mf[, irsp][, 3] == 0, 1],
              mf[, irsp][mf[, irsp][, 3] == 3, 2])
  if (tmax == tmax - control$eps)
    stop("Observations large relative to epsilon. Use a smaller time scale.")
  time <- maximalint(mf[, irsp], eps = control$eps)
  ntime <- nrow(time$int) - 1
  time$int <- time$int[1:ntime, ]
  A <- coxaalen.ineq(mm[, jadd[-1]], ntime)
  ## initial parameter values
  if (is.null(init)) init <- list()
  init.timereg <- init.timereg & !is.null(fit.timereg[[1]]$coef)
  if (init.timereg) {
    init$coef <- fit.timereg[[1]]$coef
    init$basehaz <- fit.timereg[[1]]$basehaz
  }
  if (is.null(init$coef)) init$coef <- rep(0, nprp)
  else {
    if (length(init$coef) == 1) init$coef <- rep(init$coef, nprp)
    else if (length(init$coef) != nprp)
      stop("Initial value needs ", nprp, " regression coefficients.")
  }
  if (is.null(init$basehaz)) {
    init$basehaz <- cbind(with(time, int[, 2] / int[ntime, 2]),
                          matrix(0, ntime, nadd - 1))
    basehaz <- NULL
  }
  else {
    if (ncol(init$basehaz) != nadd + 1)
      stop("Invalid initial baseline cumulative hazard data frame.")
    if (!is.data.frame(init$basehaz)) init$basehaz <- data.frame(init$basehaz)
    colnames(init$basehaz)[match(c("intercept", "Intercept"),
                                 colnames(init$basehaz))] <- "(Intercept)"
    if (all(c("time", colnames(mm)[jadd]) %in% colnames(init$basehaz)))
      init$basehaz <- init$basehaz[c("time", colnames(mm)[jadd])]
    else names(init$basehaz) <- c("time", colnames(mm)[jadd])
    init$basehaz <- init$basehaz[order(init$basehaz$time), ]
    basehaz <- init$basehaz
  }
  if (!is.null(basehaz)) {
    ## evaluate on tied right endpoints of maximal intersections
    if (nadd > 1)
      init$basehaz <-
        apply(basehaz[, -1], 2,
              function(x) linapprox(data.frame(basehaz$time, x), time$int[, 3]))
    else init$basehaz <- linapprox(basehaz[, 1:2], time$int[, 3])
    ## add to baseline hazard if linear constraints are not met
    basehaz <- matrix(A %*% as.vector(t(init$basehaz)), ncol = nadd,
                      byrow = TRUE)
    basehaz[basehaz > 0] <- 0
    basehaz <- cumsum(apply(-basehaz, 1, max))
    if (length(dim(init$basehaz)))
      init$basehaz[, "(Intercept)"] <- init$basehaz[, "(Intercept)"] + basehaz
    else init$basehaz <- init$basehaz + basehaz
  }
  fit <- .C("coxaalen",
            coef = as.double(init$coef),
            basehaz = as.double(t(init$basehaz)),
            as.integer(ntime),
            as.double(as.matrix(mm[, jprp])),
            as.integer(n),
            as.integer(nprp),
            as.double(as.matrix(mm[, jadd])),
            as.integer(nadd),
            as.integer(time$ind),
            as.double(A),
            as.integer(nrow(A)),
            as.double(control$eps),
            as.integer(control$eps.norm == "max"),
            as.integer(control$iter.max),
            as.double(control$armijo),
            as.integer(control$var.coef),
            as.double(control$coef.typ),
            as.double(control$coef.max),
            as.integer(control$trace),
            as.integer(control$thread.max),
            var = as.double(rep(0, nprp^2)),
            loglik = as.double(rep(0, control$iter.max + 1)),
            iter = as.integer(0),
            maxnorm = as.double(0),
            gradnorm = as.double(0),
            cputime = as.double(0),
            flag = as.integer(0))
  if (fit$flag == 1)
    stop("Parameter estimation failed; coefficient Hessian not invertible.")
  if (with(fit,
           any(is.na(coef), is.nan(coef), is.na(basehaz), is.nan(basehaz))))
    stop("Parameter estimation failed.")
  if (fit$flag == 2)
    stop("Variance estimation failed; profile information not invertible.")
  if (fit$flag == 3)
    stop("Variance estimation failed; profile maximizer not found.")
  if (with(fit, any(is.na(diag(var)), is.nan(diag(var)))))
    stop("Variance estimation failed.")
  if (with(fit, iter == control$iter.max & maxnorm > control$eps))
    warning("Maximum iterations reached before convergence.")
  names(fit$coef) <- names(init$coef) <- colnames(mm)[jprp]
  if (control$var.coef) {
    var <- matrix(fit$var, nprp)
    rownames(var) <- colnames(var) <- colnames(mm)[jprp]
  }
  else var <- matrix(NA, nprp, nprp)
  init$init.timereg <- init.timereg
  init$basehaz <- data.frame(time$int[, 2], init$basehaz)
  init$basehaz <- rbind(0, init$basehaz)
  basehaz <- data.frame(time$int[, 2], t(matrix(fit$basehaz, nadd)))
  basehaz <- rbind(0, basehaz)
  names(basehaz) <- names(init$basehaz) <- c("time", colnames(mm)[jadd])
  censor.rate <-
    cbind("Exact" = sum(mf[, irsp][, 3] == 1),
          "Left" = sum(mf[, irsp][, 1] == 0),
          "Interval" = sum(mf[, irsp][, 1] > 0 & mf[, irsp][, 3] == 3),
          "Right" = sum(mf[, irsp][, 3] == 0)) / n
  rownames(censor.rate) <- "Censoring rate"
  fit <- list(call = cl, n = n, p = nprp, coef = fit$coef, var = var,
              basehaz = basehaz, init = init,
              loglik = n * with(fit, loglik[1:(iter + 1)]), iter = fit$iter,
              maxnorm = fit$maxnorm, gradnorm = abs(fit$gradnorm),
              cputime = fit$cputime, timereg = fit.timereg,
              na.action = attr(mf, "na.action"), censor.rate = censor.rate,
              control = control)
  if (control$data)
    fit$data <- list(maximalint = time$int, prop = mm[, jprp], add = mm[, jadd])
  if (length(fit$timereg) == 1) fit$timereg <- fit$timereg[[1]]
  class(fit) <- c("coxaalen", "coxinterval")
  fit
}
