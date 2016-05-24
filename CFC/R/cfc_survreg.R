cfc.prepdata <- function(formul, dat) {
  vars.all <- all.vars(formul)
  var.time <- vars.all[1]
  var.status <- vars.all[2]
  vars.expl <- vars.all[-c(1,2)]
  
  causes.plus.censoring <- sort(unique(dat[, var.status]))
  causes <- setdiff(causes.plus.censoring, 0)
  K <- length(causes)
  vars.newstatus <- paste("status", causes, sep = "_")
  formul.list <- list()
  for (k in 1:K) {
    dat[, vars.newstatus[k]] <- 1 * (dat[, var.status] == causes[k])
    formul.list[[k]] <- as.formula(paste0("Surv(", var.time, ",", vars.newstatus[k], ") ~ ", paste(vars.expl, collapse = "+")))
  }
  formul.noresp <- as.formula(paste0("~", paste(vars.expl, collapse = "+")))
  return (list(K = K, dat = dat, formula.list = formul.list, formula.noresp = formul.noresp, tmax = max(dat[, var.time])))
}

cfc.survreg.survprob <- function(t, args, n) { # predicting survival probability at a given time from index using object of class "survreg"
  dist.name <- args$dist
  mydist <- unlist(unname(survreg.distributions[dist.name]))
  if (any(names(mydist) == "dist")) { # derived distribution
    mydist.base <- unlist(unname(survreg.distributions[mydist$dist]))
    p <- mydist.base$density((mydist$trans(t) - t(args$coefficients) %*% args$x[n, ]) / args$scale)[, 1]
  } else { # base distribution
    p <- mydist$density((t - t(args$coefficients) %*% args$x[n, ]) / args$scale)[, 1]
  }
  return (1 - p)
}

cfc.survreg <- function(formula, data, newdata = NULL, dist = "weibull", control = survreg.control()
  , tout, Nmax = 100L, rel.tol = 1e-5) {  # TODO: fix this to use newdata, not data during prediction
  
  # prepare data for cause-specific survival regression
  ret <- cfc.prepdata(formul = formula, dat = data)
  K <- ret$K
  formula.list <- ret$formula.list
  formula.noresp <- ret$formula.noresp
  dat <- ret$dat
  if (missing(tout)) tout <- seq(from = 0.0, to = ret$tmax, length.out = 100L)
  
  # survival regression on each cause
  if (length(dist) == 1) dist <- rep(dist, K)
  reg.list <- list()
  for (k in 1:K) reg.list[[k]] <- survreg(formula.list[[k]], dat, dist = dist[k], control = control, x = TRUE)
  
  # preparing new data
  if (!is.null(newdata)) {
    for (k in 1:K) {
      mf <- model.frame(formula.noresp, newdata)
      mm <- model.matrix(formula.noresp, mf)
      reg.list[[k]]$x <- mm[, colnames(reg.list[[k]]$x)]
    }
  }
  
  # calculating cumulative incidence
  f.list <- list(cfc.survreg.survprob, cfc.survreg.survprob)
  cscr.out <- cscr.samples.R(f.list, reg.list, tout, Nmax = Nmax, rel.tol = rel.tol, nrow(reg.list[[1]]$x))
  class(cscr.out) <- "cfc"
  
  ret.final <- list(K = K, formulas = formula.list, regs = reg.list, tout = tout, cfc = cscr.out)
  class(ret.final) <- c("cfc.survreg", class(ret.final))
  
  return (ret.final)
}

summary.cfc.survreg <- function(object, obs.idx = "all", ...) {
  if (obs.idx[1] == "all") obs.idx <- 1:nrow(object$regs[[1]]$x)
  ci.mean <- apply(X = object$cfc$ci, MARGIN = c(1,2), FUN = mean)
  s.mean <- apply(X = object$cfc$s, MARGIN = c(1,2), FUN = mean)
  
  ret <- list(tout = object$tout, ci = ci.mean, s = s.mean)
  class(ret) <- "summary.cfc.survreg"
  
  return (ret)
}

plot.summary.cfc.survreg <- function(x, which = c(1, 2), ...) {
  K <- ncol(x$ci)
  
  if (1 %in% which) {
    ylimm <- range(x$ci)
    col.vec <- 1:K
    plot(x$tout, x$ci[, 1], type = "l", ylim = ylimm, col = col.vec[1]
         , xlab = "time from index", ylab = "cumulative incidence")
    for (k in 2:K) {
      lines(x$tout, x$ci[, k], col = col.vec[k])
    }
    legend("topleft", legend = paste("cause", 1:K), col = col.vec, lty = rep(1, K))
  }
  
  if (2 %in% which) {
    for (k in 1:K) {
      plot(x$tout, x$ci[, k], type = "l", ylim = range(x$ci[, k], 1 - x$s[, k])
           , xlab = "time from index", ylab = "cumulative incidence", main = paste("cause", k))
      lines(x$tout, 1 - x$s[, k], lty = 2)
      legend("topleft", legend = c("competing-risk adjustment", "no adjustment"), col = c(1, 1), lty = c(1,2))
    }
  }
}




