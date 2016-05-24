predict.rma <-
function (object, newmods, intercept, tau2.levels, gamma2.levels, 
    addx = FALSE, level, digits, transf, targs, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- object
    if (missing(newmods)) 
        newmods <- NULL
    if (missing(intercept)) 
        intercept <- x$intercept
    if (missing(tau2.levels)) 
        tau2.levels <- NULL
    if (missing(gamma2.levels)) 
        gamma2.levels <- NULL
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (x$knha) {
        crit <- qt(alpha/2, df = x$dfs, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    if (x$int.only && !is.null(newmods)) 
        stop("Cannot specify new moderator values for models without moderators.")
    if (is.null(newmods)) {
        if (!is.element("rma.mv", class(object))) {
            if (x$int.only) {
                k.new <- 1
                X.new <- cbind(1)
            }
            else {
                k.new <- x$k.f
                X.new <- x$X.f
            }
        }
        else {
            if (x$int.only) {
                if (!x$withG) {
                  k.new <- 1
                  X.new <- cbind(1)
                }
                if (x$withG && x$withH) {
                  if (is.null(tau2.levels) && is.null(gamma2.levels)) {
                    k.new <- x$tau2s * x$gamma2s
                    X.new <- cbind(rep(1, k.new))
                    if (x$tau2s == 1) {
                      tau2.levels <- rep(1, k.new)
                    }
                    else {
                      tau2.levels <- rep(levels(x$mf.g.f$inner), 
                        each = x$gamma2s)
                    }
                    if (x$gamma2s == 1) {
                      gamma2.levels <- rep(1, k.new)
                    }
                    else {
                      gamma2.levels <- rep(levels(x$mf.h.f$inner), 
                        times = x$tau2s)
                    }
                  }
                  if ((!is.null(tau2.levels) && is.null(gamma2.levels)) || 
                    (is.null(tau2.levels) && !is.null(gamma2.levels))) 
                    stop("Either specify both of tau2.levels and gamma2.levels or neither.")
                  if (!is.null(tau2.levels) && !is.null(gamma2.levels)) {
                    if (length(tau2.levels) != length(gamma2.levels)) 
                      stop("Length of tau2.levels and gamma2.levels must be the same.")
                    k.new <- length(tau2.levels)
                    X.new <- cbind(rep(1, k.new))
                  }
                }
                if (x$withG && !x$withH) {
                  if (is.null(tau2.levels)) {
                    k.new <- x$tau2s
                    X.new <- cbind(rep(1, k.new))
                    if (x$tau2s == 1) {
                      tau2.levels <- rep(1, k.new)
                    }
                    else {
                      tau2.levels <- levels(x$mf.g.f$inner)
                    }
                  }
                  else {
                    k.new <- length(tau2.levels)
                    X.new <- cbind(rep(1, k.new))
                  }
                  gamma2.levels <- rep(1, k.new)
                }
            }
            else {
                k.new <- x$k.f
                X.new <- x$X.f
                if (!is.null(tau2.levels) || !is.null(gamma2.levels)) 
                  warning("Arguments tau2.levels and gamma2.levels ignored when obtaining fitted values.")
                tau2.levels <- as.character(x$mf.g.f$inner)
                gamma2.levels <- as.character(x$mf.h.f$inner)
            }
        }
    }
    else {
        if ((!x$int.incl && x$p == 1L) || (x$int.incl && x$p == 
            2L)) {
            k.new <- length(newmods)
            X.new <- cbind(c(newmods))
        }
        else {
            if (is.vector(newmods) || nrow(newmods) == 1L) {
                k.new <- 1
                X.new <- rbind(newmods)
            }
            else {
                k.new <- nrow(newmods)
                X.new <- cbind(newmods)
            }
        }
        if (x$int.incl) {
            if (intercept) {
                X.new <- cbind(intrcpt = rep(1, k.new), X.new)
            }
            else {
                X.new <- cbind(intrcpt = rep(0, k.new), X.new)
            }
        }
        if (ncol(X.new) != x$p) 
            stop("Dimensions of 'newmods' do not match dimensions of the model.")
    }
    if (is.element("rma.mv", class(object)) && x$withG) {
        if (x$tau2s > 1) {
            if (is.null(tau2.levels)) {
                warning("Need to specify 'tau2.levels' argument to obtain credibility intervals.")
            }
            else {
                if (!is.numeric(tau2.levels) && anyNA(pmatch(tau2.levels, 
                  x$g.levels.f[[1]], duplicates.ok = TRUE))) 
                  stop("Non-existing levels specified via 'tau2.levels' argument.")
                if (is.numeric(tau2.levels)) {
                  tau2.levels <- round(tau2.levels)
                  if (any(tau2.levels < 1) || any(tau2.levels > 
                    x$g.nlevels.f[1])) 
                    stop("Non-existing tau^2 values specified via 'tau2.levels' argument.")
                }
                if (length(tau2.levels) == 1L) 
                  tau2.levels <- rep(tau2.levels, k.new)
                if (length(tau2.levels) != k.new) 
                  stop("Length of 'tau2.levels' does not match number of predicted values.")
            }
        }
        else {
            tau2.levels <- rep(1, k.new)
        }
    }
    if (is.element("rma.mv", class(object)) && x$withH) {
        if (x$gamma2s > 1) {
            if (is.null(gamma2.levels)) {
                warning("Need to specify 'gamma2.levels' argument to obtain credibility intervals.")
            }
            else {
                if (!is.numeric(gamma2.levels) && anyNA(pmatch(gamma2.levels, 
                  x$h.levels.f[[1]], duplicates.ok = TRUE))) 
                  stop("Non-existing levels specified via 'gamma2.levels' argument.")
                if (is.numeric(gamma2.levels)) {
                  gamma2.levels <- round(gamma2.levels)
                  if (any(gamma2.levels < 1) || any(gamma2.levels > 
                    x$h.nlevels.f[1])) 
                    stop("Non-existing gamma^2 values specified via 'gamma2.levels' argument.")
                }
                if (length(gamma2.levels) == 1L) 
                  gamma2.levels <- rep(gamma2.levels, k.new)
                if (length(gamma2.levels) != k.new) 
                  stop("Length of 'gamma2.levels' does not match number of predicted values.")
            }
        }
        else {
            gamma2.levels <- rep(1, k.new)
        }
    }
    pred <- rep(NA_real_, k.new)
    vpred <- rep(NA_real_, k.new)
    for (i in seq_len(k.new)) {
        Xi.new <- X.new[i, , drop = FALSE]
        pred[i] <- Xi.new %*% x$b
        vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
    }
    se <- sqrt(vpred)
    ci.lb <- pred - crit * se
    ci.ub <- pred + crit * se
    if (!is.element("rma.mv", class(object))) {
        cr.lb <- pred - crit * sqrt(vpred + x$tau2)
        cr.ub <- pred + crit * sqrt(vpred + x$tau2)
    }
    else {
        if (!x$withG) {
            cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2))
            cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2))
        }
        if (x$withG && !x$withH) {
            if (x$tau2s == 1) {
                cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2)
                cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2)
            }
            else {
                if (is.null(tau2.levels)) {
                  cr.lb <- rep(NA, k.new)
                  cr.ub <- rep(NA, k.new)
                  tau2.levels <- rep(NA, k.new)
                }
                else {
                  if (!is.numeric(tau2.levels)) 
                    tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], 
                      duplicates.ok = TRUE)
                  cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels])
                  cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels])
                  tau2.levels <- x$g.levels.f[[1]][tau2.levels]
                }
            }
        }
        if (x$withG && x$withH) {
            if (x$tau2s == 1 && x$gamma2s == 1) {
                cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2 + x$gamma2)
                cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                  x$tau2 + x$gamma2)
            }
            else {
                if (is.null(tau2.levels) || is.null(gamma2.levels)) {
                  cr.lb <- rep(NA, k.new)
                  cr.ub <- rep(NA, k.new)
                  tau2.levels <- rep(NA, k.new)
                  gamma2.levels <- rep(NA, k.new)
                }
                else {
                  if (!is.numeric(tau2.levels)) 
                    tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], 
                      duplicates.ok = TRUE)
                  if (!is.numeric(gamma2.levels)) 
                    gamma2.levels <- pmatch(gamma2.levels, x$h.levels.f[[1]], 
                      duplicates.ok = TRUE)
                  cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels] + x$gamma2[gamma2.levels])
                  cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + 
                    x$tau2[tau2.levels] + x$gamma2[gamma2.levels])
                  tau2.levels <- x$g.levels.f[[1]][tau2.levels]
                  gamma2.levels <- x$h.levels.f[[1]][gamma2.levels]
                }
            }
        }
    }
    if (is.function(transf)) {
        if (is.null(targs)) {
            pred <- sapply(pred, transf)
            se <- rep(NA, k.new)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            cr.lb <- sapply(cr.lb, transf)
            cr.ub <- sapply(cr.ub, transf)
        }
        else {
            pred <- sapply(pred, transf, targs)
            se <- rep(NA, k.new)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            cr.lb <- sapply(cr.lb, transf, targs)
            cr.ub <- sapply(cr.ub, transf, targs)
        }
        transf <- TRUE
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    cr.bounds <- cbind(cr.lb, cr.ub)
    rev.order <- ifelse(cr.ub < cr.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    cr.bounds[rev.order, ] <- cr.bounds[rev.order, 2:1]
    cr.lb <- cr.bounds[, 1]
    cr.ub <- cr.bounds[, 2]
    if (is.null(newmods) && !x$int.only) {
        slab <- x$slab
    }
    else {
        slab <- seq_len(k.new)
    }
    if (k.new == 1L) 
        slab <- ""
    if (na.act == "na.omit") {
        not.na <- !is.na(pred)
    }
    else {
        not.na <- rep(TRUE, k.new)
    }
    out <- list(pred = pred[not.na], se = se[not.na], ci.lb = ci.lb[not.na], 
        ci.ub = ci.ub[not.na], cr.lb = cr.lb[not.na], cr.ub = cr.ub[not.na])
    if (is.element("rma.mv", class(object)) && x$withG && x$tau2s > 
        1) 
        out$tau2.level <- tau2.levels
    if (is.element("rma.mv", class(object)) && x$withH && x$gamma2s > 
        1) 
        out$gamma2.level <- gamma2.levels
    if (addx) 
        out$X <- matrix(X.new[not.na, ], ncol = x$p)
    out$slab <- slab[not.na]
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (addx) 
        colnames(out$X) <- colnames(x$X)
    if (x$method == "FE") {
        out$cr.lb <- NULL
        out$cr.ub <- NULL
    }
    out$digits <- digits
    out$method <- x$method
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
