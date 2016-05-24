superpc.cv=
function (fit, data, n.threshold = 20, n.fold = NULL, folds = NULL, 
    n.components = 3, min.features = 5, max.features = nrow(data$x), 
    compute.fullcv = TRUE, compute.preval = TRUE, xl.mode = c("regular", 
        "firsttime", "onetime", "lasttime"), xl.time = NULL, 
    xl.prevfit = NULL) 
{
    this.call <- match.call()
    xl.mode = match.arg(xl.mode)
    if (xl.mode == "regular" | xl.mode == "firsttime") {
        type <- fit$type
        if (n.components > 5) {
            cat("Max # of components is 5", fill = TRUE)
        }
        n.components <- min(5, n.components)
        if (type == "survival" & is.null(n.fold) & is.null(folds)) {
            n.fold = 2
        }
        if (type == "regression" & is.null(n.fold) & is.null(folds)) {
            n.fold = 10
        }
        n <- ncol(data$x)
        cur.tt <- fit$feature.scores
        lower <- quantile(abs(cur.tt), 1 - (max.features/nrow(data$x)))
        upper <- quantile(abs(cur.tt), 1 - (min.features/nrow(data$x)))
        if (!is.null(folds)) {
            n.fold = length(folds)
        }
        if (is.null(folds)) {
            if (type == "survival" & sum(data$censoring.status == 
                0) > 0) {
                folds = balanced.folds(data$censoring.status, 
                  nfolds = n.fold)
                folds = c(folds, balanced.folds(data$censoring.status, 
                  nfolds = n.fold))
                folds = c(folds, balanced.folds(data$censoring.status, 
                  nfolds = n.fold))
                folds = c(folds, balanced.folds(data$censoring.status, 
                  nfolds = n.fold))
                folds = c(folds, balanced.folds(data$censoring.status, 
                  nfolds = n.fold))
                n.fold = length(folds)
            }
            else {
                folds <- vector("list", n.fold)
                breaks <- round(seq(from = 1, to = (n + 1), length = (n.fold + 
                  1)))
                cv.order <- sample(1:n)
                for (j in 1:n.fold) {
                  folds[[j]] <- cv.order[(breaks[j]):(breaks[j + 
                    1] - 1)]
                }
            }
        }
        featurescores.folds <- matrix(nrow = nrow(data$x), ncol = n.fold)
        thresholds <- seq(from = lower, to = upper, length = n.threshold)
        nonzero <- rep(0, n.threshold)
        scor <- array(NA, c(n.components, n.threshold, n.fold))
        scor.preval <- matrix(NA, nrow = n.components, ncol = n.threshold)
        scor.lower = NULL
        scor.upper = NULL
        v.preval <- array(NA, c(n, n.components, n.threshold))
    }
    if (xl.mode == "onetime" | xl.mode == "lasttime") {
        type = xl.prevfit$type
        scor = xl.prevfit$scor
        scor.preval = xl.prevfit$scor.preval
        scor.lower = xl.prevfit$scor.lower
        scor.upper = xl.prevfit$scor.upper
        v.preval = xl.prevfit$v.preval
        folds = xl.prevfit$folds
        n.fold = xl.prevfit$n.fold
        nonzero = xl.prevfit$nonzero
        featurescores.folds = xl.prevfit$featurescores.folds
        n.threshold = xl.prevfit$n.threshold
        thresholds = xl.prevfit$thresholds
        compute.fullcv = compute.fullcv
        compute.preval = compute.preval
    }
    if (xl.mode == "regular") {
        first = 1
        last = n.fold
    }
    if (xl.mode == "firsttime") {
        first = 1
        last = 1
    }
    if (xl.mode == "onetime") {
        first = xl.time
        last = xl.time
    }
    if (xl.mode == "lasttime") {
        first = n.fold
        last = n.fold
    }
    for (j in first:last) {
        cat("", fill = TRUE)
        cat(c("fold=", j), fill = TRUE)
        data.temp = list(x = data$x[, -folds[[j]]], y = data$y[-folds[[j]]], 
            censoring.status = data$censoring.status[-folds[[j]]])
        cur.tt <- superpc.train(data.temp, type = type, s0.perc = fit$s0.perc)$feature.scores
        featurescores.folds[, j] <- cur.tt
        for (i in 1:n.threshold) {
            cat(i)
            cur.features <- (abs(cur.tt) > thresholds[i])
            if (sum(cur.features) > 1) {
                nonzero[i] <- nonzero[i] + sum(cur.features)/n.fold
                cur.svd <- mysvd(data$x[cur.features, -folds[[j]]], 
                  n.components = n.components)
                xtemp = data$x[cur.features, folds[[j]], drop = FALSE]
                xtemp <- t(scale(t(xtemp), center = cur.svd$feature.means, 
                  scale = F))
                cur.v.all <- scale(t(xtemp) %*% cur.svd$u, center = FALSE, 
                  scale = cur.svd$d)
                n.components.eff <- min(sum(cur.features), n.components)
                cur.v <- cur.v.all[, 1:n.components.eff,drop=FALSE]
                v.preval[folds[[j]], 1:n.components.eff, i] <- cur.v
                if (compute.fullcv) {
                  for (k in 1:ncol(cur.v)) {
                    if (type == "survival") {
                      require(survival)
                      junk <- coxph(Surv(data$y[folds[[j]]], 
                        data$censoring.status[folds[[j]]]) ~ 
                        cur.v[, 1:k], control = coxph.control(iter.max = 10))$loglik
                      scor[k, i, j] <- 2 * (junk[2] - junk[1])
                    }
                    else {
                      junk <- summary(lm(data$y[folds[[j]]] ~ 
                        cur.v[, 1:k]))
                      scor[k, i, j] <- junk$fstat[1]
                    }
                  }
                }
            }
        }
    }
    cat("\n")
    if (xl.mode == "regular" | xl.mode == "lasttime") {
        mean.na <- function(x) {
            mean(x[!is.na(x)])
        }
        se.na <- function(x) {
            val = NA
            if (sum(!is.na(x)) > 0) {
                val = sqrt(var(x[!is.na(x)])/sum(!is.na(x)))
            }
            return(val)
        }
        lscor = apply(log(scor), c(1, 2), mean.na)
        se.lscor = apply(log(scor), c(1, 2), se.na)
        scor.lower = exp(lscor - se.lscor)
        scor.upper = exp(lscor + se.lscor)
        scor <- exp(lscor)
        if (compute.preval) {
            for (i in 1:n.threshold) {
                for (j in 1:n.components) {
                  if (sum(is.na(v.preval[, 1:j, i])) == 0) {
                    if (type == "survival") {
                      require(survival)
                      junk <- coxph(Surv(data$y, data$censoring.status) ~ 
                        v.preval[, 1:j, i])$loglik
                      scor.preval[j, i] <- 2 * (junk[2] - junk[1])
                    }
                    else {
                      junk <- summary(lm(data$y ~ v.preval[, 
                        1:j, i]))
                      scor.preval[j, i] <- junk$fstat[1]
                    }
                  }
                }
            }
        }
    }
    junk <- list(thresholds = thresholds, n.threshold = n.threshold, 
        nonzero = nonzero, scor.preval = scor.preval, scor = scor, 
        scor.lower = scor.lower, scor.upper = scor.upper, folds = folds, 
        n.fold = n.fold, featurescores.folds = featurescores.folds, 
        v.preval = v.preval, compute.fullcv = compute.fullcv, 
        compute.preval = compute.preval, type = type, call = this.call)
    class(junk) <- "superpc.cv"
    return(junk)
}

