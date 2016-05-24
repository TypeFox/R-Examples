quant.boost <-
function (formula, data, mstop = NA, quantiles = NA, cv = TRUE) 
{
    if (any(is.na(quantiles)) || !is.vector(quantiles) || any(quantiles > 
        1) || any(quantiles < 0)) {
        pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        row.grid = 3
        col.grid = 4
    }
    else {
        pp <- quantiles
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np <- length(pp)
    if (any(is.na(mstop))) 
        mstop = rep(4000, np)
    else if (length(mstop) == 1) 
        mstop = rep(mstop, np)
    else if (length(mstop) < np) {
        mstop = c(mstop, rep(max(mstop), times = np - length(mstop)))
        warning("mstop doesn't match number of quantiles")
    }
    blsstr <- labels(terms(formula))
    types = list()
    bls <- list()
    x <- list()
    bnd = list()
    yy = eval(parse(text = formula[2]), envir = data, enclos = environment(formula))
    attr(yy, "name") = deparse(formula[[2]])
    if (length(blsstr) == 0) {
        types[[1]] = "parametric"
        x[[1]] = rep(1, length(yy))
        blsstr = 1
    }
    else for (i in 1:length(blsstr)) {
        types[[i]] = strsplit(blsstr[i], "(", fixed = TRUE)[[1]][1]
        if (types[[i]] == blsstr[i]) {
            types[[i]] = "parametric"
            x[[i]] = eval(parse(text = blsstr[i]), envir = data, 
                enclos = environment(formula))
            names(x[[i]]) = blsstr[i]
            bls[[i]] = NULL
        }
        else {
            bls[[i]] <- eval(parse(text = blsstr[i]), envir = data, 
                enclos = environment(formula))
            xi <- NULL
            for (k in 1:length(bls[[i]]$get_names())) xi <- cbind(xi, 
                eval(parse(text = bls[[i]]$get_names()[k]), envir = data, 
                  enclos = environment(formula)))
            x[[i]] = xi
            names(x[[i]]) = bls[[i]]$get_names()
        }
        if (types[[i]] != "bmrf") 
            bnd[[i]] = NA
        else bnd[[i]] = environment(bls[[i]]$dpp)$args$bnd
    }
    m = length(yy)
    kfld <- 10
    ntest <- floor(m/kfld)
    cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, m)), kfld - 
        1), rep(0, m * kfld - (kfld - 1) * (m + ntest))), nrow = m)
    values <- list()
    coef <- list()
    z <- list()
    helper <- list()
    for (k in 1:length(blsstr)) {
        values[[k]] = matrix(NA, nrow = m, ncol = np)
        coef[[k]] = NULL
        fitted = matrix(NA, nrow = m, ncol = np)
        helper[[k]] = NA
        if (types[[k]] == "bbs" || types[[k]] == "bmono") {
            z[[k]] = values[[k]]
            types[[k]] = "pspline"
        }
        else if (types[[k]] == "parametric") 
            z[[k]] = values[[k]]
        else if (types[[k]] == "bols") {
            z[[k]] = values[[k]]
            types[[k]] = "parametric"
        }
        else if (types[[k]] == "bmrf") {
            z[[k]] = matrix(NA, nrow = length(attr(bnd[[k]], 
                "regions")), ncol = np)
            types[[k]] = "markov"
            helper[[k]] = list(bnd[[k]], NULL)
        }
        else if (types[[k]] == "brandom") {
            z[[k]] = matrix(NA, nrow = length(unique(x[[k]])), 
                ncol = np)
            types[[k]] = "random"
        }
        else if (types[[k]] == "bspatial") {
            z[[k]] = list()
            types[[k]] = "2dspline"
        }
        else if (types[[k]] == "brad") {
            z[[k]] = list()
            types[[k]] = "radial"
        }
    }
    dummy.reg <- function(p, formula, data, mstop, pp, cv10f, 
        types, x, blsstr, bnd) {
        values <- list()
        z <- list()
        coef <- list()
        inb <- gamboost(formula = formula, data = data, control = boost_control(mstop = mstop[p], 
            nu = 0.1, risk = "inbag"), family = QuantReg(pp[p]))
        if (cv) {
            cvr = cvrisk(inb, folds = cv10f)
            cat("Quantile ", pp[p], ", Boosting iterations: ", 
                mstop(cvr), "\n")
            inb = inb[mstop(cvr)]
        }
        else cat("Quantile ", pp[p], "\n")
        fitted = fitted(inb)
        for (k in 1:length(blsstr)) {
            coef[[k]] = coef(inb, which = k)[[1]]
            if (types[[k]] == "pspline") {
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                values[[k]] = predict(inb, d.tmp, which = k) + 
                  inb$offset
                z[[k]] = values[[k]]
            }
            else if (types[[k]] == "2dspline" || types[[k]] == 
                "radial") {
                x.min = apply(x[[k]], 2, min)
                x.max = apply(x[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 20), times = 20), rep(seq(x.min[2], 
                  x.max[2], length = 20), each = 20))
                independent = c(which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[1]), which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[2]))
                d.tmp = matrix(0, nrow = 400, ncol = dim(data)[2])
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                d.tmp[, independent] = x.gitter
                dimnames(d.tmp) = list(1:400, dimnames(data)[[2]])
                for (j in 1:dim(d.val)[2][-independent]) d.val[, 
                  j] = NA
                z[[k]] <- predict(inb, data.frame(d.tmp), which = k) + 
                  inb$offset
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] = t(matrix(z[[k]], nrow = 20, ncol = 20))
            }
            else if (types[[k]] == "random") {
                districts = unique(x[[k]])
                d.tmp = matrix(0, nrow = length(districts), ncol = dim(data)[2])
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                for (j in 1:dim(d.val)[2][-independent]) d.val[, 
                  j] = NA
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] <- predict(inb, d.tmp, which = k) + inb$offset
            }
            else if (types[[k]] == "markov") {
                districts = attr(bnd[[k]], "regions")
                d.tmp = data.frame(matrix(0, nrow = length(districts), 
                  ncol = dim(data)[2]))
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                for (j in (1:dim(d.val)[2])[-independent]) d.val[, 
                  j] = NA
                values[[k]] = predict(inb, d.val, which = k) + 
                  inb$offset
                z[[k]] <- predict(inb, d.tmp, which = k) + inb$offset
            }
            else if (types[[k]] == "parametric") {
                independent = which(dimnames(data)[[2]] == blsstr[k])
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] = values[[k]]
            }
        }
        gc()
        list(values, z, fitted, inb, coef)
    }
    if (.Platform$OS.type == "unix") 
        coef.vector = mclapply(1:np, function(i) dummy.reg(i, 
            formula, data, mstop, pp, cv10f, types, x, blsstr, 
            bnd), mc.cores = max(1, min(detectCores() - 1, 2)))
    else if (.Platform$OS.type == "windows") 
        coef.vector = mclapply(1:np, function(i) dummy.reg(i, 
            formula, data, mstop, pp, cv10f, types, x, blsstr, 
            bnd), mc.cores = 1)
    boost.object = list()
    for (k in 1:length(blsstr)) {
        coef[[k]] = matrix(NA, nrow = length(coef.vector[[1]][[5]][[k]]), 
            ncol = np)
    }
    for (i in 1:np) {
        boost.object[[i]] = coef.vector[[i]][[4]]
        fitted[, i] = coef.vector[[i]][[3]]
        for (k in 1:length(blsstr)) {
            values[[k]][, i] = coef.vector[[i]][[1]][[k]]
            coef[[k]][, i] = coef.vector[[i]][[5]][[k]]
            if (types[[k]] == "2dspline" || types[[k]] == "radial") 
                z[[k]][[i]] = coef.vector[[i]][[2]][[k]]
            else z[[k]][, i] = coef.vector[[i]][[2]][[k]]
        }
    }
    result = list(values = values, response = yy, covariates = x, 
        formula = formula, asymmetries = pp, effects = types, 
        helper = helper, fitted = fitted, coefficients = coef, 
        intercepts = rep(0, np), mboost = boost.object)
    result$predict <- function(newdata = NULL) {
        values = list()
        fitted = matrix(NA, nrow = nrow(newdata), ncol = np)
        for (i in 1:np) {
            fitted[, i] = predict(boost.object[[i]], newdata = newdata)
        }
        for (k in 1:length(blsstr)) {
            values[[k]] = matrix(NA, nrow = nrow(newdata), ncol = np)
            for (i in 1:np) {
                values[[k]][, i] = predict(boost.object[[i]], 
                  which = k, newdata = newdata) + boost.object[[i]]$offset
            }
        }
        list(fitted = fitted, values = values)
    }
    result$plotpredict <- function(which = 1) {
        return(z[[which]])
    }
    class(result) = c("expectreg", "boost")
    result
}
