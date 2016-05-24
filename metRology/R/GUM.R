#
# Changes:
#
# 2014-03-04 (S Ellison) Removed explicit require and associated check:
#            require() is now deprected in code when the packages is
#            listed in Depends or Imports
GUM <- function (var.name, x.i, u.i, nu.i, measurement.fnc, correlation = diag(length(var.name)), 
    shared.u.i = var.name, cl = 0.95, cov.factor = "Student's t", 
    sig.digits.U = 2, ...) 
{
    msgs = character(15)
    meq = parse(text = measurement.fnc)
    if (!is.matrix(x.i)) {
        x.i <- matrix(x.i, nrow = 1)
        u.i <- matrix(u.i, nrow = 1)
        nu.i <- matrix(nu.i, nrow = 1)
    }
    n <- ncol(x.i)
    if (!is.numeric(x.i) | !is.numeric(u.i) | !is.numeric(nu.i) | !is.numeric(correlation))
    {
        msgs[1] <- "Found non-numeric input where numeric input is required."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    if (min(u.i)<0)
    {
        msgs[1] <- "Found one or more negative values for standard uncertainty in input."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    if (min(nu.i)<0)
    {
        msgs[1] <- "Found one or more negative values for degrees of freedom in input."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    if (measurement.fnc == "" | is.na(measurement.fnc) | is.numeric(measurement.fnc)) {
        msgs[1] <- "No valid measurement function found in input."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    meq.names <- all.vars(meq)
    if (length(meq.names) != length(var.name)) 
        if (any(!(meq.names %in% var.name))) {
            msgs[1] <- "Variable names in measurement function do not match input variable names."
            return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, 
                k = NA, U = NA, contributions = t(rep(NA, n)), 
                sensitivities = t(rep(NA, n)), msgs = msgs))
        }
        else msgs[which.max(msgs == "")] <- "Warning: Not all input variables are used in the measurement function."
    else if (any(!(sort(meq.names) == sort(var.name)))) {
        msgs[1] <- "Variable names in measurement function do not match input variable names."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    rshared.u.i <- rank(shared.u.i)
    gname <- unique(rshared.u.i)
    for (grp in gname) if (min(nu.i[, grp == rshared.u.i]) != 
        max(nu.i[, grp == rshared.u.i])) 
        msgs[which.max(msgs == "")] <- "Warning: Degrees of freedom are expected to be the same for standard uncertainties based on the same shared value."
    if (min(eigen(correlation)$values) <= .Machine$double.eps) {
        msgs[1] <- "Correlation matrix is not valid (i.e. not positive definite) and could cause negative uncertainties to be reported."
        return(list(y = NA, uc = NA, nu.eff = NA, cl = cl, k = NA, 
            U = NA, contributions = t(rep(NA, n)), sensitivities = t(rep(NA, 
                n)), msgs = msgs))
    }
    for (i in 1:n) assign(var.name[i], x.i[, i])
    cmat = try(attr(eval(deriv(meq, var.name)), "gradient"), 
        silent = TRUE)
    if (inherits(cmat, "try-error")) {
        #require() in code with listed Depends/Imports is deprecated - SLRE
        #if (!require(numDeriv, quietly = TRUE)) {
        #    cat("ERROR: R package numDeriv is not available\n")
        #    return(rep(NA, 2))
        #}
        cmat <- NULL
        for (i in 1:nrow(x.i)) {
            x <- x.i[i, ]
            f <- function(x) {
                for (j in 1:n) assign(var.name[j], x[j])
                eval(meq)
            }
            cmat <- rbind(cmat, grad(f, x, ...))
        }
    }
    cu.i <- cmat * u.i
    uc <- sqrt(diag(cu.i %*% correlation %*% t(cu.i)))
    contributions <- cu.i^2/uc^2
    if (sum(correlation == diag(n)) < n^2) 
        contributions <- t(rep(NA, n))
    group <- function(cmat, umat, dfmat, correlation, namevec, 
        u.group) {
        if (sum(u.group == namevec) != length(u.group)) {
            nreps <- nrow(umat)
            u.group <- rank(u.group)
            gname <- unique(u.group)
            gumat <- matrix(nrow = nreps, ncol = length(gname))
            gcmat <- matrix(0, nrow = nreps, ncol = length(gname))
            gdfmat <- gumat
            i <- 1
            for (grp in gname) {
                gdfmat[, i] <- matrix(dfmat[, grp == u.group], 
                  nrow = nreps)[, 1]
                gumat[, i] <- matrix(umat[, grp == u.group], 
                  nrow = nreps)[, 1]
                subumat <- matrix(umat[, grp == u.group], nrow = nreps)
                subcmat <- matrix(cmat[, grp == u.group], nrow = nreps)
                subcorr <- correlation[grp == u.group, grp == 
                  u.group]
                if (min(subumat[, 1]^2) != 0) {
                  nsubcmat <- subcmat * subumat/subumat[, 1]
                  gcmat[, i] <- sqrt(diag(nsubcmat %*% subcorr %*% 
                    t(nsubcmat)))
                }
                i <- i + 1
            }
            dfmat <- gdfmat
            cmat <- gcmat
            umat <- gumat
        }
        list(dfmat = dfmat, cumat = cmat * umat)
    }
    grouped = group(cmat, u.i, nu.i, correlation, var.name, shared.u.i, 
        ...)
    cu.i = grouped$cumat
    nu.i = grouped$dfmat
    nu.eff <- pmax(1, uc^4/diag((cu.i^2/sqrt(nu.i)) %*% t(cu.i^2/sqrt(nu.i))))
    nu.eff[is.na(nu.eff)] <- 99999 
    if (cov.factor == "Student's t") 
        k <- qt((1 + cl)/2, nu.eff)
    else {
        k <- 2
        nu.eff <- NA
        cl <- NA
        msgs[which.max(msgs == "")] <- "Effective degrees of freedom and nominal confidence level not available when k=2 is specified."
    }
    sensitivities <- cmat
    sigden = 10^floor(log10(k * uc))
    sigden = pmax(10^floor(log10(k * uc)),.Machine$double.eps)
    sigden[sigden==Inf] <- .Machine$double.eps
    uc = round(uc/sigden, sig.digits.U) * sigden
    U = round(k * uc/sigden, sig.digits.U - 1) * sigden
    y = round(eval(meq)/sigden, sig.digits.U - 1) * sigden
    if (sum(correlation == diag(n)) < n^2) 
        msgs[which.max(msgs == "")] <- "Contributions to uncertainty not available due to correlation among input values."
    if (min(abs(sensitivities)) <= .Machine$double.eps) 
        msgs[which.max(msgs == "")] <- "Warning: Results may not be accurate since one or more sensitivity coefficients are zero."
    if (msgs[1] == "") 
        msgs[1] <- "No errors or warnings."
    list(y = y, uc = uc, nu.eff = nu.eff, cl = cl, k = k, U = U, 
        contributions = contributions, sensitivities = sensitivities, 
        msgs = msgs)
}
