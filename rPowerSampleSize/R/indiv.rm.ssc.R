################################################################################
#       Function for sample size computation                                   #
#       Decision Rules: At least r endpoints significant among m               #
#       Authors: Delorme P., Lafaye de Micheaux P., Liquet B. and Riou J.      #
#       Date: Le 17/09/2015   Version 1                                        #
################################################################################

indiv.rm.ssc <- function(method, asympt = FALSE, r, m, p = m, nCovernE = 1, muC = NULL, muE = NULL, d = NULL, delta = NULL, SigmaC = NULL, SigmaE = NULL, power = 0.8, alpha = 0.05, interval = c(2, 2000), q = 1, maxpts = 25000, abseps = 0.001, releps = 0, nbcores = 1, LB = FALSE, orig.Hochberg = FALSE) {

    if (missing(method)) stop("Missing 'method' argument.")
    if (class(method) != "character") stop("The 'method' argument should be of type character.")
    if ( (method != "Bonferroni") & (method != "Hochberg") & (method != "Holm") ) stop("The 'method' argument is misspecified.")
    if (missing(r)) stop("The 'r' argument is missing.")
    if (missing(m)) stop("The 'm' argument is missing.")  
    if (!is.numeric(r)) stop("The 'r' argument should be numeric.")
    if (!is.numeric(q)) stop("The 'q' argument should be numeric.")
    if (!is.numeric(m)) stop("The 'm' argument should be numeric.")
    if (!is.numeric(power)) stop("The 'power' argument should be numeric.")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")  
    if (!is.vector(interval) || (length(interval) != 2)) stop("The 'interval' argument should be vector of length 2.")
    if (interval[1] >= interval[2]) stop("'interval' is misspecified.")
    if (interval[1] <= 1) stop("'interval' should start at least at 2.")
    if (r > m) stop("The 'r' argument should be lower than or equal to the argument 'm'.")
    if (r > p) stop("The 'r' argument should be lower than or equal to the argument 'p'.")
    if (q > m) stop("The 'q' argument should be lower than or equal to the argument 'm'.")
    if ((power < 0) | (power > 1)) stop("The 'power' argument should be between 0 and 1.")
    if ((alpha < 0) | (alpha > 1)) stop("The 'alpha' argument should be between 0 and 1.")
    if (missing(SigmaC)) stop("The 'SigmaC' argument is missing.")
    if (missing(SigmaE)) stop("The 'SigmaE' argument is missing.")
    if (!is.matrix(SigmaC)) stop("The 'SigmaC' argument should be a matrix.")
    if (!is.matrix(SigmaE)) stop("The 'SigmaE' argument should be a matrix.")
    if ((ncol(SigmaC) != m) && (nrow(SigmaC) != m)) stop("The 'SigmaC' should be a squared matrix of dimension m.")
    if ((ncol(SigmaE) != m) && (nrow(SigmaE) != m)) stop("The 'SigmaE' should be a squared matrix of dimension m.")
    
    if (is.null(delta)) {
        if (missing(d)) stop("The 'd' argument is missing.")
        if (missing(muC)) stop("The 'muC' argument is missing.")
        if (missing(muE)) stop("The 'muE' argument is missing.")
        if (!is.vector(muC)) stop("The 'muC' argument should be a vector.")
        if (!is.vector(muE)) stop("The 'muE' argument should be a vector.")
        if (!is.vector(d)) stop("The 'd' argument should be vector.")
        if (length(muC) != m) stop("The 'muC' vector should be of length 'm'.")
        if (length(muE) != m) stop("The 'muC' vector should be of length m.")
        if (length(d) != m) stop("The 'd' vector should be of length m.")
        # Article p. 5:
        delta <- muE - muC - d # muE, muC and d should be given so that delta is the value under H_1 if we want to compute the power, or is null if we want the size of the test.
    } else {
        if (!is.null(muE) && !is.null(muC) && !is.null(d)) stop("'muE', 'muC', and 'd' should all be NULL.")
    }
    if (!is.numeric(nbcores) || (nbcores < 1)) stop("'nbcores' should be an integer greater than 1.")
    nbcores <- as.integer(nbcores)
    if (!is.logical(LB)) stop("'LB' should be a logical.")
    if (orig.Hochberg && (q != 1)) stop("'orig.Hochberg' should be set to TRUE only when q=1.") 

    SigmaEC.wo.nE <- SigmaE + SigmaC / nCovernE
    Vm1over2.wo.nE <- diag(1 / sqrt(diag(SigmaEC.wo.nE)))
    R <- Vm1over2.wo.nE %*% SigmaEC.wo.nE %*% Vm1over2.wo.nE # Correlation matrix.
    # Note: If \Sigma^g = \sigma^{2,g} K_{\rho}, g = C, E, then R is in fact equal to K_{\rho}.

    exchangeable <- if ((length(unique(R[lower.tri(R)])) == 1) && (length(unique(Vm1over2.wo.nE %*% as.matrix(delta))) == 1)) TRUE else FALSE

    if (exchangeable) {
        res <- .indive.rm.ssc(method = method, asympt = asympt, r = r, m = m, p = m, nCovernE = nCovernE, 
                              delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, power = power, alpha = alpha,
                              interval = interval, q = q, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)
        
    } else {
        res <- .indivne.rm.ssc(method = method, asympt = asympt, r = r, m = m, p = m, nCovernE = nCovernE, 
                               delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, power = power, alpha = alpha,
                               interval = interval, q = q, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)
    }
    return(res)    
}



# Fonctions, une pour chaque cas de calcul de puissance (echangeable ou pas),
# car c est interessant de pouvoir calculer la puissance et pas seulement le sample size.
Psirms <- function(r, m, p = m, nE, nCovernE = 1, delta, SigmaC, SigmaE, alpha = 0.05, q = 1, asympt = FALSE, maxpts = 25000, abseps = 0.001, releps = 0, nbcores = 1, LB = FALSE) {

    SigmaEC.wo.nE <- SigmaE + SigmaC / nCovernE
    Vm1over2.wo.nE <- diag(1 / sqrt(diag(SigmaEC.wo.nE)))
    R <- Vm1over2.wo.nE %*% SigmaEC.wo.nE %*% Vm1over2.wo.nE # Correlation matrix.
    # Note: If \Sigma^g = \sigma^{2,g} K_{\rho}, g = C, E, then R is in fact equal to K_{\rho}.
    exchangeable <- if ((length(unique(R[lower.tri(R)])) == 1) && (length(unique(Vm1over2.wo.nE %*% as.matrix(delta))) == 1)) TRUE else FALSE
 
    if (exchangeable) {
        res <- .Psirmse(r = r, m = m, p = p, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
    } else {
        res <- .Psirmsne(r = r, m = m, p = p, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)
    }
    return(list(pow = res$pow, error = res$error))
}

Psirmu <- function(r, m, p = m, nE, nCovernE = 1, delta, SigmaC, SigmaE, alpha = 0.05, q = 1, asympt = FALSE, maxpts = 25000, abseps = 0.001, releps = 0, nbcores = 1, LB = FALSE, orig.Hochberg = FALSE) {

    SigmaEC.wo.nE <- SigmaE + SigmaC / nCovernE
    Vm1over2.wo.nE <- diag(1 / sqrt(diag(SigmaEC.wo.nE)))
    R <- Vm1over2.wo.nE %*% SigmaEC.wo.nE %*% Vm1over2.wo.nE # Correlation matrix.
    # Note: If \Sigma^g = \sigma^{2,g} K_{\rho}, g = C, E, then R is in fact equal to K_{\rho}.
    exchangeable <- if ((length(unique(R[lower.tri(R)])) == 1) && (length(unique(Vm1over2.wo.nE %*% as.matrix(delta))) == 1)) TRUE else FALSE
    
    if (exchangeable) {
        res <- .Psirmue(r = r, m = m, p = m, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)
    } else {
        res <- .Psirmune(r = r, m = m, p = p, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)
    }
    return(list(pow = res$pow, error = res$error))
}

Psirmd <- function(r, m, p = m, nE, nCovernE = 1, delta, SigmaC, SigmaE, alpha = 0.05, q = 1, asympt = FALSE, maxpts = 25000, abseps = 0.001, releps = 0, nbcores = 1, LB = FALSE) {

    SigmaEC.wo.nE <- SigmaE + SigmaC / nCovernE
    Vm1over2.wo.nE <- diag(1 / sqrt(diag(SigmaEC.wo.nE)))
    R <- Vm1over2.wo.nE %*% SigmaEC.wo.nE %*% Vm1over2.wo.nE # Correlation matrix.
    # Note: If \Sigma^g = \sigma^{2,g} K_{\rho}, g = C, E, then R is in fact equal to K_{\rho}.
    exchangeable <- if ((length(unique(R[lower.tri(R)])) == 1) && (length(unique(Vm1over2.wo.nE %*% as.matrix(delta))) == 1)) TRUE else FALSE
    
    if (exchangeable) {
        res <- .Psirmde(r = r, m = m, p = p, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)
    } else {
        res <- .Psirmdne(r = r, m = m, p = p, nE = nE, nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)
    }
    return(list(pow = res$pow, error = res$error))
}


matrix.type.compute <- function(SigmaE, SigmaC, display.type = FALSE) {
        # We determine the type of the variance matrices
        compoundE <- ((length(unique(diag(SigmaE))) == 1) && (length(unique(SigmaE[lower.tri(SigmaE)])) == 1))
        compoundC <- ((length(unique(diag(SigmaC))) == 1) && (length(unique(SigmaC[lower.tri(SigmaC)])) == 1))
        multvarcompE <- isTRUE(all.equal(SigmaE, diag(diag(SigmaE)))) # Is a diagonal matrix ?
        multvarcompC <- isTRUE(all.equal(SigmaC, diag(diag(SigmaC)))) # Is a diagonal matrix ?
        sphericityE <- ((length(unique(diag(SigmaE))) == 1) && multvarcompE)
        sphericityC <- ((length(unique(diag(SigmaC))) == 1) && multvarcompC)

        if (sphericityE && sphericityC) {
            matrix.type <- 1
            if (display.type) print("SigmaE and SigmaC both have a sphericity structure.")
        } else if (multvarcompE && multvarcompC) {
            matrix.type <- 2
            if (display.type) print("SigmaE and SigmaC both have a variance components structure.")
        } else if (compoundE && compoundC) {
            matrix.type <- 3
            if (display.type) print("SigmaE and SigmaC both have a compound symmetry structure.")
        } else if ((length(unique(SigmaE[lower.tri(SigmaE)])) == 1) && (length(unique(SigmaC[lower.tri(SigmaC)])) == 1)) {
            matrix.type <- 4
            if (display.type) print("SigmaE and SigmaC both have a compound symmetry structure with unequal individual variances.")
        } else {
            matrix.type <- 5
            if (display.type) print("Either SigmaE or SigmaC is unstructured.")
        }
        return(matrix.type)
}

.tr <- function(M) sum(diag(M))
        
df.compute <- function(nE, nC, SigmaE = NULL, SigmaC = NULL, matrix.type = NULL, equalSigmas = NULL, m = NULL) {
    if (is.null(matrix.type)) {
        if (!is.null(equalSigmas)) stop("'equalSigmas' should be NULL.")
        m <- nrow(SigmaE)
        equalSigmas <- isTRUE(all.equal(SigmaE, SigmaC))
        matrix.type <- matrix.type.compute(SigmaE, SigmaC, display.type = FALSE)
    } else {
        if (is.null(m)) stop("'m' should be given a value.")
    }
    switch (matrix.type,
            "1" = {
                if (equalSigmas) {
                    f <- nE + nC - 2
                } else {
                    f <- nE + nC - 2
                }
            },
            "2" = {
                if (equalSigmas) {
                    f <- nE + nC - 2
                } else {
                    f <- nE + nC - 2
                }
            },
            "3" = {
                if (equalSigmas) {
                    f <- m * (nE + nC - 2)
                } else {
                    rhoE <- SigmaE[1, 2] / sqrt(SigmaE[1, 1] * SigmaE[2, 2])
                    rhoC <- SigmaE[1, 2] / sqrt(SigmaE[1, 1] * SigmaE[2, 2])
                    if (isTRUE(all.equal(rhoE, rhoC))) {
                        rho <- rhoE
                    } else {
                        stop("Both variance matrices should have the same rho in the compound case.")
                    }
                    sigma2E <- SigmaE[1, 1]
                    sigma2C <- SigmaC[1, 1]
                    sigma2EC <- sigma2E / nE + sigma2C / nC
                    
                    f <- m * sigma2EC ^ 2 / ((sigma2E / nE) ^ 2 / (nE - 1) + (sigma2C / nC) ^ 2 / (nC - 1))
                }
            },
            "4" = {
                if (equalSigmas) {
                    f <- nE + nC - 2
                } else {
                    SigmaEC <- SigmaE / nE + SigmaC / nC
                    f <- (.tr(SigmaEC %*% SigmaEC) + (.tr(SigmaEC)) ^ 2) / ((.tr(SigmaE %*% SigmaE / (nE ^ 2)) + (.tr(SigmaE / nE)) ^ 2) / (nE - 1) + (.tr(SigmaC %*% SigmaC / (nC ^ 2)) + (.tr(SigmaC / nC)) ^ 2) / (nC - 1))
                }
            },
            "5" = {
                if (equalSigmas) {
                    f <- nE + nC - 2
                } else {
                    SigmaEC <- SigmaE / nE + SigmaC / nC
                    f <- (.tr(SigmaEC %*% SigmaEC) + (.tr(SigmaEC)) ^ 2) / ((.tr(SigmaE %*% SigmaE / (nE ^ 2)) + (.tr(SigmaE / nE)) ^ 2) / (nE - 1) + (.tr(SigmaC %*% SigmaC / (nC ^ 2)) + (.tr(SigmaC / nC)) ^ 2) / (nC - 1))
                }
            })
    return(as.integer(round(f)))
}

