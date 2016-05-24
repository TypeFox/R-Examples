fitBMAgamma0 <-
function (ensembleData, control = controlBMAgamma0(), exchangeable = NULL) 
{
    ZERO <- 1e-100
    
    powfun <- function(x, power) x^power
    powinv <- function(x, power) x^(1/power)
    if (is.null(exchangeable)) 
        exchangeable <- ensembleGroups(ensembleData)
    if (length(unique(exchangeable)) == length(exchangeable)) 
        exchangeable <- NULL
    if (!(nullX <- is.null(exchangeable))) {
        namX <- as.character(exchangeable)
        uniqueX <- unique(namX)
        nX <- length(uniqueX)
    }
    maxIter <- control$maxIter
    tol <- eps <- control$tol
    ensembleData <- ensembleData[!dataNA(ensembleData, dates = FALSE), ]
    nObs <- dataNobs(ensembleData)
    if (!nObs) 
        stop("no observations")
    obs <- dataVerifObs(ensembleData)
    ensMemNames <- ensembleMembers(ensembleData)
    nForecasts <- length(ensMemNames)
    Y0 <- obs == 0
    n0obs <- sum(Y0)
    if (sum(!Y0) < 2) 
        stop("less than 2 nonzero obs")
    ensembleData <- as.matrix(ensembleForecasts(ensembleData))
    Xvar <- ensembleData[!Y0, ,drop=F]
    ensembleData <- apply(ensembleData, 2, function(x) sapply(x, 
        powfun, power = control$power))
    miss <- as.vector(as.matrix(is.na(ensembleData)))
    logisticFunc <- function(x, y) {
        x <- as.matrix(x)
        n <- ncol(x)
        x <- as.vector(x)
        y <- rep(y, n)
        nax <- is.na(x)
        x <- x[!nax]
        y <- y[!nax]
        if (!all(x == 0) && !all(x != 0)) {
            fit <- glm(y ~ x + (x == 0), family = binomial(logit))
            coefs <- fit$coefficients
        }
        else if (!all(x == 0)) {
            fit <- glm(y ~ x, family = binomial(logit))
            coefs <- c(fit$coefficients, 0)
        }
        else {
            fit <- glm(y ~ 1, family = binomial(logit))
            coefs <- c(fit$coefficients, 0, 0)
        }
        coefs[is.na(coefs)] <- 0
        if (!all(coefs[2:3] == 0) && coefs[2] > 0 && coefs[3] < 
            0) {
           fit <- glm(y ~ 1, family = binomial(logit))
           coefs <- c(fit$coefficients[1], 0, 0)
        }
        else if (coefs[2] > 0) {
            if (all(x != 0)) {
                fit <- glm(y ~ 1, family = binomial(logit))
                coefs <- c(fit$coefficients[1], 0, 0)
            }
            else {
                fit <- glm(y ~ (x == 0), family = binomial(logit))
                coefs <- c(fit$coefficients[1], 0, fit$coefficients[2])
                if (coefs[3] < 0) {
                  fit <- glm(y ~ 1, family = binomial(logit))
                  coefs <- c(fit$coefficients[1], 0, 0)
                 }
            }
        }
        else if (coefs[3] < 0) {
            fit <- glm(y ~ x, family = binomial(logit))
            coefs <- c(fit$coefficients, 0)
            if (coefs[2] >= 0)  {
              fit <- glm(y ~ 1, family = binomial(logit))
              coefs <- c(fit$coefficients[1], 0, 0)
            }
        }
        fit$coefficients <- coefs
        fit
    }
    if (any(Y0)) {
        if (nullX) {
            prob0fit <- apply(as.matrix(ensembleData), 2, logisticFunc, 
                y = Y0)
            prob0coefs <- lapply(prob0fit, function(x) x$coefficients)
            prob0coefs <- as.matrix(data.frame(prob0coefs))
            dimnames(prob0coefs) <- NULL
            PROB0 <- matrix(NA, nObs, nForecasts)
            PROB0[!miss] <- unlist(lapply(prob0fit, function(x) x$fitted.values))
        }
        else {
            prob0coefs <- matrix(NA, 3, nForecasts)
            dimnames(prob0coefs) <- NULL
            PROB0 <- matrix(NA, nObs, nForecasts)
            for (labX in uniqueX) {
                I <- namX == labX
                fit <- logisticFunc(ensembleData[, I, drop = FALSE], 
                  Y0)
                prob0coefs[, I] <- fit$coefficients
                miss <- is.na(ensembleData[, I, drop = FALSE])
                miss <- as.vector(as.matrix(miss))
                PROB0[, I][!miss] <- fit$fitted.values
            }
            p0 <- apply(PROB0, 1, max, na.rm = TRUE) == 1
            if (any(p0 & !Y0)) 
                stop("PROB0 == 1 for a nonzero obs")
        }
        POP <- (1 - PROB0)[!Y0, , drop = FALSE]
        PROB0 <- PROB0[Y0, , drop = FALSE]
    }
    else {
        POP <- matrix(1, nObs, nForecasts)
        PROB0 <- NULL
        prob0coefs <- matrix(0, 3, nForecasts)
        dimnames(prob0coefs) <- NULL
    }
    obs <- sapply(obs, function(x) sapply(x, powfun, power = control$power))
    obs <- obs[!Y0]
    nPrecip <- length(obs)
    ensembleData <- ensembleData[!Y0, ]
    miss <- as.vector(as.matrix(is.na(ensembleData)))
    lmFunc <- function(x, y) {
        beta0 <- min(y)
        x <- as.matrix(x)
        n <- ncol(x)
        x <- as.vector(x)
        nax <- is.na(x)
        x <- x[!nax]
        y <- rep(y, n)[!nax]
        if (all(!x)) {
            fit <- list(coefficients = c(mean(y), 0), fitted.values = rep(mean(y), 
                length(y)))
        }
        else {
            fit <- lm(y ~ x)
            coefs <- fit$coefficients
            if (coefs[1] <= 0) {
                coefs[1] <- beta0
                coefs[2] <- sum((y - beta0) * x)/sum(x * x)
                fit$coefficients <- coefs
                fit$fitted.values <- cbind(1, x) %*% coefs
            }
        }
        fit
    }
    if (nullX) {
        meanFit <- apply(as.matrix(ensembleData), 2, lmFunc, y = obs)
        biasCoefs <- lapply(meanFit, function(x) x$coefficients)
        biasCoefs <- as.matrix(data.frame(biasCoefs))
        MEAN <- matrix(NA, nPrecip, nForecasts)
        MEAN[!miss] <- unlist(lapply(meanFit, function(x) x$fitted.values))
    }
    else {
        biasCoefs <- matrix(NA, 2, nForecasts)
        MEAN <- matrix(NA, nPrecip, nForecasts)
        i <- 1
        for (labX in uniqueX) {
            I <- namX == labX
            fit <- lmFunc(ensembleData[, I, drop = FALSE], obs)
            biasCoefs[, I] <- fit$coefficients
            miss <- is.na(ensembleData[, I, drop = FALSE])
            miss <- as.vector(as.matrix(miss))
            MEAN[, I][!miss] <- fit$fitted.values
        }
    }
    miss <- is.na(Xvar)
    completeDataLLmiss <- function(z, w, m, p0, p1, X, obs, Y0) {
        objective <- function(par) {
            nObs <- length(obs)
            nFor <- ncol(X)
            v <- par[1]^2 + (par[2]^2) * X
            r <- m/v
            W <- matrix(w, nObs, length(w), byrow = TRUE)
            if (any(miss <- is.na(X))) {
                W[miss] <- 0
                W <- sweep(W, MARGIN = 1, FUN = "/", STATS = apply(W, 
                  1, sum))
            }
            q <- array(NA, dim(z))
            q[!Y0, ][!miss] <- dgamma(matrix(obs, nPrecip, nForecasts)[!miss], 
                shape = (r * m)[!miss], rate = r[!miss], log = TRUE)
            Wzero <- (p1 == 0) | (W == 0)
            include <- !miss & !Wzero
            -sum(z[!Y0, ][include] * (q[!Y0, ][include] + log(p1[include] * 
                W[include])))
        }
        objective
    }
    varCoefs <- if (is.null(control$init$varCoefs)) 
        c(1, 1)
    else control$init$varCoefs
    varCoefs <- pmax(varCoefs, 1e-04)
    names(varCoefs) <- names(weights) <- NULL
    weights <- if (is.null(control$init$weights)) 
        1
    else control$init$weights
    if (length(weights) == 1) 
        weights <- rep(weights, nForecasts)
    weights <- weights/sum(weights)
    weights <- pmax(weights, 1e-04)
    weights <- weights/sum(weights)
    if (!is.null(names(weights))) 
        weights <- weights[ensMemNames]
    if (!nullX) {
        for (labX in uniqueX) {
            I <- namX == labX
            weights[I] <- mean(weights[I])
        }
    }
    nIter <- 0
    z <- matrix(1/nForecasts, ncol = nForecasts, nrow = nObs)
    objold <- 0
    newLL <- 0
    MEAN[MEAN <= 0] <- max(min(obs),1.e-4)
    while (TRUE) {
        VAR <- varCoefs[1] + varCoefs[2] * Xvar
        RATE <- MEAN/VAR
        SHAPE <- RATE * MEAN
        z[!Y0, ][!miss] <- dgamma(matrix(obs, nPrecip, nForecasts)[!miss], 
            shape = SHAPE[!miss], rate = RATE[!miss], log = TRUE)
        zmax1 <- apply(z[!Y0, ,drop=F], 1, max, na.rm = TRUE)
        z[!Y0, ] <- sweep(z[!Y0, ,drop=F], MARGIN = 1, FUN = "-", STATS = zmax1)
        z[!Y0, ] <- sweep(POP, MARGIN = 2, FUN = "*", STATS = weights) * 
            exp(z[!Y0, ,drop=F])
        if (!is.null(PROB0)) 
            z[Y0, ] <- sweep(PROB0, MARGIN = 2, FUN = "*", STATS = weights)
        oldLL <- newLL
        newLL <- sum(zmax1 + log(apply(z[!Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        newLL <- newLL + sum(log(apply(z[Y0, , drop = FALSE], 
            1, sum, na.rm = TRUE)))
        z <- z/apply(z, 1, sum, na.rm = TRUE)
        z[z < ZERO] <- 0
        wold <- weights
        zsum2 <- apply(z, 2, sum, na.rm = TRUE)
        weights <- zsum2/sum(zsum2)
        weights[weights < ZERO] <- 0
        if (!nullX) {
            weights <- sapply(split(weights, namX), mean)[namX]
        }
        weps <- max(abs(wold - weights)/(1 + abs(weights)))
        fn <- completeDataLLmiss(z, weights, MEAN, PROB0, POP, 
            Xvar, obs, Y0)
        optimResult = if (is.null(control$optim.control)) {
                       optim(sqrt(varCoefs), fn=fn, method = "BFGS")
                    }
                    else {
                       optim(sqrt(varCoefs), fn=fn, method = "BFGS",
                             control = control$optim.control)
               		         }
        if (optimResult$convergence) 
            warning("optim does not converge")
        varOld <- varCoefs
        varCoefs <- optimResult$par^2
        veps <- max(abs(varOld - varCoefs)/(1 + abs(varCoefs)))
        ERROR <- abs(objold - optimResult$value)/(1 + abs(optimResult$value))
        objold <- optimResult$value
        if (nIter > 0) {
            error <- abs(oldLL - newLL)/(1 + abs(newLL))
            if (error < eps) 
                break
        }
        nIter <- nIter + 1
        if (nIter >= maxIter) 
            break
    }
    if (nIter >= maxIter && error > eps) 
        warning("iteration limit reached")
    dimnames(biasCoefs) <- list(NULL, ensMemNames)
    dimnames(prob0coefs) <- list(NULL, ensMemNames)
    names(weights) <- ensMemNames
    structure(list(prob0coefs = prob0coefs, biasCoefs = biasCoefs, 
        varCoefs = varCoefs, weights = weights, nIter = nIter, 
        loglikelihood = newLL, power = control$power), class = "fitBMAgamma0")
}
