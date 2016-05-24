"as.bic.glm" <-
function (g, ...) 
{
    UseMethod("as.bic.glm")
}


summary.glib <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, index.phi=1, ...) 
{
   x<- object
   gobj<- as.bic.glm.glib(x, index.phi=index.phi) 
   summary(gobj, n.models=n.models, digits=digits, conditional=conditional, display.dropped=FALSE, ...)

}

"as.bic.glm.glib" <-
function (g, index.phi = 1, ...) 
{
    probne0 <- 100*(1-g$posterior$prob0[, index.phi])
    postmean <- c(NA, g$posterior$mean[, index.phi])
    postsd <- c(NA, g$posterior$sd[, index.phi])
    size <- g$bf$npar
    bic <- g$bf$twologB10[, index.phi]
    postprob <- g$bf$postprob[, index.phi]
    deviance <- g$bf$deviance
    disp <- 1
    which <- g$models == 1
    family <- g$inputs$error
    if (!is.null(g$glm.out)) 
        reduced <- g$glm.out$reduced
    else reduced <- FALSE
    if (!is.null(g$glm.out)) 
        dropped <- g$glm.out$dropped
    else dropped <- NULL
    call <- g$call
    n.models <- length(postprob)
    n.vars <- length(probne0)
    factor.type <- FALSE
    x <- g$inputs$x
    y <- g$input$y

    if (is.null(colnames(x)))
	colnames(x)<- paste("X",1:ncol(x),sep="")

    nms <- c("Intercept", colnames(x))
    assign <- list()
    for (i in 1:length(nms)) assign[[i]] <- i
    names(assign) <- nms
    output.names <- list()
    for (i in 2:length(nms)) output.names[[i - 1]] <- NA
    names(output.names) <- nms[-1]
    namesx <- colnames(x)
    if (!is.null(g$glm.out)) 
        design <- g$glm.out$design
    else design <- NULL
    label <- rep(NA, times = n.models)
    for (i in 1:n.models) label[i] <- paste(namesx[which[i, ]], 
        sep = "", collapse = ",")
    mle <- matrix(rep(0, times = n.models * (n.vars + 1)), ncol = n.vars + 
        1)
    for (i in 1:n.models) mle[i, c(FALSE, which[i, ])] <- rbind(g$posterior.bymodel$mean[[i]][-1, 
        index.phi])
    se <- matrix(rep(0, times = n.models * (n.vars + 1)), ncol = n.vars + 
        1)
    for (i in 1:n.models) se[i, c(FALSE, which[i, ])] <- rbind(g$posterior.bymodel$sd[[i]][-1, 
        index.phi])

    mle[,1]<- rep(NA, times=nrow(mle))
    se[,1]<- rep(NA, times=nrow(se))

    nests <- n.vars + 1
    CSDbi <- rep(NA, n.vars)
    CEbi <- CSDbi
    for (i in (1:ncol(x))) {
        sel <- which[, i]
        if (sum(sel) > 0) {
            cpp <- rbind(postprob[sel]/sum(postprob[sel]))
            CEbi[i + 1] <- as.numeric(cpp %*% mle[sel, i + 1])
            CSDbi[i + 1] <- sqrt(cpp %*% (se[sel, i + 1]^2) + 
                cpp %*% ((mle[sel, i + 1] - CEbi[i + 1])^2))
        }
    }
    prior.model.weights <- NULL
    prior.param <- NULL
    result <- list(postprob = postprob, label = label, deviance = deviance, 
        size = size, bic = bic, prior.param = prior.param, prior.model.weights = prior.model.weights, 
        family = family, disp = disp, which = which, probne0 = probne0, 
        postmean = postmean, postsd = postsd, condpostmean = CEbi, 
        condpostsd = CSDbi, mle = mle, se = se, namesx = namesx, 
        reduced = reduced, dropped = dropped, call = call, n.models = n.models, 
        n.vars = n.vars, nests = nests, output.names = output.names, 
        assign = assign, factor.type = factor.type, design = design, 
        x = x, y = y)
    class(result) <- "bic.glm"
    return(result)
}

"glib" <-
function (x, ...) 
UseMethod("glib")


"glib.data.frame" <-
function (x, y, n = rep(1, nrow(x)), error = "poisson", link = "log", 
    scale = 1, models = NULL, phi = c(1, 1.65, 5), psi = 1, nu = 0, 
    pmw = rep(1, nrow(models)), glimest = TRUE, glimvar = FALSE, output.priorvar = FALSE, 
    post.bymodel = TRUE, output.postvar = FALSE, priormean = NULL, priorvar = NULL, 
    nbest = 150, call = NULL, ...) 
{
    glim.pmean.A <- function(x, y, n, error, link, scale, nu, 
        prior.spec) {
        glimot <- glim(x, y, n, error = error, link = link, scale = scale)
        coef <- matrix(glimot$coef, ncol = 1)
        eta <- cbind(rep(1, nrow(x)), x) %*% coef
        if (link == "identity") {
            mu <- eta
            deta <- rep(1,length(mu))
        }
        if (link == "log") {
            mu <- exp(eta)
            deta <- 1/mu
        }
        if (link == "logit") {
            mu <- exp(eta)/(1 + exp(eta)) * n
            deta <- 1/mu + 1/(n - mu)
        }
        if (link == "probit") {
            mu <- pnorm(eta)
            deta <- 1/(dnorm(qnorm(mu)))
        }
        if (link == "sqrt") {
            mu <- eta^2
            deta <- 1/(2 * sqrt(mu))
        }
        if (link == "inverse") {
            mu <- 1/eta
            deta <- -1/(mu^2)
        }
        if (link == "loglog") {
            mu <- 1 - exp(-exp(eta))
            deta <- -1/((1 - mu) * log(1 - mu))
        }
        if (error == "gaussian") 
            v <- rep(1,length(mu))
        if (error == "poisson") 
            v <- mu
        if (error == "binomial") 
            v <- mu * (1 - mu/n)
        if (error == "gamma") 
            v <- mu^2
        if (error == "inverse gaussian") 
            v <- mu^3
        w <- 1/(deta^2 * v)
        w <- w/sum(w)
        z <- eta + (y - mu) * deta
        xbar <- t(x) %*% w
        xbar2 <- t(x^2) %*% w
        s2x <- xbar2 - xbar^2
        zbar <- as.numeric(weighted.mean(z, w))
        zbar2 <- t(z^2) %*% w
        if (error == "binomial") 
            s2z <- weighted.mean((eta - zbar)^2 + 2 * (eta - 
                zbar) * (deta * n) * (y - mu)/n + (deta * n)^2 * 
                (y - 2 * mu * y/n + mu^2/n)/n, w)
        else s2z <- zbar2 - zbar^2
        varz <- as.numeric(s2z)
        if (is.numeric(prior.spec$mean)) 
            pmean <- prior.spec$mean
        else pmean <- c((zbar + nu * sqrt(varz)), rep(0, ncol(x)))
        if (is.numeric(prior.spec$var)) 
            A <- prior.spec$var
        else {
            A <- diag(as.vector(1/sqrt(s2x)), nrow = length(s2x))
            A <- rbind(t(-xbar/sqrt(s2x)), A)
            Acol1 <- c(1, rep(0, ncol(x)))
            A <- sqrt(varz) * cbind(Acol1, A)
        }
        list(x = x, y = y, n = n, error = error, link = link, 
            pmean = pmean, A = A, scale = glimot$scale, varz = varz)
    }
    glim.pvar <- function(model, phi = 1, psi = 1, Aot, prior.spec) {
        model1 <- c(1, model + 1)
        if (is.numeric(prior.spec$var)) {
            sigma <- Aot$A
            sig11 <- sigma[model1, model1]
            sig12 <- sigma[model1, -model1]
            sig21 <- sigma[-model1, model1]
            sig22 <- sigma[-model1, -model1]
            if (length(model1) != ncol(Aot$A)) 
                pvar <- sig11 - sig12 %*% solve(sig22) %*% sig21
            else pvar <- sigma
        }
        else {
            A <- Aot$A[model1, model1]
            p <- length(model)
            B <- diag(c(psi^2, rep(phi^2, p)), nrow = (p + 1))
            pvar <- A %*% B %*% t(A)
        }
        list(model = model, phi = phi, psi = psi, pvar = pvar)
    }
    E1 <- function(model, phi, psi, glimot1, Aot, prior.spec) {
        V <- glimot1$var
        A <- solve(V)
        logdetV <- sum(log(eigen(V)$values))
        pmean <- matrix(Aot$pmean[c(1, model + 1)], ncol = 1)
        thetahat <- matrix(glimot1$coef, ncol = 1)
        nphi <- length(phi)
        app1 <- rep(0, nphi)
        npar <- length(model) + 1
        pvar <- array(rep(0, npar * npar * nphi), dim = c(npar, 
            npar, nphi))
        for (j in (1:nphi)) {
            phij <- phi[j]
            pvar[, , j] <- glim.pvar(model, phij, psi, Aot, prior.spec)$pvar
            B <- solve(pvar[, , j])
            C <- solve(A + B)
            I <- diag(rep(1, length(model) + 1), nrow = length(model) + 
                1)
            logdetW <- sum(log(eigen(pvar[, , j])$values))
            F1 <- B %*% C %*% A %*% C %*% B + (t(I - C %*% B)) %*% 
                B %*% (I - C %*% B)
            app1[j] <- -(t(thetahat - pmean)) %*% F1 %*% (thetahat - 
                pmean) - logdetW - sum(log(eigen(A + B)$values))
        }
        list(model = model, phi = phi, psi = psi, app1 = app1)
    }
    E0 <- function(psi = 1, glimot0, Aot) {
        V <- glimot0$var
        A <- 1/V
        logV <- log(V)
        pmean <- Aot$pmean[1]
        thetahat <- glimot0$coef
        pvar <- Aot$varz * psi^2
        B <- 1/pvar
        C <- 1/(A + B)
        logW <- log(pvar)
        F0 <- B^2 * C^2 * A + (1 - C * B)^2 * B
        app1 <- F0 * (thetahat - pmean)^2 - logW - log(A + B)
        list(psi = psi, app1 = app1)
    }
    glim <- function(x, y, n, error = "gaussian", link = "identity", 
        wt = 1, resid = "none", init, intercept = TRUE, scale, offset = 0, 
        sequence, eps = 1e-04, iter.max = 10) {
        error.int <- charmatch(error, c("gaussian", "poisson", 
            "binomial", "gamma", "inverse gaussian"))
        if (is.na(error.int)) 
            stop("Invalid error type")
        else if (error.int == 0) 
            stop("Ambiguous error type")
        resid.int <- charmatch(resid, c("none", "pearson", "deviance"))
        if (is.na(resid.int)) 
            stop("Invalid residual type")
        if (!missing(scale) && !is.numeric(scale)) {
            temp <- charmatch(scale, c("pearson", "deviance"))
            if (is.na(temp)) 
                stop("Invalid scale option")
        }
        if (!(is.numeric(y) && is.numeric(x))) 
            stop("Invalid y or x values")
        x <- as.matrix(x)
        nn <- nrow(x)
        nvar <- ncol(x)
        if (length(y) != nn) 
            stop("Dimensions of x and y don't match")
        if (!missing(wt)) {
            if (length(wt) != nn) 
                stop("Dimensions of x and wt don't match")
            if (any(wt < 0, na.rm = TRUE)) 
                stop("Weights must be >=0")
        }
        if (!missing(offset) && (length(offset) != nn)) 
            stop("Dimensions of x and offset don't match")
        if (error.int != 3) 
            n <- rep(1, nn)
        else {
            if (missing(n)) 
                stop("Binomial fits require the n vector")
            else if (length(n) != nn) 
                stop("Length of n vector is incorrect")
            n <- as.vector(n)
        }
        nomiss <- !(is.na(y) | is.na(wt) | is.na(n) | is.na(x %*% 
            rep(1, nvar)))
        if (sum(nomiss) < nn) {
            warning(paste(nn - sum(nomiss), "deleted due to missing values"))
            x <- x[nomiss, , drop = FALSE]
            y <- y[nomiss]
            n <- n[nomiss]
            if (!missing(wt)) 
                wt <- wt[nomiss]
            if (!missing(offset)) 
                offset <- offset[nomiss]
            nn <- sum(nomiss)
        }
        if (missing(sequence)) 
            sequence <- c(0, nvar)
        else {
            if (max(sequence) > nvar) 
                stop("Invalid sequence argument")
            if (min(sequence) < 0) 
                stop("Invalid sequence argument")
            if (any(diff(sequence) <= 0)) 
                stop("Invalid sequence argument")
        }
        xn <- dimnames(x)[[2]]
        if (is.null(xn)) 
            xn <- paste("X", 1:nvar, sep = "")
        if (intercept) {
            x <- cbind(rep(1, nn), x)
            xn <- c("Intercept", xn)
            nvar <- nvar + 1
            sequence <- sequence + 1
        }
        dimnames(x) <- list(dimnames(x)[[1]], xn)
        if (error.int != 1 && any(y < 0)) 
            stop("Illegal y values")
        if (error.int >= 4 && any(y == 0)) 
            stop("Illegal y values")
        if (error.int == 3) {
            y <- y/n
            if (any(y > 1)) 
                stop("Illegal y values")
        }
        var <- switch(error.int, function(x, n) rep(1, length(x)), 
            function(x, n) x, function(x, n) (x - x^2)/n, function(x, 
                n) x^2, function(x, n) x^3)
        dev <- switch(error.int, function(y, mu, n) (y - mu)^2, 
            function(y, mu, n) 2 * (y * log(ifelse(y == 0, 1, 
                y/mu)) - (y - mu)), function(y, mu, n) 2 * n * 
                (y * log(ifelse(y == 0, 1, y/mu)) + (1 - y) * 
                  log(ifelse(y == 1, 1, (1 - y)/(1 - mu)))), 
            function(y, mu, n) -2 * (log(y/mu) - (y - mu)/mu), 
            function(y, mu, n) (y - mu)^2/(y * mu^2))
        link.int <- charmatch(link, c("identity", "log", "logit", 
            "sqrt", "inverse", "probit", "loglog"))
        if (is.na(link.int)) 
            stop("Invalid link type")
        else if (link.int == 0) 
            stop("Ambiguous link type")
        f <- switch(link.int, function(x) x, function(x) exp(x), 
            function(x) {
                z <- exp(ifelse(x > 80, 80, ifelse(x < -80, -80, 
                  x)))
                z/(z + 1)
            }, function(x) x^2, function(x) 1/x, pnorm, function(x) {
                z <- exp(ifelse(x > 80, 80, ifelse(x < -80, -80, 
                  x)))
                1 - exp(-z)
            })
        deriv <- switch(link.int, function(x) rep(1, length(x)), 
            function(x) exp(x), function(x) exp(x)/((1 + exp(x))^2), 
            function(x) 2 * x, function(x) -1/(x^2), function(x) 0.3989422 * 
                exp(-0.5 * x^2), function(x) exp(x) * exp(-exp(x)))
        temp.y <- switch(error.int, y, y + 0.01, (n * y + 0.5)/(n + 
            1), y, y)
        eta <- switch(link.int, temp.y, log(temp.y), log(temp.y/(1 - 
            temp.y)), sqrt(temp.y), 1/temp.y, qnorm(temp.y), 
            log(-log(1 - temp.y)))
        if (sum(is.na(eta)) > 0) 
            stop("Illegal value(s) of y for this link function")
        if (missing(init)) {
            tempd <- deriv(eta)
            z <- eta + (y - temp.y)/tempd - offset
            w <- sqrt((wt * tempd^2)/var(temp.y, n))
            fit <- qr(x * w)
            if (fit$rank < nvar) 
                stop("X matrix is not full rank")
            init.qr <- fit$qr[1:nvar, ]
            init.qty <- qr.qty(fit, z * w)[1:nvar]
        }
        else if (length(init) != nvar) 
            stop("Argument 'init' is the wrong length")
        if (!missing(wt)) 
            nn <- sum(wt > 0)
        deviance <- df <- NULL
        if (sequence[1] == 1 && intercept && missing(offset)) {
            if (missing(wt)) 
                yhat <- sum(y * n)/sum(n)
            else yhat <- sum((y * n * wt)[wt > 0])/sum((n * wt)[wt > 
                0])
            if ((yhat > 1) && (link.int == 3 || link.int == 6 || 
                link.int == 7)) 
                yhat <- 1
            if ((yhat < 0) && link.int > 1 && link.int != 5) 
                yhat <- 0
            deviance <- sum(dev(y, yhat, n) * wt)
            df <- nn - 1
            sequence <- sequence[-1]
        }
        else if (sequence[1] == 0) {
            deviance <- sum(dev(y, f(offset), n) * wt)
            df <- nn
            sequence <- sequence[-1]
        }
        for (modnum in sequence) {
            model <- 1:modnum
            if (missing(init)) 
                coef <- backsolve(init.qr, init.qty, k = modnum)[model]
            else coef <- init[model]
            tempx <- as.matrix(x[, model, drop = FALSE])
            nvar <- length(model)
            eta <- c(tempx %*% coef)
            yhat <- f(eta + offset)
            olddev <- sum(dev(y, yhat, n) * wt)
            fini <- FALSE
            for (i in seq(length = iter.max)) {
                tempd <- deriv(eta + offset)
                if (any(is.na(yhat) | is.na(tempd))) {
                  warning("Coef vector is diverging, iteration truncated")
                  break
                }
                leave.in <- ((1 + tempd) != 1)
                z <- eta + (y - yhat)/tempd
                w <- sqrt((wt * tempd^2)/var(yhat, n))
                fit <- qr(tempx[leave.in, ] * w[leave.in])
                if (fit$rank < nvar) {
                  warning(paste("Weighted X matrix no longer of full rank at iteration", 
                    i))
                  break
                }
                coef <- qr.coef(fit, (z * w)[leave.in])
                eta <- c(tempx %*% coef)
                yhat <- f(eta + offset)
                newdev <- sum(dev(y, yhat, n) * wt)
                if (abs((olddev - newdev)/(newdev + 1)) < eps) {
                  fini <- TRUE
                  break
                }
                else olddev <- newdev
            }
            if (fini == FALSE) 
                warning(paste("Model with", length(model), "variables did not converge"))
            deviance <- c(deviance, newdev)
            df <- c(df, nn - nvar)
        }
        coef <- as.vector(coef)
        names(coef) <- dimnames(tempx)[[2]]
        if (missing(scale)) 
            scale <- switch(error.int, newdev/(nn - nvar), 1, 
                1, newdev/(nn - nvar), newdev/(nn - nvar))
        else if (!is.numeric(scale)) 
            scale <- switch(charmatch(scale, c("pearson", "deviance")), 
                sum((y - yhat)^2/(var(yhat) * (nn - nvar))), 
                newdev/(nn - nvar))
        varmat <- backsolve(fit$qr[1:nvar,,drop=F], diag(nvar))
        varmat <- varmat %*% t(varmat)
        temp <- list(coef = coef, var = varmat * scale, deviance = deviance, 
            df = df, scale = scale, intercept = intercept)
        if (resid.int == 1) 
            temp
        else if (resid.int == 3) {
            resid <- rep(NA, length(nomiss))
            resid.good <- sqrt(dev(y, yhat, n))
            resid[nomiss] <- ifelse(y > yhat, resid.good, -resid.good)
            c(temp, list(resid = resid))
        }
        else {
            resid <- rep(NA, length(nomiss))
            resid[nomiss] <- (y - yhat)/ifelse(y == yhat, 1, 
                sqrt(var(yhat, n)))
            c(temp, list(resid = resid))
        }
    }


    ##### start of function proper #####

    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    glm.out <- NULL
    prior.spec <- list(mean = FALSE, var = FALSE)
    if (!is.null(priormean)) 
        prior.spec$mean <- priormean
    if (!is.null(priorvar)) {
        prior.spec$var <- priorvar
        phi <- 1
    }
    if (is.null(models)) {
        glm.out <- bic.glm(x = x, y = y, glm.family = error, 
            nbest = nbest, dispersion = NULL, occam.window = FALSE, strict = FALSE, 
            factor.type = FALSE, ...)

        tempdata <- data.frame(glm.out$x, glm.out$y)
        mm <- model.matrix(terms(y ~ ., data = tempdata), data = tempdata)
        df <- data.frame(mm)
        df <- df[, -ncol(df)]
        df <- df[, -1]
        xin <- as.matrix(df)
        y <- glm.out$y
        models <- glm.out$which + 0
        whichin <- matrix(rep(NA, times = nrow(models) * ncol(xin)), 
            ncol = ncol(xin))
        for (i in 1:ncol(models)) whichin[, glm.out$assign[[i + 
            1]] - 1] <- models[, i]
        models <- whichin
        family <- glm.out$family
        if (is.character(family)) 
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family)) 
            family <- family()
        error <- family$family
        link <- family$link
	x<- xin

    }
    if (is.data.frame(x)) {
        nx <- names(x)
        x <- as.matrix(x)
        dimnames(x) <- list(NULL, nx)
    }
    glimot0 <- glim(rep(1, length(y)), y, n, error = error, link = link, 
        intercept = FALSE, scale = scale)
    Aot <- glim.pmean.A(x, y, n, error, link, scale, nu, prior.spec)
    scale <- Aot$scale
    pmean1 <- Aot$pmean
    E0ot <- E0(psi, glimot0, Aot)
    nmodel <- nrow(models)
    nphi <- length(phi)
    chi2 <- rep(0, nmodel)
    npar <- rep(0, nmodel)
    app1 <- matrix(rep(0, nmodel * nphi), ncol = nphi, nrow = nmodel)
    glim.coef <- as.list(rep(0, nmodel))
    glim.var <- as.list(rep(0, nmodel))
    glim.se <- as.list(rep(0, nmodel))
    prior.var <- as.list(rep(0, nmodel))
    postbym.mean <- as.list(rep(0, nmodel))
    postbym.sd <- as.list(rep(0, nmodel))
    postbym.var <- as.list(rep(0, nmodel))
    deviance <- rep(0, nmodel)
    df <- rep(0, nmodel)
    prob0 <- matrix(rep(0, (ncol(x) * nphi)), ncol = nphi)
    post.mean <- matrix(rep(0, (ncol(x) * nphi)), ncol = nphi)
    post.sd <- matrix(rep(0, (ncol(x) * nphi)), ncol = nphi)
    postprob <- matrix(rep(0, (nmodel * nphi)), ncol = nphi)
    prior.var <- as.list(rep(0, nmodel))
    for (i in (1:nmodel)) {
        if (sum(models[i, ]) == 0) {
            npar[i] <- 1
            prior.var[[i]] <- array(rep(0, nphi), dim = c(1, 
                1, nphi))
            for (j in (1:nphi)) prior.var[[i]][, , j] <- psi^2 * 
                Aot$varz
        }
        else {
            model <- (1:ncol(x))[models[i, ] == 1]
            npar[i] <- length(model) + 1
            prior.var[[i]] <- array(rep(0, npar[i] * npar[i] * 
                nphi), dim = c(npar[i], npar[i], nphi))
            for (j in (1:nphi)) prior.var[[i]][, , j] <- glim.pvar(model, 
                phi[j], psi, Aot, prior.spec)$pvar
        }
    }
    for (i in (1:nmodel)) 
    {
        if (sum(models[i, ]) == 0) 
	{
            chi2[i] <- 0
            npar[i] <- 1
            app1[i, ] <- 0
            betahat <- glimot0$coef
            glim.coef[[i]] <- betahat
            V <- glimot0$var
            glim.var[[i]] <- V
            glim.se[[i]] <- sqrt(diag(glim.var[[i]]))
            deviance[i] <- glimot0$deviance[2]
            df[i] <- glimot0$df[2]
        }
        else 
	{
            model <- (1:ncol(x))[models[i, ] == 1]
            glimot1 <- glim(x[, model], y, n = n, error = error, 
                link = link, scale = scale)
            chi2[i] <- (glimot0$deviance[2] - glimot1$deviance[2])/scale
            npar[i] <- sum(models[i, ]) + 1
            E1ot <- E1(model, phi, psi, glimot1, Aot, prior.spec)
            app1[i, ] <- chi2[i] + (E1ot$app1 - E0ot$app1)
            betahat <- glimot1$coef
            glim.coef[[i]] <- betahat
            V <- glimot1$var
            glim.var[[i]] <- V
            glim.se[[i]] <- sqrt(diag(glim.var[[i]]))
            deviance[i] <- glimot1$deviance[2]
            df[i] <- glimot1$df[2]
        }


        if (sum(models[i, ]) == 0) 
	{
            prior.mean <- pmean1[1]
            postbym.mean[[i]] <- matrix(rep(0, nphi), ncol = nphi)
            postbym.sd[[i]] <- matrix(rep(0, nphi), ncol = nphi)
            postbym.var[[i]] <- array(rep(0, nphi), dim = c(1, 
                1, nphi))
            for (j in (1:nphi)) 
	    {
                W <- prior.var[[i]][, , j]
                postbym.mean[[i]][, j] <- betahat - V %*% solve(V + 
                  W) %*% (betahat - prior.mean)
                postbym.var[[i]][, , j] <- solve(solve(V) + solve(W))
                postbym.sd[[i]][, j] <- sqrt(postbym.var[[i]][, 
                  , j])
            }
        }
        else 
	{
            prior.mean <- c(pmean1[1], rep(0, (npar[i] - 1)))
            postbym.mean[[i]] <- matrix(rep(0, (npar[i] * nphi)), 
                ncol = nphi)
            postbym.sd[[i]] <- matrix(rep(0, (npar[i] * nphi)), 
                ncol = nphi)
            postbym.var[[i]] <- array(rep(0, (npar[i] * npar[i] * 
                nphi)), dim = c(npar[i], npar[i], nphi))
            for (j in (1:nphi)) 
	    {
                W <- prior.var[[i]][, , j]
                postbym.mean[[i]][, j] <- betahat - V %*% solve(V + 
                  W) %*% (betahat - prior.mean)
                postbym.var[[i]][, , j] <- solve(solve(V) + solve(W))
                postbym.sd[[i]][, j] <- sqrt(diag(postbym.var[[i]][, 
                  , j]))
            }
        }
    }

    for (j in (1:nphi)) 
    {
        tmp <- app1[, j] - max(app1[, j])
        tmp <- pmw * exp(0.5 * tmp)
        postprob[, j] <- tmp/sum(tmp)
        for (k in (1:ncol(x))) prob0[k, j] <- sum(postprob[models[, k] == 0, j])
 
        for (k in (1:ncol(x))) 
	{
            tmp1 <- 0
            tmp2 <- 0
            for (i in (1:nmodel)) 
	    {
                if (models[i, k] == 1) 
		{
                    pos <- sum(models[i, (1:k)]) + 1
                    tmp1 <- tmp1 + postbym.mean[[i]][pos, j] * postprob[i, j]
                    tmp2 <- tmp2 + (postbym.sd[[i]][pos, j]^2 + 
                            postbym.mean[[i]][pos, j]^2) * postprob[i, j]
                }
                post.mean[k, j] <- tmp1/(1 - prob0[k, j])
                post.sd[k, j] <- sqrt(tmp2/(1 - prob0[k, j]) - post.mean[k, j]^2)
            }
        }
    }
    inputs <- list(x = x, y = y, n = n, error = error, link = link, 
        models = models, phi = phi, psi = psi, nu = nu)
    if (is.numeric(prior.spec$var)) 
        inputs <- inputs[-(7:8)]
    if (is.numeric(prior.spec$mean)) 
        inputs <- inputs[-match("nu", names(inputs))]
    bf <- list(twologB10 = app1, postprob = postprob, deviance = deviance, 
        df = df, chi2 = chi2, npar = npar, scale = scale)
    posterior <- list(prob0 = prob0, mean = post.mean, sd = post.sd)
    if (glimest == TRUE) {
        if (glimvar == FALSE) 
            glim.est <- list(coef = glim.coef, se = glim.se)
        if (glimvar == TRUE) 
            glim.est <- list(coef = glim.coef, se = glim.se, 
                var = glim.var)
    }
    if (glimest == FALSE) 
        glim.est <- NULL
    if (post.bymodel == TRUE) {
        if (output.postvar == FALSE) 
            posterior.bymodel <- list(mean = postbym.mean, sd = postbym.sd)
        if (output.postvar == TRUE) 
            posterior.bymodel <- list(mean = postbym.mean, sd = postbym.sd, 
                var = postbym.var)
    }
    if (post.bymodel == FALSE) 
        posterior.bymodel <- NULL
    if (output.priorvar == FALSE) 
        prior <- list(mean = pmean1)
    if (output.priorvar == TRUE) 
        prior <- list(mean = pmean1, var = prior.var)
    result <- list(inputs = inputs, bf = bf, posterior = posterior, 
        glim.est = glim.est, posterior.bymodel = posterior.bymodel, 
        prior = prior, x = x, models = models, glm.out = glm.out, 
        call = cl)
    class(result) <- "glib"
    result
}




glib.matrix <- glib.data.frame

glib.bic.glm<- function (x, scale = 1, phi = 1, psi = 1, nu = 0, 
    glimest = TRUE, glimvar = FALSE, output.priorvar = FALSE, post.bymodel = TRUE, 
    output.postvar = FALSE, priormean = NULL, priorvar = NULL, call = NULL, ...) 
{
    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    tempdata <- data.frame(x$x, x$y)
    mm <- model.matrix(terms(y ~ ., data = tempdata), data = tempdata)
    df <- data.frame(mm)
    df <- df[, -ncol(df)]
    df <- df[, -1]
    xin <- as.matrix(df)
    y <- x$y
    models <- x$which + 0
    whichin <- matrix(rep(NA, times = nrow(models) * ncol(xin)), 
        ncol = ncol(xin))
    for (i in 1:ncol(models)) whichin[, x$assign[[i + 1]] - 1] <- models[, 
        i]
    models <- whichin
    family <- x$family
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    error <- family$family
    link <- family$link
    glib(x = xin, y = y, error = error, link = link, models = models, 
        scale = scale, phi = phi, psi = psi, nu = nu, glimest = glimest, 
        glimvar = glimvar, output.priorvar = output.priorvar, 
        post.bymodel = post.bymodel, output.postvar = output.postvar, 
        priormean = priormean, priorvar = priorvar, call = cl, ...)
}
