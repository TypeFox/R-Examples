spreml <-
function (formula, data, index = NULL, w, w2=w, lag = FALSE,
          errors = c("semsrre", "semsr", "srre", "semre",
          "re", "sr", "sem","ols", "sem2srre", "sem2re"),
          pvar = FALSE, hess = FALSE, quiet = TRUE,
          initval = c("zeros", "estimate"),
          x.tol = 1.5e-18, rel.tol = 1e-15, ...)
{
    ## mod from spreml5.R to include experimental versions of functions in /optim

    #require(maxLik)

    trace <- as.numeric(!quiet)
    if (pvar)
        print("<implement pvar>")
    if (!is.null(index)) {
        #require(plm)
        data <- plm.data(data, index)
    }
    index <- data[, 1]
    tindex <- data[, 2]
    cl <- match.call()
    #require(nlme)
    if (!is.matrix(w)) {
        if ("listw" %in% class(w)) {
 #           require(spdep)
 #           w <- listw2mat(w)
        }
        else {
            stop("w has to be either a 'matrix' or a 'listw' object")
        }
    }
    if (dim(data)[[1]] != length(index))
        stop("Non conformable arguments")
#    X <- model.matrix(formula, data = data)
#    y <- model.response(model.frame(formula, data = data))
    pmod <- plm(formula, data, model="pooling")
    X <- model.matrix(pmod)
    y <- pmodel.response(pmod)
    names(index) <- row.names(data)
    ind <- index[which(names(index) %in% row.names(X))]
    tind <- tindex[which(names(index) %in% row.names(X))]
    oo <- order(tind, ind)
    X <- X[oo, , drop=FALSE]
    y <- y[oo]
    ind <- ind[oo]
    tind <- tind[oo]
    n <- length(unique(ind))
    k <- dim(X)[[2]]
    t <- max(tapply(X[, 1], ind, length))
    nT <- length(ind)
#    if (dim(w)[[1]] != n)
#        stop("Non conformable spatial weights") # temporary: adapt to listw or mat
    balanced <- n * t == nT
    if (!balanced)
        stop("Estimation method unavailable for unbalanced panels")
    sv.length <- switch(match.arg(errors), semsrre = 3, semsr = 2,
          srre = 2, semre = 2, re = 1, sr = 1, sem = 1, ols = 0,
          sem2srre = 3, sem2re = 2)
    errors. <- match.arg(errors)
    if (is.numeric(initval)) {
        if (length(initval) != sv.length) {
            stop("Incorrect number of initial values supplied for error vcov parms")
        }
        coef0 <- initval
    }
    else {
        switch(match.arg(initval), zeros = {
            coef0 <- rep(0, sv.length)
        }, estimate = {
            if (nchar(errors.) < 4) {
                stop("Pre-estimation of unique vcov parm is meaningless: \n please select (default) option 'zeros' or supply a scalar")
            }
            coef0 <- NULL
            if (grepl("re", errors.)) {
                REmodel <- REmod(X, y, ind, tind, n, k, t, nT,
                  w, coef0 = 0, hess = FALSE, trace = trace,
                  x.tol = 1.5e-18, rel.tol = 1e-15, ...)
                coef0 <- c(coef0, REmodel$errcomp)
            }
            if (grepl("sr", errors.)) {
                ARmodel <- ssrmod(X, y, ind, tind, n, k, t, nT,
                  w, coef0 = 0, hess = FALSE, trace = trace,
                  x.tol = 1.5e-18, rel.tol = 1e-15, ...)
                coef0 <- c(coef0, ARmodel$errcomp)
            }
            if (grepl("sem", errors.)) {
                SEMmodel <- semmod(X, y, ind, tind, n, k, t,
                  nT, w, coef0 = 0, hess = FALSE, trace = trace,
                  x.tol = 1.5e-18, rel.tol = 1e-15, ...)
                coef0 <- c(coef0, SEMmodel$errcomp)
            }
        })
    }
    if (lag) {
        est.fun <- switch(match.arg(errors), semsrre = {
           saremsrREmod
        }, sem2srre = {
           sarem2srREmod
        }, semsr = {
            saremsrmod
        }, srre = {
            sarsrREmod
        }, semre = {
            saremREmod
        }, re = {
            sarREmod
        }, sr = {
            sarsrmod
        }, sem = {
            saremmod
        }, ols = {
            sarmod
        }, sem2re = {
            sarem2REmod
        })
	  coef0 <- c(coef0, 0)
    } else {
        est.fun <- switch(match.arg(errors), semsrre = {
            semsrREmod
        }, sem2srre = {
           sem2srREmod
        }, semsr = {
            semsrmod
        }, srre = {
            ssrREmod
        }, semre = {
            semREmod
        }, re = {
            REmod
        }, sr = {
            ssrmod
        }, sem = {
            semmod
        }, ols = {
            olsmod
            #stop("No lag and no covariance parameters selected: use lm()")
        }, sem2re = {
            sem2REmod
        })
        arcoef <- NULL
    }
    RES <- est.fun(X, y, ind, tind, n, k, t, nT, w = w, w2 = w2,
                   coef0 = coef0, hess = hess, trace = trace,
                   x.tol = x.tol, rel.tol = rel.tol, ...)
    y.hat <- as.vector(X %*% RES$betas)
    res <- y - y.hat
    nam.rows <- dimnames(X)[[1]]
    names(y.hat) <- nam.rows
    names(res) <- nam.rows
    model.data <- data.frame(cbind(y, X[, -1]))
    dimnames(model.data)[[1]] <- nam.rows
    type <- "random effects ML"
    sigma2v <- RES$sigma2
    sigma2mu <- if(is.null(RES$errcomp["phi"])) {0} else {
      as.numeric(sigma2v*RES$errcomp["phi"])
    }
    sigma2.1 <- sigma2mu + sigma2v
    sigma2 <- list(one = sigma2.1, idios = sigma2v, id = sigma2mu)
    spmod <- list(coefficients = RES$betas, arcoef = RES$arcoef,
        errcomp = RES$errcomp, vcov = RES$covB, vcov.arcoef = RES$covAR,
        vcov.errcomp = RES$covPRL, residuals = res, fitted.values = y.hat,
        sigma2 = sigma2, model = model.data, type = type, call = cl,
        errors = errors, logLik = RES$ll)
    class(spmod) <- "splm"
    return(spmod)
}
