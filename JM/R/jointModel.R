jointModel <- function (lmeObject, survObject, timeVar, parameterization = c("value", "slope", "both"), 
        method = c("weibull-PH-aGH", "weibull-PH-GH", "weibull-AFT-aGH", "weibull-AFT-GH", 
        "piecewise-PH-aGH", "piecewise-PH-GH", "Cox-PH-aGH", "Cox-PH-GH", "spline-PH-aGH", 
        "spline-PH-GH", "ch-Laplace"), interFact = NULL, derivForm = NULL, lag = 0, 
        scaleWB = NULL, CompRisk = FALSE, init = NULL, control = list(), ...) {
    cl <- match.call()
    if (!inherits(lmeObject, "lme"))
        stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1)
        stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct))
        warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!is.null(lmeObject$modelStruct$varStruct))
        warning("variance structure in 'lmeObject' is ignored.\n")        
    if (!inherits(survObject, "coxph") && !inherits(survObject, "survreg"))
        stop("\n'survObject' must inherit from class coxph or class survreg.")
    if (!is.matrix(survObject$x))
        stop("\nuse argument 'x = TRUE' in ", 
            if (inherits(survObject, "coxph")) "'coxph()'." else "'survreg()'.")
    if (length(timeVar) != 1 || !is.character(timeVar))
        stop("\n'timeVar' must be a character string.")
    method. <- match.arg(method)
    method <- switch(method., "weibull-AFT-GH" = , "weibull-AFT-aGH" = "weibull-AFT-GH",
        "weibull-PH-GH" = , "weibull-PH-aGH" = "weibull-PH-GH", 
        "piecewise-PH-GH" = , "piecewise-PH-aGH" = "piecewise-PH-GH",
        "Cox-PH-GH" = , "Cox-PH-aGH" = "Cox-PH-GH", "spline-PH-GH" =, 
        "spline-PH-aGH" = "spline-PH-GH", "ch-Laplace" = "ch-Laplace")
    parameterization <- match.arg(parameterization)
    if (method == "Cox-PH-GH" && !inherits(survObject, "coxph"))
        stop("\nfor 'method = Cox-PH-GH', 'survObject' must inherit from class coxph.")
    if (parameterization %in% c("slope", "both") && method %in% c("Cox-PH-GH", "ch-Laplace"))
        stop("\nthe slope parameterization is not currently available for methods 'Cox-PH-GH' & 'ch-Laplace'.")
    if (parameterization %in% c("slope", "both") && is.null(derivForm)) {
        stop("\nwhen parameterization is 'slope' or 'both' you need to specify the 'derivForm' argument.")
    }
    if (parameterization %in% c("slope", "both") && !is.list(derivForm)) {
        stop("\nthe 'derivForm' argument must be a list with components 'fixed' (a formula),\n\t'indFixed'", 
            "(a numeric vector), 'random' (a formula) and 'indRandom' (a numeric vector).")
    }
    if (!is.null(interFact) && !is.list(interFact)) {
        stop("\nthe 'interFact' argument must be a list -- check the help file for more info.")
    }
    if (!is.null(interFact) && method %in% c("Cox-PH-GH", "ch-Laplace")) {
        stop("\nincluding interaction terms is not currently available for methods 'Cox-PH-GH' & 'ch-Laplace'.")
    }
    if (CompRisk && (method != "spline-PH-GH" || is.null(survObject$strata))) {
        stop("\nto fit a competing risks joint model you must choose as method 'spline-PH-GH'",
            " and include a strata() in the specification of the coxph().")
    }
    # survival process
    formT <- formula(survObject)
    if (inherits(survObject, "coxph")) {
        W <- survObject$x
        keepW <- suppressWarnings(!is.na(survObject$coefficients))
        W <- W[, keepW, drop = FALSE]
        if (CompRisk) {
            nRisks <- length(unique(survObject$strata))  
        } else {
            nRisks <- 1
        }
        surv <- survObject$y
        if (attr(surv, "type") == "right") {
            LongFormat <- FALSE
            Time <- survObject$y[, 1]
            d <- survObject$y[, 2]
        } else if (attr(surv, "type") == "counting") {
            LongFormat <- TRUE
            if (is.null(survObject$model))
                stop("\nplease refit the Cox model including in the ", 
                    "call to coxph() the argument 'model = TRUE'.")
            Time <- survObject$y[, 2]
            d <- survObject$y[, 3]
        }
        idT <- if (!is.null(survObject$model$cluster)) {
            as.vector(unclass(survObject$model$cluster))
        } else {
            if (!CompRisk) seq_along(Time)
            else rep(seq_len(length(Time)/nRisks), each = nRisks)
        }
        idT <- match(idT, unique(idT))
    } else {
        W <- survObject$x[, -1, drop = FALSE]
        Time <- exp(survObject$y[, 1])
        d <- survObject$y[, 2]
        idT <- seq_along(Time)
        LongFormat <- FALSE
        nRisks <- 1
    }
    nT <- length(unique(idT))
    if (LongFormat && is.null(survObject$model$cluster))
        stop("\nuse argument 'model = TRUE' and cluster() in coxph().")
    if (!length(W))
        W <- NULL
    if (sum(d) < 5)
        warning("\nmore than 5 events are required.")
    WintF.vl <- WintF.sl <- as.matrix(rep(1, length(Time)))
    if (!is.null(interFact)) {
        if (!is.null(interFact$value))
            WintF.vl <- if (is.null(survObject$model)) {
                model.matrix(interFact$value, data = interFact$data)
            } else {
                model.matrix(interFact$value, data = survObject$model)
            }
        if (!is.null(interFact$slope))
            WintF.sl <- if (is.null(survObject$model)) {
                model.matrix(interFact$slope, data = interFact$data)
            } else {
                model.matrix(interFact$slope, data = survObject$model)
            }
    }
    # longitudinal process
    id <- as.vector(unclass(lmeObject$groups[[1]]))
    b <- data.matrix(ranef(lmeObject))
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT)
        stop("sample sizes in the longitudinal and event processes differ; ", 
            "maybe you forgot the cluster() argument.\n")
    TermsX <- lmeObject$terms
    data <- lmeObject$data[all.vars(TermsX)]
    data <- data[complete.cases(data), ]
    formYx <- formula(lmeObject)
    mfX <- model.frame(TermsX, data = data)
    X <- model.matrix(formYx, mfX)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])    
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ)
    y.long <- model.response(mfX, "numeric")
    data.id <- data[!duplicated(id), ]
    data.id <- data.id[idT, ]
    if (!timeVar %in% names(data))
        stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
    # check if there are any longitudinal measurements after the event times
    max.timeY <- tapply(data[[timeVar]], id, max)
    max.timeT <- tapply(Time, idT, max)
    if (!all(max.timeT >= max.timeY)) {
        idnams <- factor(lmeObject$groups[[1]])
        stop("\nit seems that there are longitudinal measurements taken after the event times for some subjects ",
            "(i.e., check subject(s): ", paste(levels(idnams)[(max.timeT < max.timeY)], collapse = ", "), ").")
    }
    # extra design matrices for the longitudinal part
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    if (parameterization %in% c("value", "both")) {
        mfX.id <- model.frame(TermsX, data = data.id)
        mfZ.id <- model.frame(TermsZ, data = data.id)
        Xtime <- model.matrix(formYx, mfX.id)
        Ztime <- model.matrix(formYz, mfZ.id)
        long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ])
    }
    if (parameterization %in% c("slope", "both")) {
        mfX.deriv <- model.frame(terms(derivForm$fixed), data = data)
        TermsX.deriv <- attr(mfX.deriv, "terms")
        mfZ.deriv <- model.frame(terms(derivForm$random), data = data)
        TermsZ.deriv <- attr(mfZ.deriv, "terms")
        mfX.deriv.id <- model.frame(TermsX.deriv, data = data.id)
        mfZ.deriv.id <- model.frame(TermsZ.deriv, data = data.id)      
        Xtime.deriv <- model.matrix(derivForm$fixed, mfX.deriv.id)
        Ztime.deriv <- model.matrix(derivForm$random, mfZ.deriv.id)
        Xderiv <- model.matrix(derivForm$fixed, mfX.deriv)
        Zderiv <- model.matrix(derivForm$random, mfZ.deriv)        
        long.deriv <- as.vector(c(Xderiv %*% fixef(lmeObject)[derivForm$indFixed]) + 
            if (length(derivForm$indRandom) > 1 || derivForm$indRandom) 
                rowSums(Zderiv * b[id, derivForm$indRandom, drop = FALSE])
            else
                rep(0, nrow(Zderiv)))
    }
    if (parameterization == "value")
        long.deriv <- NULL
    if (parameterization == "slope")
        long <- NULL        
    # response vectors and design matrices
    y <- list(y = y.long, logT = log(Time), d = d, lag = lag)
    x <- list(X = X, Z = Z, W = W, WintF.vl = WintF.vl, 
        WintF.sl = WintF.sl, idT = idT, nRisks = nRisks)
    x <- switch(parameterization, 
        "value" = c(x, list(Xtime = Xtime, Ztime = Ztime)),
        "slope" = c(x, list(Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
        "both" = c(x, list(Xtime = Xtime, Ztime = Ztime, Xtime.deriv = Xtime.deriv, 
            Ztime.deriv = Ztime.deriv)))
    # control values
    ind.noadapt <- method. %in% c("weibull-AFT-GH", "weibull-PH-GH", "piecewise-PH-GH", 
        "Cox-PH-GH", "spline-PH-GH")
    con <- list(only.EM = FALSE, iter.EM = if (method == "spline-PH-GH") 120 else 50, 
        iter.qN = 350, optimizer = "optim", tol1 = 1e-03, tol2 = 1e-04, 
        tol3 = if (!CompRisk) sqrt(.Machine$double.eps) else 1e-09, numeriDeriv = "fd", eps.Hes = 1e-06, 
        parscale = NULL, step.max = 0.1, backtrackSteps = 2, 
        knots = NULL, ObsTimes.knots = TRUE,
        lng.in.kn = if (method == "piecewise-PH-GH") 6 else 5, ord = 4, 
        equal.strata.knots = TRUE, typeGH = if (ind.noadapt) "simple" else "adaptive", 
        GHk = if (ncol(Z) < 3 && nrow(Z) < 2000) 15 else 9, 
        GKk = if (method == "piecewise-PH-GH" || length(Time) > nRisks*nT) 7 else 15, verbose = FALSE)
    if (method == "Cox-PH-GH") {
        con$only.EM <- TRUE
        con$iter.EM <- 200
        con$GHk <- if (ncol(Z) == 1) 15 else if (ncol(Z) == 2) 11 else 9
    }
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (con$typeGH != "simple" && !"GHk" %in% namc) {
        con$GHk <- if (ncol(Z) <= 3 && nrow(Z) < 2000) 5 else 3
    }
    if (length(noNms <- namc[!namc %in% namC]) > 0) 
        warning("unknown names in 'control': ", paste(noNms, collapse = ", "))
    if (method == "Cox-PH-GH" && !con$only.EM)
        stop("with method 'Cox-PH-GH' only the EM algorithm is used.\n")
    if (method == "Cox-PH-GH" && any(!is.na(match(c("iter.qN", "optimizer"), namc))))
        warning("method 'Cox-PH-GH' uses only the EM algorithm.\n")
    # extra design matrices for 'method = "weibull-AFT-GH"' and 'method = "weibull-PH-GH"'
    # extra design matrices for 'method = "spline-PH-GH"' and 'method = "spline-PH-Laplace"'
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", "spline-PH-GH", "spline-PH-Laplace")) {
        wk <- gaussKronrod(con$GKk)$wk
        sk <- gaussKronrod(con$GKk)$sk
        if (LongFormat) {
            Time0 <- survObject$y[, 1]
            P <- (Time - Time0) / 2
            P1 <- (Time + Time0) / 2
            st <- outer(P, sk) + P1
        } else {
            P <- as.vector(Time)/2
            st <- outer(P, sk + 1)
        }
        dimnames(st) <- names(P) <- NULL
        id.GK <- rep(seq_along(Time), each = con$GKk)
        data.id2 <- data.id[id.GK, , drop = FALSE]
        data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
        if (parameterization %in% c("value", "both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            Xs <- model.matrix(formYx, mfX)
            Zs <- model.matrix(formYz, mfZ)
        }
        if (parameterization %in% c("slope", "both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
        }
        Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        x <- c(x, list(P = P, st = c(t(st)), wk = wk, Ws.intF.vl = Ws.intF.vl, 
            Ws.intF.sl = Ws.intF.sl))
        x <- switch(parameterization,
            "value" = c(x, list(Xs = Xs, Zs = Zs)),
            "slope" = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv)),
            "both" = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv, Xs = Xs, Zs = Zs)))
        if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
            strt <- if (is.null(survObject$strata)) gl(1, length(Time)) else survObject$strata
            nstrt <- length(levels(strt))
            split.Time <- split(Time, strt)
            ind.t <- if (LongFormat) {
                unlist(tapply(idT, idT, 
                    function (x) c(rep(FALSE, length(x) - 1), TRUE)))
            } else {
                rep(TRUE, length(Time))
            }
            kn <- if (con$equal.strata.knots) {
                kk <- if (is.null(con$knots)) {
                    pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
                    pp <- tail(head(pp, -1), -1)
                    quantile(Time[ind.t], pp, names = FALSE)
                } else {
                    con$knots
                }
                kk <- kk[kk < max(Time)]
                rr <- rep(list(sort(c(rep(range(Time, st), con$ord), kk))), nstrt)
                names(rr) <- names(split.Time)
                rr
            } else {
                spt <- if (length(Time) > nT & !CompRisk) 
                    mapply(function (x, y){
                        x[unlist(tapply(y, y, 
                            function (z) c(rep(FALSE, length(z) - 1), TRUE)))]
                    }, split.Time, split(idT, strt), SIMPLIFY = FALSE)
                else split.Time
                lapply(spt, function (t) {
                    kk <- if (is.null(con$knots)) {
                        pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
                        pp <- tail(head(pp, -1), -1)
                        quantile(t, pp, names = FALSE)
                    } else {
                        con$knots
                    }
                    kk <- kk[kk < max(t)]
                    sort(c(rep(range(Time, st), con$ord), kk))
                })
            }
            con$knots <- kn
            W2 <- mapply(function (k, t) splineDesign(k, t, ord = con$ord), kn, 
                            split.Time, SIMPLIFY = FALSE)
            if (any(sapply(W2, colSums) == 0))
                stop("\nsome of the knots of the B-splines basis are set outside the range",
                    "\n   of the observed event times for one of the strata; refit the model", 
                    "\n   setting the control argument 'equal.strata.knots' to FALSE.")
            W2 <- mapply(function (w2, ind) {
                out <- matrix(0, length(Time), ncol(w2))
                out[strt == ind, ] <- w2
                out
            }, W2, levels(strt), SIMPLIFY = FALSE)
            W2 <- do.call(cbind, W2)
            strt.s <- rep(strt, each = con$GKk)
            split.Time <- split(c(t(st)), strt.s)
            W2s <- mapply(function (k, t) splineDesign(k, t, ord = con$ord), 
                kn, split.Time, SIMPLIFY = FALSE)
            W2s <- mapply(function (w2s, ind) {
                out <- matrix(0, length(Time) * con$GKk, ncol(w2s))
                out[strt.s == ind, ] <- w2s
                out
            }, W2s, levels(strt), SIMPLIFY = FALSE)
            W2s <- do.call(cbind, W2s)
            y <- c(y, list(strata = strt))
            x <- c(x, list(W2 = W2, W2s = W2s))
        }
    }
    # extra design matrices for 'method = "piecewise-PH-GH"'
    if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(con$GKk)$wk
        sk <- gaussKronrod(con$GKk)$sk
        nk <- length(sk)
        if (is.null(con$knots) || !is.numeric(con$knots)) {
            Q <- con$lng.in.kn + 1
            qs <- if (con$ObsTimes.knots) {
                unique(quantile(Time, seq(0, 1, len = Q + 1), 
                    names = FALSE)[-c(1, Q + 1)])
            } else {
                unique(quantile(Time[d == 1], seq(0, 1, len = Q - 1), 
                    names = FALSE))                
            }
            qs <- qs + 1e-06
            if (max(qs) > max(Time))
                qs[which.max(qs)] <- max(Time) - 1e-06
            con$knots <- qs
            qs <- c(0, qs, max(Time) + 1)
            Q <- length(qs) - 1
        } else {
            qs <- c(0, sort(con$knots), max(Time) + 1)
            Q <- length(qs) - 1
        }
        ind <- findInterval(Time, qs, rightmost.closed = TRUE)
        D <- matrix(0, length(ind), Q)
        D[cbind(seq_along(ind), ind)] <- 1
        D <- D * d
        Tiq <- outer(Time, qs, pmin)
        Lo <- Tiq[, 1:Q]
        Up <- Tiq[, 2:(Q+1)]
        T <- Up - Lo
        P <- T / 2
        P[P < con$tol3] <- as.numeric(NA)
        P1 <- (Up + Lo) / 2
        st <- matrix(0, nY, nk*Q)
        skQ <- rep(sk, Q)
        for (i in seq_len(nY)) {
            st[i, ] <- rep(P[i, ], each = nk) * skQ + rep(P1[i, ], each = nk)
        }
        y <- c(y, list(ind.D = ind))
        id.GK <- rep(seq_len(nY), rowSums(!is.na(st)))
        P <- c(t(P))
        data.id2 <- data.id[rep(seq_len(nY), each = nk*Q), ]
        data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
        data.id2 <- data.id2[!is.na(data.id2[[timeVar]]), ]
        if (parameterization %in% c("value", "both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            Xs <- model.matrix(formYx, mfX)
            Zs <- model.matrix(formYz, mfZ)
        }
        if (parameterization %in% c("slope", "both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
        }
        Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        x <- c(x, list(P = P[!is.na(P)], st = st[!is.na(st)], wk = wk, 
            id.GK = id.GK, Q = Q, Ws.intF.vl = Ws.intF.vl, Ws.intF.sl = Ws.intF.sl))
        x <- switch(parameterization,
            "value" = c(x, list(Xs = Xs, Zs = Zs)),
            "slope" = c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv)),
            "both" = c(x, list(Xs.deriv = Xs.deriv, 
                Zs.deriv = Zs.deriv, Xs = Xs, Zs = Zs))
        )
    }
    # extra design matrices for 'method = "Cox-PH-GH"' with event times prior to observed time for the ith subject
    if (method == "Cox-PH-GH") {
        unqT <- sort(unique(Time[d == 1]))
        times <- lapply(Time, function (t) unqT[t >= unqT])
        ind.len <- sapply(times, length)
        indT <- rep(1:nrow(data.id), ind.len)
        data.id2 <- data.id[indT, ]
        data.id2[timeVar] <- pmax(unlist(times, use.names = FALSE) - lag, 0)
        if (parameterization %in% c("value", "both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            Xtime2 <- model.matrix(formYx, mfX)
            Ztime2 <- model.matrix(formYz, mfZ)
        }
        if (parameterization %in% c("slope", "both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            Xtime2.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            Ztime2.deriv <- model.matrix(derivForm$random, mfZ.deriv)
        }
        x <- c(x, list(indT = indT))
        x <- switch(parameterization,
            "value" = c(x, list(Xtime2 = Xtime2, Ztime2 = Ztime2)),
            "slope" = c(x, list(Xtime2.deriv = Xtime2.deriv, 
                Ztime2.deriv = Ztime2.deriv)),
            "both" = c(x, list(Xtime2.deriv = Xtime2.deriv, 
                Ztime2.deriv = Ztime2.deriv, Xtime2 = Xtime2, Ztime2 = Ztime2))
        )
    }
    # initial values
    VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", 
        lmeObject$sigma^2)[[1]]
    if (con$typeGH != "simple") {
        Vs <- vector("list", nY)
        inv.VC <- solve(VC)
        for (i in 1:nY) {
            Z.i <- Z[id == i, , drop = FALSE]
            Vs[[i]] <- solve(crossprod(Z.i) / lmeObject$sigma^2 + inv.VC)        
        }
        con$inv.chol.VCs <- lapply(Vs, function (x) solve(chol(solve(x))))
        con$det.inv.chol.VCs <- sapply(con$inv.chol.VCs, det)
    }
    con$inv.chol.VC <- solve(chol(solve(VC)))
    con$det.inv.chol.VC <- det(con$inv.chol.VC)
    con$ranef <- b
    if (all(VC[upper.tri(VC)] == 0))
        VC <- diag(VC)
    init.surv <- initial.surv(Time, d, W, WintF.vl, WintF.sl, id, 
        times = data[[timeVar]], method, parameterization, long = long, 
        long.deriv = long.deriv, 
        extra = list(W2 = x$W2, control = con, ii = idT, strata = survObject$strata),
        LongFormat = CompRisk | length(Time) > nT)
    if (method == "Cox-PH-GH" && length(init.surv$lambda0) < length(unqT))
        init.surv$lambda0 <- basehaz(survObject)$hazard
    initial.values <- c(list(betas = fixef(lmeObject), sigma = lmeObject$sigma, D = VC), init.surv)
    if (!is.null(init)) {
        nams1 <- names(init)
        nams2 <- names(initial.values)
        if (!is.list(init) || length(noNms <- nams1[!nams1 %in% nams2])) {
            warning("unknown names in 'init': ", paste(noNms, collapse = ", "))
        } else {
            initial.values[nams1] <- init
        }
    }
    # remove objects
    rmObjs <- c(names(x), "y.long", "mfX", "mfZ", "data.id2")
    rm(list = rmObjs); gc()
    # joint model fit
    out <- switch(method,
        "Cox-PH-GH" = phGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
        "weibull-AFT-GH" = weibullAFTGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
        "weibull-PH-GH" = weibullPHGH.fit(x, y, id, initial.values, scaleWB, parameterization, derivForm, con),
        "piecewise-PH-GH" = piecewisePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
        "spline-PH-GH" = splinePHGH.fit(x, y, id, initial.values, parameterization, derivForm, con),
        "ch-Laplace" = chLaplace.fit(x, y, id, initial.values, b, parameterization, derivForm, con))
    # check for problems with the Hessian at convergence
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite.\n")
    }
    out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
    out$x <- x
    out$y <- y
    out$times <- data[[timeVar]]
    out$data <- data
    out$data.id <- data.id
    out$method <- method
    out$termsYx <- TermsX
    out$termsYz <- TermsZ
    if (parameterization %in% c("slope", "both")) {
        out$termsYx.deriv <- TermsX.deriv
        out$termsYz.deriv <- TermsZ.deriv
    }
    out$termsT <- survObject$terms
    out$formYx <- formYx
    out$formYz <- formYz
    out$formT <- formT
    out$timeVar <- timeVar
    out$control <- con
    out$parameterization <- parameterization
    out$derivForm <- derivForm
    out$interFact <- interFact
    out$CompRisk <- CompRisk
    out$LongFormat <- LongFormat
    out$assignY <- attr(lmeObject$fixDF, "assign")[-1]
    out$assignT <- survObject$assign
    out$call <- cl
    class(out) <- "jointModel"
    out
}
