simulateJM <-
function (nsim, nsub, thetas, times, formulas, Data = NULL, method = c("weibull-PH", 
        "weibull-AFT", "piecewise-PH", "spline-PH"), lag = 0, censoring = "uniform", 
        max.FUtime = NULL, seed = NULL, return.ranef = FALSE) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (!is.list(thetas))
        stop("'thetas' must be a list of parameter values.\n")
    if (!is.numeric(times))
        stop("'times' must be a numeric vector denoting the time points ", 
            "at which the longitudinal measurements are collected.\n")
    if (!is.list(formulas))
        stop("'formulas' must be a list of formulas.\n")
    if (any(!names(formulas) %in% c("Yfixed", "Yrandom", "Tfixed", "timeVar")))
        stop("'formulas' should contain these components: 'Yfixed', 'Yrandom', ", 
            "'timeVar', 'Tfixed' - check the help file.\n")
    method <- match.arg(method)
    if (is.null(max.FUtime))
        max.FUtime <- max(times) + 2 * IQR(times)
    # function to compute the inverse survival function
    invS <- function (t, u, i) {
        TD <- function (v) {
            # function to compute the time-dependent 
            # part for patient i at time v
            dd <- Data[rep(i, length(v)), , drop = FALSE]
            dd[[timeVar]] <- pmax(v - lag, 0)
            if (parameterization %in% c("value", "both")) {
                XX <- model.matrix(formYx, data = dd)
                ZZ <- model.matrix(formYz, data = dd)
                Y <- as.vector(XX %*% betas + 
                    rowSums(ZZ * b[rep(i, nrow(ZZ)), , drop = FALSE]))
                out <- alpha * Y
            }
            if (parameterization %in% c("slope", "both")) {
                XX.deriv <- model.matrix(derivForm$fixed, data = dd)
                ZZ.deriv <- model.matrix(derivForm$random, data = dd)
                Y.deriv <- as.vector(XX.deriv %*% betas[derivForm$indFixed] + 
                    rowSums(ZZ.deriv * b[rep(i, nrow(ZZ.deriv)), 
                        derivForm$indRandom, drop = FALSE]))
                out <- if (parameterization == "both") out + 
                    Dalpha * Y.deriv else Dalpha * Y.deriv
            }
            out
        }
        h <- function (s) {
            TD.i <- TD(s)
            switch(method, 
                "weibull-PH" = exp(log(sigma.t) + (sigma.t - 1) * log(s) + 
                    eta.t[i] + TD.i),
                "weibull-AFT" = {
                    ff <- function (v) exp(TD(v))
                    V <- exp(eta.t[i]) * integrate(ff, lower = 0, upper = s)$value
                    exp(log(sigma.t) + (sigma.t - 1) * log(V) + eta.t[i] + TD.i)
                },
                "piecewise-PH" = {
                    ind <- findInterval(s, object$control$knots, 
                        rightmost.closed = TRUE)
                    xi[ind] * exp(eta.t[i] + TD.i)
                },
                "spline-PH" = {
                    W2 <- splineDesign(object$control$knots, s, 
                        ord = object$control$ord, outer.ok = TRUE)
                    exp(c(W2 %*% gammas.bs) + eta.t[i] + TD.i)
                }
            )
        }
        integrate(h, lower = 0, upper = t)$value + log(u)
    }
    # coefficients
    betas <- thetas$betas
    sigma <- thetas$sigma
    D <- thetas$D
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    sigma.t <- thetas$sigma.t
    xi <- thetas$xi
    gammas.bs <- thetas$gammas.bs
    # design matrices
    formYx <- formulas$Yfixed
    formYz <- formulas$Yrandom
    formT <- formulas$Tfixed    
    timeVar <- formulas$timeVar
    id <- rep(seq_len(nsub), each = length(times))
    times <- rep(times, nsub)
    parameterization <- "value"
    if (is.null(Data))
        Data <- data.frame()
    DD <- Data[id, , drop = FALSE]
    DD[[timeVar]] <- times
    X <- model.matrix(formYx, data = DD)
    Z <- model.matrix(formYz, data = DD)
    W <- model.matrix(formT, Data)
    ncz <- ncol(Z)
    if (method %in% c("spline-PH", "piecewise-PH"))
        W <- W[, -1, drop = FALSE]
    # simulation
    val <- vector("list", nsim)
    ranef <- vector("list", nsim)
    for (ii in seq_len(nsim)) {
        # simulate random effects
        ranef[[ii]] <- b <- mvrnorm(nsub, rep(0, ncz), D)
        # simulate event times
        eta.t <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, nsub)
        u <- runif(nsub)
        trueTimes <- numeric(nsub)    
        for (i in 1:nsub) {
            Root <- try(uniroot(invS, interval = c(1e-05, max.FUtime), u = u[i], 
                i = i)$root, TRUE)
            #tries <- 5
            while(inherits(Root, "try-error")) {
                b[i, ] <- c(mvrnorm(1, rep(0, ncz), D))
                Root <- try(uniroot(invS, interval = c(1e-05, max.FUtime), 
                    u = u[i], i = i)$root, TRUE)
                #tries <- tries - 1
            }
            trueTimes[i] <- Root
        }
        # simulate longitudinal responses
        eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
        y <- rnorm(nrow(DD), eta.y, sigma)
        if (!is.null(censoring)) {
            if (censoring == "uniform") {
                Ctimes <- runif(nsub, 0, max.FUtime)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            } else if (is.numeric(censoring)) {
                Ctimes <- rep(censoring, length.out = nsub)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            } else {
                warning("not appropriate value for argument 'censoring';", 
                    " uniform censoring is used instead.\n")
                Ctimes <- runif(nsub, 0, max.time)
                Time <- pmin(trueTimes, Ctimes)
                event <- as.numeric(trueTimes <= Ctimes)
            }
        } else {
            Time <- trueTimes
            event <- rep(1, length(Time))
        }
        ni <- tapply(id, id, length)
        DD$y <- y
        DD$Time <- rep(Time, ni)
        DD$Event <- rep(event, ni)
        DD$id <- id
        DD <- DD[DD[[timeVar]] <= DD$Time, ]
        row.names(DD) <- seq_len(nrow(DD))
        val[[ii]] <- DD
    }
    if (return.ranef)
        attr(val, "ranef") <- ranef
    attr(val, "seed") <- RNGstate
    val
}
