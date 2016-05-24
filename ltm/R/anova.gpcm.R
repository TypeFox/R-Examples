anova.gpcm <-
function (object, object2, simulate.p.value = FALSE, B = 200, verbose = getOption("verbose"), seed = NULL, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    if (missing(object2))
        stop("anova.gpcm() computes LRTs between two fitted Generalized Partial Credit Models.")
    if (!inherits(object2, "gpcm"))
        stop(deparse(substitute(object2)), " must inherit from class 'grm'.")
    if (!object$constraint %in% c("1PL", "rasch"))
        stop(deparse(substitute(object)), " must be a constrained GPCM in order to be nested in ", 
                deparse(substitute(object2)))
    if (!isTRUE(all.equal(object$X, object2$X)))
        warning("it seems that the two objects represent models fitted in different data sets.")
    L0 <- logLik(object)
    L1 <- logLik(object2)
    nb0 <- attr(L0, "df")
    nb1 <- attr(L1, "df")
    df <- nb1 - nb0
    if (df < 0)
        stop("'object' is not nested in 'object2'.\n")
    LRT <- - 2 * (L0 - L1)
    attributes(LRT) <- NULL
    if (LRT < 0)
        warning("it seems that the two models are not nested.\n")
    p.value <- if (!simulate.p.value) {
        pchisq(LRT, df, lower.tail = FALSE)
    } else {
        null.thetas <- object$coefficients
        ncatg <- sapply(null.thetas, length)
        p <- length(null.thetas)
        n <- nrow(object$X)
        vec.thetas <- if (object$constraint == "1PL") {
            null.thetas[seq(1, p - 1)] <- lapply(null.thetas[seq(1, p - 1)], function (x) x[-length(x)])
            unlist(null.thetas, use.names = FALSE)
        } else {
            null.thetas <- lapply(null.thetas, function (x) x[-length(x)])
            unlist(null.thetas, use.names = FALSE)
        }
        Var.thetas <- vcov(object)
        Lout <- numeric(B)
        if (!is.null(seed))
            set.seed(seed)
        for (b in 1:B) {
            stat <- try({
                new.thetas <- mvrnorm(1, vec.thetas, Var.thetas)
                new.thetas <- betas.gpcm(new.thetas, p, ncatg, object$constraint)
                newData <- rmvordlogis(n, new.thetas, IRT = object$IRT.param, model = "gpcm")
                fit0 <- gpcm(newData, constraint = object$constraint, start.val = new.thetas, control = object$control)
                fit1 <- gpcm(newData, constraint = object2$constraint, start.val = new.thetas, control = object2$control)
                - 2 * (fit0$log.Lik - fit1$log.Lik)
            }, TRUE)
            Lout[b] <- if (!inherits(stat, "try-error")) stat else as.numeric(NA)
            if (verbose)
                cat("\nDataset:", b, "finished.")
        }
        good <- !is.na(Lout) & Lout > 0
        if ((newB <- sum(good)) < B) {
            warning("the fit failed in ", B - newB, " datasets.\n")
            B <- newB
        }
        (sum(Lout[good] >= LRT) + 1) / (B + 1)
    }
    out <- list(nam0 = deparse(substitute(object)), L0 = L0, nb0 = nb0, aic0 = AIC(object), 
                bic0 = AIC(object, k = log(attr(L0, "n"))), nam1 = deparse(substitute(object2)), L1 = L1, 
                nb1 = nb1, aic1 = AIC(object2), bic1 = AIC(object2, k = log(attr(L1, "n"))), LRT = LRT, df = df, 
                p.value = p.value, simulate.p.value = simulate.p.value, call = object$call)
    if (simulate.p.value) {
        out$B <- B
        out$LRTvals <- Lout[good]
    }
    class(out) <- "aov.gpcm"
    out
}
