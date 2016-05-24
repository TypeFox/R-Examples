GoF.gpcm <-
function (object, simulate.p.value = TRUE, B = 99, seed = NULL, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    nas <- any(na.ind <- is.na(object$X))
    pearson.chi <- function (object) {
        R <- resid(object)
        res <- (R[, "Resid"])^2
        if (nas)
            sum(res, na.rm = TRUE)
        else
            sum(res, na.rm = TRUE) + sum(R[, "Obs"]) - sum(R[, "Exp"])
    }
    Tobs <- pearson.chi(object)
    betas <- object$coefficients
    ncatg <- sapply(betas, length)
    p <- length(betas)
    df <- prod(ncatg) - attr(logLik(object), "df") - 1
    p.val <- if (!simulate.p.value) {
        pchisq(Tobs, df, lower.tail = FALSE)
    } else {
        constraint <- object$constraint
        IRT.param <- object$IRT.param
        vec.betas <- if (constraint == "gpcm") {
            unlist(betas, use.names = FALSE)
        } else if (constraint == "1PL") {
            betas[seq(1, p - 1)] <- lapply(betas[seq(1, p - 1)], function (x) x[-length(x)])
            unlist(betas, use.names = FALSE)
        } else {
            betas <- lapply(betas, function (x) x[-length(x)])
            unlist(betas, use.names = FALSE)
        }
        Var.betas <- vcov(object)
        n <- nrow(object$X)
        Ts <- numeric(B)
        if (!is.null(seed))
            set.seed(seed)
        old <- options(warn = (-1))
        on.exit(options(old))
        for (i in 1:B) {
            tstat <- try({
                new.betas <- mvrnorm(1, vec.betas, Var.betas)
                new.betas <- betas.gpcm(new.betas, p, ncatg, constraint)
                newData <- rmvordlogis(n, new.betas, IRT.param, "gpcm")
                if (nas)
                    newData[na.ind] <- NA
                fit <- gpcm(newData, constraint, start.val = new.betas, control = object$control)
                pearson.chi(fit)
            }, TRUE)
            Ts[i] <- if (!inherits(tstat, "try-error")) tstat else as.numeric(NA)
        }
        good <- !is.na(Ts)
        if ((newB <- sum(good)) < B) {
            warning("the fit failed in ", B - newB, " datasets.\n")
            B <- newB
        }
        (1 + sum(Ts[good] >= Tobs)) / (B + 1)
    }
    out <- list(Tobs = Tobs, p.value = p.val, df = df, simulate.p.value = simulate.p.value, B = B, call = object$call)
    class(out) <- "GoF.gpcm"
    out
}
