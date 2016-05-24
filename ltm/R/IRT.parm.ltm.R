IRT.parm.ltm <-
function (object, standard.errors = FALSE, robust = FALSE, digits.abbrv = 6, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    if (object$ltst$factors > 1 || object$ltst$quad.z1)
        return(list(parms = object$coef, se = if (standard.errors) sqrt(diag(vcov(object, robust = robust))) else NULL))
    thetas <- object$coef
    parms <- cbind(-thetas[, 1] / thetas[, 2], thetas[, 2])
    dimnames(parms) <- list(rownames(thetas), abbreviate(c("Difficulty", "Discrimination"), digits.abbrv))
    out <- list(parms = parms)
    out$se <- if (standard.errors) {
        p <- nrow(thetas)
        constraint <- object$constraint
        Var <- if (!is.null(constraint)) {
            ind <- (constraint[, 2] - 1) * p + constraint[, 1]
            V <- matrix(NA, 2 * p, 2 * p)
            seq.ind <- seq(1, 2 * p)[-ind]
            V[seq.ind, seq.ind] <- vcov(object, robust = robust)
            V
        } else
            vcov(object, robust = robust)
        ses <- rep(NA, p)
        for (i in seq(along = ses)) {
            Vi <- Var[c(i, p + i), c(i, p + i)]
            ses[i] <- if (is.na(Vi[2, 2])) {
                    Vi[1, 1] / thetas[i, 2]
                } else {
                    deltamethod(~ -x1 / x2, c(thetas[i, 1], thetas[i, 2]), Vi)
                }
        }
        ses <- c(ses, sqrt(diag(Var[seq(p + 1, 2 * p), seq(p + 1, 2 * p)])))
        if (!is.null(constraint))
            ses[-ind]
        else
            ses
    } else
        NULL
    out
}
