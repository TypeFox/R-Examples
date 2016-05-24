IRT.parm.rasch <-
function (object, standard.errors = FALSE, robust = FALSE, digits.abbrv = 6, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    thetas <- object$coef
    parms <- c(-thetas[, 1] / thetas[1, 2], thetas[1, 2])
    difc <- abbreviate(c("Difficulty", "Discrimination"), digits.abbrv)
    names(parms) <- c(paste(difc[1], rownames(thetas), sep = "."), difc[2])
    out <- list(parms = parms)
    out$se <- if (standard.errors) {
        p <- nrow(thetas)
        constraint <- object$constraint
        if (!is.null(constraint) && any(constraint[, 1] == p + 1)) {
            Var <- vcov(object, robust = robust)
            ses <- sqrt(diag(Var)) / constraint[which(constraint[, 1] == p + 1), 2]
            return(list(parms = parms, se = ses))
        }
        ind <- if (!is.null(constraint)) (1:p)[-sort(constraint[, 1])] else 1:p
        n.ind <- length(ind)
        Var <- vcov(object, robust = robust)
        ses <- numeric(n.ind)
        for (i in seq(along = ind))
            ses[i] <- deltamethod(~ -x1 / x2, c(thetas[ind[i], 1], thetas[1, 2]), Var[c(i, n.ind + 1), c(i, n.ind + 1)])
        c(ses, sqrt(Var[n.ind + 1, n.ind + 1]))
    } else
        NULL
    out
}
