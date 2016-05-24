IRT.parm.tpm <-
function (object, standard.errors = FALSE, digits.abbrv = 6, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    thetas <- object$coef
    parms <- cbind(thetas[, 1], -thetas[, 2] / thetas[, 3], thetas[, 3])
    dimnames(parms) <- list(rownames(thetas), abbreviate(c("Guessing", "Difficulty", "Discrimination"), digits.abbrv))
    out <- list(parms = parms)
    out$se <- if (standard.errors) {
        opt <- options(warn = (-1))
        on.exit(options(opt))
        p <- nrow(thetas)
        constraint <- object$constraint
        type <- object$type
        Var <- if (!is.null(constraint)) {
            ind <- if (type == "rasch" && any(ind.rasch <- constraint[, 2] == 3)) {
                c((constraint[!ind.rasch, 2] - 1) * p + constraint[!ind.rasch, 1], 2 * p + 1)
            } else {
                (constraint[, 2] - 1) * p + constraint[, 1]
            }
            dimV <- if (type == "rasch") 2 * p + 1 else 3 * p
            V <- matrix(NA, dimV, dimV)
            seq.ind <- seq(1, dimV)[-ind]
            V[seq.ind, seq.ind] <- vcov(object)
            V
        } else
            vcov(object)
        ses <- rep(NA, p)
        for (i in seq(along = ses)) {
            di <- if (type == "rasch") 2 * p + 1 else 2 * p + i
            Vi <- Var[c(p + i, di), c(p + i, di)]
            ses[i] <- if (is.na(Vi[2, 2])) {
                Vi[1, 1] / thetas[i, 3]
            } else {
                deltamethod(~ - x1 / x2, c(thetas[i, 2], thetas[i, 3]), Vi)
            }
        }
        c(sqrt(diag(Var[seq(1, p), seq(1, p)])), ses, sqrt(diag(as.matrix(Var[seq(2 * p + 1, di), seq(2 * p + 1, di)]))))
    } else
        NULL
    out
}
