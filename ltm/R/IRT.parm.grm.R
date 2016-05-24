IRT.parm.grm <-
function (object, standard.errors = FALSE, digits.abbrv = 6, ...) {
    if(!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    thetas <- object$coef
    difc <- abbreviate(c("Extremity", "Discrimination"), digits.abbrv)
    parms <- lapply(thetas, function (x) {
        nx <- length(x)
        out <- c(x[-nx] / x[nx], x[nx])
        names(out) <- c(paste(difc[1], seq(1, nx - 1), sep = ""), difc[2])
        out
    })
    out <- list(parms = parms)
    out$se <- if (standard.errors && !is.null(object$hessian)) {
        thets <- unlist(thetas)
        ncatg <- sapply(thetas, length)
        Var <- vcov(object)
        ind1 <- if (object$constrained) sum(ncatg) - length(ncatg) + 1 else cumsum(ncatg)
        ind2 <- if (object$constrained) seq(1, sum(ncatg - 1)) else (1:length(thets))[-ind1]
        if (object$constrained) {
            ll <- cumsum(ncatg)
            thets <- thets[-ll[-length(ll)]]
        }
        ses <- numeric(nrow(Var))
        ses[ind1] <- sqrt(diag(as.matrix(Var[ind1, ind1])))
        ind1 <- if (object$constrained) rep(ind1, length(ind2)) else rep(ind1, ncatg - 1)        
        for (i in seq(along = ind2)) {
            ses[ind2[i]] <- deltamethod(~ x1 / x2,
                c(thets[ind2[i]], thets[ind1[i]]), Var[c(ind2[i], ind1[i]), c(ind2[i], ind1[i])])
        }
        names(ses) <- colnames(Var)
        ses
    } else
        NULL
    out
}
