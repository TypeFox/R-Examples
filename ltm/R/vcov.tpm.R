vcov.tpm <-
function (object, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    inv.hes <- solve(object$hessian)
    p <- nrow(object$coef)
    thetas <- object$coef 
    onams <- nams <- as.vector(t(outer(c("c.", "beta.1", "beta.2"), as.character(1:p), paste, sep = "")))
    if (!is.null(constraint <- object$constraint)) {
        nams <- nams[-((constraint[, 2] - 1) * p + constraint[, 1])]
        thetas <- thetas[-((constraint[, 2] - 1) * p + constraint[, 1])]
        cp <- p - sum(constraint[, 2] == 1)
        dfp <- p - sum(constraint[, 2] == 2)
        dip <- if ((type <- object$type) == "rasch") {
            if (any(constraint[, 2] == 3)) 0 else 1
        } else 
            p - sum(constraint[, 2] == 3) 
    } else {
        cp <- dfp <- p
        dip <- if ((type <- object$type) == "rasch") 1 else p
    }
    nams <- nams[seq(1, cp + dfp + dip)]
    thetas <- unique(c(thetas))
    if (type == "rasch" && any(constraint[, 2] == 3))
        thetas <- thetas[-length(thetas)]
    formlLis.c <- if (cp == 0) NULL else lapply(paste("x", 1:cp, sep = ""), function (strg) {
        as.formula(paste("~ 1 / (1 + exp(-", strg, "))", sep = ""))
    })
    formlLis.b <- lapply(paste("x", seq(cp + 1, length(thetas)), sep = ""), function (strg) {
        as.formula(paste("~ ", strg, sep = ""))
    })    
    formlLis <- c(formlLis.c, formlLis.b)
    inv.hes[] <- deltamethod(formlLis, thetas, inv.hes, FALSE)
    if (cp > 0)
        inv.hes[1:cp, 1:cp] <- inv.hes[1:cp, 1:cp] * object$max.guessing^2
    inv.hes[seq(cp + 1, nrow(inv.hes)), 1:cp] <- inv.hes[seq(cp + 1, nrow(inv.hes)), 1:cp] * object$max.guessing
    inv.hes[1:cp, seq(cp + 1, nrow(inv.hes))] <- inv.hes[1:cp, seq(cp + 1, nrow(inv.hes))] * object$max.guessing
    dimnames(inv.hes) <- list(nams, nams)
    attr(inv.hes, "drop.ind") <- onams %in% nams
    inv.hes
}
