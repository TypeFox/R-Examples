H.longPC <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Ys <- as.vector(Xs %*% betas) + Zsb
        Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
        eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
        eta.s <- if (parameterization == "both")
            eta.s + Ws.intF.sl.alph * Ys.deriv 
        else
            Ws.intF.sl.alph * Ys.deriv
    }
    exp.eta.tw <- exp(eta.tw)
    H1 <- XtX / sigma^2
    Int <- xi[ind.K] * wkP * exp(eta.s)
    H2 <- H1
    H2 <- matrix(0, ncx, ncx)
    for (i in 1:ncx) {
        for (j in i:ncx) {
            XX <- if (parameterization == "value") {
                Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
            } else if (parameterization == "slope") {
                if (i %in% indFixed && j %in% indFixed) {
                    ii <- match(i, indFixed)
                    jj <- match(j, indFixed)
                    Ws.intF.sl.alph^2 * Xs.deriv[, ii] * Xs.deriv[, jj]
                } else
                    0
            } else {
                if (i %in% indFixed && j %in% indFixed) {
                    ii <- match(i, indFixed)
                    jj <- match(j, indFixed)
                    (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * Xs.deriv[, ii]) * 
                        (Ws.intF.vl.alph * Xs[, j] + Ws.intF.sl.alph * Xs.deriv[, jj])
                } else if (i %in% indFixed && !j %in% indFixed) {
                    ii <- match(i, indFixed)
                    (Ws.intF.vl.alph * Xs[, i] + Ws.intF.sl.alph * Xs.deriv[, ii]) * 
                        (Ws.intF.vl.alph * Xs[, j])
                } else if (!i %in% indFixed && j %in% indFixed) {
                    jj <- match(j, indFixed)
                    (Ws.intF.vl.alph * Xs[, i]) * (Ws.intF.vl.alph * Xs[, j] + 
                        Ws.intF.sl.alph * Xs.deriv[, jj])
                } else {
                    Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
                }
            }
            ki <- exp.eta.tw * rowsum(Int * XX, id.GK, reorder = FALSE)
            kii <- c((p.byt * ki) %*% wGH)
            H2[i, j] <- sum(kii, na.rm = TRUE)
        }
    }
    H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
    H1 + H2
}
