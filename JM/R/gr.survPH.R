gr.survPH <-
function (thetas) {
    gammas <- thetas[seq_len(ncww)]
    alpha <- thetas[ncww + 1]
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Y2
    exp.eta.s <- exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    Int <- lambda0[ind.L1] * exp.eta.s
    sc.gammas <- if (!is.null(WW)) {
        S1 <- matrix(0, n, k)
        S1[unq.indT, ] <- rowsum(Int, indT, reorder = FALSE)
        - colSums(WW * (d - c((p.byt * (exp.eta.tw * S1)) %*% wGH)), na.rm = TRUE)
    } else 
        NULL
    S2 <- matrix(0, n, k)
    S2[unq.indT, ] <- rowsum(Int * Y2, indT, reorder = FALSE)
    sc.alpha <- - sum((p.byt * (d * Y - exp.eta.tw * S2)) %*% wGH, na.rm = TRUE)
    c(sc.gammas, sc.alpha)
}
