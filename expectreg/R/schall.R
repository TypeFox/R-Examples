schall <-
function (yy, B, pp, DD, nb, lala, constmat, center, types) 
{
    glatterms = which(types != "parametric")
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    dc = 1
    dw = 1
    w <- rep(1, times = m)
    it = 1
    while ((dc >= 0.01 || dw != 0) && it < 100) {
        aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, constmat)
        vector.a.ma.schall <- aa$a
        diag.hat = aa$diag.hat.ma
        w0 <- w
        l0 <- lala
        for (i in glatterms) {
            partbasis = (sum(nb[0:(i - 1)]) + 1):(sum(nb[0:i]))
            if (center) {
                partB = B[, -1, drop = FALSE][, partbasis, drop = FALSE]
                partDD = DD[, -1, drop = FALSE][-1, , drop = FALSE][, 
                  partbasis, drop = FALSE]
                partaa = aa$a[-1][partbasis]
            }
            else {
                partB = B[, partbasis, drop = FALSE]
                partDD = DD[, partbasis, drop = FALSE]
                partaa = aa$a[partbasis]
            }
            if (any(partDD != 0)) {
                v <- partDD %*% partaa
                z <- aa$fitted
                w <- aa$weight
                H = solve(t(partB) %*% (w * partB) + lala[i] * 
                  t(partDD) %*% partDD)
                H = apply(sqrt(w) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig2 <- sum(w * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma, na.rm = TRUE))
                tau2 <- sum(v^2, na.rm = TRUE)/sum(H, na.rm = TRUE) + 
                  1e-06
                lala[i] = max(sig2/tau2, 1e-10, na.rm = TRUE)
            }
        }
        dc <- max(abs(log10(l0 + 1e-06) - log10(lala + 1e-06)))
        dw <- sum(w != w0, na.rm = TRUE)
        it = it + 1
    }
    if (it == 100) 
        warning("Schall algorithm did not converge. Stopping after 100 iterations.")
    list(aa$a, lala, aa$diag.hat.ma, aa$fitted)
}
