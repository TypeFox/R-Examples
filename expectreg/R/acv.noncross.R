acv.noncross <-
function (penalty, yy, Bx, Bp, P, Amat, bvec, a_p, pp) 
{
    lambda = abs(penalty)
    difw <- 1
    B <- Bp %x% Bx
    W_temp <- matrix(rep(0.5), ncol = ncol(B), nrow = nrow(B), 
        byrow = FALSE)
    z <- B %*% a_p
    bdegp = 1
    mp = length(pp)
    n = length(yy)/mp
    pw <- c()
    for (k in 1:mp) {
        pw <- c(pw, rep(pp[k], n))
    }
    while (difw > 0) {
        W_temp[, 1:(ncol(Bx) * ncol(Bp))] <- pw * (yy > z) + 
            (1 - pw) * (yy <= z)
        W_tempB <- W_temp * B
        Dmat <- t(B) %*% W_tempB + 2 * lambda * (diag(1, mp - 
            1 + bdegp, mp - 1 + bdegp) %x% P)
        dvec <- t(yy) %*% W_tempB
        a_pq <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
            bvec = bvec)
        a_sol <- as.vector(unlist(a_pq$solution))
        z <- B %*% a_sol
        W1 <- W_temp
        W_temp[, 1:(ncol(Bx) * ncol(Bp))] <- pw * (yy > z) + 
            (1 - pw) * (yy <= z)
        difw <- sum(abs(W_temp - W1))
    }
    BW1B <- t(B) %*% (W1 * B)
    H <- solve(BW1B + 2 * lambda * (diag(1, mp - 1 + bdegp, mp - 
        1 + bdegp) %x% P)) %*% BW1B
    score = length(yy) * sum(pw * (yy - z)^2)/sum(1 - diag(diag(H)))^2
    score
}
