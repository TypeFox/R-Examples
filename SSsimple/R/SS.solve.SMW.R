SS.solve.SMW <-
function (Z, F, H, Q, inv.R, length.out, P0, beta0 = 0) {
	R <- inv.R
    rm(inv.R)
    d <- ncol(H)
    n <- nrow(H)
    T <- length.out
    params <- .internal.chk.mod.params(F, H, Q, R, P0 = P0, beta0 = beta0, d, n) ### R IS INVERSE !!!
    F <- params$F
    H <- params$H
    Q <- params$Q
    R <- params$R ############################## R IS INVERSE !!!
    P0 <- params$P0
    beta0 <- params$beta0
    B.apri <- matrix(NA, T, d)
    B.apos <- matrix(NA, T, d)
    P <- P0
    I <- diag(1, d)
    tHR <- crossprod(H, R)
    RH <- R %*% H
    tHRH <- tHR %*% H
    for (j in 1:T) {
        if (j == 1) {
            B.apri[j, ] = F %*% beta0
        }
        else {
            B.apri[j, ] = F %*% B.apos[j - 1, ]
        }
        P <- F %*% tcrossprod(P, F) + Q
        PtH <- tcrossprod(P, H)
		
		SMW.HPHR <- R - RH %*% solve( solve( P, tol=0 ) + tHRH,  tol=0 ) %*% tHR ###### R IS INVERSE !!!
		
        K <- PtH %*% SMW.HPHR
        B.apos[j, ] <- B.apri[j, ] + K %*% (Z[j, ] - H %*% B.apri[j, ])
        P <- (I - K %*% H) %*% P
    }
    return(list("B.apri" = B.apri, "B.apos" = B.apos))
}
