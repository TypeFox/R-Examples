kalsmo <- function(object)
    UseMethod("kalsmo")

"kalsmo.car" <- function(object)
{
  par <- object$phi
  n <- object$n.used
  p <- object$order
  k <- object$scale
    TRAN <- array(0, dim = c(p, p, n))
    SSCV <- array(0, dim = c(p, p, n))
    PSCV <- array(0, dim = c(p, p, n))
    FSCV <- array(0, dim = c(p, p, n))
    FSES <- matrix(0, n, p)
    PSES <- matrix(0, n, p)
    SSES <- matrix(0, n, p)
    SSER <- matrix(0, n, p)
    SVAR <- rep(0, n)
    filcomp <- matrix(0, n, p)
    dl <- as.complex(rep(0, p))
#    delta <- diff(series[, 1])
  delta <- diff(object$tim)
    delta <- c(0, delta)
#    z.obs <- series[, 2] - par[p + 1]
  z.obs <- object$ser - object$x.mean# detrend mean
    K <- c(rep(0, p - 1), 1)
    H <- choose(p - 1, k = 0:(p - 1))/k^(0:(p - 1))
#    w <- sort(polyroot(rev(c(1, par[1:p]))))
    w <- sort(round(polyroot(rev(c(1, par[1:p]))),4)) 
#   sortind <- sort.list(abs(object$rootr), decreasing = TRUE)
#    carroot <- complex(real = object$rootr, imag = object$rooti)
#    w <- carroot <- rev(sort(carroot, partial = sortind))
    r <- -k * (1 - w)/(1 + w)
    sortind <- sort.list(abs(Re(r)), decreasing = TRUE)
    r <- rev(r[sortind]) #changed  
    ###r <- rev(sort(r, partial = sortind))
    R <- matrix(0, p, p)
    E <- matrix(0, p, p)
    P <- matrix(0, p, p)
    for (i in 1:p) R[i, ] <- r^(i - 1)
    R <- round(R, 4)
    D <- diag(r)
   J <- solve(R)[, p]
    G <- H %*% R
    state <- matrix(0, n + 1, p)
    innovation <- v <- rep(0, n)
    z.est <- rep(0, n)
    sse <- 0
    P0 <- matrix(0, p, p)
    for (i in 1:p) for (j in 1:p) P0[i, j] <- -J[i] * Conj(J[j])/(r[i] + 
        Conj(r[j]))
    P <- P0
    i <- 1
    while (i < n + 1) {
        dl <- exp(r * delta[i])
        E <- diag(dl)
        TRAN[, , i] <- E
        state[i, ] <- dl * state[i, ]
        PSES[i, ] <- state[i, ]
        P <- E %*% (P - P0) %*% Conj(E) + P0
        PSCV[, , i] <- P
        z.est[i] <- Re(G %*% state[i, ])
        innovation[i] <- z.obs[i] - z.est[i]
        v[i] <- Re(G %*% P %*% t(Conj(G)))
        gain <- P %*% Conj(t(G))/v[i]
        state[i + 1, ] <- state[i, ] + gain * innovation[i]
        FSES[i, ] <- state[i + 1, ]
        filcomp[i, ] <- G * state[i + 1, ]
        sse <- sse + innovation[i]^2/v[i]
        P <- P - gain %*% Conj(t(gain)) * v[i]
        FSCV[, , i] <- P
        i <- i + 1
    }
    i <- n
    while (i > 0) {
        if (i < n) {
            P.star <- FSCV[, , i] %*% Conj(TRAN[, , (i + 1)]) %*% 
                solve(PSCV[, , (i + 1)])
        }
        SSES[i, ] <- FSES[i, ]
        SSCV[, , i] <- FSCV[, , i]
        if (i < n) {
            SSES[i, ] <- FSES[i, ] + P.star %*% (SSES[i + 1, 
                ] - PSES[i + 1, ])
            SSCV[, , i] <- FSCV[, , i] + P.star %*% (SSCV[, , 
                (i + 1)] - PSCV[, , (i + 1)]) %*% Conj(P.star)
        }
        SSER[i, ] <- G * SSES[i, ]
        SVAR[i] <- Re(G %*% SSCV[, , i] %*% t(Conj(G)))
        i <- i - 1
    }
    innovation <- round(innovation, 4)
    v <- round(v, 4)
    state <- round(state, 4)
    vwhite <- sse/(n - object$np)
    RET <- list(tim=object$tim,ser=object$ser,sse = Re(sse), vwhite = Re(vwhite), state = state, innovation = innovation, 
        filcomp = filcomp, v = v, roots = r, sses = SSES, sser = SSER, 
        svar = SVAR)
}
