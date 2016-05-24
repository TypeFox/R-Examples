##
##  q p s p e c i a l . R
##


qpspecial <- function(G, x, maxit = 100) {
    stopifnot(is.numeric(G), is.matrix(G))

    m <- nrow(G); n <- ncol(G)
    if (m*n <= 0) {
        warning("qpspecial: Matrix 'G' is empty; nothing can be done.")
        return(list(x = c(), d = c(), q = Inf, niter = 0, info = 2))
    }
    maxit <- max(floor(maxit), 5)

    e <- matrix(1, n, 1)
    if (missing(x)) {
        x <- e
    } else {
        x <- as.matrix(c(x))
        nx <- length(x)
        if (any(x < 0) || nx != n)
            x <- e
    }

    idx <- seq(1, (n*n), by = n+1)
    Q <- t(G) %*% G
    z <- x
    y <- 0
    eta <- 0.9995
    delta <- 3
    mu0 <- sum(x*z)/n
  
    tolmu <- 1e-5
    tolrs <- 1e-5
    kmu <- tolmu * mu0
  
    nQ <- norm(Q, "I") + 2
    krs <- tolrs * nQ

    ap <- 0; ad <- 0

    for (k in 1:maxit) {

        r1 <- -Q %*% x + e*y + z
        r2 <- -1 + sum(x)
        r3 <- -x * z
        rs <- norm(rbind(r1, r2), "I")
        mu <- -sum(r3)/n
    
        if (mu < kmu) {
            if (rs < krs) {
                niter <- k-1; info <- 0
                break
            }
        }

        zdx <- z / x
        QD <- Q
        QD[idx] <- QD[idx] + zdx
        C <- chol(QD)
        KT <- solve(t(C), e)
        M <- sum(KT * KT)

        r4 <- r1 + r3/x
        r5 <- sum(KT * solve(t(C), r4))
        r6 <- r2 + r5
        dy <- -r6/M
        r7 <- r4 + e*dy
        dx <- solve(C, solve(t(C), r7))
        dz <- (r3 - z*dx)/x

        p <- -x / dx
        p0 <- p[p > 0]
        if (length(p0) > 0) { ap <- min(p0, 1)
        } else              { ap <- 1 }
        p <- -z / dz
        p0 <- p[p > 0]
        if (length(p0) > 0) { ad <- min(p0, 1)
        } else              { ad <- 1 }

        mauff <- sum((x + ap*dx) * (z + ad*dz)) / n
        sig <- (mauff/mu)^delta
    
        r3 <- r3 + sig*mu
        r3 <- r3 - dx*dz
        r4 <- r1 + r3/x
        r5 <- sum(KT * solve(t(C), r4))
        r6 <- r2 + r5
        dy <- -r6/M
        r7 <- r4 + e*dy
        dx <- solve(C, solve(t(C), r7))
        dz <- (r3 - z*dx)/x

        p <- -x / dx
        p0 <- p[p > 0]
        if (length(p0) > 0) { ap <- min(p0, 1)
        } else              { ap <- 1 }
        p <- -z/dz
        p0 <- p[p > 0]
        if (length(p0) > 0) { ad <- min(p0, 1)
        } else              { ad <- 1 }

        x <- x + eta * ap * dx
        y <- y + eta * ad * dy
        z <- z + eta * ad * dz
    }

    if (k == maxit) info <- 1
    x <- pmax(x,0)
    x <- x/sum(x)
    d <- G %*% x
    q <- sum(d * d)
    list(x = x, d = d, q = q, niter = k, info = info)
}
