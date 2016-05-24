##
##  g m r e s . R  Generalized Linear equation solver
##


gmres <- function(A, b, x0 = rep(0, length(b)), 
                  errtol = 1e-6, kmax = length(b)+1, reorth = 1) {
    stopifnot(is.numeric(A), is.numeric(b), is.matrix(A))
    b <- as.matrix(b)
    n <- length(b)
    if (nrow(A) != n || ncol(A) != n)
        stop("Matrix 'A' must be square and compatible with 'b'.")

    # initialization
    h <- zeros(kmax)
    v <- zeros(n, kmax)
    c <- zeros(kmax+1, 1)
    s <- zeros(kmax+1, 1)

    normF <- function(x) norm(as.matrix(x), type = 'F')

    x <- as.matrix(x0)
    if (norm(x, 'F') != 0) {
        r <- b - A %*% x
    } else {
        r <- b
    }

    rho <- norm(r, 'F')
    g <- rho*eye(kmax+1, 1)
    errtol <- errtol * norm(b, 'F')
    error <- c()

    # test for termination on entry
    error <- c(error, rho)
    niter <- 0
    if(rho < errtol)
        return(list(x = x, error = error, niter = niter))

    v[, 1] <- r/rho
    beta <- rho

    # GMRES iteration
    k <- 0
    while (rho > errtol && k < kmax) {
        k <- k + 1
        v[, k+1] <- A %*% v[, k]
        normav <- normF(v[, k+1])

        # modified Gram-Schmidt
        for (j in 1:k) {
            h[j, k] <- t(v[, j]) %*% v[, k+1]
            v[, k+1] <- v[, k+1] - h[j,k] * v[, j]
        }
        h[k+1, k] <- normF(v[, k+1])
        normav2 <- h[k+1, k]

        # reorthogonalize
        if ((reorth == 1 && normav + 0.001*normav2 == normav) || reorth == 3) {
            for (j in 1:k) {
                hr <- t(v[, j]) %*% v[, k+1]
                h[j, k] <- h[j, k] + hr
                v[, k+1] = v[, k+1] - hr*v[, j]
            }
            h[k+1, k] <- normF(v[, k+1])
        }

        # watch out for happy breakdown
        if (h[k+1, k] != 0)
            v[, k+1] <- v[, k+1] / h[k+1, k]

        # form and store the information for the new Givens rotation
        if (k > 1)
            h[1:k, k] <- .givapp(c[1:(k-1)], s[1:(k-1)], h[1:k, k], k-1)

        nu <- normF(h[k:(k+1), k])
        if (nu != 0) {
            c[k] <- Conj(h[k, k] / nu)
            s[k] <- -h[k+1, k] / nu
            h[k, k] <- c[k] * h[k, k] - s[k] * h[k+1, k]
            h[k+1, k] <- 0
            g[k:(k+1)] <- .givapp(c[k], s[k], g[k:(k+1)], 1)
        }

        # update the residual norm
        rho <- abs(g[k+1])
        error <- c(error,rho)
    }

    # at this point either k > kmax or rho < errtol
    y <- qr.solve(h[1:k, 1:k], g[1:k])
    x <- x0 + v[1:n, 1:k] %*% y

    return(list(x = x, error = error, niter = k))
}


.givapp <- function(c, s, v_in, k) {
    v_rot <- v_in
    for (i in 1:k) {
        w1 <- c[i] * v_rot[i] - s[i] * v_rot[i+1]
        w2 <- s[i] * v_rot[i] + Conj(c[i]) * v_rot[i+1]
        v_rot[i:(i+1)] <- c(w1, w2)
    }
    v_rot
}
