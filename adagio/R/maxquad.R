initMaxquad <- function(n, m) {
    matr_A <- matrix(NA, n*m, n)
    matr_B <- matrix(NA, n, m)
    for (l in 1:m) {
        for (i in 1:n) {
            if (i < n) {
            for (j in (i+1):n) {
                matr_A[i+(l-1)*n, j] <- exp(i/j) * cos(i*j) * sin(l) 
                matr_A[j+(l-1)*n, i] <- matr_A[i+(l-1)*n, j]
            }
            }
            matr_B[i, l] <- exp(i/l) * sin(i*l)
        }
        for (i in 1:n) {
            j <- which((1:n) != i)
            matr_A[i+(l-1)*n, i] <- i*abs(sin(l))/n + 
                                      sum(abs(matr_A[i+(l-1)*n, j]))
        }
    }
    return(list(A = matr_A, B = matr_B))
}


maxquad <- function(n, m) {
    matr <- initMaxquad(n, m)
    matr_A <- matr$A
    matr_B <- matr$B

    maxquadf <- function(x) {
        n <- length(x)
        f <- t(matr_A[1:n, ] %*% x) %*% x - c(t(matr_B[, 1]) %*% x)
        for (l in 2:m) {
            d <- t(matr_A[((l-1)*n+1):(l*n), ] %*% x) %*% x - t(matr_B[, l]) %*% x
            if (d > f) f <- d
        }
        drop(f)
    }

    maxquadg <- function(x) {
        n <- length(x)
        f <- t(matr_A[1:n, ] %*% x) %*% x - t(matr_B[, 1]) %*% x
        k <- 1
        for (l in 2:m) {
            d <- t(matr_A[((l-1)*n+1):(l*n), ] %*% x) %*% x - t(matr_B[, l]) %*% x
            if (d > f) {
                f <- d
                k <- l
            }
        }
        g <- 2 * matr_A[((k-1)*n+1):(k*n), ] %*% x - matr_B[, k]
        drop(g)
    }


    list(fn = maxquadf, gr = maxquadg)
}
