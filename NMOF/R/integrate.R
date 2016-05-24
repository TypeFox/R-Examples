changeInterval <- function(nodes, weights,
                           oldmin, oldmax, newmin, newmax) {
    newrange <- newmax - newmin
    oldrange <- oldmax - oldmin
    nodes <- (newrange*nodes + newmin*oldmax - newmax*oldmin)/oldrange
    weights <- weights * newrange/oldrange
    list(nodes = nodes, weights = weights)
}

xwGauss <- function(n, method = "legendre") {
    n <- makeInteger(n, "'n'", 1L)
    method <- match.arg(tolower(method),
                        choices = c("legendre", "laguerre", "hermite"))
    switch(method,
           legendre = {
               ##if (n == 1L)
               ##    return(list(nodes = 0, weights = 2))
               ind <- seq_len(n - 1L)
               eta <- 1 / sqrt(4 - (ind)^(-2))
               A <- array(0, dim = c(n, n))
               A[cbind(ind,    ind+1L)] <- eta
               A[cbind(ind+1L, ind)] <- eta
               W <- 2
           },
           laguerre = {
               ind <- seq_len(n)
               delta <- 2 * ind - 1
               ind <- ind[-n]
               eta <- ind
               A <- array(0, dim = c(n, n))
               diag(A) <- delta
               A[cbind(ind,    ind+1L)] <- eta
               A[cbind(ind+1L, ind)] <- eta
               W <- 1
           },
           hermite = {
               ind <- seq_len(n - 1L)
               eta <- sqrt(ind/2)
               A <- array(0, dim = c(n, n))
               A[cbind(ind,    ind+1L)] <- eta
               A[cbind(ind+1L, ind)] <- eta
               W <- sqrt(pi)
           },
           stop("unknown method")
       )
    eig <- eigen(A, symmetric = TRUE)
    x <- eig$values
    i <- order(x)
    x <- x[i]
    w <- W *eig$vectors[1L, i] * eig$vectors[1L, i]
    list(nodes = x, weights = w)
}
