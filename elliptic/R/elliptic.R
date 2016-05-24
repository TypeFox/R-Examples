"amn" <-
function (u) 
{
    "asub" <- function(m, n) {
        if ((m < 0) | (n < 0)) {
            return(0)
        }
        else {
            return(a[n + 1, m + 1])
        }
    }
    index <- cbind(unlist(sapply(1:u, function(i) {
        1:i
    })), unlist(sapply(1:u, function(i) {
        i:1
    })))
    a <- matrix(0, u, u)
    a[1, 1] <- 1
    for (i in 2:nrow(index)) {
        m <- index[i, 2] - 1
        n <- index[i, 1] - 1
        a[index[i, , drop = FALSE]] <- 3 * (m + 1) * asub(m + 
            1, n - 1) + 16/3 * (n + 1) * asub(m - 2, n + 1) - 
            1/3 * (2 * m + 3 * n - 1) * (4 * m + 6 * n - 1) * 
                asub(m - 1, n)
    }
    return(a)
}
"as.primitive" <-
function (p, n = 3, tol = 1e-05, give.answers = FALSE) 
{
    if (class(p) == "primitive") {
        return(p)
    }
    tau <- p[2]/p[1]
    if (Im(tau) == 0) {
        stop("period ratio real")
    }
    jj <- as.matrix(expand.grid(0:n, n:-n))[-n * (n + 1) - 1, 
        ]
    magnitudes <- abs(jj %*% p)
    o <- order(magnitudes, abs(jj[, 1] - 1), abs(jj[, 2]))
    first.row <- jj[o[1], ]
    found.second.row <- FALSE
    for (i in 2:nrow(jj)) {
        second.row.candidate <- jj[o[i], ]
        jj.matrix <- cbind(first.row, second.row.candidate)
        if (abs(det(jj.matrix)) > tol) {
            second.row <- second.row.candidate
            found.second.row <- TRUE
            break()
        }
    }
    if (!found.second.row) {
        stop("unimodular transformation out of range.  Try increaing n")
    }
    M <- rbind(first.row, second.row)
    p.prim <- M %*% p
    if (Re(p.prim[1]) < 0) {
        p.prim[1] <- -p.prim[1]
        M[1, ] <- -M[1, ]
    }
    if (Im(p.prim[2]/p.prim[1]) < -tol) {
        M[2, ] <- -M[2, ]
    }
    out <- as.vector(M %*% p)
    class(out) <- "primitive"
    if (give.answers) {
        return(list(M = M, p = out, mags = abs(out)))
    }
    else {
        return(out)
    }
}
"cc" <-
function (u, m, ...) 
{
    1
}
"cd" <-
function (u, m, ...) 
{
    theta.c(u, m = m, ...)/theta.d(u, m = m, ...)
}
"ck" <-
function (g, n = 20) 
{
    if (n < 3) {
        stop("error: n must be >3")
    }
    out <- rep(0, n)
    g2 <- g[1]
    g3 <- g[2]
    out[1] <- 0
    out[2] <- g2/20
    out[3] <- g3/28
    for (k in 4:n) {
        out[k] <- 3/((2 * k + 1) * (k - 3)) * sum(out[2:(k - 
            2)] * out[(k - 2):2])
    }
    return(out)
}
"cn" <-
function (u, m, ...) 
{
    theta.c(u, m = m, ...)/theta.n(u, m = m, ...)
}
"congruence" <-
function (a, l = 1) 
{
    l <- as.integer(l)
    m <- as.integer(a[1])
    n <- as.integer(a[2])
    zero <- as.integer(0)
    one <- as.integer(1)
    if (m == zero & n == one) {
        return(NULL)
    }
    if (m == one & n == zero) {
        return(c(NA, 1))
    }
    if (m == 1) {
        return(rbind(c(a),c(1, n + l), c(0, l)))
    }
    if (n == 1) {
        return(rbind(c(a),c(m - l, 1), c(l, 0)))
    }
    q1 <- which((+l + (1:m) * n)%%m == zero)
    if (!any(q1)) {
        q1 <- NA
    }
    q2 <- which((-l + (1:n) * m)%%n == zero)
    if (!any(q2)) {
        q2 <- NA
    }
    out <- rbind(a, cbind(q1, q2))
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(out)
}
"coqueraux" <-
function (z, g, N = 5, use.fpp = FALSE, give = FALSE) 
{
    if (use.fpp) {
        if (class(g) != "parameters") {
            g <- parameters(g = g)
            g2 <- g$g[1]
            g3 <- g$g[2]
        }
        z <- fpp(z, 2 * g$Omega)
    }
    if (class(g) != "parameters") {
        g2 <- g[1]
        g3 <- g[2]
    }
    else {
        g2 <- g$g[1]
        g3 <- g$g[2]
    }
    z0 <- z/2^N
    z0.sq <- z0^2
    out <- 1/z0.sq + z0.sq * (g2/20 + z0.sq * g3/28)
    for (i in 1:N) {
        out <- -2 * out + (6 * out^2 - g2/2)^2/(4 * (4 * out^3 - 
            g2 * out - g3))
    }
    if (give) {
        error <- abs(g2^2)/2450/2^(8 * N) * abs(z)^9 * sqrt(abs(4 * 
            out^3 - g2 * out - g3))
        return(list(out = out, error = error))
    }
    else {
        return(out)
    }
}
"cs" <-
function (u, m, ...) 
{
    theta.c(u, m = m, ...)/theta.s(u, m = m, ...)
}
"dc" <-
function (u, m, ...) 
{
    theta.d(u, m = m, ...)/theta.c(u, m = m, ...)
}
"dd" <-
function (u, m, ...) 
{
    theta.d(u, m = m, ...)/theta.d(u, m = m, ...)
}
"divisor" <-
function (n, k = 1) 
{
    if (length(n) > 1) {
        return(sapply(n, match.fun(sys.call()[[1]]), k = k))
    }
    if (n == 1) {
        return(1)
    }
    x <- factorize(n)
    xx <- unique(x)
    jj <- rbind(vals = xx, cnts = tabulate(match(x, xx)))
    f <- function(u) {
        p <- u[1]
        alpha <- u[2]
        (p^((alpha + 1) * k) - 1)/(p^k - 1)
    }
    if (k > 0) {
        return(as.integer(prod(apply(jj, 2, f))))
    }
    else {
        return(as.integer(prod(1 + jj[2, ])))
    }
}

"liouville" <- function(n){
out <- ifelse(sapply(factorize(n), length)%%2, -1L, 1L)
out[n==1] <- 1L
return(out)
}
  
"dn" <-
function (u, m, ...) 
{
    theta.d(u, m = m, ...)/theta.n(u, m = m, ...)
}
"ds" <-
function (u, m, ...) 
{
    theta.d(u, m = m, ...)/theta.s(u, m = m, ...)
}
"e16.28.1" <-
function (z, m, ...) 
{
    theta1(z, m = m, ...)^2 * theta4(0, m = m, ...)^2 - theta3(z, 
        m = m, ...)^2 * theta2(0, m = m, ...)^2 + theta2(z, m = m, 
        ...)^2 * theta3(0, m = m, ...)^2
}
"e16.28.2" <-
function (z, m, ...) 
{
    theta2(z, m = m, ...)^2 * theta4(0, m = m, ...)^2 - theta4(z, 
        m = m, ...)^2 * theta2(0, m = m, ...)^2 + theta1(z, m = m, 
        ...)^2 * theta3(0, m = m, ...)^2
}
"e16.28.3" <-
function (z, m, ...) 
{
    theta3(z, m = m, ...)^2 * theta4(0, m = m, ...)^2 - theta4(z, 
        m = m, ...)^2 * theta3(0, m = m, ...)^2 + theta1(z, m = m, 
        ...)^2 * theta2(0, m = m, ...)^2
}
"e16.28.4" <-
function (z, m, ...) 
{
    theta4(z, m = m, ...)^2 * theta4(0, m = m, ...)^2 - theta3(z, 
        m = m, ...)^2 * theta3(0, m = m, ...)^2 + theta2(z, m = m, 
        ...)^2 * theta2(0, m = m, ...)^2
}
"e16.28.5" <-
function (m, ...) 
{
    theta2(0, m = m, ...)^4 + theta4(0, m = m, ...)^4 - theta3(0, 
        m = m, ...)^4
}
"e16.36.6a" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    2 * K * theta1(v, m = m, ...)/(pi * theta2(0, m = m, ...) * 
        theta3(0, m = m, ...) * theta4(0, m = m, ...))
}
"e16.36.6b" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    theta2(v, m = m, ...)/theta2(0, m = m, ...)
}
"e16.36.7a" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    theta3(v, m = m, ...)/theta3(0, m = m, ...)
}
"e16.36.7b" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    theta4(v, m = m, ...)/theta4(0, m = m, ...)
}
"e16.37.1" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 1
    for (n in 1:maxiter) {
        out.new <- out * (1 - 2 * q^(2 * n) * cos(2 * v) + q^(4 * 
            n))
        if (near.match(out, out.new)) {
            return((16 * q/(m * (1 - m)))^(1/6) * sin(v) * out)
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.37.2" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 1
    for (n in 1:maxiter) {
        out.new <- out * (1 + 2 * q^(2 * n) * cos(2 * v) + q^(4 * 
            n))
        if (near.match(out, out.new)) {
            return((16 * q * sqrti(1 - m)/m)^(1/6) * cos(v) * 
                out)
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.37.3" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 1
    for (n in 1:maxiter) {
        out.new <- out * (1 + 2 * q^(2 * n - 1) * cos(2 * v) + 
            q^(4 * n - 2))
        if (near.match(out, out.new)) {
            return((m * (1 - m)/(16 * q))^(1/12) * out)
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.37.4" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 1
    for (n in 1:maxiter) {
        out.new <- out * (1 - 2 * q^(2 * n - 1) * cos(2 * v) + 
            q^(4 * n - 2))
        if (near.match(out, out.new)) {
            return((m/(16 * q * (1 - m)^2))^(1/12) * out)
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.38.1" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + (-1)^n * q^(n * (n + 1)) * sin((2 * 
            n + 1) * v)
        if (near.match(out, out.new)) {
            return(out.new * sqrt(2 * pi * sqrt(q)/(sqrti(m) * 
                sqrti(1 - m) * K)))
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.38.2" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + q^(n * (n + 1)) * cos((2 * n + 1) * 
            v)
        if (near.match(out, out.new)) {
            return(out.new * sqrt(2 * pi * sqrti(q)/(sqrti(m) * 
                K)))
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.38.3" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + q^(n * n) * cos(2 * n * v)
        if (near.match(out, out.new)) {
            return((1 + 2 * out.new) * sqrt(pi/(2 * K)))
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e16.38.4" <-
function (u, m, maxiter = 30) 
{
    q <- nome(m)
    K <- K.fun(m)
    v <- pi * u/(2 * K)
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + (-1)^n * q^(n * n) * cos(2 * n * v)
        if (near.match(out, out.new)) {
            return((1 + 2 * out.new) * sqrt(pi/(2 * sqrti(1 - 
                m) * K)))
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"e18.10.9" <-
function (parameters) 
{
    Omega <- parameters$Omega
    e <- parameters$e
    out <- rep(NA, 3)
    q <- exp((1i) * pi * Omega[2]/Omega[1])
    out[1] <- 12 * Omega[1]^2 * e[1] - pi^2 * (theta3(z = 0, 
        q = q)^4 + theta4(z = 0, q = q)^4)
    out[2] <- 12 * Omega[1]^2 * e[2] - pi^2 * (theta2(z = 0, 
        q = q)^4 - theta4(z = 0, q = q)^4)
    out[3] <- 12 * Omega[1]^2 * e[3] + pi^2 * (theta2(z = 0, 
        q = q)^4 + theta3(z = 0, q = q)^4)
    return(out)
}
"e1e2e3" <-
function (g, use.laurent = TRUE, AnS = is.double(g), Omega = NULL, 
    tol = 1e-06) 
{
    g2 <- g[1]
    g3 <- g[2]
    e <- polyroot(c(-g3, -g2, 0, 4))
    if (AnS) {
        if (!is.double(g)) {
            stop("AnS only consider real g")
        }
        Delta <- g2^3 - 27 * g3^2
        if (Delta < 0) {
            pos.im <- which.max(Im(e))
            neg.im <- which.min(Im(e))
            e <- c(e[pos.im], e[-c(pos.im, neg.im)], e[neg.im])
            e <- massage(e)
        }
        else {
            e <- sort(e, decreasing = TRUE)
            e <- Re(e)
        }
        names(e) <- c("e1", "e2", "e3")
        return(e)
    }
    out.accurate <- e
    if (use.laurent) {
        e.approx <- rep(NA, 3)
        out <- e.approx
        if (is.null(Omega)) {
            Omega <- half.periods(e = e)
        }
        if (AnS) {
            e.approx[1] <- P.laurent(z = Omega[1], g = g, tol = tol)
            e.approx[2] <- P.laurent(z = Omega[2], g = g, tol = tol)
            e.approx[3] <- P.laurent(z = fpp(sum(Omega), Omega * 
                2), g = g, tol = tol)
        }
        else {
            e.approx[1] <- coqueraux(z = Omega[1], g = g)
            e.approx[3] <- coqueraux(z = Omega[2], g = g)
            e.approx[2] <- -e.approx[1] - e.approx[3]
        }
        out[1] <- out.accurate[which.min(abs(e - e.approx[1]))]
        out[2] <- out.accurate[which.min(abs(e - e.approx[2]))]
        out[3] <- out.accurate[which.min(abs(e - e.approx[3]))]
    }
    names(out) <- c("e1", "e2", "e3")
    return(out)
}
"eee.cardano" <-
function (g) 
{
    g2 <- g[1]
    g3 <- g[2]
    Delta <- g3^2 - g2^3/27
    epsilon <- exp(pi * (0+2i)/3)
    if (g2 == 0) {
        e1 <- (g3/4)^(1/3)
        return(c(e1 = e1, e2 = epsilon * e1, e3 = e1/epsilon))
    }
    if (g3 == 0) {
        e1 <- g2^(1/2)/2
        return(c(e1 = e1, e2 = -e1, e3 = 0))
    }
    alpha <- (g2/12)^(1/2)
    beta <- (g3 + sqrt(Delta))^(1/3)/(2 * alpha)
    gamma <- 1/beta
    e1 <- alpha * (beta + gamma)
    e2 <- alpha * (epsilon * beta + gamma/epsilon)
    e3 <- alpha * (1/epsilon * beta + gamma * epsilon)
    return(c(e1, e2, e3))
}
"equianharmonic" <-
function (...) 
{
    jj <- gamma(1/3)^3/(4 * pi)
    omega1 <- jj/2 - (1i) * jj * sqrt(3)/2
    omega2 <- Conj(omega1)
    Omega <- c(omega1, omega2)
    epsilon <- exp(pi * (1i)/3)
    e <- c(4^(-1/3) * epsilon^2, 4^(-1/3), 4^(-1/3) * epsilon^(-2))
    names(e) <- c("e1", "e2", "e3")
    eta <- epsilon * pi/(2 * omega2 * sqrt(3))
    etadash <- -epsilon^(-1) * pi/(2 * omega2 * sqrt(3))
    Eta <- c(etadash, eta - etadash, -eta)
    out <- list(Omega = Omega, q = exp(pi * (1i) * omega2/omega1), 
        e = e, g = c(g2 = 0, g3 = 1), Delta = -27, Eta = Eta, 
        is.AnS = TRUE, given = "d")
    class(out) <- "parameters"
    return(out)
}
"eta" <-
function (z, ...) 
{
    f <- function(u) {
        q <- exp(pi * (1i) * u * 3)
        return(theta3((1/2 + u/2) * pi, q = q, ...))
    }
    out <- sapply(z, f) * exp(pi * (1i) * z/12)
    attributes(out) <- attributes(z)
    return(out)
}
"eta.series" <-
function (z, maxiter = 300) 
{
    jj <- 1
    for (n in 1:maxiter) {
        jj <- jj * (1 - exp(2 * pi * (1i) * n * z))
    }
    return(exp(pi * (1i) * z/12) * jj)
}
"factorize" <-
function (n) 
{
    if (!is.numeric(n)) 
        stop("cannot factorize non-numeric arguments")
    if (length(n) > 1) {
        l <- list()
        for (i in seq(along = n)) l[[i]] <- Recall(n[i])
        return(l)
    }
    if (n != round(n) || n < 2) 
        return(n)
    tab <- primes(n)
    fac <- numeric(0)
    while (length(tab <- tab[n%%tab == 0]) > 0) {
        n <- n/prod(tab)
        fac <- c(fac, tab)
    }
    as.integer(sort(fac))
}
"farey" <-
function (n, print = FALSE, give.series = FALSE) 
{
    a <- outer(0:n, 0:n, "/")
    a <- as.vector(a[is.finite(a) & a <= 1])
    a <- MASS::fractions(unique(MASS::rational(sort(a))))
    a <- attributes(a)$fracs
    a[1] <- "0/1"
    a[length(a)] <- "1/1"
    if (print) {
        print(a)
        return(invisible(a))
    }
    a <- strsplit(a, "/")
    a <- do.call("rbind", a)
    mode(a) <- "integer"
    if (give.series) {
        colnames(a) <- c("num", "den")
        return(a)
    }
    if (n == 0) {
        return(c(0, 1, 1, 0))
    }
    if (n == 1) {
        return(c(0, 1, 1, 1))
    }
    r <- nrow(a)
    a <- a[c(1, rep(2:(r - 1), each = 2), r), ]
    a <- as.vector(t(a))
    dim(a) <- c(4, r - 1)
    return(a)
}
"fpp" <-
function (z, p, give = FALSE) 
{
    attr <- attributes(z)
    z <- as.vector(z)
    jj.mn <- round(mn(z, p))
    z <- z - jj.mn %*% p
    attributes(z) <- attr
    if (give) {
        return(list(answer = z, mn = jj.mn))
    }
    else {
        return(z)
    }
}
"g2.fun" <-
function (b, use.first = TRUE, ...) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    jj2 <- theta2(0, q = q, ...)
    jj3 <- theta3(0, q = q, ...)
    jj4 <- theta4(0, q = q, ...)
    if (use.first) {
        return((1/12) * (pi/p1)^4 * (jj2^8 - jj3^4 * jj2^4 + 
            jj3^8))
    }
    else {
        return((2/3) * (pi/(2 * p1))^4 * (jj2^8 + jj3^8 + jj4^8))
    }
}
"g2.fun.direct" <-
function (b, nmax = 50, tol = 1e-10) 
{
    if (length(b) == 1) {
        b <- c(1, b)
    }
    jj <- expand.grid(-nmax:nmax, -nmax:nmax)
    jj <- jj[-(2 * nmax * nmax + 2 * nmax + 1), ]
    jj <- as.matrix(jj)
    return(massage(60 * sum(1/(2 * jj %*% b)^4), tol = tol))
}
"g2.fun.divisor" <-
function (b, nmax = 50, strict = TRUE, tol = 1e-10) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    s <- 0
    for (n in 1:nmax) {
        s.new <- s + divisor(n, 3) * q^(2 * n)
        if (near.match(s, s.new)) {
            return(massage((pi/p1)^4 * (1/12 + 20 * s), tol = tol))
        }
        s <- s.new
    }
    if (strict) {
        stop("series not converged")
    }
    else {
        warning("series not converged.  Partial sum returned")
        return(massage((2 * pi/p1)^4 * (1/12 + 20 * s), tol = tol))
    }
}
"g2.fun.fixed" <-
function (b, nmax = 50, tol=1e-10, give = FALSE) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    jj <- 1:nmax
    ee <- q^(2 * jj)
    out <- jj^3 * ee/(1 - ee)
    if (give) {
        out <- cumsum(out)
    }
    else {
        out <- sum(out)
    }
    return(massage((pi/p1)^4 * (1/12 + 20 * out), tol = tol))
}
"g2.fun.lambert" <-
function (b, nmax = 50, tol = 1e-10, strict = TRUE) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q.sq <- exp(2 * pi * (1i) * tau)
    s <- 0
    q.sq.power.n <- q.sq
    for (n in 1:nmax) {
        s.new <- s + n^3 * q.sq.power.n/(1 - q.sq.power.n)
        if (near.match(s, s.new)) {
            return(massage((pi/p1)^4 * (1/12 + 20 * s), tol = tol))
        }
        s <- s.new
        q.sq.power.n <- q.sq.power.n * q.sq
    }
    if (strict) {
        stop("series not converged")
    }
    else {
        warning("series not converged.  Partial sum returned")
        return(massage((pi/p1)^4 * (1/12 + 20 * s), tol = tol))
    }
}
"g2.fun.vectorized" <-
function (b, nmax = 50, tol = 1e-10, give = FALSE) 
{
    if (is.vector(b)) {
        p1 <- rep(1, length(b))
        tau <- b
    }
    else {
        p1 <- b[, 1]
        tau <- b[, 2]/p1
    }
    if (any(Im(tau) < 0)) {
        stop("Im(tau)<0")
    }
    q2 <- exp(pi * (0+2i) * tau)
    jj <- 1:nmax
    out <- outer(jj, q2, function(n, q) {
        n^3 * q^n/(1 - q^n)
    })
    if (give) {
        out <- apply(out, 2, cumsum)
    }
    else {
        out <- apply(out, 2, sum)
    }
    return(massage((pi/p1)^4 * (1/12 + 20 * out), tol = tol))
}
"g3.fun" <-
function (b, use.first = TRUE, ...) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    jj2 <- theta2(0, q = q, ...)^4
    jj3 <- theta3(0, q = q, ...)^4
    jj4 <- theta4(0, q = q, ...)^4
    if (use.first) {
        return((pi/(2 * p1))^6 * ((8/27) * (jj2^3 + jj3^3) - 
            (4/9) * (jj2 + jj3) * jj2 * jj3))
    }
    else {
        return((4/27) * (pi/(2 * p1))^6 * (jj2 + jj3) * (jj3 + 
            jj4) * (jj4 - jj2))
    }
}
"g3.fun.direct" <-
function (b, nmax = 50, tol = 1e-10) 
{
    if (length(b) == 1) {
        b <- c(1, b)
    }
    jj <- expand.grid(-nmax:nmax, -nmax:nmax)
    jj <- jj[-(2 * nmax * nmax + 2 * nmax + 1), ]
    jj <- as.matrix(jj)
    return(massage(140 * sum(1/(2 * jj %*% b)^6), tol = tol))
}
"g3.fun.divisor" <-
function (b, nmax = 50, strict = TRUE, tol = 1e-10) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    s <- 0
    for (n in 1:nmax) {
        s.new <- s + divisor(n, 5) * q^(2 * n)
        if (near.match(s, s.new)) {
            return(massage((pi/p1)^6 * (1/216 - 7/3 * s), tol = tol))
        }
        s <- s.new
    }
    if (strict) {
        stop("series not converged")
    }
    else {
        warning("series not converged.  Partial sum returned")
        return(massage((2 * pi/p1)^6 * (1/216 - 7/3 * s), tol = tol))
    }
}
"g3.fun.fixed" <-
function (b, nmax = 50, tol = 1e-10, give = FALSE) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q <- exp(pi * (1i) * tau)
    jj <- 1:nmax
    ee <- q^(2 * jj)
    out <- jj^5 * ee/(1 - ee)
    if (give) {
        out <- cumsum(out)
    }
    else {
        out <- sum(out)
    }
    return(massage((pi/p1)^6 * (1/216 - 7/3 * out), tol = tol))
}
"g3.fun.lambert" <-
function (b, nmax = 50, tol = 1e-10, strict = TRUE) 
{
    jj <- p1.tau(b)
    p1 <- jj$p1
    tau <- jj$tau
    q.sq <- exp(2 * pi * (1i) * tau)
    s <- 0
    q.sq.power.n <- q.sq
    for (n in 1:nmax) {
        s.new <- s + n^5 * q.sq.power.n/(1 - q.sq.power.n)
        if (near.match(s, s.new)) {
            return(massage((pi/p1)^6 * (1/216 - 7/3 * s), tol = tol))
        }
        s <- s.new
        q.sq.power.n <- q.sq.power.n * q.sq
    }
    if (strict) {
        stop("series not converged")
    }
    else {
        warning("series not converged.  Partial sum returned")
        return(massage((pi/p1)^4 * (1/12 + 20 * s), tol = tol))
    }
}
"g3.fun.vectorized" <-
function (b, nmax = 50, tol = 1e-10, give = FALSE) 
{
    if (is.vector(b)) {
        p1 <- rep(1, length(b))
        tau <- b
    }
    else {
        p1 <- b[, 1]
        tau <- b[, 2]/p1
    }
    if (any(Im(tau) < 0)) {
        stop("Im(tau)<0")
    }
    q2 <- exp(pi * (0+2i) * tau)
    jj <- 1:nmax
    out <- outer(jj, q2, function(n, q) {
        n^5 * q^n/(1 - q^n)
    })
    if (give) {
        out <- apply(out, 2, cumsum)
    }
    else {
        out <- apply(out, 2, sum)
    }
    return(massage((pi/p1)^6 * (1/216 - 7/3 * out), tol = tol))
}
"g.fun" <-
function (b, ...) 
{
    c(g2 = g2.fun(b, ...), g3 = g3.fun(b, ...))
}
"H" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v = pi * u/(2 * K)
    return(theta1(v, m = m, ...))
}
"H1" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v = pi * u/(2 * K)
    return(theta2(v, m = m, ...))
}
"half.periods" <-
function (ignore = NULL, e = NULL, g = NULL, primitive = TRUE) 
{
    if (!xor(is.null(e), is.null(g))) {
        stop("supply exactly one of e, g")
    }
    if (is.null(e)) {
        e <- e1e2e3(g)
    }
    omega1 <- K.fun((e[2] - e[3])/(e[1] - e[3]))/sqrti(e[1] - 
        e[3])
    omega2 <- (1i)/sqrti(e[1] - e[3]) * K.fun(1 - (e[2] - e[3])/(e[1] - 
        e[3]))
    if (primitive) {
        return(as.primitive(c(omega1, omega2)))
    }
    else {
        return(c(omega1, omega2))
    }
}
"Im<-" <-
function (x, value) 
{
    if (is.complex(value)) {
        stop("RHS must be pure real")
    }
    if (all(value == 0)) {
        return(Re(x))
    }
    else {
        return(Re(x) + (1i) * value)
    }
}
"integrate.contour" <-
function (f, u, udash, ...) 
{
    myintegrate(function(x, ...) {
        f(u(x), ...) * udash(x)
    }, lower = 0, upper = 1, ...)
}
"integrate.segments" <-
function (f, points, close = TRUE, ...) 
{
    if (isTRUE(close)) {
        points <- c(points, points[1])
    }
    out <- 0
    for (i in 1:(length(points) - 1)) {
        u <- function(z) {
            points[i] + (points[i + 1] - points[i]) * z
        }
        udash <- function(z) {
            points[i + 1] - points[i]
        }
        out <- out + integrate.contour(f, u, udash, ...)
    }
    return(out)
}
"is.primitive" <-
function (p, n = 3, tol = 1e-05) 
{
    all(abs(p - as.primitive(p, n = n, tol = tol)) < tol)
}
"J" <-
function (tau, use.theta = TRUE, ...) 
{
    if (use.theta) {
        q <- exp(pi * (1i) * tau)
        jj.2 <- theta2(z = 0, q = q, ...)
        jj.3 <- theta3(z = 0, q = q, ...)
        jj.4 <- theta4(z = 0, q = q, ...)
        return((jj.2^8 + jj.3^8 + jj.4^8)^3/(jj.2 * jj.3 * jj.4)^8/54)
    }
    else {
        return(1/(1 - 27 * g3.fun(tau, ...)^2/g2.fun(tau, ...)^3))
    }
}
"K.fun" <-
function (m, strict = TRUE, maxiter = 7) 
{
    a.old <- 1
    b.old <- sqrti(1 - m)
    for (i in 1:maxiter) {
        a.new <- 0.5 * (a.old + b.old)
        b.new <- sqrti(a.old * b.old)
        if (near.match(a.new, a.old)) {
            return(pi/(2 * a.new))
        }
        a.old <- a.new
        b.old <- b.new
    }
    if (strict) {
        stop("iteration not stable")
    }
    else {
        warning("iteration not stable: partial result returned")
        return(pi/(2 * a.new))
    }
}
"lambda" <-
function (tau, ...) 
{
    q <- exp(pi * (1i) * tau)
    (theta2(z = 0, q = q, ...)/theta3(z = 0, q = q, ...))^4
}
"latplot" <-
function (p, n = 10, do.lines = TRUE, ...) 
{
    p1 <- p[1]
    p2 <- p[2]
    plot(lattice(p, n), xaxt = "n", yaxt = "n", bty = "n", pch = 16, 
        ...)
    axis(1, pos = 0, lwd = 2)
    axis(2, pos = 0, lwd = 2)
    slope1 <- Im(p1)/Re(p1)
    slope2 <- Im(p2)/Re(p2)
    int1 <- Im(p2) - Re(p2) * (slope1)
    int2 <- Im(p1) - Re(p1) * (slope2)
    if (do.lines) {
        for (u in -n:n) {
            if (is.finite(slope1)) {
                abline(u * int1, slope1, col = "gray")
            }
            else {
                abline(v = u * Re(p2), col = "gray")
            }
            if (is.finite(slope2)) {
                abline(u * int2, slope2, col = "gray")
            }
            else {
                abline(v = u * Re(p1), col = "gray")
            }
        }
    }
    points(Re(p1), Im(p1), pch = 16, cex = 3, col = "red")
    points(Re(p2), Im(p2), pch = 16, cex = 3, col = "green")
}
"lattice" <-
function (p, n) 
{
    outer(p[1] * (-n:n), p[2] * (-n:n), "+")
}
"lemniscatic" <-
function (...) 
{
    omega1 <- gamma(1/4)^2/(4 * sqrt(pi))
    omega2 <- (1i) * omega1
    Omega <- c(omega1, omega2)
    e <- c(1/2, 0, -1/2)
    names(e) <- c("e1", "e2", "e3")
    jj <- pi/4/omega1
    Eta <- c(jj, -jj * (1i), jj * (1i - 1))
    out <- list(Omega = Omega, q = exp(pi * (1i) * omega2/omega1), 
        e = e, g = c(g2 = 1, g3 = 0), Delta = 1, Eta = Eta, is.AnS = TRUE, 
        given = "d")
    class(out) <- "parameters"
    return(out)
}
"limit" <-
function (x, upper = quantile(Re(x), 0.99, na.rm = TRUE), lower = quantile(Re(x), 
    0.01, na.rm = TRUE), na = FALSE) 
{
    if (is.complex(x)) {
        return(Recall(Re(x), upper = upper, lower = lower, na = na) + 
            (1i) * Recall(Im(x), upper = upper, lower = lower, 
                na = na))
    }
    if (na) {
        x[x < lower] <- NA
        x[x > upper] <- NA
    }
    else {
        x <- pmax(x, lower)
        x <- pmin(x, upper)
    }
    return(x)
}
"massage" <-
function (z, tol = 1e-10) 
{
    if (length(z) == 1) {
        if (abs(Im(z)) < tol) {
            return(Re(z))
        }
        else {
            if (abs(Re(z)) < tol) {
                return((1i) * Im(z))
            }
            else {
                return(z)
            }
        }
    }
    Im(z[abs(Im(z)) < tol]) <- 0
    Re(z[abs(Re(z)) < tol]) <- 0
    if (all(Im(z) == 0)) {
        z <- Re(z)
    }
    return(z)
}
"mn" <-
function (z, p) 
{
    p1 <- p[1]
    p2 <- p[2]
    m <- (Re(z) * Im(p2) - Im(z) * Re(p2))/(Re(p1) * Im(p2) - 
        Im(p1) * Re(p2))
    n <- (Re(z) * Im(p1) - Im(z) * Re(p1))/(Re(p2) * Im(p1) - 
        Im(p2) * Re(p1))
    cbind(m, n)
}
"mob" <-
function (M, x) 
{
    (M[1] * x + M[3])/(M[2] * x + M[4])
}
"%mob%" <-
function (M, x) 
{
    mob(M, x)
}
"mobius" <-
function (n) 
{
    if (length(n) > 1) {
        return(sapply(n, mobius))
    }
    if (n == 1) {
        return(1)
    }
    jj <- table(factorize(n))
    if (any(jj > 1)) {
        return(as.integer(0))
    }
    else {
        return(as.integer((-1)^length(jj)))
    }
}
"myintegrate" <-
function (f, lower, upper, ...) 
{
    f.real <- function(x, ...) {
        Re(f(x, ...))
    }
    f.imag <- function(x, ...) {
        Im(f(x, ...))
    }
    jj.1 <- integrate(f.real, lower = lower, upper = upper, ...)
    jj.2 <- integrate(f.imag, lower = lower, upper = upper, ...)
    jj.1$value + (1i) * jj.2$value
}
"nc" <-
function (u, m, ...) 
{
    theta.n(u, m = m, ...)/theta.c(u, m = m, ...)
}
"nd" <-
function (u, m, ...) 
{
    theta.n(u, m = m, ...)/theta.d(u, m = m, ...)
}
"near.match" <-
function (x, y, tol = NULL) 
{
    if (is.null(tol)) {
        tol <- .Machine$double.eps * 2
    }
    return(isTRUE(all.equal(x, y, tol = tol)))
}
"newton_raphson" <-
function (initial, f, fdash, maxiter, give=TRUE, tol = .Machine$double.eps) 
{
    old.guess <- initial
    for (i in seq_len(maxiter)) {
        new.guess <- old.guess - f(old.guess)/fdash(old.guess)
        jj <- f(new.guess)
        if(is.na(jj) | is.infinite(jj)){
          break
        }
        if (near.match(new.guess, old.guess) | abs(jj) < tol) {
          if(give){
            return(list(root=new.guess,
                        f.root=jj,
                        iter=i))
          } else {
            return(new.guess)
          }
        }
        old.guess <- new.guess
    }
    stop("did not converge")
}

"nn" <-
function (u, m, ...) 
{
    theta.n(u, m = m, ...)/theta.n(u, m = m, ...)
}
"nome" <-
function (m) 
{
    K <- K.fun(m)
    Kdash <- K.fun(1 - m)
    return((exp(-pi * Kdash/K)))
}
"nome.k" <-
function (k) 
{
    K <- K.fun(m = sqrt(k))
    Kdash <- K.fun(sqrt(1 - k^2))
    return((exp(-pi * Kdash/K)))
}
"ns" <-
function (u, m, ...) 
{
    theta.n(u, m = m, ...)/theta.s(u, m = m, ...)
}
"P" <-
function (z, g = NULL, Omega = NULL, params = NULL, use.fpp = TRUE, 
    give.all.3 = FALSE, ...) 
{
    if (is.null(params)) {
        params <- parameters(g = g, Omega = Omega)
    }
    if (use.fpp) {
        z <- fpp(z, p = 2 * params$Omega)
    }
    e <- params$e
    q <- params$q
    omega <- params$Omega[1]
    v <- pi * z/(2 * omega)
    out1 <- e[1] + pi^2/(4 * omega^2) * (theta1.dash.zero.q(q) * 
        theta2(v, q = q, ...)/theta2(0, q = q, ...)/theta1(v, 
        q = q, ...))^2
    if (give.all.3) {
        out2 <- e[2] + pi^2/(4 * omega^2) * (theta1.dash.zero.q(q) * 
            theta3(v, q = q, ...)/theta3(0, q = q, ...)/theta1(v, 
            q = q, ...))^2
        out3 <- e[3] + pi^2/(4 * omega^2) * (theta1.dash.zero.q(q) * 
            theta4(v, q = q, ...)/theta4(0, q = q, ...)/theta1(v, 
            q = q, ...))^2
        return(drop(cbind(out1, out2, out3)))
    }
    else {
        attributes(out1) <- attributes(z)
        return(out1)
    }
}
"p1.tau" <-
function (b) 
{
    if (length(b) == 2) {
        p1 <- b[1]
        tau <- b[2]/b[1]
    }
    else {
        if (identical(ncol(b), 2)) {
            p1 <- b[1, ]
            tau <- b[2, ]/b[1, ]
        }
        else {
            p1 <- 1
            tau <- b
        }
    }
    if (any(Im(tau) < 0)) {
        warning("g2 and g3 not defined where Im(p2/p1)")
    }
    return(list(p1 = p1, tau = tau))
}
"parameters" <-
function (Omega = NULL, g = NULL, description = NULL) 
{
    jj <- c(!is.null(Omega), !is.null(g), !is.null(description))
    if (sum(jj) != 1) {
        stop("function must be supplied with exactly one argument")
    }
    if (!is.null(description)) {
        return(switch(description, equianharmonic = equianharmonic(), 
            lemniscatic = lemniscatic(), pseudolemniscatic = pseudolemniscatic(), 
            ))
    }
    if (is.null(Omega)) {
        given <- "g"
        g2 <- g[1]
        g3 <- g[2]
        Omega <- half.periods(g = g)
    }
    else {
        given <- "o"
        if (!is.primitive(Omega)) {
            warning("Omega supplied not a primitive pair of half periods.  Function converting Omega to a primitive pair ")
            Omega <- as.primitive(Omega)
        }
        g <- g.fun(Omega)
        g2 <- g[1]
        g3 <- g[2]
    }
    e <- e1e2e3(g, Omega = Omega)
    Delta <- g2^3 - 27 * g3^2
    omega1 <- Omega[1]
    omega2 <- Omega[2]
    omega3 <- -omega1 - omega2
    p1 <- 2 * omega1
    p2 <- 2 * omega2
    jj.q <- exp(pi * (1i) * omega2/omega1)
    eta1 <- -pi^2 * theta1dashdashdash(0, q = jj.q)/(12 * omega1 * 
        theta1.dash.zero.q(q = jj.q))
    eta2 <- omega2/omega1 * eta1 - pi * (1i)/(2 * omega1)
    eta3 <- -eta2 - eta1
    Eta <- c(eta1, eta2, eta3)
    out <- list(Omega = Omega, q = exp(pi * (1i) * (omega2)/omega1), 
        e = e, g = g, Delta = Delta, Eta = Eta, 
        is.AnS = FALSE, given = given)
    class(out) <- "parameters"
    return(out)
}
"Pdash" <-
function (z, g = NULL, Omega = NULL, params = NULL, use.fpp = TRUE, 
    ...) 
{
    if (is.null(params)) {
        params <- parameters(g = g, Omega = Omega)
    }
    if (use.fpp) {
        z <- fpp(z, p = 2 * params$Omega)
    }
    q <- params$q
    omega <- params$Omega[1]
    v <- pi * z/(2 * omega)
    out <- -pi^3/(4 * omega^3) * theta2(v, q = q, ...) * theta3(v, 
        q = q, ...) * theta4(v, q = q, ...) * theta1dash(0, q = q, 
        ...)^3/(theta2(0, q = q, ...) * theta3(0, q = q, ...) * 
        theta4(0, q = q, ...) * theta1(v, q = q, ...)^3)
    return(out)
}
"Pdash.laurent" <-
function (z, g = NULL, nmax = 80) 
{
    g2 <- g[1]
    g3 <- g[2]
    ckn <- ck(g = c(g2, g3), n = nmax)
    psum <- z * 0
    z.squared <- z * z
    zz <- z
    for (k in 2:nmax) {
        psum.new <- psum + zz * ckn[k] * (2 * k - 2)
        if (near.match(psum, psum.new) & ckn[k] > 0) {
            return(-2/z^3 + psum)
        }
        psum <- psum.new
        zz <- zz * z.squared
    }
    warning("series not converged.  See p636 for radius of convergence")
    return(-2/z^3 + psum)
}
"P.laurent" <-
function (z, g = NULL, tol = 0, nmax = 80) 
{
    g2 <- g[1]
    g3 <- g[2]
    ckn <- ck(g = c(g2, g3), n = nmax)
    psum <- z * 0
    z.squared <- z * z
    zz <- z.squared
    for (k in 2:nmax) {
        psum.new <- psum + zz * ckn[k]
        if (near.match(psum, psum.new, tol = tol) & abs(ckn[k]) > 
            0) {
            return(1/z^2 + psum)
        }
        psum <- psum.new
        zz <- zz * z.squared
    }
    warning("series not converged; partial sum returned.  See p636 for radius of convergence")
    return(1/z^2 + psum)
}
"P.pari" <-
function (z, Omega, pari.fun = "ellwp", numerical = TRUE) 
{
    attr <- attributes(z)
    z <- as.vector(z)
    a <- cbind(z, Omega[1], Omega[2])
    out <- NULL
    pari.complex <- function(x) {
        gsub("i", "*I", x)
    }
    for (i in 1:nrow(a)) {
        string <- paste("echo '", pari.fun, "([", pari.complex(2 * 
            a[i, 2]), ",", pari.complex(2 * a[i, 3]), "],", pari.complex(a[i, 
            1]), ")' | gp -q")
        jj <- gsub(" ", "", sub("\\*I", "i", system(string, intern = TRUE)))
        if (numerical) {
            jj <- as.complex(jj)
        }
        out <- c(out, jj)
    }
    attributes(out) <- attr
    return(out)
}
"primes" <-
function (n) 
{
    if ((M2 <- max(n)) <= 1) 
        return(numeric(0))
    x <- 1:M2
    x[1] <- 0
    p <- 1
    M <- floor(sqrt(M2))
    while ((p <- p + 1) <= M) if (x[p] != 0) 
        x[seq(p^2, n, p)] <- 0
    as.integer(x[x > 0])
}
"pseudolemniscatic" <-
function (...) 
{
    jj <- gamma(1/4)^2/(4 * sqrt(2 * pi))
    Omega <- c(jj * (1 - (1i)), jj * (1 + (1i)))
    e <- c(1/2, 0, -1/2) * (1i)
    names(e) <- c("e1", "e2", "e3")
    jj <- pi/4/Omega[1]
    Eta <- c(jj, -jj * (1i), jj * (1i - 1))
    out <- list(Omega = Omega, q = exp(pi * (1i) * Omega[2]/Omega[1]), 
        e = e, g = c(g2 = -1, g3 = 0), Delta = 1, Eta = Eta, 
        is.AnS = TRUE, given = "d")
    class(out) <- "parameters"
    return(out)
}
"residue" <- function(f, z0, r, O=z0, ...){
    if(r <= abs(z0-O)){
        warning("contour does not wrap round z0.  Either increase r or move O closer to z0")
    }

    if(is.complex(r)){
        warning('imaginary part of r discarded')
        r <- Re(r)
    }
        
    u <- function(x){O+r*exp(pi*2i*x)} # 0 <= x <= 1 
    udash <- function(x){r*pi*2i*exp(pi*2i*x)}
    integrate.contour(function(z,...){f(z,...)/(z-z0)},u,udash,...)/(pi*2i)
}
"Re<-" <-
function (x, value) 
{
    if (is.complex(value)) {
        stop("RHS must be pure real")
    }
    return((1i) * Im(x) + value)
}
"sc" <-
function (u, m, ...) 
{
    theta.s(u, m = m, ...)/theta.c(u, m = m, ...)
}
"sd" <-
function (u, m, ...) 
{
    theta.s(u, m = m, ...)/theta.d(u, m = m, ...)
}
"sigma" <-
function (z, g = NULL, Omega = NULL, params = NULL, use.theta = TRUE, 
    ...) 
{
    if (is.null(params)) {
        params <- parameters(g = g, Omega = Omega)
    }
    Omega <- params$Omega
    Eta <- params$Eta
    if (use.theta) {
        o <- Omega[1]
        q <- exp(pi * (1i) * Omega[2]/Omega[1])
        return(2 * o/pi * exp(Eta[1] * z^2/2/o) * theta1(pi * 
            z/2/o, q = q, ...)/theta1dash(0, q = q, ...))
    }
    jj <- fpp(z, 2 * Omega, give = TRUE)
    z <- jj$answer
    M <- jj$mn[, 1]
    N <- jj$mn[, 2]
    (-1)^(M + N + M * N) * Recall(z, params = params, use.theta = TRUE, 
        ...) * exp((z + M * Omega[1] + N * Omega[2]) * (2 * M * 
        Eta[1] + 2 * N * Eta[2]))
}
"sigmadash.laurent" <-
function (z, g = NULL, nmax = 8, give.error = FALSE) 
{
    g2 <- g[1]
    g3 <- g[2]
    attr <- attributes(z)
    z <- as.vector(z)
    jj <- amn(nmax)
    if (give.error) {
        minor.diag <- row(jj) == 1:nmax & col(jj) == nmax:1
        jj[!minor.diag] <- 0
    }
    m <- col(jj) - 1
    n <- row(jj) - 1
    non.z <- as.vector(jj * (g2/2)^m * (2 * g3)^n/factorial(4 * 
        m + 6 * n))
    power.z <- as.vector(4 * m + 6 * n)
    out <- outer(z, power.z, "^")
    out <- sweep(out, 2, non.z, "*")
    if (give.error) {
        out <- apply(out, 1, function(x) {
            sum(abs(x))
        })
    }
    else {
        out <- apply(out, 1, sum)
    }
    attributes(out) <- attr
    return(out)
}
"sigma.laurent" <-
function (z, g = NULL, nmax = 8, give.error = FALSE) 
{
    g2 <- g[1]
    g3 <- g[2]
    attr <- attributes(z)
    z <- as.vector(z)
    jj <- amn(nmax)
    if (give.error) {
        minor.diag <- row(jj) == 1:nmax & col(jj) == nmax:1
        jj[!minor.diag] <- 0
    }
    m <- col(jj) - 1
    n <- row(jj) - 1
    non.z <- as.vector(jj * (g2/2)^m * (2 * g3)^n/factorial(4 * 
        m + 6 * n + 1))
    power.z <- as.vector(4 * m + 6 * n + 1)
    out <- outer(z, power.z, "^")
    out <- sweep(out, 2, non.z, "*")
    if (give.error) {
        out <- apply(out, 1, function(x) {
            sum(abs(x))
        })
    }
    else {
        out <- apply(out, 1, sum)
    }
    attributes(out) <- attr
    return(out)
}
"sn" <-
function (u, m, ...) 
{
    theta.s(u, m = m, ...)/theta.n(u, m = m, ...)
}
"sqrti" <-
function (x) 
{
    if (is.complex(x)) {
        return(sqrt(x))
    }
    if (any(x < 0)) {
        return(sqrt(x + 0i))
    }
    else {
        return(sqrt(x))
    }
}
"ss" <-
function (u, m, ...) 
{
    1
}
"Theta" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v = pi * u/(2 * K)
    return(theta4(v, m = m, ...))
}
"theta.00" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + q^(n^2) * cos(2 * z * n)
        if (near.match(out, out.new)) {
            ans <- 1 + 2 * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta.01" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + (-1)^n * q^(n^2) * cos(2 * z * n)
        if (near.match(out, out.new)) {
            ans <- 1 + 2 * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta1" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + (-1)^n * q^(n * (n + 1)) * sin((2 * 
            n + 1) * z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"Theta1" <-
function (u, m, ...) 
{
    K <- K.fun(m)
    v = pi * u/(2 * K)
    return(theta3(v, m = m, ...))
}
"theta.10" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + q^(n * (n + 1)) * cos((2 * n + 1) * 
            z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta.11" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + (-1)^n * q^(n * (n + 1)) * sin((2 * 
            n + 1) * z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta1dash" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + (-1)^n * q^(n * (n + 1)) * (2 * n + 
            1) * cos((2 * n + 1) * z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta1dashdash" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out - (-1)^n * q^(n * (n + 1)) * (2 * n + 
            1)^2 * sin((2 * n + 1) * z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta1dashdashdash" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out - (-1)^n * q^(n * (n + 1)) * (2 * n + 
            1)^3 * cos((2 * n + 1) * z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta1.dash.zero" <-
function (m, ...) 
{
    theta2(0, m = m, ...) * theta3(0, m = m, ...) * theta4(0, 
        m = m, ...)
}
"theta1.dash.zero.q" <-
function (q, ...) 
{
    theta2(0, q = q, ...) * theta3(0, q = q, ...) * theta4(0, 
        q = q, ...)
}
"theta2" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 0:maxiter) {
        out.new <- out + q^(n * (n + 1)) * cos((2 * n + 1) * 
            z)
        if (near.match(out, out.new)) {
            ans <- 2 * q^(1/4) * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta3" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + q^(n^2) * cos(2 * z * n)
        if (near.match(out, out.new)) {
            ans <- 1 + 2 * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta4" <-
function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE, 
    maxiter = 30) 
{
    if (!xor(is.null(m), is.null(q))) {
        stop("supply exactly one of m, q")
    }
    if (is.null(q)) {
        q <- nome(m)
    }
    out <- 0
    for (n in 1:maxiter) {
        out.new <- out + (-1)^n * q^(n^2) * cos(2 * z * n)
        if (near.match(out, out.new)) {
            ans <- 1 + 2 * out
            if (give.n) {
                return(list(iterations = n, ans = ans))
            }
            else {
                return(ans)
            }
        }
        out <- out.new
    }
    stop("maximum iterations reached")
}
"theta.c" <-
function (u, m, method = "16.36.6", ...) 
{
    switch(method, "16.36.6" = e16.36.6b(u, m, ...), "16.36.6b" = e16.36.6b(u, 
        m, ...), "16.36.1b" = e16.36.6b(u, m, ...), "16.37.2" = e16.37.2(u, 
        m, ...), "16.38.2" = e16.38.2(u, m, ...), stop("method not recognized"))
}
"theta.d" <-
function (u, m, method = "16.36.7", ...) 
{
    switch(method, "16.36.7" = e16.36.7a(u, m, ...), "16.36.7a" = e16.36.7a(u, 
        m, ...), "16.36.2a" = e16.36.7a(u, m, ...), "16.37.3" = e16.37.3(u, 
        m, ...), "16.38.3" = e16.38.3(u, m, ...), stop("method not recognized"))
}
"theta.n" <-
function (u, m, method = "16.36.7", ...) 
{
    switch(method, "16.36.7" = e16.36.7b(u, m, ...), "16.36.7b" = e16.36.7b(u, 
        m, ...), "16.36.2b" = e16.36.7b(u, m, ...), "16.37.4" = e16.37.4(u, 
        m, ...), "16.38.4" = e16.38.4(u, m, ...), stop("method not recognized"))
}
"theta.s" <-
function (u, m, method = "16.36.6", ...) 
{
    switch(method, "16.36.6" = e16.36.6a(u, m, ...), "16.36.6a" = e16.36.6a(u, 
        m, ...), "16.36.1a" = e16.36.6a(u, m, ...), "16.37.1" = e16.37.1(u, 
        m, ...), "16.38.1" = e16.38.1(u, m, ...), stop("method not recognized"))
}
"totient" <-
function (n) 
{
    if (length(n) > 1) {
        return(sapply(n, match.fun(sys.call()[[1]])))
    }
    as.integer(n * prod(1 - 1/unique(factorize(n))))
}
"unimodular" <-
function (n) 
{
    if (n == 1) {
        return(array(diag(2), c(2, 2, 1)))
    }
    out <- do.call("cbind", sapply(0:n, farey))
    out <- unique(out, MARGIN = 2)
    dim(out) <- c(2, 2, length(out)/4)
    return(out[2:1, , ])
}
"unimodularity" <-
function (n, o, FUN, ...) 
{
    u <- unimodular(n)
    N <- dim(u)[3]
    out <- rep(0, N)
    for (i in 1:N) {
        out[i] <- FUN(drop(u[, , i] %*% o), ...)
    }
    return(out)
}
"view" <-
function (x, y, z, scheme = 0, real.contour = TRUE, imag.contour = real.contour, 
    default = 0, col = "black", r0 = 1, power = 1, show.scheme = FALSE, 
    ...) 
{
    if (is.numeric(scheme)) {
        f <- switch(as.character(scheme), "0" = function(z) {
            u <- 2/pi * atan(Mod(z)/r0)
            s0 <- Re(u * 0 + 1)
            s0[u > 0.5] <- (2 * (1 - u[u > 0.5]))^power
            v0 <- Re(u * 0 + 1)
            v0[u < 0.5] <- (2 * u[u < 0.5])^power
            return(hsv(h = scale(Arg(z)), s = s0, v = v0))
        }, "1" = function(z) {
            hsv(h = scale(Arg(z)), s = scale(Im(z)), v = 1)
        }, "2" = function(z) {
            hsv(h = g(Arg(z)), s = scale(abs(z)), v = 1)
        }, "3" = function(z) {
            hsv(h = scale(Re(z)), s = 1, v = scale(Mod(z))^power)
        }, "4" = function(z) {
            hsv(h = 0.4, s = 1, v = scale(Arg(z))^power)
        }, "5" = function(z) {
            hsv(h = 0.4, s = 0, v = 0.5 + 0.5 * (Im(z) > 0))
        }, "6" = function(z) {
            hsv(h = 0.4, s = 1, v = 0.5 + 0.5 * (Im(z) > 0))
        }, "7" = function(z) {
            hsv(h = scale(Re(z))^power, s = 1, v = scale(Mod(z))^power)
        }, "8" = function(z) {
            hsv(h = wrap(Arg(z)))
        }, "9" = function(z) {
            hsv(h = wrap(Arg(z)), v = scale(Mod(z))^power)
        }, "10" = function(z) {
            hsv(h = wrap(Arg(z)), v = scale(exp(-Mod(z))))
        }, "11" = function(z) {
            hsv(h = wrap(Arg(z)), s = scale(Mod(z))^power)
        }, "12" = function(z) {
            hsv(h = wrap(Arg(z)), s = scale(exp(-Mod(z))))
        }, "13" = function(z) {
            hsv(h = 0.3, s = 1, v = (floor(Re(z)) + floor(Im(z)))%%2)
        }, "14" = function(z) {
            hsv(h = wrap(Arg(z)), s = 1, v = (floor(Re(z)) + 
                floor(Im(z)))%%2)
        }, "15" = function(z) {
            hsv(h = wrap(Arg(z)), s = 1, v = 0.4 + 0.4 * (floor(Re(z)) + 
                floor(Im(z)))%%2)
        }, "16" = function(z) {
            hcl(h = 360 * wrap(Arg(z)), l = 100 * scale(Mod(z))^power)
        }, "17" = function(z) {
            hcl(h = 360 * wrap(Arg(z)), c = 100 * scale(Mod(z))^power)
        }, "18" = function(z) {
            rgb(red = scale(Re(z)), green = 1 - scale(Re(z))^power, 
                blue = scale(Im(z))^power)
        }, "19" = function(z) {
            rgb(red = scale(Re(z)), green = scale(Im(z))^power, blue = 0)
        }, function(z) {
            hsv(s = 0, v = 1)
        })
    }
    else {
        f <- scheme
        environment(f) <- environment()
    }
    if (show.scheme) {
        return(f)
    }
    jj <- z
    jj[] <- (1:length(z))/length(z)
    jj <- Re(jj)
    breakup <- function(x) {
        ifelse(x > 1/2, 3/2 - x, 1/2 - x)
    }
    g <- function(x) {
        0.5 + atan(x)/pi
    }
    scale <- function(x) {
        (x - min(x))/(max(x) - min(x))
    }
    wrap <- function(x) {
        1/2 + x/(2 * pi)
    }
    if (!is.na(default)) {
        z[is.na(z)] <- default
        z[is.infinite(z)] <- default
    }
    suppressWarnings(image(x, y, z = jj, col = f(z), asp = 1, ...))
    if (real.contour) {
      suppressWarnings(contour(x, y, Re(z), add = TRUE, lty = 1, col = col, 
                               ...))
    }
    if (imag.contour) {
      suppressWarnings(contour(x, y, Im(z), add = TRUE, lty = 2, col = col, ...))
    }
}
"zeta" <-
function (z, g = NULL, Omega = NULL, params = NULL, use.fpp = TRUE, 
    ...) 
{
    if (is.null(params)) {
        params <- parameters(g = g, Omega = Omega)
    }
    Omega <- params$Omega
    Eta <- params$Eta
    if (use.fpp) {
        jj <- fpp(z, 2 * Omega, give = TRUE)
        z <- jj$answer
        M <- jj$mn[, 1]
        N <- jj$mn[, 2]
        return(Recall(z, params = params, use.fpp = FALSE, ...) + 
            2 * M * Eta[1] + 2 * N * Eta[2])
    }
    else {
        o <- Omega[1]
        q <- exp(pi * (1i) * Omega[2]/Omega[1])
        jj <- pi * z/(2 * o)
        return(z * Eta[1]/o + pi * theta1dash(jj, q = q, ...)/(2 * 
            o * theta1(jj, q = q, ...)))
    }
}
"zeta.laurent" <-
function (z, g = NULL, nmax = 80) 
{
    g2 <- g[1]
    g3 <- g[2]
    ckn <- ck(g = c(g2, g3), n = nmax)
    psum <- z * 0
    z.squared <- z * z
    zz <- z * z.squared
    for (k in 2:nmax) {
        psum.new <- psum + zz * ckn[k]/(2 * k - 1)
        if (near.match(psum, psum.new) & abs(ckn[k]) > 0) {
            return(1/z - psum)
        }
        psum <- psum.new
        zz <- zz * z.squared
    }
    warning("series not converged.  See p636 for radius of convergence")
    return(1/z - psum)
}

# Following code for sigma1 sigma2 and sigma3 commented out
# because output does not agree with Mathematica.  I'll investigate
# before adding them to the package. rksh.

#sigma1 <- function(u,params){
#  o1 <- params$Omega[1]
#  eta1 <- params$Eta[1]
#exp(-eta1*u)*
#  sigma(o1+u,params=params)/
#    sigma(o1,params=params)
#}
#
#sigma2 <- function(u,params){
#  o2 <-  -sum(params$Omega)
#  eta2 <- params$Eta[3]
#exp(-eta2*u)*
#  sigma(o2+u,params=params)/
#    sigma(o2,params=params)
#}
#
#sigma3 <- function(u,params){
#  o3 <- params$Omega[2]
#  eta3 <- params$Eta[2]
#  exp(-eta3*u)*
#    sigma(o3+u,params=params)/
#      sigma(o3,params=params)
#}
#
#

