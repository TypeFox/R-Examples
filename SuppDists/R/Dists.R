dFriedman <-
function (x, r, N, log = FALSE) 
{
    M <- max(length(x), length(r), length(N))
    x <- rep(x, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    value <- .C("dFriedmanR", as.double(x), as.integer(r), as.integer(N), 
        as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dghyper <-
function (x, a, k, N, log = FALSE) 
{
    M <- max(length(x), length(a), length(k), length(N))
    x <- rep(x, length.out = M)
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("dghyperR", as.integer(x), as.double(a), as.double(k), 
        as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dinvGauss <-
function (x, nu, lambda, log = FALSE) 
{
    N <- max(length(x), length(nu), length(lambda))
    x <- rep(x, length.out = N)
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    value <- .C("dinvGaussR", as.double(x), as.double(nu), as.double(lambda), 
        as.integer(N), lambda = double(N),PACKAGE="SuppDists")$lambda
    if (log == TRUE) 
        value <- log(value)
    value
}
dJohnson <-
function (x, parms, log = FALSE) 
{
    tfun <- function(x) if (x == "SN") 
        1
    else if (x == "SL") 
        2
    else if (x == "SU") 
        3
    else 4
    vecFromList <- function(item, aList) {
        if (!is.list(aList[[1]])) 
            return(aList[[item]])
        else {
            tVec <- vector(length = 0)
            for (i in 1:length(aList)) {
                tVec <- append(tVec, (aList[[i]])[[item]])
            }
        }
        tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- max(length(gamma), length(x))
    x <- rep(x, length.out = N)
    gamma <- rep(gamma, length.out = N)
    delta <- rep(delta, length.out = N)
    xi <- rep(xi, length.out = N)
    lambda <- rep(lambda, length.out = N)
    type <- rep(type, length.out = N)
    value <- .C("dJohnsonR", as.double(x), as.double(gamma), 
        as.double(delta), as.double(xi), as.double(lambda), as.integer(type), 
        as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dKendall <-
function (x, N, log = FALSE) 
{
    M <- max(length(x), length(N))
    x <- rep(x, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("dKendallR", as.integer(N), as.double(x), as.integer(M), 
        val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dKruskalWallis <-
function (x, c, N, U, log = FALSE) 
{
    M <- max(length(x), length(c), length(N), length(U))
    x <- rep(x, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    value <- .C("dKruskalWallisR", as.double(x), as.integer(c), 
        as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
        val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dmaxFratio <-
function (x, df, k, log = FALSE) 
{
    if (log == TRUE) 
        p <- exp(p)
    N <- max(length(x), length(df), length(k))
    x <- rep(x, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    .C("dmaxFratioR", as.double(x), as.integer(df), as.integer(k), 
        as.integer(N), val = double(N),PACKAGE="SuppDists")$val
}
dNormScore <-
function (x, c, N, U, log = FALSE) 
{
    M <- max(length(x), length(c), length(N), length(U))
    x <- rep(x, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    value <- .C("dKruskalWallisR", as.double(x), as.integer(c), 
        as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
        val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dPearson <-
function (x, N, rho = 0, log = FALSE) 
{
    M <- max(length(x), length(rho), length(N))
    x <- rep(x, length.out = M)
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("dcorrR", as.double(x), as.double(rho), as.integer(N), 
        as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
dSpearman <-
function (x, r, log = FALSE) 
{
    M <- max(length(x), length(r))
    x <- rep(x, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    value <- .C("dFriedmanR", as.double(x), as.integer(r), as.integer(N), 
        as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    if (log == TRUE) 
        value <- log(value)
    value
}
JohnsonFit <-
function (t, moment = "quant")
{
    firstChar=substring(moment,1,1)
    if (firstChar=="f") {
        mom <- moments(t)
        mu <- mom[[1]]
        sigma <- mom[[2]]
        skew <- mom[[3]]
        kurt <- mom[[4]]
        value <- .C("JohnsonMomentFitR", as.double(mu), as.double(sigma), 
            as.double(skew), as.double(kurt), gamma = double(1), 
            delta = double(1), xi = double(1), lambda = double(1), 
            type = integer(1),PACKAGE="SuppDists")
    }
    else if (firstChar=="u") {
      mu<-t[1]
      sigma<-sqrt(t[2])
      skew<-t[3]/sigma^3
      kurt<-(t[4]/t[2]^2)-3
      value <- .C("JohnsonMomentFitR", as.double(mu), as.double(sigma), 
          as.double(skew), as.double(kurt), gamma = double(1), 
          delta = double(1), xi = double(1), lambda = double(1), 
          type = integer(1),PACKAGE="SuppDists")
     }
    else if (firstChar=="q") {
        input <- quantile(t, probs = c(0.05, 0.206, 0.5, 0.794, 
            0.95), names = FALSE)
        x5 <- input[[1]]
        x20.6 <- input[[2]]
        x50 <- input[[3]]
        x79.4 <- input[[4]]
        x95 <- input[[5]]
        value <- .C("JohnsonFitR", as.double(x95), as.double(x79.4), 
            as.double(x50), as.double(x20.6), as.double(x5), 
            gamma = double(1), delta = double(1), xi = double(1), 
            lambda = double(1), type = integer(1),PACKAGE="SuppDists")
    }
    else return(NA)
    types <- c("SN", "SL", "SU", "SB")
    list(gamma = value$gamma, delta = value$delta, xi = value$xi, 
        lambda = value$lambda, type = types[value$type])
}
makeStatList <-
function (head, mn, med, var, mod, third, fourth, dig) 
{
    sd <- sqrt(var)
    skew <- sign(third) * abs(third)/sd^3
    kurt <- -3 + fourth/var^2
    pskew <- (mn - mod)/sd
    if (dig > 0) {
        mn <- round(mn, digits = dig)
        med <- round(med, digits = dig)
        mod <- round(mod, digits = dig)
        var <- round(var, digits = dig)
        sd <- round(sd, digits = dig)
        third <- round(third, digits = dig)
        fourth <- round(fourth, digits = dig)
        pskew <- round(pskew, digits = dig)
        skew <- round(skew, digits = dig)
        kurt <- round(kurt, digits = dig)
    }
    theList <- list(Mean = mn, Median = med, Mode = mod, Variance = var, 
        SD = sd, ThirdCentralMoment = third, FourthCentralMoment = fourth, 
        PearsonsSkewness...mean.minus.mode.div.SD = pskew, Skewness...sqrtB1 = skew, 
        Kurtosis...B2.minus.3 = kurt)
    c(head, theList)
}
moments <-
function (x) 
{
    N <- length(x)
    v <- ((N - 1)/N) * var(x)
    sigma <- sqrt(v)
    m3 <- (sum((x - mean(x))^3))/N
    skew <- m3/sigma^3
    m4 <- (sum((x - mean(x))^4))/N
    kurt <- (m4/v^2) - 3
    c(mean = mean(x), sigma = sigma, skew = skew, kurt = kurt)
}
normOrder <-
function (N) 
{
    N <- if (length(N) > 1) 
        length(N)
    else N
    M <- N%/%2
    value <- .C("normOrdR", val = double(M), as.integer(N), as.integer(M),PACKAGE="SuppDists")$val
    if (0 == N%%2) 
        c(-value, rev(value))
    else c(-value, 0, rev(value))
}
pFriedman <-
function (q, r, N, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(r), length(N))
    q <- rep(q, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pFriedmanR", as.double(q), as.integer(r), 
            as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uFriedmanR", as.double(q), as.integer(r), 
            as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pghyper <-
function (q, a, k, N, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(a), length(k), length(N))
    q <- rep(q, length.out = M)
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pghyperR", as.integer(q), as.double(a), 
            as.double(k), as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("ughyperR", as.integer(q), as.double(a), 
            as.double(k), as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pinvGauss <-
function (q, nu, lambda, lower.tail = TRUE, log.p = FALSE) 
{
    N <- max(length(q), length(nu), length(lambda))
    q <- rep(q, length.out = N)
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    if (lower.tail == TRUE) {
        value <- .C("pinvGaussR", as.double(q), as.double(nu), 
            as.double(lambda), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uinvGaussR", as.double(q), as.double(nu), 
            as.double(lambda), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pJohnson <-
function (q, parms, lower.tail = TRUE, log.p = FALSE) 
{
    tfun <- function(x) if (x == "SN") 
        1
    else if (x == "SL") 
        2
    else if (x == "SU") 
        3
    else 4
    vecFromList <- function(item, aList) {
        if (!is.list(aList[[1]])) 
            return(aList[[item]])
        else {
            tVec <- vector(length = 0)
            for (i in 1:length(aList)) {
                tVec <- append(tVec, (aList[[i]])[[item]])
            }
        }
        tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- max(length(gamma), length(q))
    q <- rep(q, length.out = N)
    gamma <- rep(gamma, length.out = N)
    delta <- rep(delta, length.out = N)
    xi <- rep(xi, length.out = N)
    lambda <- rep(lambda, length.out = N)
    type <- rep(type, length.out = N)
    if (lower.tail == TRUE) {
        value <- .C("pJohnsonR", as.double(q), as.double(gamma), 
            as.double(delta), as.double(xi), as.double(lambda), 
            as.integer(type), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uJohnsonR", as.double(q), as.double(gamma), 
            as.double(delta), as.double(xi), as.double(lambda), 
            as.integer(type), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pKendall <-
function (q, N, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(N))
    q <- rep(q, length.out = M)
    N <- rep(N, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pKendallR", as.integer(N), as.double(q), 
            as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uKendallR", as.integer(N), as.double(q), 
            as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pKruskalWallis <-
function (q, c, N, U, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(c), length(N), length(U))
    q <- rep(q, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pKruskalWallisR", as.double(q), as.integer(c), 
            as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
            val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uKruskalWallisR", as.double(q), as.integer(c), 
            as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
            val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pmaxFratio <-
function (q, df, k, lower.tail = TRUE, log.p = FALSE) 
{
    N <- max(length(q), length(df), length(k))
    q <- rep(q, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    if (lower.tail == TRUE) {
        value <- .C("pmaxFratioR", as.double(q), as.integer(df), 
            as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("umaxFratioR", as.double(q), as.integer(df), 
            as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pNormScore <-
function (q, c, N, U, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(c), length(N), length(U))
    q <- rep(q, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pKruskalWallisR", as.double(q), as.integer(c), 
            as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
            val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uKruskalWallisR", as.double(q), as.integer(c), 
            as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
            val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pPearson <-
function (q, N, rho = 0, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(rho), length(N))
    q <- rep(q, length.out = M)
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pcorrR", as.double(q), as.double(rho), as.integer(N), 
            as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("ucorrR", as.double(q), as.double(rho), as.integer(N), 
            as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
pSpearman <-
function (q, r, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(r))
    q <- rep(q, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C("pFriedmanR", as.double(q), as.integer(r), 
            as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
        value <- .C("uFriedmanR", as.double(q), as.integer(r), 
            as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}
qFriedman <-
function (p, r, N, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(r), length(N))
    p <- rep(p, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    .C("qFriedmanR", as.double(p), as.integer(r), as.integer(N), 
        as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
}
qghyper <-
function (p, a, k, N, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(a), length(k), length(N))
    p <- rep(p, length.out = M)
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("qghyperR", as.double(p), as.double(a), as.double(k), 
        as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    value
}
qinvGauss <-
function (p, nu, lambda, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    N <- max(length(p), length(nu), length(lambda))
    p <- rep(p, length.out = N)
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    .C("qinvGaussR", as.double(p), as.double(nu), as.double(lambda), 
        as.integer(N), value = double(N),PACKAGE="SuppDists")$value
}
qJohnson <-
function (p, parms, lower.tail = TRUE, log.p = FALSE) 
{
    tfun <- function(x) if (x == "SN") 
        1
    else if (x == "SL") 
        2
    else if (x == "SU") 
        3
    else 4
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    vecFromList <- function(item, aList) {
        if (!is.list(aList[[1]])) 
            return(aList[[item]])
        else {
            tVec <- vector(length = 0)
            for (i in 1:length(aList)) {
                tVec <- append(tVec, (aList[[i]])[[item]])
            }
        }
        tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- max(length(gamma), length(p))
    p <- rep(p, length.out = N)
    gamma <- rep(gamma, length.out = N)
    delta <- rep(delta, length.out = N)
    xi <- rep(xi, length.out = N)
    lambda <- rep(lambda, length.out = N)
    type <- rep(type, length.out = N)
    .C("qJohnsonR", as.double(p), as.double(gamma), as.double(delta), 
        as.double(xi), as.double(lambda), as.integer(type), as.integer(N), 
        val = double(N),PACKAGE="SuppDists")$val
}
qKendall <-
function (p, N, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(N))
    p <- rep(p, length.out = M)
    N <- rep(N, length.out = M)
    .C("qKendallR", as.integer(N), as.double(p), as.integer(M), 
        val = double(M),PACKAGE="SuppDists")$val
}
qKruskalWallis <-
function (p, c, N, U, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(c), length(N), length(U))
    p <- rep(p, length.out = M)
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    .C("qKruskalWallisR", as.double(p), as.integer(c), as.integer(N), 
        as.double(U), as.integer(Ns), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
}
qmaxFratio <-
function (p, df, k, lower.tail = TRUE, log.p = FALSE) 
{
    if (lower.tail == FALSE) 
        p <- 1 - p
    if (log.p == TRUE) 
        p <- exp(p)
    N <- max(length(p), length(df), length(k))
    p <- rep(p, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    .C("qmaxFratioR", as.double(p), as.integer(df), as.integer(k), 
        as.integer(N), val = double(N),PACKAGE="SuppDists")$val
}
qNormScore <-
function (p, c, N, U, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(c), length(N), length(U))
    p <- rep(p, length.out = M)
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    .C("qKruskalWallisR", as.double(p), as.integer(c), as.integer(N), 
        as.double(U), as.integer(Ns), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
}
qPearson <-
function (p, N, rho = 0, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(rho), length(N))
    p <- rep(p, length.out = M)
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    .C("qcorrR", as.double(p), as.double(rho), as.integer(N), 
        as.integer(M), val = double(M),PACKAGE="SuppDists")$val
}
qSpearman <-
function (p, r, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p == TRUE) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    M <- max(length(p), length(r))
    p <- rep(p, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    .C("qFriedmanR", as.double(p), as.integer(r), as.integer(N), 
        as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
}
rFriedman <-
function (n, r, N) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- max(length(r), length(N))
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    .C("rFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
        as.integer(n), as.integer(M), value = double(n),PACKAGE="SuppDists")$value
}
rghyper <-
function (n, a, k, N) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    K <- max(length(a), length(k), length(N))
    a <- rep(a, length.out = K)
    k <- rep(k, length.out = K)
    N <- rep(N, length.out = K)
    .C("rghyperR", as.double(a), as.double(k), as.double(N), 
        as.integer(n), as.integer(K), value = double(n),PACKAGE="SuppDists")$value
}
rinvGauss <-
function (n, nu, lambda) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    N <- max(length(nu), length(lambda))
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    .C("rinvGaussR", as.double(nu), as.double(lambda), as.integer(n), 
        as.integer(N), value = double(n),PACKAGE="SuppDists")$value
}
rJohnson <-
function (n, parms) 
{
    tfun <- function(x) if (x == "SN") 
        1
    else if (x == "SL") 
        2
    else if (x == "SU") 
        3
    else 4
    vecFromList <- function(item, aList) {
        if (!is.list(aList[[1]])) 
            return(aList[[item]])
        else {
            tVec <- vector(length = 0)
            for (i in 1:length(aList)) {
                tVec <- append(tVec, (aList[[i]])[[item]])
            }
        }
        tVec
    }
    n <- if (length(n) > 1) 
        length(n)
    else n
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    M <- length(gamma)
    .C("rJohnsonR", as.double(gamma), as.double(delta), as.double(xi), 
        as.double(lambda), as.integer(type), as.integer(n), as.integer(M), 
        val = double(n),PACKAGE="SuppDists")$val
}
rKendall <-
function (n, N) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- length(N)
    .C("rKendallR", as.integer(N), as.integer(n), as.integer(M), 
        val = double(n),PACKAGE="SuppDists")$val
}
rKruskalWallis <-
function (n, c, N, U) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    .C("rKruskalWallisR", randArray = double(n), as.integer(n), 
        as.integer(M), as.integer(c), as.integer(N), as.double(U), 
        as.integer(Ns),PACKAGE="SuppDists" )$randArray
}
rmaxFratio <-
function (n, df, k) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- max(length(df), length(k))
    df <- rep(df, length.out = M)
    k <- rep(k, length.out = M)
    .C("rmaxFratioR", as.integer(df), as.integer(k), as.integer(n), 
        as.integer(M), value = double(n),PACKAGE="SuppDists")$value
}
rMWC1019 <-
function (n, new.start = FALSE, seed = 556677) 
{
    n <- if (length(n) == 1) 
        n
    else length(n)
    .C("MWC1019R", val = double(n), as.integer(n), as.integer(new.start), 
        as.integer(seed),PACKAGE="SuppDists")$val
}
rNormScore <-
function (n, c, N, U) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    .C("rKruskalWallisR", randArray = double(n), as.integer(n), 
        as.integer(M), as.integer(c), as.integer(N), as.double(U), 
        as.integer(Ns),PACKAGE="SuppDists" )$randArray
}
rPearson <-
function (n, N, rho = 0) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- max(length(rho), length(N))
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    .C("rcorrR", as.double(rho), as.integer(N), as.integer(n), 
        as.integer(M), val = double(n),PACKAGE="SuppDists")$val
}
rSpearman <-
function (n, r) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    M <- length(r)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    .C("rFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
        as.integer(n), as.integer(M), value = double(n),PACKAGE="SuppDists")$value
}
rziggurat <-
function (n, normal = TRUE, new.start = FALSE, seed = 556677) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
   .C("ziggR", val = double(n), as.integer(n), as.integer(normal), 
        as.integer(new.start), as.integer(seed),PACKAGE="SuppDists")$val
}
sFriedman <-
function (r, N) 
{
    M <- max(length(r), length(N))
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    value <- .C("sFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
        as.integer(M), mn = double(M), med = double(M), mod = double(M), 
        var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Friedman's chi-square", r = r, N = N)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, -1)
}
sghyper <-
function (a, k, N) 
{
    M <- max(length(a), length(k), length(N))
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("sghyperR", as.double(a), as.double(k), as.double(N), 
        as.integer(M), mn = double(M), med = double(M), mod = double(M), 
        var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Generalized Hypergeometric", a = a, 
        k = k, N = N)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, -1)
}
sinvGauss <-
function (nu, lambda) 
{
    N <- max(length(nu), length(lambda))
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    med <- qinvGauss(0.5, nu, lambda)
	nu[nu<=0]<-NA
	lambda[lambda<=0]<-NA
    factor <- (nu^2)/lambda
    var <- nu * factor
    k3 <- 3 * var * factor
    k4 <- 5 * k3 * factor
    mod <- -1.5 * factor + nu * sqrt(1 + 2.25 * (nu/lambda)^2)
    third <- k3
    fourth <- k4 + 3 * var^2
    aList <- list(title = "Inverse Gaussian", nu = nu, lambda = lambda)
    makeStatList(aList, nu, med, var, mod, third, fourth, -1)
}
sJohnson <-
function (parms) 
{
    tfun <- function(x) if (x == "SN") 
        1
    else if (x == "SL") 
        2
    else if (x == "SU") 
        3
    else 4
    vecFromList <- function(item, aList) {
        if (!is.list(aList[[1]])) 
            return(aList[[item]])
        else {
            tVec <- vector(length = 0)
            for (i in 1:length(aList)) {
                tVec <- append(tVec, (aList[[i]])[[item]])
            }
        }
        tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- length(gamma)
    value <- .C("sJohnsonR", as.double(gamma), as.double(delta), 
        as.double(xi), as.double(lambda), as.integer(type), as.integer(N), 
        mn = double(N), med = double(N), mod = double(N), var = double(N), 
        third = double(N), fourth = double(N),PACKAGE="SuppDists")
    aList <- list(title = "Johnson Distribution", gamma = gamma, 
        delta = delta, xi = xi, lambda = lambda, type = type)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, -1)
}
sKendall <-
function (N) 
{
    M <- length(N)
    mn <- rep(0, length.out = M)
    med <- rep(0, length.out = M)
    mod <- rep(0, length.out = M)
    third <- rep(0, length.out = M)
    var <- (4 * N + 10)/(9 * N * (N - 1))
    fourth <- .C("fourthKendallR", as.integer(N), as.integer(M), 
        val = double(M),PACKAGE="SuppDists")$val
    aList <- list(title = "Kendall's Tau", N = N)
    makeStatList(aList, mn, med, var, mod, third, fourth, -1)
}
sKruskalWallis <-
function (c, N, U) 
{
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    value <- .C("sKruskalWallisR", as.integer(c), as.integer(n), 
        as.double(U), as.integer(Ns), as.integer(M), var = double(M), 
        mod = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    mn <- (c - 1)
    aList <- list(title = "Kruskal Wallis", c = c, N = n, U = U)
    median <- qKruskalWallis(0.5, c, n, U, Ns)
    makeStatList(aList, mn, median, value$var, value$mod, value$third, 
        value$fourth, -1)
}
smaxFratio <-
function (df, k) 
{
    N <- max(length(df), length(k))
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    value <- .C("smaxFratioR", as.integer(df), as.integer(k), 
        as.integer(N), mn = double(N), med = double(N), mod = double(N), 
        var = double(N), third = double(N), fourth = double(N),PACKAGE="SuppDists")
    aList <- list(title = "Maximum F ratio", df = df, k = k)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, 2)
}
sNormScore <-
function (c, N, U) 
{
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    value <- .C("sKruskalWallisR", as.integer(c), as.integer(n), 
        as.double(U), as.integer(Ns), as.integer(M), var = double(M), 
        mod = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    mn <- (c - 1)
    aList <- list(title = "Normal Scores", c = c, N = n, U = U)
    median <- qNormScore(0.5, c, n, U)
    makeStatList(aList, mn, median, value$var, value$mod, value$third, 
        value$fourth, -1)
}
sPearson <-
function (N, rho = 0) 
{
    M <- max(length(rho), length(N))
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("scorrR", as.double(rho), as.integer(N), as.integer(M), 
        mn = double(M), med = double(M), mod = double(M), var = double(M), 
        third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Correlation coefficient", rho = rho, 
        N = N)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, -1)
}
sSpearman <-
function (r) 
{
    M <- length(r)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    value <- .C("sFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
        as.integer(M), mn = double(M), med = double(M), mod = double(M), 
        var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Spearman's rho", r = r)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
        value$third, value$fourth, -1)
}
tghyper <-
function (a, k, N) 
{
    value <- .C("tghyperR", as.double(a), as.double(k), as.double(N), 
        strn =paste(rep(" ", 128), collapse=""),PACKAGE="SuppDists"  )
	value$strn
}
