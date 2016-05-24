"davies.moment" <-
function (n = 1, i = 1, order = 1, params) 
{
    C <- params[1]
    l1 <- params[2]
    l2 <- -params[3]  # sign change
    determinant <- 1 + n + order * (l1 + l2)
    if (is.na(determinant)) {
        return(determinant)
    }
    if (determinant > 0) {
        return(exp(order * log(C) + lgamma(i + order * l1) - 
            lgamma(i) + lgamma(n + 1) - lgamma(n + 1 - i) + lgamma(n + 
            1 - i + order * l2) - lgamma(n + 1 + order * l1 + 
            order * l2)))
    }
    else {
        return(NaN)
    }
}

"davies.start" <-
function (x, threeps = c(0.1, 0.5, 0.9), small = 0.01) 
{
    x <- sort(x)
    force.in.01 <- function(i) {
        min(max(i, 0), 1)
    }
    plow <- force.in.01(threeps[1])
    pmed <- force.in.01(threeps[2])
    phigh <- force.in.01(threeps[3])
    number.of.methods <- 3
    params <- matrix(NA, number.of.methods, 3)
    med <- abs(quantile(x, pmed))
    left <- abs(quantile(x, plow))
    right <- abs(quantile(x, phigh))
    C <- med
    l1 <- log(left/med)/log(plow)
    l2 <- log(right/med)/log(plow)
    params[1, ] <- as.vector(c(C, l1, l2))
    l1 <- log(med/right)/log(pmed/phigh)
    C <- right/(phigh^l1)
    l2 <- -small
    params[2, ] <- as.vector(c(C, l1, l2))
    l1 <- small
    l2 <- log(med/right)/log(pmed/plow)
    C <- right/(plow^l2)
    params[3, ] <- as.vector(c(C, l1, l2))

    params[,3] <- -params[,3] #sign change

    
    f <- function(n) {
        objective(params[n, ], x)
    }
    obj <- sapply(1:number.of.methods, f)
    return(params[which.min(obj), ])
}

"ddavies" <-
function (x, params) 
{
    ddavies.p(pdavies(x, params), params)
}

"ddavies.p" <-
function (x, params) 
{
    C <- params[1]
    l1 <- params[2]
    l2 <- -params[3]  #sign change
    1/(C * (l1 * x^(l1 - 1) * (1 - x)^l2 - l2 * x^l1 * (1 - x)^(l2 - 
        1)))
}
"dgld" <-
function (x, params) 
{
    l1 <- params[1]
    l2 <- params[2]
    l3 <- params[3]
    l4 <- params[4]
    dgld.p(pgld(x, params), params)
}
"dgld.p" <-
function (x, params) 
{
    l1 <- params[1]
    l2 <- params[2]
    l3 <- params[3]
    l4 <- params[4]
    1/((l4 * (1 - x)^(l4 - 1) + l3 * x^(l3 - 1))/l2)
}

"expected.gld" <-
function (n = 1, i = 1, params) 
{
    l1 <- params[1]
    l2 <- params[2]
    l3 <- params[3]
    l4 <- params[4]
    return(l1 + davies.moment(n = n, i = i, params = c(1/l2, 
        l3, 0)) - davies.moment(n = n, i = i, params = c(1/l2, 
        0, l4)))
}

"expected.gld.approx" <-
function (n = 1, i = 1, params) 
{
    qgld(i/(n + 1), params)
}

"expected.value" <-
function (n, i, params) 
{
    return(davies.moment(n = n, i = i, order = 1, params))
}

"expected.value.approx" <-
function (n, i, params) 
{
    C <- params[1]
    l1 <- params[2]
    l2 <- -params[3]  #sign convention
    p.i <- i/(n + 1)
    return(C * (p.i)^l1 * (1 - p.i)^l2)
}

"fit.davies.p" <-
function (x, print.fit = FALSE, use.q = TRUE, params = NULL, 
    small = 1e-05, ...) 
{
    x <- sort(x)
    if (is.null(params)) {
        if (use.q == TRUE) {
            jj <- least.squares(x, do.print = TRUE)
        }
        else {
            jj <- maximum.likelihood(x, do.print = TRUE)
        }
    }
    else {
        jj <- list(parameters = params, error = objective(params, 
            dataset = x))
    }
    params <- jj$parameters
    if (print.fit) {
        print(jj)
        print(params)
    }
    p <- seq(from = small, to = 1 - small, length = 1000)
    x.fit <- qdavies(p, params)
    y.fit <- ddavies.p(p, params)
    plot(x.fit, y.fit, ylim = c(0, max(y.fit)), type = "l", ...)
    points(x, x * 0, pch = 3)
}

"fit.davies.q" <-
function (x, print.fit = FALSE, use.q = TRUE, params = NULL, ...) 
{
    x <- sort(x)
    if (is.null(params)) {
        if (use.q == TRUE) {
            jj <- least.squares(x, do.print = TRUE)
        }
        else {
            jj <- maximum.likelihood(x, do.print = TRUE)
        }
        params <- jj$parameters
    }
    else {
        jj <- list(parameters = params, error = objective(params, 
            dataset = x))
    }
    if (print.fit) {
        print(jj)
    }
    n <- length(x)
    plot(sort(x), pch = 16, ...)
    points(1:n, expected.value(n, 1:n, params), type = "l")
}

"kurtosis" <-
function (params) 
{
    fourth.moment <- M(4, params) - 4 * M(3, params) * M(1, params) + 
        6 * M(2, params) * M(1, params)^2 - 3 * M(1, params)^4
    return(fourth.moment/variance(params)^2)
}

"least.squares" <-
function (data, do.print = FALSE, start.v = NULL) 
{
    data <- sort(data)
    if (is.null(start.v)) {
        start.v <- davies.start(data)
    }
    jj <- optim(start.v, objective, dataset = data)
    if (do.print != TRUE) {
        return(jj$par)
    }
    else {
        return(list(parameters = jj$par, error = jj$value))
    }
}

"likelihood" <-
function (params, data) 
{
    prod(ddavies(data, params))
}

"M" <-
function (order, params) 
{
    davies.moment(n = 1, i = 1, order = order, params = params)
}

"maximum.likelihood" <-
function (data, do.print = FALSE, start.v = NULL) 
{
    data <- sort(data)
    if (is.null(start.v)) {
        start.v <- davies.start(data)
    }
    f <- function(params, data) {
        neg.log.likelihood(params, data)
    }
    jj <- optim(start.v, f, data = data)
    if (do.print != TRUE) {
        return(jj$par)
    }
    else {
        return(list(parameters = jj$par, error = jj$value))
    }
}

"mu" <-
function (params) 
{
    M(1, params)
}

"neg.log.likelihood" <-
function (params, data) 
{
    -sum(log(ddavies(data, params)))
}

"objective" <-
function (params, dataset) 
{
    dataset <- sort(dataset)
    n <- length(dataset)
    error <- sum((dataset - expected.value(n, 1:n, params))^2)
    penalty <- function(x) {
        abs(x) * (x < 0) * 1000
    }
    return(error + penalty(params[1]) + penalty(params[2]) + 
        penalty(-params[3]))
}

"objective.approx" <-
function (params, dataset) 
{
    dataset <- sort(dataset)
    n <- length(dataset)
    sum((dataset - expected.value.approx(n, 1:n, params))^2)
}

"pdavies" <-
function (x, params) 
{
    if (any(is.na(params))) {
        return(params[is.na(params)])
    }
    f <- function(p, a) {
        qdavies(p = p, params = a[1:3]) - a[4]
    }
    options(warn = -1)
    if (length(x) <= 1) {
        if (length(x) == 0) {
            return(x)
        }
        if (is.na(x)) {
            return(x)
        }
        if (x <= 0) {
            return(0)
        }
        return(uniroot(f, c(0, 1), tol = 1e-06, a = c(params, 
            x))$root)
    }
    else {
        val <- sapply(x, pdavies, params)
        attributes(val) <- attributes(x)
        return(val)
    }
}

"pgld" <-
function (q, params) 
{
    f <- function(p, a) {
        qgld(p = p, params = a[1:4]) - a[5]
    }
    options(warn = -1)
    if (length(q) == 1) {
        return(uniroot(f, c(0, 1), tol = 1e-06, a = c(params, 
            q))$root)
    }
    else {
        val <- sapply(q, pgld, params)
        attributes(val) <- attributes(q)
        return(val)
    }
}

"plotcf" <-
function (y, q = 0.05) 
{
    y <- sort(y)
    x <- seq(along = y)
    is.lt <- y < q
    plot(x, y, type = "n")
    points(x[is.lt], y[is.lt], pch = 16)
    lines(x[!is.lt], y[!is.lt])
}

"qdavies" <-
function (p, params) 
{
    C <- params[1]
    l1 <- params[2]
    l2 <- -params[3]
    ans <- C * p^l1 * (1 - p)^l2
    ans[p <= 0] <- 0
    ans[p >= 1] <- Inf
    ans
}

"qgld" <-
function (p, params) 
{
    out <- p
    l1 <- params[1]
    l2 <- params[2]
    l3 <- params[3]
    l4 <- params[4]

    out[p <= 0] <- -Inf
    out[p >= 1] <- Inf
    
    l1 + (p^l3 - (1 - p)^l4)/l2
}

"rdavies" <-
function (n = 1, params) 
{
    qdavies(runif(n), params)
}

"rgld" <-
function (n = 1, params) 
{
    l1 <- params[1]
    l2 <- params[2]
    l3 <- params[3]
    l4 <- params[4]
    return(qgld(runif(n), params))
}

"rstupid" <-
function (n, a = 1, b = 2, c = 3, d = 4) 
{
    if (!((a < b) & (b < c) & (c < d))) {
        stop("should have a<b<c<d")
    }
    jj <- rbinom(n, size = 1, prob = (b - a)/((b - a) + (d - c)))
    return(jj * runif(n, a, b) + (1 - jj) * runif(n, c, d))
}

"skewness" <-
function (params) 
{
    third.moment <- M(3, params) - 3 * M(1, params) * M(2, params) + 
        2 * M(1, params)^3
    return(third.moment/(variance(params)^(3/2)))
}

"twolines.vert" <-
function (p, y1, y2, ...) 
{
    n <- length(p)
    plot(p, pmin(y1, y2), type = "n")
    plot(p, pmax(y1, y2), type = "n")
    plot(p, y1, type = "l")
    lines(p, y2, type = "l")
    segments(p, y1, p, y2, ...)
}

"variance" <-
function (params) 
{
    M(2, params) - M(1, params)^2
}
