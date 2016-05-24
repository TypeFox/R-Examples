piecewiseExp.ph <-
function (coxObject, knots = NULL, length.knots = 6) {
    Time <- coxObject$y[, 1]
    d <- coxObject$y[, 2]
    n <- length(Time)
    if (is.null(knots)) {
        Q <- length.knots + 1
        knots <- unique(quantile(Time, seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
        knots <- knots + 1e-06
        if (max(knots) > max(Time))
            knots[which.max(knots)] <- max(Time) - 1e-06
    }
    knots <- c(0, sort(knots), max(Time) + 1)
    Q <- length(knots) - 1
    ind <- findInterval(Time, knots, rightmost.closed = TRUE)
    D <- matrix(0, n, Q)
    D[cbind(seq_along(ind), ind)] <- 1
    D <- c(D * d)
    Tiq <- outer(Time, knots, pmin)
    T <- c(Tiq[, 2:(Q+1)] - Tiq[, 1:Q])
    X <- coxObject$x[rep(seq_len(n), Q), ]
    ND <- suppressWarnings(data.frame(Time = T, D = D, X, xi = gl(Q, n),
        check.names = FALSE)[T > 0, ])
    glm(D ~ . + offset(log(Time)) - Time - 1,
        family = poisson, data = ND)
}
