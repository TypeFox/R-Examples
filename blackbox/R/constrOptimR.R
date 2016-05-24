## version of constrOptim that tests within hullness en rational
constrOptimR <- function (theta, f, grad, ui, ci, mu = 1e-04, control = list(),
    method = if (is.null(grad)) "Nelder-Mead" else "BFGS", outer.iterations = 100,
    outer.eps = 1e-05, ...) {
    if (!is.null(control$fnscale) && control$fnscale < 0)
        mu <- -mu
    R <- function(theta, q.theta.old, ...) {
        q.theta <- d2q(theta)
        q.ui.theta <- qmatmult(ui, matrix(q.theta))
        gi <- q2d(qmq(q.ui.theta, ci))
        if (any(gi < 0))
            return(NaN)
        gi.old <- q2d(qmq(qmatmult(ui, matrix(q.theta.old)),
            ci))
        bar <- sum(gi.old * log(gi) - q2d(q.ui.theta))
        if (!is.finite(bar))
            bar <- -Inf
        f(theta, ...) - mu * bar
    }
    dR <- function(theta, q.theta.old, ...) {
        q.theta <- d2q(theta)
        q.ui.theta <- qmatmult(ui, matrix(q.theta))
        gi <- drop(qmq(q.ui.theta, ci))
        gi.old <- drop(qmq(qmatmult(ui, matrix(q.theta.old)),
            ci))
        pif <- qmatmult(t(matrix(qdq(gi.old, gi))), ui)
        paf <- qmatmult(t(matrix(rep("1", nrow(ui)))), ui)
        dbar <- qmq(pif, paf)
        grad(theta, ...) - mu * q2d(dbar)
    }
    q.theta <- d2q(theta)
    if (any(q2d(qmq(qmatmult(ui, matrix(q.theta)), ci)) <= 0))
        stop.redef("initial value not feasible")
    obj <- f(theta, ...)
    r <- R(theta, q.theta, ...)
    s.mu <- sign(mu)
    for (i in 1L:outer.iterations) {
        obj.old <- obj
        r.old <- r
        theta.old <- theta
        q.theta.old <- q.theta
        fun <- function(theta, ...) {
            R(theta, q.theta.old, ...)
        }
        gradient <- function(theta, ...) {
            dR(theta, q.theta.old, ...)
        }
        a <- optim(theta.old, fun, gradient, control = control,
            method = method, ...)
        r <- a$value
        if (is.finite(r) && is.finite(r.old) && abs(r - r.old) <
            (0.001 + abs(r)) * outer.eps)
            break
        theta <- a$par
        q.theta <- d2q(theta)
        obj <- f(theta, ...)
        if (s.mu * obj > s.mu * obj.old)
            break
    }
    if (i == outer.iterations) {
        a$convergence <- 7
        a$message <- "Barrier algorithm ran out of iterations and did not converge"
    }
    if (mu > 0 && obj > obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function increased at outer iteration",
            i)
    }
    if (mu < 0 && obj < obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function decreased at outer iteration",
            i)
    }
    a$outer.iterations <- i
    a$barrier.value <- a$value
    a$value <- f(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a
}
