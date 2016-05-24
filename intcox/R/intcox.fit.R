intcox.fit <-
function (formula, data, covar, eps, itermax)                       # Iterated Convex Minorant Algorithm
{
    folge <- sort.list(c(data$left, data$right[data$cens == 3]))    # sorting the pooled ends
    rang <- rank(c(data$left, data$right[data$cens == 3]))          # ranks of the pooled ends
    lamb.beta <- intcox.breslow(formula, data, covar)
    lambda0v <- intcox.hazard0(data$right[data$cens == 3], data,
        lamb.beta)                                                  # first estimation of the baseline-hazard for the interval ends
    beta <- lamb.beta$fit                                           # first estimation for beta
    lambda0u <- intcox.hazard0.beg(data$left, data, lamb.beta)      # derived baseline-hazard for the interval beginnings
    lambda0g <- c(lambda0u, lambda0v)
    likeli.vec <- NULL
    null.vector <- rep(0, length(lambda0g))
    itmax <- itermax
    it <- 0                                                         # iteration counter
    abbruch <- 0                                                    # break condition not fulfilled
    while (abbruch == 0) {                                          # until the break condition is fulfilled
        it <- it + 1
        if (it > itmax)
            abbruch <- 3
        derivs.wert <- intcox.derivs(data, covar, lambda0u, lambda0v,
            beta)                                                   # first and second derivatives and likelihood
        if (any(derivs.wert$g1 <= 0)) {
            abbruch <- 4
        }
        if (any(!is.finite(derivs.wert$l1))) {
            abbruch <- 4
        }
        if (any(derivs.wert$g2 <= 0)) {
            abbruch <- 4
        }
        if (any(!is.finite(derivs.wert$l2))) {
            abbruch <- 4
        }
        likeli.vec <- c(likeli.vec, derivs.wert$likeli)
        g1.inv <- 1/(derivs.wert$g1)
        g2.inv <- 1/(derivs.wert$g2)
        lambda0g <- c(lambda0u, lambda0v)
        alpha <- 1                                                  # step size
        if (any(!is.finite(lambda0g + alpha * g1.inv * derivs.wert$l1))) {
            abbruch <- 4
        }
        ind <- rep(1, length(lambda0g))
        if (abbruch == 0) {
            repeat {                                                # step size trials
                lambda0g.new <- intcox.pavaC(derivs.wert$g1[folge],
                  pmax((lambda0g + alpha * g1.inv * derivs.wert$l1)[folge],
                    null.vector))[rang]                             # PAVA with C
                beta.neu <- beta + alpha * g2.inv * derivs.wert$l2
                lambda0u.neu <- lambda0g.new[1:length(data$left)]
                lambda0v.neu <- lambda0g.new[-(1:length(data$left))]
                likeli.neu <- intcox.derivs(data, covar, lambda0u.neu,
                  lambda0v.neu, beta.neu)$likeli
                likeli.diff <- likeli.neu - derivs.wert$likeli
                if (is.na(likeli.diff)) {
                  abbruch <- 4
                  break
                }
                if (likeli.diff > 0) {
                  break
                }
                else {
                  alpha <- alpha * 1/2                              # try the half step size
                }
                if (alpha < 0.5^40) {
                  abbruch <- 2
                  break
                }
            }
        }
        if (abbruch == 0) {
            alpha.u <- alpha * 1/2
            lambda0g.u <- intcox.pavaC(derivs.wert$g1[folge], pmax((lambda0g +
                alpha.u * g1.inv * derivs.wert$l1)[folge], null.vector))[rang]  # PAVA with C
            beta.u <- beta + alpha.u * g2.inv * derivs.wert$l2
            lambda0u.u <- lambda0g.u[1:length(data$left)]
            lambda0v.u <- lambda0g.u[-(1:length(data$left))]
            likeli.u <- intcox.derivs(data, covar, lambda0u.u, lambda0v.u,
                beta.u)$likeli
            likeli.diff.u <- likeli.u - derivs.wert$likeli
        }
        if (abbruch == 0) {
            alpha.o <- alpha * 3/2
            lambda0g.o <- intcox.pavaC(derivs.wert$g1[folge], pmax((lambda0g +
                alpha.o * g1.inv * derivs.wert$l1)[folge], null.vector))[rang]  # PAVA with C
            beta.o <- beta + alpha.o * g2.inv * derivs.wert$l2
            lambda0u.o <- lambda0g.o[1:length(data$left)]
            lambda0v.o <- lambda0g.o[-(1:length(data$left))]
            likeli.o <- intcox.derivs(data, covar, lambda0u.o, lambda0v.o,
                beta.o)$likeli
            likeli.diff.o <- likeli.o - derivs.wert$likeli
        }
        zweig <- 0
        if (abbruch == 0) {
            if (likeli.neu >= likeli.o && likeli.neu >= likeli.u) {
                lambda0u <- lambda0u.neu
                lambda0v <- lambda0v.neu
                beta <- beta.neu
                zweig <- 2
            }
            if (likeli.u >= likeli.o && likeli.u >= likeli.neu) {
                lambda0u <- lambda0u.u
                lambda0v <- lambda0v.u
                likeli.diff <- likeli.diff.u
                beta <- beta.u
                zweig <- 1
            }
            if (likeli.o >= likeli.u && likeli.o >= likeli.neu) {
                lambda0u <- lambda0u.o
                lambda0v <- lambda0v.o
                likeli.diff <- likeli.diff.o
                beta <- beta.o
                zweig <- 3
            }
            likeli.rel <- abs(likeli.diff/derivs.wert$likeli)
            if (is.na(likeli.rel)) {
                abbruch <- 4
            }
            else {
                if (likeli.rel < eps)
                  abbruch <- 1                                  # if no break
            }
        }
    }
    if (abbruch == 2)
        likeli.vec <- likeli.vec
    else likeli.vec <- c(likeli.vec, likeli.neu)
    time.point <- c(data$left, data$right[data$cens == 3])
    if (abbruch == 2)
        lambda0 <- lambda0g
    else lambda0 <- lambda0g.new
    ordnung <- sort.list(time.point)
    time.point <- time.point[ordnung]
    lambda0 <- lambda0[ordnung]
    time.point <- time.point[!duplicated(lambda0)]
    lambda0 <- lambda0[!duplicated(lambda0)]
    intcox.ret <- list(coefficients = beta, lambda0 = lambda0, time.point = time.point,
        likeli.vec = likeli.vec, iter = it, termination = abbruch)
    return(intcox.ret)
}
