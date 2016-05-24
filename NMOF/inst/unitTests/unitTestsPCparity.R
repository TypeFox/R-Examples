## -*- truncate-lines: t; -*-
test.putCallParity <- function() {
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; vol <- 0.3
    call <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    put <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEqualsNumeric(putCallParity("call", call, put, S, X, tau, r, q), call)
    checkEqualsNumeric(putCallParity("put", call, put, S, X, tau, r, q), put)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    call <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    put <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEqualsNumeric(putCallParity("call", call, put, S, X, tau, r, q), call)
    checkEqualsNumeric(putCallParity("put", call, put, S, X, tau, r, q), put)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.04; vol <- 0.3
    call <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    put <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEqualsNumeric(putCallParity("call", call, put, S, X, tau, r, q), call)
    checkEqualsNumeric(putCallParity("put", call, put, S, X, tau, r, q), put)

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.04; vol <- 0.3
    call <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    put <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEqualsNumeric(putCallParity("call", call, put, S, X, tau, r, q), call)
    checkEqualsNumeric(putCallParity("put", call, put, S, X, tau, r, q), put)

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.0; vol <- 0.3
    D <- 20; tauD <- 0.5 ## dividends
    call <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, tauD = tauD, D = D, type = "call")$value
    put <- vanillaOptionEuropean(S, X, tau, r, q, vol^2,  tauD = tauD, D = D, type = "put")$value
    checkEqualsNumeric(putCallParity("call", call, put, S, X, tau, r, q, tauD, D), call)
    checkEqualsNumeric(putCallParity("put", call, put, S, X, tau, r, q, tauD, D), put)
}
