coarseLine <-
function (fun, x, dx, beta = 3, maximise = FALSE, ...) 
{
    sign = 1 - 2 * maximise
    if (missing(dx)) {
        dx = numeric(length(x))
        for (i in seq_along(x)) {
            x2 = x
            x2[i] = x[i] + 1e-08
            dx[i] = -sign * 1e+08 * (fun(x2, ...) - fun(x, ...))
        }
    }
    if (beta <= 1) 
        stop("adjustment factor beta must be larger than 1")
    adj = 1
    tol = .Machine$double.eps
    best = sign * fun(x, ...)
    nok = TRUE
    while (nok) {
        best = sign * fun(x + adj * dx, ...)
        if (!is.na(best)) 
            nok = FALSE
        else adj = adj/beta
        if (adj < tol) 
            stop("Failed to find valid move")
    }
    try = sign * fun(x + dx * adj * beta, ...)
    if (is.na(try) || try >= best) 
        down = TRUE
    else down = FALSE
    if (down) 
        fact = 1/beta
    else fact = beta
    nok = TRUE
    while (nok) {
        try = sign * fun(x + dx * adj * fact, ...)
        if (is.na(try) || try >= best) {
            nok = FALSE
            if (adj < tol || 1/adj < tol) 
                return(list(x = x + dx * adj, adj = adj, move = dx * 
                  adj, best = sign * best, code = 2))
        }
        else {
            adj = adj * fact
            best = try
            if (adj < tol || 1/adj < tol) 
                return(list(x = x + dx * adj, adj = adj/fact, 
                  move = dx * adj/fact, best = sign * best, code = 3))
        }
    }
    return(list(x = x + adj * dx, adj = adj, move = dx * adj, 
        best = sign * best, code = 1))
}
