armijo <-
function (fun, x, dx, beta = 3, sigma = 0.5, grad, maximise = FALSE, 
    searchup = TRUE, adj.start = 1, ...) 
{
    sign = 1 - 2 * maximise
    if (beta <= 1) 
        stop("adjustment factor beta must be larger than 1")
    adj = adj.start
    tol = .Machine$double.eps
    cur = sign * fun(x, ...)
    nok = TRUE
    neg = TRUE
    if (missing(grad)) {
        grad = numeric(length(x))
        for (i in seq_along(x)) {
            x2 = x
            x2[i] = x[i] + 1e-08
            grad[i] = 1e+08 * (sign * fun(x2, ...) - cur)
        }
    }
    else grad = sign * grad
    if (missing(dx)) 
        dx = -grad
    s = -sigma * sum(dx * grad)
    moddx = sqrt(sum(dx^2))
    if (s < 0) {
        return(list(x = x, adj = 0, move = dx * 0, best = sign * 
            cur, code = 2))
    }
    while (nok) {
        try = sign * fun(x + dx * adj, ...)
        if (!is.na(try) && (cur - try) > adj * s) {
            nok = FALSE
        }
        else {
            adj = adj/beta
            searchup = FALSE
        }
        if (neg && !is.na(try) && (cur - try) > -adj * s) 
            neg = FALSE
        if (adj * moddx < tol) 
            return(list(x = x, adj = adj, move = dx * adj, best = sign * 
                cur, code = 2 * neg))
    }
    if (searchup) {
        ok = TRUE
        while (ok) {
            adj = adj * beta
            try2 = sign * fun(x + dx * adj, ...)
            if (is.na(try2)) {
                adj = adj/beta
                break
            }
            if (try - try2 < .Machine$double.eps) {
                adj = adj/beta
                break
            }
            else {
                try = try2
            }
            if (adj > 1/tol) {
                if (sqrt(sum(adj * dx^2)) < 10 * tol) 
                  return(list(x = x + dx * adj, adj = adj, move = dx * 
                    adj, best = sign * try, code = 3))
                else break
            }
        }
    }
    return(list(x = x + dx * adj, adj = adj, move = dx * adj, 
        best = sign * try, code = 1))
}
