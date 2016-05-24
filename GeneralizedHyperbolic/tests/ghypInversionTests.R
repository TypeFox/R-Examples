ghypInversionTestpq <- function(ps, intTol = .Machine$double.eps^0.25,
                                uniTol = .Machine$double.eps^0.25,
                                method = "spline", ...)
{
    qpps <- qghyp(pghyp(ps, intTol = intTol, ...), uniTol = uniTol,
                 method = method, ...)
    diffs <- qpps - ps
    result <- list(qpps = qpps, ps = ps, diffs = diffs)

    return(result)
}

ghypInversionTestqp <- function(qs, intTol = .Machine$double.eps^0.25,
                                uniTol = .Machine$double.eps^0.25,
                                method = "spline", ...)
{
    pqqs <- pghyp(qghyp(qs, uniTol = uniTol, method = method, ...),
                  intTol = intTol, ...)
    diffs <- pqqs - qs
    result <- list(pqqs = pqqs, qs = qs, diffs = diffs)

    return(result)
}
