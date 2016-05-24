dT <-
function (x, df = stop("no 'df' argument"), ncp = 0, log = FALSE, 
    max.terms = 1e+05, num.terms.inc = 20, rel.tol = .Machine$double.eps) 
{
    if (missing(ncp)) {
        y <- stats::dt(x = x, df = df, log = log)
    }
    else if (all(ncp == 0 | is.na(ncp))) {
        y <- stats::dt(x = x, df = df, ncp = ncp, log = log)
    }
    else {
        names.x <- names(x)
        arg.mat <- cbind.no.warn(x = as.vector(x), df = as.vector(df), 
            ncp = as.vector(ncp))
        na.index <- is.na.matrix(arg.mat)
        if (all(na.index)) 
            y <- rep(NA, nrow(arg.mat))
        else {
            y <- numeric(nrow(arg.mat))
            y[na.index] <- NA
            y.no.na <- y[!na.index]
            for (i in c("x", "df", "ncp")) assign(i, arg.mat[!na.index, 
                i])
            if (any(df < .Machine$double.eps)) 
                stop("All non-missing values of 'df' must be positive.")
            ncp.eq.0 <- ncp == 0
            if (any(ncp.eq.0)) {
                y.no.na[ncp.eq.0] <- stats::dt(x = x[index], 
                  df = df[index], log = log)
            }
            if (!all(ncp.eq.0)) {
                finite.index <- is.finite(x) & is.finite(df) & 
                  is.finite(ncp)
                index <- !finite.index & !ncp.eq.0
                y.no.na[index] <- stats::dt(x = x[index], df = df[index], 
                  ncp = ncp[index], log = FALSE)
                index <- finite.index & !ncp.eq.0
                if (any(index)) {
                  if (max.terms < 1 || trunc(max.terms) != max.terms) 
                    stop("'max.terms' must be a positive integer.")
                  if (num.terms.inc < 1 || num.terms.inc > max.terms || 
                    trunc(num.terms.inc) != num.terms.inc) 
                    stop(paste("'num.terms.inc' must be a positive integer,", 
                      "and must be less than or equal to 'max.terms'."))
                  if (rel.tol < .Machine$double.eps) 
                    stop(paste("'rel.tol' must be at least", 
                      .Machine$double.eps, "."))
                  y.no.na[index] <- derivative(pT, x = x[index], 
                    df = df[index], ncp = ncp[index], lower.tail = TRUE, 
                    log.p = FALSE, max.terms = max.terms, num.terms.inc = num.terms.inc, 
                    rel.tol = rel.tol)
                }
            }
            if (log) 
                y.no.na <- log(y.no.na)
            y[!na.index] <- y.no.na
        }
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    y
}
