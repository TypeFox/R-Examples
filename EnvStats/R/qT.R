qT <-
function (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, max.terms = 1e+05, 
    num.terms.inc = 20, rel.tol = .Machine$double.eps) 
{
    if (missing(ncp)) {
        q <- stats::qt(p = p, df = df, lower.tail = lower.tail, 
            log.p = log.p)
    }
    else if (all(ncp == 0 | is.na(ncp))) {
        q <- stats::qt(p = p, df = df, ncp = ncp, lower.tail = lower.tail, 
            log.p = log.p)
    }
    else {
        names.p <- names(p)
        arg.mat <- cbind.no.warn(p = as.vector(p), df = as.vector(df), 
            ncp = as.vector(ncp))
        na.index <- is.na.matrix(arg.mat)
        if (all(na.index)) 
            q <- rep(NA, nrow(arg.mat))
        else {
            q <- numeric(nrow(arg.mat))
            q[na.index] <- NA
            q.no.na <- q[!na.index]
            for (i in c("p", "df", "ncp")) assign(i, arg.mat[!na.index, 
                i])
            if (log.p) {
                if (any(p > 0)) 
                  stop("When log.p=TRUE all values of p must be less than or equal to 0")
                p <- exp(p)
            }
            else if (any(p < 0) || any(p > 1)) 
                stop("All non-missing values of 'p' must be between 0 and 1.")
            if (!lower.tail) 
                p <- 1 - p
            if (any(df < .Machine$double.eps)) 
                stop("All non-missing values of 'df' must be positive.")
            index <- ncp == 0
            if (any(index)) {
                q.no.na[index] <- stats::qt(p = p[index], df = df[index], 
                  lower.tail = TRUE, log.p = FALSE)
            }
            index <- ncp != 0
            if (any(index)) {
                q.no.na[index & p == 0] <- -Inf
                q.no.na[index & p == 1] <- Inf
                index <- (1:length(q.no.na))[index & (0 < p & 
                  p < 1)]
                if (any(index)) {
                  func.to.min <- function(q.weird, p.weird, df.weird, 
                    ncp.weird, max.terms.weird, num.terms.inc.weird, 
                    rel.tol.weird) {
                    (pT(q = q.weird, df = df.weird, ncp = ncp.weird, 
                      lower.tail = TRUE, log.p = FALSE, max.terms = max.terms.weird, 
                      num.terms.inc = num.terms.inc.weird, rel.tol = rel.tol.weird) - 
                      p.weird)^2
                  }
                  o.warn <- options("warn")
                  for (i in index) {
                    options(warn = -1)
                    start <- stats::qt(p = p[i], df = df[i], 
                      ncp = ncp[i], lower.tail = TRUE, log.p = FALSE)
                    options(o.warn)
                    q.no.na[i] <- nlminb(start = start, objective = func.to.min, 
                      p.weird = p[i], df.weird = df[i], ncp.weird = ncp[i], 
                      max.terms.weird = max.terms, num.terms.inc.weird = num.terms.inc, 
                      rel.tol.weird = rel.tol)$par
                  }
                }
            }
            q[!na.index] <- q.no.na
        }
        if (!is.null(names.p)) 
            names(q) <- rep(names.p, length = length(q))
    }
    q
}
