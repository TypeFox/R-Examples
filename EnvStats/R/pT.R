pT <-
function (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, max.terms = 1e+05, 
    num.terms.inc = 20, rel.tol = .Machine$double.eps) 
{
    if (missing(ncp)) {
        p <- stats::pt(q = q, df = df, lower.tail = lower.tail, 
            log.p = log.p)
    }
    else if (all(ncp == 0 | is.na(ncp))) {
        p <- stats::pt(q = q, df = df, ncp = ncp, lower.tail = lower.tail, 
            log.p = log.p)
    }
    else {
        names.q <- names(q)
        arg.mat <- cbind.no.warn(q = as.vector(q), df = as.vector(df), 
            ncp = as.vector(ncp))
        na.index <- is.na.matrix(arg.mat)
        if (all(na.index)) 
            p <- rep(NA, nrow(arg.mat))
        else {
            n <- nrow(arg.mat)
            p <- numeric(n)
            p[na.index] <- NA
            p.no.na <- p[!na.index]
            for (i in c("q", "df", "ncp")) assign(i, arg.mat[!na.index, 
                i])
            if (any(df < .Machine$double.eps)) 
                stop("All non-missing values of 'df' must be positive.")
            ncp.eq.0 <- ncp == 0
            if (any(ncp.eq.0)) {
                p.no.na[ncp.eq.0] <- stats::pt(q = q[ncp.eq.0], 
                  df = df[ncp.eq.0], lower.tail = TRUE, log.p = FALSE)
            }
            if (!all(ncp.eq.0)) {
                finite.index <- is.finite(q) & is.finite(df) & 
                  is.finite(ncp)
                index <- !finite.index & !ncp.eq.0
                p.no.na[index] <- stats::pt(q = q[index], df = df[index], 
                  ncp = ncp[index], lower.tail = TRUE, log.p = FALSE)
                if (any(finite.index)) {
                  if (max.terms < 1 || trunc(max.terms) != max.terms) 
                    stop("'max.terms' must be a positive integer.")
                  if (num.terms.inc < 1 || num.terms.inc > max.terms || 
                    trunc(num.terms.inc) != num.terms.inc) 
                    stop(paste("'num.terms.inc' must be a positive integer,", 
                      "and must be less than or equal to 'max.terms'."))
                  if (rel.tol < .Machine$double.eps) 
                    stop(paste("'rel.tol' must be at least", 
                      .Machine$double.eps, "."))
                  warn.no.na <- logical(length(p.no.na))
                  index <- finite.index & ncp > 0
                  if (any(index)) {
                    temp.list <- pnct.pos.ncp(q = q[index], df = df[index], 
                      ncp = ncp[index], max.terms = max.terms, 
                      num.terms.inc = num.terms.inc, rel.tol = rel.tol)
                    p.no.na[index] <- temp.list$p
                    warn.no.na[index] <- temp.list$warn
                  }
                  index <- finite.index & ncp < 0
                  if (any(index)) {
                    temp.list <- pnct.pos.ncp(q = -q[index], 
                      df = df[index], ncp = -ncp[index], max.terms = max.terms, 
                      num.terms.inc = num.terms.inc, rel.tol = rel.tol)
                    p.no.na[index] <- 1 - temp.list$p
                    warn.no.na[index] <- temp.list$warn
                  }
                  if (any(warn.no.na)) {
                    if (n == 1) {
                      p.no.na <- NA
                      warning(paste("No convergence for the finite", 
                        "sum approximation to 'p'.  Change the values", 
                        "of 'max.terms' and/or 'rel.tol'."))
                    }
                    else {
                      p.no.na[warn.no.na] <- NA
                      warning(paste("No convergence for the finite ", 
                        "sum approximation to 'p' for the following ", 
                        "element index or indices:", "\n\n\t", 
                        paste((1:n)[!na.index][warn.no.na], collapse = " "), 
                        "\n\n  ", "Change the value(s) of 'max.terms' and/or ", 
                        "'rel.tol'.", sep = ""))
                    }
                  }
                }
            }
            p[!na.index] <- p.no.na
        }
        lower.tail <- as.logical(lower.tail)
        if (length(lower.tail) != 1) {
            lower.tail <- lower.tail[1]
            warning("Only first element of 'lower.tail' used")
        }
        if (!lower.tail) 
            p <- 1 - p
        log.p <- as.logical(log.p)
        if (length(log.p) != 1) {
            log.p <- log.p[1]
            warning("Only first element of 'log.p' used")
        }
        if (log.p) 
            p <- log(p)
        if (!is.null(names.q)) 
            names(p) <- rep(names.q, length = length(p))
    }
    p
}
