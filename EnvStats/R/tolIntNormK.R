tolIntNormK <-
function (n, df = n - 1, coverage = 0.95, cov.type = "content", 
    ti.type = "two-sided", conf.level = 0.95, method = "exact", 
    rel.tol = 1e-07, abs.tol = rel.tol) 
{
    if (!is.numeric(n) || length(n) > 1 || n < 2) 
        stop("'n' must be a scalar larger than 1.")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    if (cov.type == "content") {
        if (!is.numeric(conf.level) || length(conf.level) > 1 || 
            conf.level <= 0 || conf.level >= 1) 
            stop("'conf.level' must be a scalar greater than 0 and less than 1.")
        alpha <- 1 - conf.level
        if (ti.type == "two-sided") {
            method <- match.arg(method, c("exact", "wald.wolfowitz"))
            r.fcn <- function(r, x, p) {
                ((pnorm(x + r) - pnorm(x - r)) - p)^2
            }
            r <- nlminb(start = -qnorm(0.5 - coverage/2), objective = r.fcn, 
                lower = 0, x = 1/sqrt(n), p = coverage)$par
            u <- sqrt(df/qchisq(alpha, df))
            K <- r * u
            if (method == "exact") {
                F.fcn.vec <- function(x, k, df, p) {
                  R <- nlminb(start = 0, objective = r.fcn, lower = 0, 
                    x = x, p = p)$par
                  pchisq(q = df * (R/k)^2, df = df, lower.tail = FALSE)
                }
                F.fcn <- Vectorize(F.fcn.vec)
                fcn.to.integrate <- function(x, k, df, p, n) {
                  F.fcn(x = x, k = k, df = df, p = p) * exp(-n * 
                    x^2/2)
                }
                fcn.to.minimize <- function(k, df, p, n, alpha) {
                  (sqrt(n/(2 * pi)) * integrate(fcn.to.integrate, 
                    lower = -Inf, upper = Inf, k = k, df = df, 
                    p = p, n = n, rel.tol = rel.tol, abs.tol = abs.tol)$value - 
                    1 + alpha)^2
                }
                K <- nlminb(start = K, objective = fcn.to.minimize, 
                  lower = 0, df = df, p = coverage, n = n, alpha = 1 - 
                    conf.level)$par
            }
        }
        else {
            K <- qT(p = conf.level, df = df, ncp = qnorm(coverage) * 
                sqrt(n))/sqrt(n)
        }
    }
    else {
        if (ti.type == "two-sided") {
            K <- sqrt(1 + 1/n) * qt(1 - (1 - coverage)/2, df)
        }
        else {
            K <- sqrt(1 + 1/n) * qt(coverage, df)
        }
    }
    K
}
