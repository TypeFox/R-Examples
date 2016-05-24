pnct.pos.ncp.pos.q <-
function (i, q, df, ncp, max.terms = 1e+05, num.terms.inc = 20, 
    rel.tol = .Machine$double.eps) 
{
    if (any(q < .Machine$double.eps)) 
        stop("All values of 'q' must be positive.")
    q <- q[i]
    df <- df[i]
    ncp <- ncp[i]
    num.terms <- num.terms.inc
    p <- pnct.pos.ncp.pos.q.finite.sum(q, df, ncp, start = 0, 
        end = num.terms - 1)
    num.terms <- 2 * num.terms.inc
    while (num.terms <= max.terms) {
        p.inc <- pnct.pos.ncp.pos.q.finite.sum(q, df, ncp, start = num.terms - 
            num.terms.inc, end = num.terms - 1)
        if (p > 0 && (p.inc/p) < rel.tol) 
            break
        p <- p + p.inc
        num.terms <- num.terms + num.terms.inc
    }
    warn <- ifelse(num.terms > max.terms, 1, 0)
    c(p = p, warn = warn)
}
