pnct.pos.ncp <-
function (q, df, ncp, max.terms = 1e+05, num.terms.inc = 20, 
    rel.tol = .Machine$double.eps) 
{
    n <- length(q)
    if (any(length.list(q, df, ncp) != n)) 
        stop("'q', 'df', and 'ncp' must all have the same length.")
    if (any(ncp < .Machine$double.eps)) 
        stop("All values of 'ncp' must be positive.")
    p <- numeric(n)
    warn <- logical(n)
    P0 <- pnorm(-ncp)
    index <- q == 0
    if (any(index)) 
        p[index] <- P0[index]
    index <- q > 0
    if (any(index)) {
        PF <- pf(q = q^2, df1 = 1, df2 = df, ncp = ncp^2)
        new.index <- index & ((P0 + PF) == 0)
        if (any(new.index)) {
            p[new.index] <- 0
        }
        new.index <- index & (PF == 1)
        if (any(new.index)) {
            p[new.index] <- 1
        }
        new.index <- index & ((P0 + PF) > 0) & (PF < 1)
        if (any(new.index)) {
            temp.mat <- sapply(1:sum(new.index), pnct.pos.ncp.pos.q, 
                q[new.index], df[new.index], ncp[new.index], 
                max.terms, num.terms.inc, rel.tol)
            p[new.index] <- P0[new.index] + temp.mat[1, ]
            warn[new.index] <- as.logical(temp.mat[2, ])
        }
    }
    index <- q < 0
    if (any(index)) {
        new.index <- index & P0 == 0
        if (any(new.index)) {
            p[new.index] <- 0
            new.index <- index & P0 > 0
            if (any(new.index)) {
                temp.mat <- sapply(1:sum(new.index), pnct.pos.ncp.pos.q, 
                  -q[new.index], df[new.index], ncp[new.index], 
                  max.terms, num.terms.inc, rel.tol)
                p[new.index] <- P0[new.index] + temp.mat[1, ] - 
                  pf(q[new.index]^2, 1, df[new.index], ncp[new.index]^2)
                warn[new.index] <- as.logical(temp.mat[2, ])
            }
        }
        else {
            temp.mat <- sapply(1:sum(index), pnct.pos.ncp.pos.q, 
                -q[index], df[index], ncp[index], max.terms, 
                num.terms.inc, rel.tol)
            p[index] <- P0[index] + temp.mat[1, ] - pf(q[index]^2, 
                1, df[index], ncp[index]^2)
            warn[index] <- as.logical(temp.mat[2, ])
        }
    }
    list(p = p, warn = warn)
}
