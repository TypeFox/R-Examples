convest=
function (p, niter = 100, doplot = FALSE, doreport = FALSE) 
{
    if (!length(p)) 
        return(NA)
    if (any(is.na(p))) 
        stop("Missing values in p not allowed")
    if (any(p < 0 | p > 1)) 
        stop("All p-values must be between 0 and 1")
    k <- niter
    ny <- 1e-06
    ord.p=order(p)
    p <- sort(p)
    m <- length(p)
    p.c <- ceiling(100 * p)/100
    p.f <- floor(100 * p)/100
    t.grid <- (1:100)/100
    x.grid <- (0:100)/100
    t.grid.mat <- matrix(t.grid, ncol = 1)
    f.hat <- rep(1, 101)
    f.hat.p <- rep(1, m)
    theta.hat <- 0.01 * which.max(apply(t.grid.mat, 1, function(theta) sum((2 * 
        (theta - p) * (p < theta)/theta^2))))
    f.theta.hat <- 2 * (theta.hat - x.grid) * (x.grid < theta.hat)/theta.hat^2
    f.theta.hat.p <- 2 * (theta.hat - p) * (p < theta.hat)/theta.hat^2
    i <- 1
    j <- 0
    thetas <- numeric()
    for (j in 1:k) {
        if (sum((f.hat.p - f.theta.hat.p)/f.hat.p) > 0) 
            eps <- 0
        else {
            l <- 0
            u <- 1
            while (abs(u - l) > ny) {
                eps <- (l + u)/2
                if (sum(((f.hat.p - f.theta.hat.p)/((1 - eps) * 
                  f.hat.p + eps * f.theta.hat.p))[f.hat.p > 0]) < 
                  0) 
                  l <- eps
                else u <- eps
            }
        }
        f.hat <- (1 - eps) * f.hat + eps * f.theta.hat
        pi.0.hat <- f.hat[101]
        d <- -sum((f.theta.hat.p - f.hat.p)/f.hat.p)
        if (doreport == TRUE) {
            cat("j:", j, "\tpi0:", pi.0.hat, "\ttheta.hat:", 
                theta.hat, "\t\tepsilon:", eps, "\tD:", d, "\n")
        }
        f.hat.p <- 100 * (f.hat[100 * p.f + 1] - f.hat[100 * 
            p.c + 1]) * (p.c - p) + f.hat[100 * p.c + 1]
        theta.hat <- 0.01 * which.max(apply(t.grid.mat, 1, function(theta) sum((2 * 
            (theta - p) * (p < theta)/theta^2)/f.hat.p)))
        f.theta.hat <- 2 * (theta.hat - x.grid) * (x.grid < theta.hat)/theta.hat^2
        f.theta.hat.p <- 2 * (theta.hat - p) * (p < theta.hat)/theta.hat^2
        if (sum(f.theta.hat.p/f.hat.p) < sum(1/f.hat.p)) {
            theta.hat <- 0
            f.theta.hat <- rep(1, 101)
            f.theta.hat.p <- rep(1, m)
        }
        if (sum(thetas == theta.hat) == 0) {
            thetas[i] <- theta.hat
            thetas <- sort(thetas)
            i <- i + 1
        }
        pi.0.hat <- f.hat[101]
        if (doplot == TRUE) {
            plot(x.grid, f.hat, type = "l", main = paste(format(round(pi.0.hat, 
                5), digits = 5)), ylim = c(0, 1.2))
            points(thetas, f.hat[100 * thetas + 1], pch = 20, 
                col = "blue")
        }
    }
    ans=pi.0.hat
    fp=approxfun(x.grid, f.hat, rule=2)
    lfdr=ans/fp(p)
    lfdr[ord.p]=lfdr
    attr(ans, 'lfdr')=lfdr
    class(ans)='convest'
    return(ans)
}
