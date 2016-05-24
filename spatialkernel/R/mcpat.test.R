##Monte Carlo test of change of pattern over time (marks)
mcpat.test <- function(pts, marks, t, h, ntest=100, proc=TRUE)
{
    p2k <- function(p, r) { ## find k, sum_{j=1}^{k-1} p_j < r <= sum_{j=1}^k
        k <- 1
        pa <- p[1]
        while(r>pa) {
            k <- k+1
            pa <- pa+p[k]
        }
        k
    }
    
    tpfun <- function(p){
        sum(apply(p, 3, function(x) sum((x-apply(p, c(1,2), mean))^2)))
    }
        
    n <- length(marks)
    mtypes <- sort(unique(marks))
    m <- length(mtypes)
    tt <- sort(unique(t))
    ntt <- length(tt)
    ps <- array(NA, dim=c(n, m, ntt), dimnames = list(NULL, mtypes, NULL)) 
    tp <- NULL
    nh <- length(h)
    if(nh > 1) {
      lcp1 <- cvloglp(pts, marks, t, h)$cv
      hopt <- h[which.max(lcp1)]
    } else hopt <- h
    for(i in 1:ntest){
        if(proc) cat("\rProcessing No.", i, "out of", ntest)
        if(i==1) {
            y1 <- marks
        } else {
            runifn <- runif(n, min=0, max=1)
            for(j in 1:n) {
                y1[j] <- mtypes[p2k(p0[j,], runifn[j])]
            }
        }
        for(j in 1:ntt) {
            ndx <- which(t==tt[j])
            ps[,,j] <- phat(pts, pts[ndx,], y1[ndx], hopt)$p[, mtypes]
        }
        if(i==1) {
            p0 <- apply(ps, 1:2, mean) ## mean of p_j(x, t) over t
        }
        tp <- c(tp, tpfun(ps))
    }
    pv <- (ntest+1-rank(tp)[1])/ntest
    ##pv1 <- (ntest-rank(tp[2:ntest])[1])/(ntest-1)
    invisible(list(pvalue=pv, pts=pts, marks=marks, t=t, h=h, ntest=ntest))
}
