## Calculate Kendall's tau and Spearman's rho for 2-D copula density estimate
summary.sscopu <- function(object,...)
{
    ## Check input
    if (class(object)!="sscopu") stop("gss error in summary.sscopu: not a sscopu object")
    if (dim(object$mdsty)[2]!=2) stop("gss error in summary.sscopu: not a 2-D copula")
    ## Set up quadrature
    hsz <- 40
    qdsz <- 2*hsz
    qd <- gauss.quad(qdsz,c(0,1))
    gap <- diff(qd$pt)
    g.wk <- gap[hsz]/2
    for (i in 1:(hsz-2)) g.wk <- c(g.wk,gap[hsz+i]-g.wk[i])
    g.wk <- 2*g.wk
    pp <- qd$pt[1]/(1/2-sum(g.wk))
    adj <- c(pp,rep(.5,qdsz-2),1-pp)
    qd.pt <- cbind(rep(qd$pt,qdsz),rep(qd$pt,rep(qdsz,qdsz)))
    ## Calculate cdf
    d.qd <- dsscopu(object,qd.pt)
    d.qd.wk <- matrix(d.qd,qdsz,qdsz)
    f.qd <- NULL
    for (i in 1:qdsz) {
        for (j in 1:qdsz) {
            wt1 <- qd$wt[1:i]
            wt1[i] <- wt1[i]*adj[i]
            wt2 <- qd$wt[1:j]
            wt2[j] <- wt2[j]*adj[j]
            f.qd <- c(f.qd,sum(d.qd.wk[1:i,1:j]*outer(wt1,wt2)))
        }
    }
    ## Calculate tau and rho
    tau <- 4*sum(f.qd*d.qd*outer(qd$wt,qd$wt))-1
    rho <- 12*sum(d.qd*outer(qd$pt*qd$wt,qd$pt*qd$wt))-3
    ## return
    list(tau=tau,rho=rho)
}
