qprodnormalMC <-function(p, mu.x, mu.y, se.x, se.y, rho=0, lower.tail=TRUE, n.mc=1e7){
    if (!lower.tail)
        p <- 1-p
    mean.v <- c(mu.x,mu.y)
    var.mat <- matrix(c(se.x^2,se.x*se.y*rho,se.x*se.y*rho,se.y^2),2)
    a_b <- matrix(rnorm(2*n.mc),ncol=n.mc)
    a_b <- crossprod(chol(var.mat),a_b)+mean.v
    a_b <- t(a_b)
    ab <- a_b[,1]*a_b[,2]
    q <- quantile(ab,c(p))
    names(q) <- NULL
    x <- ab<q
    mean.x <- sum(x)/n.mc
    error.mean.x <- sd(x)/n.mc
    return(list(q=q, error=error.mean.x))
}
