tpval <- function(x,mu=0) {
    xbar <- mean(x) 
    n <- length(x)
    s <- sd(x)
    t <- (xbar - mu) / (s / sqrt(n))
    return ( 2 * pt(-abs(t), df=n-1) )
}
ddd <- data.frame(
    pval = c(
        replicate( 2000, tpval(rt(5,df=3)) ),
        replicate( 2000, tpval(rt(10,df=3)) ),
        replicate( 2000, tpval(rt(20,df=3)) ),
        replicate( 2000, tpval(rt(30,df=3)) ),
        replicate( 2000, tpval(rexp(5,1), mu=1) ),
        replicate( 2000, tpval(rexp(10,1), mu=1) ),
        replicate( 2000, tpval(rexp(20,1), mu=1) ),
        replicate( 2000, tpval(rexp(30,1), mu=1) )
        ),
    dist = c(
        rep("T(3); n=05",2000), rep("T(3); n=10",2000),
        rep("T(3); n=20",2000), rep("T(3); n=30",2000),
        rep("Exp(1); n=05",2000), rep("Exp(1); n=10",2000),
        rep("Exp(1); n=20",2000), rep("Exp(1); n=30",2000)
        ))
myplot1 <- xqqmath(~pval|dist,ddd,dist=qunif,idline=TRUE,cex=0.4)
myplot2 <- xqqmath(~pval|dist,ddd,dist=qunif,idline=TRUE,cex=0.4,
    xlim=c(0,0.2), ylim=c(0,0.2))
