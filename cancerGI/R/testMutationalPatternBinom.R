testMutationalPatternBinom <- function (x, y, PRINT=FALSE) {
    xy.tmp <- cbind (x, y)
    xy.tmp.nona <- xy.tmp[which (!is.na (x) & !is.na (y)),]
    
    x <- xy.tmp.nona[,1]
    y <- xy.tmp.nona[,2]
    
    p.x <- c(sum (x==0), sum(x==1), sum(x==2)) / length (x)
    p.y <- c(sum (y==0), sum(y==1), sum(y==2)) / length (y)
    
    xy <- c(sum (x==0 & y==0), sum (x==1 & y==0), sum(x==2 & y==0), sum (x==0 & y==1), sum (x==1 & y==1), sum(x==2 & y==1), sum (x==0 & y==2), sum (x==1 & y==2), sum(x==2 & y==2))
    
    if (PRINT) {
        print (matrix (xy, ncol=3))
        print (p.x %*% t(p.y) * length (x))
    }
    xy.LL.binom <- binom.test (x=xy[1], n=length (x), p=p.x[1]*p.y[1], alternative="greater")
    xy.GL.binom <- binom.test (x=xy[3], n=length (x), p=p.x[3]*p.y[1], alternative="greater")
    xy.LG.binom <- binom.test (x=xy[7], n=length (x), p=p.x[1]*p.y[3], alternative="greater")
    xy.GG.binom <- binom.test (x=xy[9], n=length (x), p=p.x[3]*p.y[3], alternative="greater")
    
    xy.ME.binom <- binom.test (x=sum (xy[c(2,4,6,8)], na.rm=TRUE), n=length (x), p=(p.x[1]+p.x[3])*(1-p.y[1]-p.y[3]) + (p.y[1]+p.y[3])*(1-p.x[1]-p.x[3]), alternative="greater")
    
    return (c(xy.LL.binom$p.value, xy.GL.binom$p.value, xy.LG.binom$p.value, xy.GG.binom$p.value, xy.ME.binom$p.value))
}

