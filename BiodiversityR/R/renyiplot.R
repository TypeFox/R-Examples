`renyiplot` <- function(xr, addit=F, pch=1,
    xlab="alpha", ylab="H-alpha", ylim=c(0,m),
    labelit=T, legend=T, col=1, cex=1, rainbow=T, evenness=F,...) 
{
    x <- xr
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    m <- max(x,na.rm=T)
    names <- colnames(x)
    names <- as.factor(names)
    pos <- -1000
    ylab <- "H-alpha"
    if(evenness==T) {
        pos <- 1000
        if (ylab== "H-alpha") {ylab <- "E-alpha"}
        x[,] <- x[,]-x[,1]
        m <- min(x,na.rm=T)
        ylim <- c(m,0)
    }
    if(addit==F) {
        graphics::plot(names, rep(pos,p), xlab=xlab, ylab=ylab, ylim=ylim, bty="l",...)
    }
    if (n > 25) {
        warning("Symbol size was kept constant as there were more than 25 profiles (> number of symbols that are currently used in R)")
        rainbow <- T
    }
    if (rainbow==T && n > 1) {
        grDevices::palette(rainbow(n))
        for (i in 1:n) {
            if (n<26) {graphics::points(c(1:p), x[i,], pch=i, col=i, cex=cex, type="o")}
            if (n>25) {graphics::points(c(1:p), x[i,], pch=19, col=i, cex=cex, type="o")}
            if (labelit==T) {
                graphics::text(1, x[i,1], labels=rownames(x)[i], pos=2, col=i, cex=cex)
                graphics::text(p, x[i,p], labels=rownames(x)[i], pos=4, col=i, cex=cex)
            }
        }
        if (legend==T && n<26) {legend(graphics::locator(1), legend=rownames(x), pch=c(1:n), col=c(1:n))}
        if (legend==T && n>25) {legend(graphics::locator(1), legend=rownames(x), pch=rep(19,n), col=c(1:n))}
    }else{
        for (i in 1:n) {
            graphics::points(c(1:p), x[i,], pch=pch, col=col, cex=cex, type="o")
            if (labelit==T) {
                graphics::text(1, x[i,1], labels=rownames(x)[i], pos=2, col=col, cex=cex)
                graphics::text(p, x[i,p], labels=rownames(x)[i],pos=4, col=col, cex=cex)
            }
        }
        if (legend==T) {legend(graphics::locator(1), legend=rownames(x), pch=c(1:n))}
    }
    grDevices::palette("default")
}






