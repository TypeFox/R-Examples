plot.BiogTest <-
function (x, max.poly = 50,...) {
    if(class(x)[1] != "BiogTest")
        stop("x is not a 'BiogTest' class object")

    m <- matrix(c(1,2,1,3),2,2)
    layout(m, widths=c(1, 1), heights=c(1.2, 2))

    par(mar=c(4,4,4,2))
    m <- min(x$Zsim)
    M <- max(x$Zsim)
    dif <- M-m
    if((x$Z < M*1.1) && (x$Z > m - M*0.1)) {
        h1 <- hist(x$Zsim, xlab=expression(paste("Distribution of simulated Z"[sim], "values")), breaks=25, xlim=c(m - M*0.1, M*1.1), col="grey90", border="grey", main="")
        title("Biogeographical null model test", line=1.5) 
        title(paste("P(Obs Z <= null) =", x$pval[2]), line=0.5, font.main=1, cex.main=1)
        abline(v=x$Z, col="red", lty=3, lwd=3)
        abline(v=mean(x$Zsim), lwd=2)
        abline(v=quantile(x$Zsim, 0.975), lwd=2, lty=3)
        abline(v=quantile(x$Zsim, 0.025), lwd=2, lty=3)
    } else
    if(x$Z > M + M*0.1) {
        h1 <- hist(x$Zsim, xlab=expression(paste("Distribution of simulated Z"[sim], "values")), breaks=25, xlim=c(m, M+dif*0.5), col="grey90", border="grey", xaxt = "n", main="")
        title("Ecological null model test", line=1.5) 
        title(paste("P(Obs Z <= null) =", x$pval[2]), line=0.5, font.main=1, cex.main=1)
        locs <- as.integer(c(m, (m+M+dif*0.05)/2, M+dif*0.05, M+dif*0.4))
        axis(1, at=locs[1:3])
        axis(1, labels=FALSE, tck=FALSE)
        axis(1, at=locs[4], labels=as.integer(x$Z))
        abline(v=locs[4], col="red", lty=3, lwd=3)
        polygon(c(M+dif*0.19, M+dif*0.21, M+dif*0.21, M+dif*0.19), rep(c(-max(h1$counts)*0.1, max(h1$counts)*0.1), each=2), border="white", col="white")
        segments(x0=M+dif*0.19, y0=-max(h1$counts)*0.1, y1=max(h1$counts)*0.1)
        segments(x0=M+dif*0.21, y0=-max(h1$counts)*0.1, y1=max(h1$counts)*0.1)
        abline(v=mean(x$Zsim), lwd=2)
        abline(v=quantile(x$Zsim, 0.975), lwd=2, lty=3)
        abline(v=quantile(x$Zsim, 0.025), lwd=2, lty=3)
    } else {
        h1 <- hist(x$Zsim, xlab=expression(paste("Distribution of simulated Z"[sim], "values")), breaks=25, xlim=c(m-dif*0.4, M), col="grey90", border="grey", main="", xaxt = "n")
        title("Ecological null model test", line=1.5)
        title(paste("P(Obs Z <= null) =", x$pval[2]), line=0.5, font.main=1, cex.main=1)
        locs <- as.integer(c(m-dif*0.05, (m+M-dif*0.05)/2, M, m-dif*0.4))
        axis(1, at=locs[1:3])
        axis(1, labels=FALSE, tck=FALSE)
        axis(1, at=locs[4], labels=as.integer(x$Z))
        abline(v=locs[4], col="red", lty=3, lwd=3)
        polygon(c(m-dif*0.19, m-dif*0.21, m-dif*0.21, m-dif*0.19), rep(c(-max(h1$counts)*0.1, max(h1$counts)*0.1), each=2), border="white", col="white")
        segments(x0=m-dif*0.19, y0=-max(h1$counts)*0.1, y1=max(h1$counts)*0.1)
        segments(x0=m-dif*0.21, y0=-max(h1$counts)*0.1, y1=max(h1$counts)*0.1)
        abline(v=mean(x$Zsim), lwd=2)
        abline(v=quantile(x$Zsim, 0.975), lwd=2, lty=3)
        abline(v=quantile(x$Zsim, 0.025), lwd=2, lty=3)
    }

    par(mar=c(5,4,4,2))
    max.n <- max(do.call("rbind", lapply(x$obs, function(i) max(i[,1]))))
    max.x <- ifelse(x$type=="Coverage measure", 1, max.n+0.2*max.n)
    max.y <- max(do.call("rbind", lapply(x$obs, function(i) max(i[,2]))))
    xlab <- ifelse(colnames(x$obs[[1]])[1]=="sample-size", ifelse(x$subclass == "Abundance based method", "Number of individuals", "Number of samples"), "Expected coverage")
    ylab <- ifelse(colnames(x$obs[[1]])[2]=="Hill (q=0)", "Expected species richness", colnames(x$obs[[1]])[2])
    lrc <- do.call("c", lapply(x$obs, function(t) max(t[,1])))
    plot(x$obs[[1]][,2]~x$obs[[1]][,1], type="n", xlab=xlab, ylab=ylab, xlim=c(0, max.x), ylim=c(0, max.y))
    nL <- length(x$obs)
    for(i in 1:nL) {
        lines(x$obs[[i]][,1], x$obs[[i]][,2])
    }
    p <- combn(1:length(x$obs), m = 2)
    for(i in 1:ifelse(dim(p)[2]<=max.poly, dim(p)[2], max.poly)) {
        max.cov <- min(lrc[p[,i]])
        nxax1 <- findInterval(max.cov, x$obs[[p[1,i]]][,1])
        nxax2 <- findInterval(max.cov, x$obs[[p[2,i]]][,1])
        xax1 <- x$obs[[p[1,i]]][1:nxax1,1]
        xax2 <- x$obs[[p[2,i]]][1:nxax2,1]
        yval1 <- x$obs[[p[1,i]]][1:nxax1,2]
        yval2 <- x$obs[[p[2,i]]][1:nxax2,2]
        polygon(cbind(c(xax2, rev(xax1)), c(yval2, rev(yval1))), density=seq(10, 100, 90/(nL-1)), border="black", col="grey")
        }
    title("Observed rarefaction\n curves")
    pos <- ifelse(x$type=="Coverage measure", "topleft", "bottomright")
    legend(pos, fill="grey", density=40, legend=expression('Area A'[i]))

    par(mar=c(5,4,4,2))
    max.n <- max(do.call("rbind", lapply(x$sim, function(i) max(i[,1]))))
    max.x <- ifelse(x$type=="Coverage measure", 1, max.n+0.2*max.n)
    max.y <- max(do.call("rbind", lapply(x$sim, function(i) max(i[,2]))))
    lrc <- do.call("c", lapply(x$sim, function(t) max(t[,1])))
    plot(x$sim[[1]][,2]~x$sim[[1]][,1], type="n", xlab=xlab, ylab=ylab, xlim=c(0, max.x), ylim=c(0, max.y))
    nL <- length(x$sim)
    for(i in 1:nL) {
        lines(x$sim[[i]][,1], x$sim[[i]][,2])
    }
    p <- combn(1:length(x$sim), m = 2)
    for(i in 1:ifelse(dim(p)[2]<=max.poly, dim(p)[2], max.poly)) {
        max.cov <- min(lrc[p[,i]])
        nxax1 <- findInterval(max.cov, x$sim[[p[1,i]]][,1])
        nxax2 <- findInterval(max.cov, x$sim[[p[2,i]]][,1])
        xax1 <- x$sim[[p[1,i]]][1:nxax1,1]
        xax2 <- x$sim[[p[2,i]]][1:nxax2,1]
        yval1 <- x$sim[[p[1,i]]][1:nxax1,2]
        yval2 <- x$sim[[p[2,i]]][1:nxax2,2]
        polygon(cbind(c(xax2, rev(xax1)), c(yval2, rev(yval1))), density=seq(10, 100, 90/(nL-1)), border="black", col="grey")
        }
    title("One randomized set of\nrarefaction curves")
    pos <- ifelse(x$type=="Coverage measure", "topleft", "bottomright")
    legend(pos, fill="grey", density=40, legend=expression('Area A'[i]))

}
