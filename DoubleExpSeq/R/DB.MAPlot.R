DB.MAPlot <-
function (y,m, groups , contrast=c(1,2) , de.tags = NULL, col = "lightgrey" , deCol = "red" , deCex = 0.2 , xlab = "Average Over Groups of log2 Mean Total Count" ,
    ylab = "logFC of Odds Ratio", pch = 19, cex = 0.2, panel.last = grid(col="red",lwd=.2), ylim=c(-15,15) , ... )
{
        cols1 <- groups == unique(groups)[contrast[1]]
        cols2 <- groups == unique(groups)[contrast[2]]
        wh.cols <- c(which(cols1),which(cols2))
        y <- y[,wh.cols]
        m <- m[,wh.cols]
        groups <- as.character(groups[ wh.cols ])
        groups <- as.factor(groups)
        cols1 <- groups == unique(groups)[1]
        cols2 <- groups == unique(groups)[2]
        g1 <- rowMeans(y[,cols1,drop=FALSE],na.rm=TRUE)/rowMeans(m[,cols1,drop=FALSE],na.rm=TRUE)
        g2 <- rowMeans(y[,cols2,drop=FALSE],na.rm=TRUE)/rowMeans(m[,cols2,drop=FALSE],na.rm=TRUE)
        m1 <- rowMeans(m[,cols1,drop=FALSE],na.rm=TRUE)
        m2 <- rowMeans(m[,cols2,drop=FALSE],na.rm=TRUE)

        ylab <- paste(ylab, ":", paste(unique(groups)[c(2,1)], collapse = "-"), paste = "")

        i <- match(de.tags, rownames(y))
        i <- i[!is.na(i)]

        A <- (log2(m1) + log2(m2))/2
        M <- log2(g2/(1-g2)) - log2(g1/(1-g1))
         M[M==Inf] <- 15
         M[M==-Inf] <- -15

        plot(A, M, col = col, ylim=ylim,cex=cex, ...)
                 points(A[abs(M)==15], M[abs(M)==15], col = "orange", cex=cex )
         points(A[de.tags], M[de.tags], col = deCol, cex=deCex , pch=pch )
         abline( h=0 , col="black" , lty="dashed" , lwd = 1.25 , pch=pch )
        invisible( list("M"=M,"A"=A) )
}
