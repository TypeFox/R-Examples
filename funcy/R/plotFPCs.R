#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

plotFPCs <- function(plotParams, legendPlace){
    t_unique <- plotParams$t_unique
    mean.timeIn <- plotParams$mean.timeIn
    mu.smoothed <- plotParams$mu.smoothed
    mean <- plotParams$mean.raw
    cov.timeIn <- plotParams$cov.timeIn
    cov.raw <- plotParams$cov.raw
    cov.timeOut <- plotParams$cov.timeOut
    cov.smoothed <- plotParams$cov.smoothed
    tildeG <- plotParams$tildeG
    hatV <- plotParams$hatV
    base <- plotParams$base
    nt <- plotParams$nt
    
    par(mfrow=c(2,3),oma=c(0,0,2,0))
    
    ##Plot mu and smoothed curve
    minY <- min(mu.smoothed, mean)
    maxY <- max(mu.smoothed, mean)
    plot(t_unique, mu.smoothed, col="red", type='l' ,xlab="time", ylab="",
         main="Smoothed mean function", ylim=c(minY, maxY))
    points(mean.timeIn, as.numeric(mean))
    legend("topright",legend=c("smoothed mu","original mu"),col=c("red","black"), box.lty=1, lwd=3)
    ##Plot Cov
    if(!requireNamespace("scatterplot3d"))
        stop("Please install scatterplot3d to use this plot.")
    
    scatterplot3d::scatterplot3d(cov.timeIn[,1], cov.timeIn[,2], cov.raw, xlab="x",
                                 ylab="y", zlab="cov(x,y)", main="Raw covariance")
    ##Plot Cov smoothed
    scatterplot3d::scatterplot3d(cov.timeOut[,1],cov.timeOut[,2],  cov.smoothed, xlab="x",
                                 ylab="y", zlab="cov(x,y)", main="Smoothed reduced covariance")
    ##Plot Sigma
    plot(t_unique, tildeG, main="Sigma", xlab="time",
         xlim=range(t_unique), ylim=c(min(tildeG,hatV),max(tildeG,hatV)), type="b")
    lines(t_unique, hatV, col=2, lwd=3)
    suppressWarnings(arrows(t_unique, tildeG,t_unique,
                            tildeG+(hatV-tildeG), length = 0.05, col=3))
    legend("topleft", legend=c("smoothed cov on diagonal","smoothed var on diagonal"),col=c("red","black"), box.lty=1, lwd=3)
    
    matplot(base[,1:min(dim(base)[2],3)],type='l',
            main="First 3 eigenfunctions", xlab="time")
    title("Generating functional principal components", outer=TRUE)
    par(mfrow=c(1,1))
    
}


