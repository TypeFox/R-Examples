plot.bacr <- function(x, ...){
    vv = colMeans(x$MY)
    vv = vv[1:(length(vv)-1)]
    #    par(mar=c(5.1, 7.1, 4.1, 2.1),cex.lab=2,cex.axis=1.5)
    plot(vv,xlab="", ylab="Posterior Inclusion Probability", xaxt="n", ...)
    axis(1,at=1:length(vv), x$predictorsY[1:(length(x$predictorsY)-1)],las=3)
}