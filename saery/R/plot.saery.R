plot.saery <-
function(direct, eblup, mse.eblup, sigma2edi){
  resid <- direct-eblup
  resid <- resid[order(sigma2edi)]
  devAskNewPage(TRUE)
  par(bty="n")
    plot(resid, type = "l", col=1, ylab = "", xaxt = "n", pch = 1)
    abline(h=0, lty = 2, lwd = 2, col = "mediumseagreen")
    title(main = expression(paste("Residuals sorted by direct estimator variances")), line = 3)
    title(ylab = "residuals")
    axis(1, line = 1, col = "gray70")
    axis(2, col = "gray70")
  
    plot(direct, eblup, type = "p", col=1, xlab = "direct", ylab = "eblup", xaxt = "n", pch = 1)
    abline(a=0, b=1, lty = 2, lwd = 2, col = "mediumseagreen")
    title(main=expression(paste("Estimates")), line = 3)
    axis(1, line = 1, col = "gray70")
    axis(2, col = "gray70")
  
    sqrt.mse.eblup <- sqrt(mse.eblup[order(sigma2edi)])
    sqrt.sigma2edi <- sqrt(sort(sigma2edi))
    u <- max(sqrt.mse.eblup, sqrt(sigma2edi))
    l <- min(sqrt.mse.eblup, sqrt(sigma2edi))
    plot(sqrt.sigma2edi, type="b", col = 2, ylim=c(l,u), ylab="", xaxt = "n", pch=1)
    lines(sqrt.mse.eblup, type="b", col = 4, pch=4)
    title(main="Root-MSE estimates", line = 3)
    axis(1, line = 1, col = "gray70")
    axis(2, col = "gray70")
    text <- c("direct", "eblup")
    legend(1, u/1.05, text, col = c(2,4), pch = c(1, 4), bty="n", lwd = 2)
  
}
