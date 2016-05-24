# plot method for Permutation smacof

plot.smacofPerm <- function(x, alpha = 0.05, main, xlab, ylab, ...)
{
# x ... object of class smacofPerm
    
    if (missing(main)) main <- "ECDF SMACOF Permutation Test" else main <- main
    if (missing(xlab)) xlab <- "Stress Values" else xlab <- xlab
    if (missing(ylab)) ylab <- "Probability" else ylab <- ylab
    
    Ecdf(x$stressvec, main = main, xlab = xlab, ylab = ylab, subtitles = FALSE, ...)
    abline(v = x$stress.obs, col = "gray", lty = 1)
    text(x$stress.obs, y = 1, labels = paste("Stress: ",round(x$stress.obs, 3), sep = ""), col = "gray", pos = 4, cex = 0.7)
    abline(h = x$pval, col = "gray", lty = 1)  
    text(max(x$stressvec)-0.01*max(x$stressvec), y = x$pval, labels = paste("p-value: ",round(x$pval, 3), sep = ""), col = "gray", pos = 3, cex = 0.7)
    abline(h = alpha, col = "gray", lty = 2) 
}


  


  

 
