plotConv <-
function (bmao, include.legend = TRUE, add.grid = TRUE, ...) 
{
    if (!is.bma(bmao)) 
        stop("submit an object of class bma")
    mat = pmp.bma(bmao, oldstyle = TRUE)
    norm_const = sum(mat[, 1])/sum(mat[, 2])
    mat = cbind(mat[, 2] * norm_const, mat[, 1])
    if (length(bmao$topmod$lik()) == 0L) {
        stop("plotConv needs at least one model stored in topmod in order to produce a plot")
    }
    cor.pmp = format(round(.cor.topmod(bmao$topmod), 4), nsmall = 4)
    dotargs = match.call(plot, expand.dots = FALSE)$...
    dotargs = .adjustdots(dotargs, lwd = 2, main = paste("Posterior Model Probabilities\n(Corr: ", 
        cor.pmp, ")", sep = ""), lty = 1, col = c("steelblue3", 
        "tomato"), cex.main = 0.8, xlab = "Index of Models", 
        ylab = "", type = "l")
    eval(as.call(c(list(as.name("matplot"), as.name("mat")), 
        as.list(dotargs))))
    if (as.logical(add.grid)) 
        grid()
    if (as.logical(include.legend)) 
        legend("topright", lty = eval(dotargs$lty), legend = c("PMP (MCMC)", 
            "PMP (Exact)"), col = eval(dotargs$col), ncol = 2, 
            bty = "n", cex = 1, lwd = eval(dotargs$lwd))
}
