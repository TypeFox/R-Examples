# plot percentile CI for matrix obtained from simSigmaP()
plotCIsigmaP <- 
function(msim, conf.level = 0.95, shade.col = "orange", 
	ordered = TRUE, xlim = NULL, xlab = expression(sigma[P]),
	las = 1, mar = c(4.5, 6.5, 2, 1), ...)
{  
    if (!inherits(msim, "simSigmaP"))
       stop("'msim' must be an object of class simSigmaP")
 
    # stats  
    med <- apply(msim, 2, mean)
    cv <- function(x) 100 * sd(x) / mean(x)
    cvs <- apply(msim, 2, cv)
    alpha = 1 - conf.level
    ll <- apply(msim, 2, quantile, alpha/2)
    ul <- apply(msim, 2, quantile, 1 - alpha/2)

    out <- data.frame(means = med, CV = cvs, 
       LCL = ll, UCL = ul, methods = names(med))
    if (ordered) {
       o <- order(out[, 1])
       out <- out[o, ]
    }

    # graph
    par(mar = mar, las = las, ...)
    plot(out[, 1], seq(1, nrow(out)), 
       yaxt = "n", ylab = "", 
       ylim = c(0, ncol(msim)+1),
       xlim = if (is.null(xlim)) c(min(ll), max(ul)) else xlim,
       xlab = xlab, ...)
    axis(2, at = 1:nrow(out), labels = out[, 5])
    abline(h = seq(1, nrow(out)), col = "gray")
    for(i in 1:nrow(out)) {
       polygon(y = c(-100, -100, i, i), 
          x = c(out[i, 4], out[i, 3], out[i, 3], out[i, 4]), 
          border = "gray", lty = 3,
          col = adjustcolor(shade.col, alpha.f = 0.1))
       arrows(out[i, 3], i, out[i, 4], i, 
          length = 0.05, angle = 90, code = 3)
    }
    invisible(out)
}
