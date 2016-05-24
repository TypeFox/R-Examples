crm.plot.new <-
function (xx, xname = deparse(substitute(xx)), crm.mean = NULL, 
    crm.sd = NULL, n.sd = 2, crm.new = 0, ylim = NULL, ...) 
{
    temp.x <- remove.na(xx)
    x <- temp.x$x[1:temp.x$n]
    n <- temp.x$n
    xmean <- mean(x)
    xsd <- sqrt(var(x))
    xrsd <- round(100 * xsd/xmean, 2)
    xmed <- median(x)
    xmad <- mad(x)
    rrsd <- round(100 * xmad/xmed, 2)
    x.range <- range(x)
    cat("  Control Reference Material:", xname, "\n\t\t\t\t  Mean =", 
        signif(xmean, 3), "\t SD =", signif(xsd, 3), "\t RSD% =", xrsd,
        "\t  N =", n)
    cat("\n\t\t\t\tMedian =", signif(xmed, 3), "\tMAD =", signif(xmad, 3),
        "\trRSD% =", rrsd, "\n")
#
    if(!is.null(crm.mean)) cat("  CRM mean provided =", crm.mean, "\n")
    else crm.mean <- xmean
    if(!is.null(crm.sd)) cat("  CRM SD provided =", crm.sd, "\n")
    else crm.sd <- xsd
    cat("  Number of new CRM determinations to be plotted =", crm.new, "\n")
    cat("  Number of SDs to be used for tolerance bounds =", n.sd, "\n")
#
    hitol <- crm.mean + n.sd * crm.sd
    lotol <- crm.mean - n.sd * crm.sd
    plot.min <- min(x.range[1], lotol, ylim[1])
    plot.max <- max(x.range[2], hitol, ylim[2])
    ylabel <- paste(xname, ", red lines indicate \261", n.sd, " SD limits",
        sep = "")
#    
    plot(seq(1:n), x, xlim = c(1, n+crm.new), ylim = c(plot.min, plot.max), 
         ylab = ylabel, xlab = paste("Ordered determinations of", xname), ...)
    abline(h = crm.mean, col = 3, lty = 3)
    abline(h = hitol, col = 2, lty = 2)
    abline(h = lotol, col = 2, lty = 2)
#
    invisible()
}
