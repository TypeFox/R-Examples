crm.plot <-
function (xx, xname = deparse(substitute(xx)), crm.mean = NULL, 
    crm.sd = NULL, n.sd = 2, crm.tol = NULL, ...) 
{
    # Function to display CRM results for QA/QC appraisal.  The expected
    # mean must be provided, with either an associated standard deviation
    # or a percentage acceptable tolerance.  A default of 2 SD is provided
    # for calculating tolerance bounds around the mean.
    #
    if (is.null(crm.mean)) 
        stop("The accepted CRM mean must be provided\n")
    if (is.null(crm.sd) && is.null(crm.tol))
        stop("Either the crm.sd or crm.tol must be provided\n")
    temp.x <- remove.na(xx)
    x <- temp.x$x[1:temp.x$n]
    n <- temp.x$n
    xmean <- mean(x)
    xsd <- sqrt(var(x))
    xrsd <- round(100 * xsd/xmean, 2)
    xmed <- median(x)
    xmad <- mad(x)
    rrsd <- round(100 * xmad/xmed, 2)
    if (is.null(crm.sd)) 
        n.sd = NULL
    cat(" Control Reference Material:", xname, "\n\t\t\t\t  Mean =", 
        crm.mean, "\t SD =", crm.sd, "\t N.SD =", n.sd, "\tPercent tolerance =", 
        crm.tol, "\n")
    cat(" This project:\t N =", n, "\t  Mean =", signif(xmean, 3), "\t SD =", 
        signif(xsd, 3), "\t RSD% =", xrsd, "\n\t\t\t\tMedian =", signif(xmed, 3),
        "\tMAD =", signif(xmad, 3), "\trRSD% =", rrsd, "\n")
    x.range <- range(x)
    if (!is.null(crm.sd)) {
        hitol <- crm.mean + n.sd * crm.sd
        lotol <- crm.mean - n.sd * crm.sd
        plot.min <- min(x.range[1], lotol)
        plot.max <- max(x.range[2], hitol)
        ylabel <- paste(xname, ", red lines indicate +/-", n.sd, 
            " SD limits", sep = "")
        plotid <- 1
    }
    if (!is.null(crm.tol)) {
        x <- 100 * (abs(x - crm.mean)/crm.mean)
        plot.min <- 0
        plot.max <- max(max(x), crm.tol)
        ylabel <- "Percentage absolute difference relative to accepted value"
        plotid <- 2
    }
    plot(seq(1:n), x, ylim = c(plot.min, plot.max), ylab = ylabel,
        xlab = paste("Ordered determinations of", xname), ...)
    if (plotid == 1) {
        abline(h = crm.mean, col = 3, lty = 3)
        abline(h = hitol, col = 2, lty = 2)
        abline(h = lotol, col = 2, lty = 2)
    }
    else {
        abline(h = crm.tol, col = 2, lty = 2)
    }
    invisible()
}
