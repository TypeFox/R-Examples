ad.plot3 <-
function (x1, x2, xname = deparse(substitute(x1)), if.order = FALSE, 
    ad.tol = NULL, ldl = NULL, maxrat = NULL, if.text = FALSE, 
    if.cpp = FALSE, ...) 
{
    if (length(x1) != length(x2)) 
        stop("The lengths of the vectors are not the same\n")
    temp.x <- remove.na(cbind(x1, x2))
    a1 <- temp.x$x[1:temp.x$n, 1]
    a2 <- temp.x$x[1:temp.x$n, 2]
    n <- temp.x$n
    frame()
    xrat <- a1/a2
    xbar <- (a1 + a2)/2
    rat.med <- median(xrat)
    rat.mad <- mad(xrat)
    rat.rrsd <- 100 * rat.mad /rat.med
    rat.bar <- mean(xrat)
    rat.sd <- sqrt(var(xrat))
    rat.rsd <- 100 * rat.sd / rat.bar
    rat.range <- range(xrat)
    rat.sem <- rat.sd/sqrt(n)
    t05 <- qt(0.975, n - 1)
    sem.limit <- rat.sem * t05
    t05.limit <- signif(t05 * rat.sd, 4)
    fact.hi <- 1 + t05.limit
    text1 <- signif(1/fact.hi, 2)
    text2 <- signif(fact.hi, 3)
    cat("\n  Range of ratios of the", n, "duplicates:", signif(rat.range[1], 
        3), "to", signif(rat.range[2], 3))
    cat("\n  Median ratio =", round(rat.med, 4), "\tMAD of ratios =", 
        signif(rat.mad, 4), "\trRSD % =", round(rat.rrsd, 1))
    cat("\n  Mean ratio   =", round(rat.bar, 4), "\t SD of ratios =", 
        signif(rat.sd, 4), "\t RSD % =", round(rat.rsd, 1),
        "\t95% CI =", t05.limit)
    delta <- abs(rat.bar - 1)
    cat("\n\n  Absolute ratio difference from 1 =", signif(delta, 4),"\tSE of Mean of ratios =",
        signif(rat.sem, 3), "\t95% CI =", signif(sem.limit, 3))
    if (delta <= sem.limit) 
        cat("\n  Mean ratio is not different from 1 at the 95% level, no bias\n")
    else cat("\n  Mean ratio is different from 1 at the 95% level, duplicates biased\n")
    if (is.null(maxrat)) 
        maxrat <- max(1/rat.range[1], rat.range[2])
    if (!if.cpp) {
        if (if.order) {
            plot(seq(1, n), xrat, ylim = c(1/maxrat, maxrat), 
                xlab = paste("Ordered determinations of", xname), 
                ylab = "Ratio of duplicates", log = "y")
            abline(h = 1, lty = 2)
            if (!is.null(ad.tol)) {
                cat(paste("\n  \261 ", ad.tol, "% tolerance lines on the ratio plotted\n", 
                  sep = ""))
                abline(h = (1 + ad.tol/100), lty = 3, col = 2)
                abline(h = (1 - ad.tol/100), lty = 3, col = 2)
            }
        }
        else {
            plot(xbar, xrat, ylim = c(1/maxrat, maxrat), xlab = paste("Mean of duplicates for", 
                xname), ylab = "Ratio of duplicates", log = "xy")
            abline(h = 1, lty = 2)
            abline(h = rat.bar, lty = 2, ... )
            abline(h = rat.bar * text2, lty = 3, ...)
            abline(h = rat.bar * text1, lty = 3, ...)
            cat("\n  95% of duplicates will fall between factors of", 
                text2, "and", text1, "times a value")
            t05.limit <- signif(qt(0.975, n - 1) * rat.mad, 4)
            fact.hi <- 1 + t05.limit
            text3 <- signif(1/fact.hi, 2)
            text4 <- signif(fact.hi, 3)
            cat("\n    Robust factor estimates based on the MAD are", 
                text4, "and", text3, "\n")
            if (!is.null(ldl)) 
                abline(v = ldl, lty = 3, col = 2)
            if (if.text) {
                text <- paste("95% of duplicates fall between factors of", 
                  text2, "and", text1, "times a value,\n                               ", 
                  "robust estimates are", text4, "and", text3)
                text(locator(1), text, adj = 0, cex = 0.9)
            }
        }
    }
    else {
        oldpar <- par()
        on.exit(par(oldpar))
        par(mfrow = c(1, 2), pty = "s")
        plot(xbar, xrat, ylim = c(1/maxrat, maxrat), xlab = paste("Mean of duplicates for", 
            xname), ylab = "Ratio of duplicates", log = "xy", 
            ...)
        abline(h = 1, lty = 2)
        abline(h = rat.bar, lty = 2, ... )
        abline(h = rat.bar * text2, lty = 3, ...)
        abline(h = rat.bar * text1, lty = 3, ...)
        if (!is.null(ldl)) 
            abline(v = ldl, lty = 3, col = 2)
        cnpplt(xrat, xlab = "Ratio of duplicates", xlim = c(1/maxrat, 
            maxrat), log = T, ifshape = T, ...)
    }
    cat("\n")
    invisible()
}
