gx.md.plt0 <-
function (md, n, p, trim = trim, ptrim = -1, proc = proc, main = main, 
    ifadd = ifadd, cexf = cexf, cex = cex, ...) 
{
    sorted.md <- sort(md)
    nmd <- length(md)
    fractile <- numeric(nmd)
    for (i in 1:nmd) fractile[i] <- (i - 0.5)/nmd
    chi2 <- qchisq(fractile, p)
    if (proc == " " | proc == "") 
        xlabel <- "Mahalanobis Distance"
    else xlabel <- paste("Robust (", proc, ") Mahalanobis Distance", 
        sep = "")
    ylabel <- paste("Chi-square (df = ", p, ")", sep = "")
    plot(sorted.md, chi2, type = "n", xlab = xlabel, ylab = ylabel, 
        main = main, ...)
    limits <- par("usr")
    xpos <- limits[1] + (limits[2] - limits[1]) * 0.05
    ypos <- limits[4] - (limits[4] - limits[3]) * 0.05
    if (trim <= 0) {
        points(sorted.md[1:nmd], chi2[1:nmd], pch = 3, cex = cex, 
            ...)
        text(xpos, ypos, paste("N =", n), adj = 0, cex = cex, 
            ...)
        if (!is.null(ifadd)) 
            gx.add.chisq(ifadd, df = p, cex = cexf)
    }
    else {
        nc <- n - trim
        points(sorted.md[1:nc], chi2[1:nc], pch = 3, cex = cex, 
            ...)
        points(sorted.md[nc + 1:nmd], chi2[nc + 1:nmd], pch = 1, 
            cex = cex, ...)
        ifflip = TRUE
        xpos <- limits[2] - (limits[2] - limits[1]) * 0.05
        ypos <- limits[3] + (limits[4] - limits[3]) * 0.23
        text(xpos, ypos, paste("N =", n), adj = 1, cex = cex, 
            ...)
        ypos <- limits[3] + (limits[4] - limits[3]) * 0.17
        text(xpos, ypos, paste("Trimmed =", trim), adj = 1, cex = cex, 
            ...)
        atrim <- round((100 * trim)/n)
        ypos <- limits[3] + (limits[4] - limits[3]) * 0.11
        text(xpos, ypos, paste("Trim % =", atrim), adj = 1, cex = cex, 
            ...)
        if (proc == "mvt") {
            ypos <- limits[3] + (limits[4] - limits[3]) * 0.05
            text(xpos, ypos, paste("Requested % =", ptrim * 100), 
                adj = 1, cex = cex, ...)
            ifflip = FALSE
        }
        if (!is.null(ifadd)) 
            gx.add.chisq(0.98, df = p, ifflip = ifflip, cex = cexf)
    }
    invisible()
}
