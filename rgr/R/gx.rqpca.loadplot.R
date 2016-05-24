gx.rqpca.loadplot <-
function (save, main = "", crit = 0.3, cex = 0.8, cex.axis = 0.7, 
    cex.main = 0.8) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(pty = "m", las = 1)
    if (main == "") 
        banner <- paste("PC loadings >", crit, "for", deparse(substitute(save)), 
            "\ndata source:", save$input)
    else banner <- main
    l <- save$rload
    range.l <- range(l)
    ly <- -1; uy <- 1 
    if (range.l[1] < -0.975) ly <- range.l[1] - 0.1
    if (range.l[2] > 0.975) uy <- range.l[2] + 0.1
    k <- dim(l)[2]
    p <- dim(l)[1]
    lnam <- save$matnames[[2]]
    plot(cbind(c(0, 1, 1, 0, 0), c(uy, uy, ly, ly, uy)), type = "l", 
        axes = FALSE, xlab = "", ylab = "")
    segments(0, 0, 1, 0)
    if (uy > 1) segments(0, 1, 1, 1, lty = 2)
    segments(0, 0.5, 1, 0.5, lty = 2)
    segments(0, -0.5, 1, -0.5, lty = 2)
    if (ly < -1) segments(0, -1, 1, -1, lty = 2)
    tplace1 = -0.3
    mtext("-1", side = 2, at = -1, line = tplace1, cex = cex.axis)
    mtext("-0.5", side = 2, at = -0.5, line = tplace1, cex = cex.axis)
    mtext("0", side = 2, at = 0, line = tplace1, cex = cex.axis)
    mtext("+0.5", side = 2, at = 0.5, line = tplace1, cex = cex.axis)
    mtext("+1", side = 2, at = 1, line = tplace1, cex = cex.axis)
    title(banner, cex.main = cex.main)
    bb <- apply(l^2, 2, sum)/sum(l^2)
    bb1 <- cumsum(bb)
    cumpct <- cumsum(save$econtrib)
    mtext("0%", side = 3, at = 0, line = tplace1, cex = 0.7)
    tplace2 = -0.5
    for (i in 1:k) {
        segments(bb1[i], uy, bb1[i], ly)
        lplot <- abs(l[, i]) > crit
        lsel <- l[lplot, i]
        names(lsel) <- lnam[lplot]
        if (i == 1) {
            mtext(paste(round(cumpct[i]), "%", sep = ""), side = 3, 
                at = bb1[i], line = tplace1, cex = cex.axis)
            chardist <- bb[1]/(length(lsel) + 1)
            text(seq(from = chardist, by = chardist, length = length(lsel)), 
                lsel, names(lsel), cex = cex)
            mtext(paste("PC-", round(i), sep = ""), side = 1, 
                at = bb1[i]/2, line = tplace2, cex = cex.axis)
        }
        else {
            if (length(lsel) >= 1) {
                mtext(paste(round(cumpct[i]), "%", sep = ""), 
                  side = 3, at = bb1[i], line = tplace1, cex = cex.axis)
                chardist <- (bb1[i] - bb1[i - 1])/(length(lsel) + 
                  1)
                text(seq(from = bb1[i - 1] + chardist, by = chardist, 
                  length = length(lsel)), lsel, names(lsel), 
                  cex = cex)
                mtext(paste("PC-", round(i), sep = ""), side = 1, 
                  at = bb[i]/2 + bb1[i - 1], line = tplace2, 
                  cex = cex.axis)
            }
        }
    }
    invisible()
}
