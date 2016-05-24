asplot <- function (locus, map, genes, flanking=1e3, best.pval=NULL, sf=c(4,4), logpmax=10, pch=21)
{
    NAME <- locus$NAME
    PVAL <- locus$PVAL
    snp <- NAME[PVAL==min(PVAL)]
    chr <- locus$CHR[1]
    hit <- locus[NAME==snp, ]
    if (is.null(best.pval)) best.pval <- hit$PVAL
    lb <- min(locus$POS) - flanking
    ub <- max(locus$POS) + flanking
    lu <- ub - lb
    center <- lb + lu / 2
    center.100kb <- round(center / 1e5) * 1e5
    offset.100kb <- round(lu / 5 / 1e5) * 1e5
    offset <- logpmax / sf[1]
    ylim <- logpmax + offset
    ylim4 <- ylim / 4
    yadj <- -offset + ylim / sf[2]
    keep <- subset(map, map$POS > lb & map$POS < ub)
    genes <- with(genes, subset(genes, (START > lb & START < ub) | (STOP > lb & STOP < ub)))
    par(mar = c(4, 4, 3, 4))
    xy <- xy.coords(keep$POS, yadj + (keep$THETA/60) * (3 * ylim4))
    plot(xy$x, xy$y, type = "l", col = "lightblue", lwd = 1, ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F)
    box()
    if (offset.100kb!=0)  center5 <- center.100kb + (-2:2) * offset.100kb
    else
    {
       p1 <- min(xy$x)
       p5 <- max(xy$x)
       p3 <- (p1+p5)/2
       p2 <- (p1+p3)/2
       p4 <- (p3+p5)/2
       center5 <- c(p1,p2,p3,p4,p5)
    }
    axis(1, at = center5, labels = round(center5 / 1e3, 2), las = 1)
    mtext(paste("Chromosome", chr, "position (kb)", sep = " "), side = 1, line = 2.5)
    axis(2, at = seq(0, logpmax, 2), labels = seq(0, logpmax, 2), las = 1)
    mtext("-log10(Observed p)", side = 2, at = logpmax / 2, line = 2)
    axis(4, at = yadj + (0:4) * ylim4, labels = paste((0:4)*20), las = 1, col.axis="blue")
    mtext("Recombination rate (cM/Mb)", side = 4, at = -offset + 2 * ylim4, line = 2, col="blue")
    lines(c(lb, ub), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    points(hit$POS, -log10(hit$PVAL), pch = pch, cex = 2.5, bg = "red")
    text(hit$POS, -log10(hit$PVAL), labels = hit$NAME, pos = 3, offset = 1)
    if (-log10(best.pval) < logpmax)
    {
       points(hit$POS, -log10(best.pval), pch = pch, cex = 2.5, bg = "blue")
       text(hit$POS, -log10(best.pval), labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 2)
    } else {
       points(hit$POS, logpmax, pch = pch, cex = 2.5, bg = "blue")
       text(hit$POS, logpmax, labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 1)
    }
    RSQR <- locus$RSQR
    strong.ld <- subset(locus, NAME != snp & RSQR >= 0.8)
    moderate.ld <- subset(locus, NAME != snp & RSQR >= 0.5 & RSQR < 0.8)
    weak.ld <- subset(locus, NAME != snp & RSQR >= 0.2 & RSQR < 0.5)
    not.in.ld <- subset(locus, NAME != snp & RSQR < 0.2)
    na.ld <- subset(locus, NAME != snp & is.na(RSQR))
    colors <- c("red","orange","yellow","grey","white")
    points(strong.ld$POS, -log10(strong.ld$PVAL), pch = pch, cex = 1.25, bg = colors[1])
    points(moderate.ld$POS, -log10(moderate.ld$PVAL), pch = pch, cex = 1.25, bg = colors[2])
    points(weak.ld$POS, -log10(weak.ld$PVAL), pch = pch, cex = 1.25, bg = colors[3])
    points(not.in.ld$POS, -log10(not.in.ld$PVAL), pch = pch, cex = 1, bg = colors[4])
    points(na.ld$POS, -log10(na.ld$PVAL), pch = pch, cex = 1, bg = colors[5])
    for (i in 1:nrow(genes))
    {
        gi <- genes[i, ]
        GENE <- gi$GENE
        cat("-",GENE,"\n")
        start <- gi$START
        finish <- gi$STOP
        center <- (start + finish)/2
        lb <- min(xy$x)
        ub <- max(xy$x)
        adj <- -offset+2*(i%%4)/3
        if (!is.na(GENE))
        {  
           if (start<lb) text(lb, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
           else if (finish>ub) text(ub, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
           else text(center, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
        }
        STRAND <- gi$STRAND
        if (STRAND == "+") arrows(start, adj, finish, adj, length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
        else arrows(start, adj, finish, adj, length = 0.05, lwd = 2, code = 1, lty = "solid", col = "darkgreen")
    }
    ltext <- rbind("0.8-","0.5-","0.2-","0.0-","NA")
    legend(min(xy$x),logpmax,legend=ltext,title=expression(r^2),fill=colors)
}
