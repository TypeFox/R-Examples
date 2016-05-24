fcrosChrPlot <- function(chrData, thr = 0.05, deb = 100, fin = 1e7) {
    a1 <- 0.5*thr
    a2 <- 1 - 0.5*thr
    positions <- chrData$Position
    ndata <- length(positions)
    mega <- 1e6
    if ((deb < 0) || (fin < 0) || (fin < deb)) {
       chr <- unique(chrData$Chromosome)
       xax <- chrData$Position/mega
       yax <- chrData$f.L2R
       ymin <- min(yax); ymin <- min(ymin, -1.1)
       ymax <- max(yax); ymax <- max(ymax, 1.1)
       fax <- chrData$f.value
       dax <- (fax <= a1); dp = c(rep(-1, length(fax[dax])))
       gax <- (fax >= a2); gp = c(rep(1, length(fax[gax])))
    } else {
       if (deb < positions[1]) { deb <- positions[1] }
       if (fin > positions[ndata]) { fin <- positions[ndata] }
       i <- 1;        t1 <- i;
       while (deb > positions[i]) { t1 <- i; i <- i+1 }
       t2 <- i
       while (positions[i] < fin) { t2 <- i; i <- i+1 }
       chr <- unique(chrData$Chromosome[t1:t2])
       xax <- chrData$Position[t1:t2]/mega
       yax <- chrData$f.L2R[t1:t2]
       ymin <- min(yax); ymin <- min(ymin, -1.1)
       ymax <- max(yax); ymax <- max(ymax, 1.1)
       fax <- chrData$f.value[t1:t2]
       dax <- (fax <= a1); dp = c(rep(-1, length(fax[dax])))
       gax <- (fax >= a2); gp = c(rep(1, length(fax[gax])))
    }

    plot(xax, yax, ylim = c(ymin, ymax), cex = .8, xlab = "x 10^6", ylab = "Log2 Ratio", main = chr)
    points(xax[dax], dp, col = "blue", type = "p", cex = 0.6, pch = 24)
    points(xax[gax], gp, col = "red", type = "p", cex = 0.6, pch = 25)
    abline(h = 0, col = "black")
}
