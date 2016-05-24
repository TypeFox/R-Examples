fcrosChrPlot2 <- function(chrData, chrSeg, deb = 100, fin = 1e7) {
    positions <- chrData$Position
    ndata <- length(positions)
    mega <- 1e6
    lBounds <- chrSeg[,3]
    uBounds <- chrSeg[,4]
    if ((deb < 0) || (fin < 0) || (fin < deb)) {
       chr <- unique(chrData$Chromosome)
       xax <- positions/mega
       yax <- chrData$f.L2R
       ymin <- min(yax); ymin <- min(ymin, -1.1)
       ymax <- max(yax); ymax <- max(ymax, 1.1)
       iseg_s <- chrSeg[,1]
       iseg_e <- chrSeg[,2]
       segVal <- chrSeg[,5]
       nseg <- length(iseg_s)
    } else {          
       if (deb < positions[1]) { deb <- positions[1] }
       if (fin > positions[ndata]) { fin <- positions[ndata] }
       i <- 1;       t1 <- i;
       while (deb > positions[i]) { t1 <- i; i <- i+1 }       
       t2 <- i
       while (positions[i] < fin) { t2 <- i; i <- i+1 }
       chr <- unique(chrData$Chromosome[t1:t2])
       xax <- positions[t1:t2]/mega
       yax <- chrData$f.L2R[t1:t2]
       ymin <- min(yax); ymin <- min(ymin, -1.1)
       ymax <- max(yax); ymax <- max(ymax, 1.1)
       iseg_s <- chrSeg[,1]
       iseg_e <- chrSeg[,2]
       segVal <- chrSeg[,5]
       nseg <- length(iseg_s)
       if ((iseg_e[nseg] < t1) || (iseg_s[1] > t2)) {
          t3 <- 0
       } else if ((iseg_s[1] > t1) && (iseg_e[nseg] < t2)) {
          t3 <- 1; t4 <- nseg
       } else if ((iseg_s[1] < t1) && (iseg_e[nseg]<t2)) {
          i <- 1;  t3 <- i
          while (iseg_s[i] < t1) { t3 <- i+1; i <- i+1}
          t4 <- nseg
       } else if ((iseg_s[1] > t1) && (iseg_e[nseg] > t2)) {
          t3 <- 1
          i <- nseg; t4 <- i
          while (iseg_s[i] > t2) { t4 <- i-1; i <- i-1}
       } else {
          i <- 1;  t3 <- i
          while (iseg_s[i] < t1) { t3 <- i+1; i <- i+1}
          i <- nseg; t4 <- i
          while (iseg_s[i] > t2) { t4 <- i-1; i <- i-1}
       }
       if (t3) {
          iseg_s <- iseg_s[t3:t4]
          iseg_e <- iseg_e[t3:t4]
          segVal <- segVal[t3:t4]
          nseg <- length(iseg_s)
       } else {
          nseg <- 0
       }
    }

    plot(xax, yax, ylim = c(ymin, ymax), cex = .8, xlab = "x 10^6", ylab = "Log2 Ratio", main = chr)
    abline(h = 0, col = "black")
    if (nseg) {
       for (i in 1:nseg) {
           fax <- segVal[i]
           np <- iseg_e[i]-iseg_s[i]+1
           xax <- positions[iseg_s[i]:iseg_e[i]]/mega
           if (fax <0) {
              yax <- c(rep(-1, np))
              lines(xax, yax, col = "blue", type = "b", cex = 0.6, pch = 24)
           } else {
             yax <- c(rep(1, np))
             lines(xax, yax, col = "red", type = "b", cex = 0.6, pch = 25)
           }
       }
    } else {
    }
}
