`plotMissing` <-
function (x, print.labels.SNPs = TRUE, 
    main = "Genotype missing data", ...) 
{ 
   if(!inherits(x,"setupSNP"))
     stop("x must be an object of class 'setupSNP'")

   colSNPs<-attr(x,"colSNPs")
   data.SNPs <- t(x[colSNPs, drop=FALSE])
   label.SNPs<- attr(x,"label.SNPs")
   genInfo<-attr(x,"gen.info")

   data.Missing <- is.na(data.SNPs)
   old.xpd <- par("xpd")
   old.las <- par("las")
   par(xpd = TRUE)
   on.exit(par(xpd = old.xpd, las = old.las))
   image(1:nrow(data.Missing), 1:ncol(data.Missing), data.Missing, 
        col = c("white", "black"), ylab = "Individuals", xlab = ifelse(print.labels.SNPs, 
            "", "SNPs"), axes = !print.labels.SNPs)
   if (print.labels.SNPs) {
        axis(1, at = c(1:length(label.SNPs)), labels = label.SNPs, 
            las = 3, cex.axis = 0.7)
        axis(2)
    }
    title(main, line = 3)
    if (!is.null(genInfo)) 
        n.snps <- table(genInfo[, 2])
    else n.snps <- length(label.SNPs)
    a <- c(0.5, cumsum(n.snps) + 0.5)
    b <- par("usr")
    if (!is.null(genInfo)) 
        col.ok <- c("black", rep("red", length(a) - 1))
    else col.ok <- c("black", rep("black", length(a) - 1))
    segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col = col.ok)
    abline(h = 0.5 + c(0, ncol(data.Missing)), xpd = FALSE)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.snps))
    if (!is.null(genInfo)) {
        segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col = c("black", 
            rep("red", length(a) - 1)))
        names.geno <- unique(genInfo[, 2])
        n.gen <- length(names.geno)
        for (i in 1:n.gen) text(mean(wh[i + c(0, 1)]), a[4] + 
            (a[4] - a[3]) * 0.025, names.geno[i], srt = 45, cex = 0.8, 
            adj = 0.2)
    }
}

