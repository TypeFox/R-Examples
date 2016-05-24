fcrosChrSummary <- function(af, xinfo, chromosomes = c(1:22,"X","Y"), alpha = 0.05) {
   tmp <- xinfo
   n <- nrow(tmp)
   tmp$Position <- 0.5*(as.numeric(xinfo$Start) + as.numeric(xinfo$End))
   tmp$f.L2R <- log2(af$FC2)
   tmp$f.value <- af$f.value
   a1 <- 0.5*alpha
   a2 <- 1-0.5*alpha
   seg <- c(rep(0, n))
   seg[af$f.value <= a1] <- -1
   seg[af$f.value >= a2] <- 1
   tmp$f.call <- seg
   chrSumm <- matrix(c(rep(0,3*length(chromosomes))), ncol = 3)
   colnames(chrSumm) <- c("Chr", "nLoss", "nGain")
   idx <- which(tmp$Chromosome == chromosomes[1])
   xda <- tmp[idx,]
   xdr <- order(xda$Position)
   xinfo.s <- xda[xdr,]
   chrSumm[1,1] <- chromosomes[1]
   chrSumm[1,2] <- sum(xda$f.call == -1)
   chrSumm[1,3] <- sum(xda$f.call == 1)
   for (i in 2:length(chromosomes)) {
       idx <- which(tmp$Chromosome == chromosomes[i])
       xda <- tmp[idx,]
       xdr <- order(xda$Position)
       xinfo.s <- rbind(xinfo.s, xda[xdr,])
       chrSumm[i,1] <- chromosomes[i]
       chrSumm[i,2] <- sum(xda$f.call == -1)
       chrSumm[i,3] <- sum(xda$f.call == 1)
   }
   list(xinfo.s = xinfo.s, chrSumm = chrSumm)
}
