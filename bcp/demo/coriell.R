##### Coriell chromosome 11 #####
data(coriell)
chrom11 <- as.vector(na.omit(coriell$Coriell.05296[coriell$Chromosome==11]))
bcp.11 <- bcp(chrom11)
plot(bcp.11, main="Coriell chromosome 11 (bcp)")

# To see bcp and Circular Binary Segmentation results, using
# base graphics (see plot.bcp.legacy for more examples):
if(require("DNAcopy")) {
  n <- length(chrom11)
  cbs <- segment(CNA(chrom11, rep(1, n), 1:n), verbose = 0)
  cbs.ests <- rep(unlist(cbs$output[6]), unlist(cbs$output[5]))
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:n, bcp.11$data[,2], col="grey", pch=20, xlab="Location",
       ylab="Posterior Mean",
       main="Coriell chromosome 11 (DNAcopy)")
  lines(cbs.ests, col="red")
  lines(bcp.11$posterior.mean, lwd=2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(1:n, bcp.11$posterior.prob, type="l", ylim=c(0,1),
       xlab="Location", ylab="Posterior Probability", main="")
  for (i in 1:(dim(cbs$output)[1]-1)) {
    abline(v=cbs$output$loc.end[i], col="red")
  }
  par(op3)
  par(op)
} else {
  cat("DNAcopy is not loaded")
}