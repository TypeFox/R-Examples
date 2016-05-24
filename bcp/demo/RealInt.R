##### U.S. Ex-Post Interest Rate #####
data("RealInt")
bcp.ri <- bcp(as.vector(RealInt))
plot(bcp.ri, main="U.S. Ex-Post Interest Rate (bcp)", xlab="Time")

# To see bcp and Bai and Perron results:
if (require("strucchange")) {
  bp <- breakpoints(RealInt ~ 1, h = 2)$breakpoints
  rho <- rep(0, length(RealInt))
  rho[bp] <- 1
  b.num<-1 + c(0,cumsum(rho[1:(length(rho)-1)]))
  bp.mean <- unlist(lapply(split(RealInt,b.num),mean))
  bp.ri <- rep(0,length(RealInt))
  for (i in 1:length(bp.ri)) bp.ri[i] <- bp.mean[b.num[i]]
  xax <- seq(1961, 1987, length=103)
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:length(bcp.ri$data[,2]), bcp.ri$data[,2], col="grey", pch=20,
       xlab="", ylab="Posterior Mean", main="U.S. Ex-Post Interest Rate (strucchange)")
  lines(bcp.ri$posterior.mean, lwd=2)
  lines(bp.ri, col="blue")
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(xax, bcp.ri$posterior.prob, yaxt="n", type="l", ylim=c(0,1),
       xlab="Year", ylab="Posterior Probability", main="")
  for (i in 1:length(bp.ri)) abline(v=xax[bp[i]], col="blue")
  axis(2, yaxp=c(0, 0.9, 3))
  par(op3)
  par(op)
} else {
  cat("strucchange is not loaded")
}
