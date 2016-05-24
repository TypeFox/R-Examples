## produces plots used in KojYan 2010 (IME)

load("frankRho.rda")
library(copula, lib.loc="../../copula.Rcheck")
source("../../copula/R/debye.R")

dRhoFrankCopula <- function(copula) {
  alpha <- copula@parameters
  return( 12 / (alpha * (exp(alpha) - 1)) - 36 / alpha^2 * debye2(alpha) + 24 / alpha^2 * debye1(alpha) )
}

thetaGrid <- seq(-.999, .999, by = .001)
alphaGrid <- 40 * atanh(thetaGrid)
rhoTrue <- sapply(alphaGrid, function(x) rho(frankCopula(x)))
dRhoTrue <- sapply(alphaGrid, function(x) dRhoFrankCopula(frankCopula(x)))

pdf("frank-rho.pdf", height=3, width=6, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0.5))
plot(thetaGrid, rhoTrue, type="l", xlab=expression(alpha), ylab=expression(rho))
curve(frankRhoFun(atanh(x) * 40), add=TRUE, col="blue")
legend("topleft", legend=c("true", "numerical"), col=c("black", "blue"), lty=c(1,2), cex=0.7)

curve(frankdRho(atanh(x) * 40), col="blue", xlab=expression(alpha), ylab=bquote(rho~"'"))
lines(thetaGrid, dRhoTrue)
legend("topleft", legend=c("true", "numerical"), col=c("black", "blue"), lty=c(1,2), cex=0.7)
dev.off()


pdf("frank-rho-err.pdf", height=3, width=6, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0.5))
plot(thetaGrid, frankRhoFun(atanh(thetaGrid) * 40) - rhoTrue, type="l", xlab=expression(alpha), ylab=expression(rho))

plot(thetaGrid, frankdRho(atanh(thetaGrid) * 40) - dRhoTrue, type="l", xlab=expression(alpha), ylab=bquote(rho~"'"))
dev.off()

load("t4Rho.rda")
rhoTrue <- sapply(thetaGrid, function(x) rho(tCopula(x)))
dRhoTrue <-  6 / (pi * sqrt(4 - thetaGrid^2))

pdf("t4-rho.pdf", height=3, width=6, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0.5))
plot(thetaGrid, rhoTrue, type="l", xlab=expression(theta), ylab=expression(rho))
curve(t4RhoFun(x), add=TRUE, col="blue", lty=2)
legend("bottomright", legend=c("normal", "numerical"), col=c("black", "blue"), lty=c(1,2), cex=0.7)

curve(t4dRho(x), col="blue", xlab=expression(theta), ylab=bquote(rho~"'"), lty = 2)
lines(thetaGrid, dRhoTrue)
legend("bottomright", legend=c("normal", "numerical"), col=c("black", "blue"), lty=c(1,2), cex=0.7)
dev.off()

pdf("t4-rho-err.pdf", height=3, width=6, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0.5))
plot(thetaGrid, t4RhoFun(thetaGrid) - rhoTrue, type="l", xlab=expression(theta), ylab=bquote("difference in "~rho), ylim=c(-.016, .016))

plot(thetaGrid, t4dRho(thetaGrid) - dRhoTrue, type="l", xlab=expression(theta), ylab=bquote("difference in "~rho~"'"))
dev.off()
