### R code from vignette source 'mixtools.Rnw'

###################################################
### code chunk number 1: faithful
###################################################
library(mixtools)
data(faithful)
attach(faithful)


###################################################
### code chunk number 2: geyser
###################################################
hist(waiting, main="Time between Old Faithful eruptions",
     xlab="Minutes",  ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.4)


###################################################
### code chunk number 3: normmixEM
###################################################
wait1 <- normalmixEM(waiting, lambda = .5, mu = c(55, 80), sigma = 5)


###################################################
### code chunk number 4: geyserEM (eval = FALSE)
###################################################
## plot(wait1, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
##      main2="Time between Old Faithful eruptions", xlab2="Minutes")


###################################################
### code chunk number 5: geyserEM
###################################################
for(i in 1:2){
  file=paste("geyserEM", i, ".pdf", sep="")
  pdf(file=file, paper="special", width=6, height=6)
  plot(wait1, whichplots=i, cex.axis = 1.4, cex.lab = 1.4, cex.main =
       1.8,   main2 = "Time between Old Faithful eruptions", xlab2 =
       "Minutes")
  dev.off()
  cat("\\includegraphics{", file, "}\n", sep="")
}


###################################################
### code chunk number 6: geyserestimates
###################################################
wait1[c("lambda", "mu", "sigma")]


###################################################
### code chunk number 7: geysersummary
###################################################
summary(wait1)


###################################################
### code chunk number 8: cutpoint
###################################################
data("Waterdata")
cutpts <- 10.5*(-6:6)
watermult <- makemultdata(Waterdata, cuts = cutpts)


###################################################
### code chunk number 9: multmixEM
###################################################
set.seed(15)
theta4 <- matrix(runif(56), ncol = 14)
theta3 <- theta4[1:3,]
mult3 <- multmixEM(watermult, lambda = rep(1, 3)/3, theta = theta3)
mult4 <- multmixEM (watermult, lambda = rep (1, 4) / 4, theta = theta4)


###################################################
### code chunk number 10: mixtools.Rnw:575-581 (eval = FALSE)
###################################################
## cdf3 <- compCDF(Waterdata, mult3$posterior, lwd=2, lab=c(7, 5, 7),
##                 xlab="Angle in degrees", ylab="Component CDFs",
##                 main="Three-Component Solution")
## cdf4 <- compCDF(Waterdata, mult4$posterior, lwd=2, lab=c(7, 5, 7),
##                 xlab="Angle in degrees", ylab="Component CDFs",
##                 main="Four-Component Solution")


###################################################
### code chunk number 11: cutpointplots
###################################################
pdf(file="WDcutpoint3comp.pdf", paper="special", width=8, height=8)
cdf3 <- compCDF(Waterdata, mult3$posterior, lwd=3,
                xlab="Angle in degrees", lab=c(7, 5, 7), ylab="Component CDFs",
                main="Three-Component Solution", cex.axis=1.4, cex.lab=1.5,
                cex.main=1.5)
ltext <- paste(round(mult3$lam*100, 1), "%", sep="")
legend("bottomright", legend=ltext, pch=15:17, cex=1.5, pt.cex=1.35)
y <- compCDF(Waterdata, mult3$posterior, x=cutpts, makeplot=F)
for(i in 1:3) points(cutpts, y[i,], pch=14+i, cex=1.35)
dev.off()
pdf(file="WDcutpoint4comp.pdf", paper="special", width=8, height=8)
cdf4 <- compCDF(Waterdata, mult4$posterior, lwd=3,
                xlab="Angle in degrees", lab=c(7, 5, 7),
                ylab="Component CDFs", main="Four-Component Solution",
                cex.axis=1.4, cex.lab=1.5, cex.main=1.5)
ltext <- paste(round(mult4$lam*100,1), "%", sep="")
legend("bottomright", legend=ltext, pch=15:18, cex=1.5, pt.cex=1.35)
y <- compCDF(Waterdata, mult4$posterior, x=cutpts, makeplot=F)
for(i in 1:4) points(cutpts, y[i,], pch=14+i, cex=1.35)
dev.off()


###################################################
### code chunk number 12: summarymult4
###################################################
summary(mult4)


###################################################
### code chunk number 13: spsymmplots
###################################################
pdf(file="spsymmfig1.pdf", paper="special", width=8, height=8)
par(mar=0.1+c(5,4.2,4,1.8))
plot(wait1, which = 2, cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5,
     main2 = "Time between Old Faithful eruptions", xlab2 = "Minutes")
wait2 <- spEMsymloc(waiting, mu0 = c(55, 80))
plot(wait2, lty = 2, newplot = FALSE, addlegend = FALSE)
dev.off()
pdf(file="spsymmfig2.pdf", paper="special", width=8, height=8)
par(mar=0.1+c(5,4.2,4,1.8))
wait2a <- spEMsymloc(waiting, mu0 = c(55, 80), bw = 1)
wait2b <- spEMsymloc(waiting, mu0 = c(55, 80), bw = 6)
plot(wait2a, lty = 1, addlegend = FALSE, cex.axis = 1.4, cex.lab = 1.5,
     cex.main = 1.5, title = "Time between Old Faithful eruptions",
     xlab = "Minutes")
plot(wait2b, lty = 2, newplot = FALSE, addlegend = FALSE)
dev.off()


###################################################
### code chunk number 14: plotspsymm (eval = FALSE)
###################################################
## plot(wait1, which = 2, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8,
##      main2 = "Time between Old Faithful eruptions", xlab2 = "Minutes")
## wait2 <- spEMsymloc(waiting, mu0 = c(55, 80))
## plot(wait2, lty = 2, newplot = FALSE, addlegend = FALSE)


###################################################
### code chunk number 15: bandwidth
###################################################
bw.nrd0(waiting)


###################################################
### code chunk number 16: plotbweffect (eval = FALSE)
###################################################
## wait2a <- spEMsymloc(waiting, mu0 = c(55, 80), bw = 1)
## wait2b <- spEMsymloc(waiting, mu0 = c(55, 80), bw = 6)
## plot(wait2a, lty = 1, addlegend = FALSE, cex.axis = 1.4,
##      cex.lab = 1.4, cex.main = 1.8, xlab = "Minutes",
##      title = "Time between Old Faithful eruptions")
## plot(wait2b, lty = 2, newplot = FALSE, addlegend = FALSE)


###################################################
### code chunk number 17: gaussexample
###################################################
m <- 2; r <- 3; n <- 300; S <- 100
lambda <- c(0.4, 0.6)
mu <- matrix(c(0, 0, 0, 3, 4, 5), m, r, byrow = TRUE)
sigma <- matrix(rep(1, 6), m, r, byrow = TRUE)


###################################################
### code chunk number 18: gaussinitial
###################################################
centers <- matrix(c(0, 0, 0, 4, 4, 4), 2, 3, byrow = TRUE)
ISE <- matrix(0, m, r, dimnames = list(Components = 1:m, Blocks = 1:r))
nblabsw <- 0


###################################################
### code chunk number 19: sqMISE
###################################################
set.seed(1000)
for (mc in 1:S) {
  x <- rmvnormmix(n, lambda, mu, sigma)
  a <- npEM(x, centers, verb = FALSE, samebw = FALSE)
  if (a$lambda[1] > a$lambda[2]) nblabsw <- nblabsw + 1
  for (j in 1:m) {
     for (k in 1:r) {
     ISE[j, k] <- ISE[j, k] + ise.npEM(a, j, k, dnorm,
         lower = mu[j, k] - 5, upper = mu[j, k] + 5, plots = FALSE,
         mean = mu[j, k], sd = sigma[j, k])$value #$
    }
  }
}
MISE <- ISE/S
print(sqMISE <- sqrt(MISE))


###################################################
### code chunk number 20: summarygauss
###################################################
summary(a)


###################################################
### code chunk number 21: plotgauss3rm (eval = FALSE)
###################################################
## plot(a)


###################################################
### code chunk number 22: gauss3rm
###################################################
pdf("gauss3rm.pdf", paper="special", width=10, height=5)
par(mfrow=c(1,3), ask=F)
plot(a)
dev.off()


###################################################
### code chunk number 23: true5rm
###################################################
pdf("truepdf5rm_block1.pdf")
par(mar=0.1+c(5,4.2,4,1.5))
x <- seq(-10, 25, len=250)
plot(x, .4* dt(x, 2, 0) + .6 * dt(x, 10, 8), type="l", lwd=3, col=2,
     cex.axis=1.4, cex.lab=1.5, cex.main=1.5, main="Block 1", xlab="",
     ylab="Density")
lines (x, .4*dt(x, 2, 0), lwd=4, lty=2)
lines (x, .6*dt(x, 10, 8), lwd=4, lty=2)
dev.off()
pdf("truepdf5rm_block2.pdf")
par(mar=0.1+c(5,4.2,4,1.5))
x <- seq(0, 1, len=250)
plot(x, .4 + .6 * dbeta(x, 1, 5), type="l", lwd=3, col=2, cex.axis=1.4,
   cex.lab=1.5, cex.main=1.5, main="Block 2", xlab="", ylab="Density",
   ylim= c(0, 3.4))
lines (x, rep(.4, 250), lwd=4, lty=2)
lines (x, .6*dbeta(x, 1, 5), lwd=4, lty=2)
dev.off()


###################################################
### code chunk number 24: parameters5rm
###################################################
m <- 2; r <- 5
lambda <- c(0.4, 0.6)
df <- c(2, 10); ncp <- c(0, 8)
sh1 <- c(1, 1) ; sh2 <- c(1, 5)


###################################################
### code chunk number 25: generate5rm
###################################################
n <- 300; z <- sample(m, n, rep = TRUE, prob = lambda)
r1 <- 3; z2 <- rep(z, r1)
x1 <- matrix(rt(n * r1, df[z2], ncp[z2]), n, r1)
r2 <- 2; z2 <- rep(z, r2)
x2 <- matrix(rbeta(n * r2, sh1[z2], sh2[z2]), n, r2)
x <- cbind(x1, x2)


###################################################
### code chunk number 26: npEM5rm
###################################################
id <- c(rep(1, r1), rep(2, r2))
centers <- matrix(c(0, 0, 0, 1/2, 1/2, 4, 4, 4, 1/2, 1/2), m, r,
                  byrow = TRUE)
b <- npEM(x, centers, id, eps = 1e-8, verb = FALSE, samebw = FALSE)


###################################################
### code chunk number 27: plot5rm (eval = FALSE)
###################################################
## plot(b, breaks = 15)


###################################################
### code chunk number 28: plot5rmcommands
###################################################
pdf("npEM5rm.pdf", width=8, height=5)
par(mfrow=c(1,2))
plot(b, breaks = 15)

dev.off()


###################################################
### code chunk number 29: ISEnpEM5rm (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## for (j in 1:2){
##   ise.npEM(b, j, 1, truepdf = dt, lower = ncp[j] - 10,
##            upper = ncp[j] + 10, df = df[j], ncp = ncp[j])
##   ise.npEM(b, j, 2, truepdf = dbeta, lower = -0.5,
##            upper = 1.5,  shape1 = sh1[j], shape2 = sh2[j])
## }


###################################################
### code chunk number 30: plotISEnpEM5rm
###################################################
options(warn=-1)
pdf("ISEnpEM5rm.pdf", width=8, height=8)
par(mfrow = c(2, 2))
for (j in 1:2){
  ise.npEM(b, j, 1, truepdf = dt, lower = ncp[j] - 10,
           upper = ncp[j] + 10, df = df[j], ncp = ncp[j])
  ise.npEM(b, j, 2, truepdf = dbeta, lower = -0.5,
           upper = 1.5,  shape1 = sh1[j], shape2 = sh2[j])
}
dev.off()


###################################################
### code chunk number 31: gnpdata
###################################################
data("CO2data")
attach(CO2data)
pdf("gnpdata.pdf")
par(mar=0.1+c(5,4.2,4,1.5))
plot(GNP, CO2, xlab="Gross National Product",
     ylab=expression(paste(CO[2]," per Capita")),
     cex.lab=1.5, cex.main=1.5, cex.axis=1.4,
     main="1996 GNP and Emissions Data")
text(GNP, CO2, country, adj=c(.5,-.5))
dev.off()


###################################################
### code chunk number 32: mixtools.Rnw:1282-1284
###################################################
data("CO2data")
attach(CO2data)


###################################################
### code chunk number 33: CO2reg
###################################################
CO2reg <- regmixEM(CO2, GNP, lambda = c(1, 3) / 4,
                   beta = matrix(c(8, -1, 1, 1), 2, 2), sigma = c(2, 1))


###################################################
### code chunk number 34: summaryCO2reg
###################################################
summary(CO2reg)


###################################################
### code chunk number 35: plotCO2reg (eval = FALSE)
###################################################
## plot(CO2reg, density = TRUE, alpha = 0.01, cex.main = 1.5, cex.lab = 1.5,
##      cex.axis = 1.4)


###################################################
### code chunk number 36: trueplotCO2reg
###################################################
for(i in 1:2){
    file=paste("CO2reg", i, ".pdf", sep="")
    pdf(file=file, paper="special", width=6, height=6)
    plot(CO2reg, whichplots=i, alpha = 0.01, cex.main = 1.5, cex.lab = 1.5,
         cex.axis = 1.4)
    dev.off()
    cat("\\includegraphics{", file, "}\n", sep="")
}


###################################################
### code chunk number 37: CO2igle
###################################################
CO2igle <- regmixEM.loc(CO2, GNP, beta = CO2reg$beta, sigma = CO2reg$sigma,
                        lambda = CO2reg$posterior, kern.l = "Beta",
                        kernl.h = 20, kernl.g = 3)


###################################################
### code chunk number 38: CO2iglesummary
###################################################
summary(CO2igle)


###################################################
### code chunk number 39: lamplot (eval = FALSE)
###################################################
## plot(GNP, CO2igle$post[,1], xlab = "GNP", cex.axis = 1.4, cex.lab = 1.5,
##   ylab = "Final posterior probabilities")
## lines(sort(GNP), CO2igle$lambda[order(GNP), 1], col=2)
## abline(h = CO2igle$lambda[1], lty = 2)


###################################################
### code chunk number 40: truelamplot
###################################################
pdf("lamplot.pdf")
plot(GNP, CO2igle$post[,1], xlab = "GNP", cex.axis = 1.4, cex.lab = 1.5,
  ylab = "Final posterior probabilities")
lines(sort(GNP), CO2igle$lambda[order(GNP), 1], col=2, lwd=2)
abline(h = CO2igle$lambda[1], lty = 2, lwd=2)
dev.off()


###################################################
### code chunk number 41: CO2boot
###################################################
set.seed(123)
CO2boot <- boot.se(CO2reg, B = 100)


###################################################
### code chunk number 42: bootresults
###################################################
rbind(range(CO2boot$beta[1,]), range(CO2boot$beta[2,]))


###################################################
### code chunk number 43: CO2bootse
###################################################
CO2boot[c("lambda.se", "beta.se", "sigma.se")]


###################################################
### code chunk number 44: modelsel
###################################################
data("Waterdata")
cutpts <- 10.5*(-6:6)
watermult <- makemultdata(Waterdata, cuts = cutpts)
set.seed(10)
multmixmodel.sel(watermult, comps = 1:4, epsilon = 0.001)


