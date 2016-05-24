### R code from vignette source 'FMEdyna.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("FME")
options(prompt = "> ")
options(width = 80)


###################################################
### code chunk number 2: FMEdyna.Rnw:137-139
###################################################
pars <- list(gmax = 0.5, eff = 0.5,
              ks = 0.5, rB = 0.01, dB = 0.01)


###################################################
### code chunk number 3: FMEdyna.Rnw:154-167
###################################################
solveBact <- function(pars, times=seq(0,50,by=0.5)) {
  derivs <- function(t, state, pars) { # returns rate of change
    with(as.list(c(state, pars)), {

      dBact <-  gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
      dSub  <- -gmax    *Sub/(Sub+ks)*Bact + dB*Bact
      return(list(c(dBact, dSub), TOC = Bact + Sub))
    })
 }
 state   <- c(Bact = 0.1, Sub = 100)
 ## ode solves the model by integration...
 return(ode(y = state, times = times, func = derivs, parms = pars))
}


###################################################
### code chunk number 4: FMEdyna.Rnw:173-174
###################################################
out <- solveBact(pars)


###################################################
### code chunk number 5: ode
###################################################
matplot(out[,1], out[,-1], type = "l", lty = 1:3, lwd = c(2, 2, 1),
   col = "black", xlab = "time, hour", ylab = "mol C/m3")

legend("topright", c("Bacteria", "Glucose", "TOC"),
       lty = 1:3, lwd = c(2, 2, 1))


###################################################
### code chunk number 6: odefig
###################################################
matplot(out[,1], out[,-1], type = "l", lty = 1:3, lwd = c(2, 2, 1),
   col = "black", xlab = "time, hour", ylab = "mol C/m3")

legend("topright", c("Bacteria", "Glucose", "TOC"),
       lty = 1:3, lwd = c(2, 2, 1))


###################################################
### code chunk number 7: FMEdyna.Rnw:208-211
###################################################
parRanges <- data.frame(min = c(0.4, 0.4, 0.0), max = c(0.6, 0.6, 0.02))
rownames(parRanges) <- c("gmax", "eff", "rB")
parRanges


###################################################
### code chunk number 8: FMEdyna.Rnw:221-227
###################################################
tout    <- 0:50
print(system.time(
sR <- sensRange(func = solveBact, parms = pars, dist = "grid",
       sensvar = c("Bact", "Sub"), parRange = parRanges[3,], num = 50)
))
head(summary(sR))


###################################################
### code chunk number 9: sens
###################################################
summ.sR <- summary(sR)
par(mfrow=c(2, 2))
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3",
    legpos = "topright", mfrow = NULL)
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3", mfrow = NULL,
     quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, "Sensitivity to rB", cex = 1.25)
par(mfrow = c(1, 1))


###################################################
### code chunk number 10: sensfig
###################################################
summ.sR <- summary(sR)
par(mfrow=c(2, 2))
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3",
    legpos = "topright", mfrow = NULL)
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3", mfrow = NULL,
     quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, "Sensitivity to rB", cex = 1.25)
par(mfrow = c(1, 1))


###################################################
### code chunk number 11: FMEdyna.Rnw:261-263
###################################################
Sens2 <- summary(sensRange(func = solveBact, parms = pars,
   dist = "latin", sensvar = "Bact", parRange = parRanges, num = 100))


###################################################
### code chunk number 12: sens2
###################################################
plot(Sens2, main = "Sensitivity gmax,eff,rB", xlab = "time, hour",
   ylab = "molC/m3")


###################################################
### code chunk number 13: sensfig2
###################################################
plot(Sens2, main = "Sensitivity gmax,eff,rB", xlab = "time, hour",
   ylab = "molC/m3")


###################################################
### code chunk number 14: FMEdyna.Rnw:296-299
###################################################
SnsBact<- sensFun(func = solveBact, parms = pars,
                 sensvar = "Bact", varscale = 1)
head(SnsBact)


###################################################
### code chunk number 15: sfun
###################################################
plot(SnsBact)


###################################################
### code chunk number 16: sfunfig
###################################################
plot(SnsBact)


###################################################
### code chunk number 17: FMEdyna.Rnw:324-325
###################################################
summary(SnsBact)


###################################################
### code chunk number 18: FMEdyna.Rnw:338-339
###################################################
summary(sensFun(solveBact, pars, varscale = 1), var = TRUE)


###################################################
### code chunk number 19: FMEdyna.Rnw:348-349
###################################################
cor(SnsBact[ ,-(1:2)])


###################################################
### code chunk number 20: pairs
###################################################
pairs(SnsBact)


###################################################
### code chunk number 21: pairsfig
###################################################
pairs(SnsBact)


###################################################
### code chunk number 22: FMEdyna.Rnw:377-382
###################################################
SF <- function (pars) {
  out <- solveBact(pars)
  return(out[nrow(out), 2:3])
}
CRL <- modCRL(func = SF, parms = pars, parRange = parRanges[1,])


###################################################
### code chunk number 23: crl
###################################################
plot(CRL)


###################################################
### code chunk number 24: crlfig
###################################################
plot(CRL)


###################################################
### code chunk number 25: FMEdyna.Rnw:409-412
###################################################
CRL2 <- modCRL(func = SF, parms = pars, parMean = c(gmax = 0.5, eff = 0.7),
               parCovar = matrix(nr = 2, data = c(0.02, 0.02, 0.02, 0.05)),
               dist = "norm", sensvar = "Bact", num = 150)


###################################################
### code chunk number 26: crl2
###################################################
pairs(CRL2)


###################################################
### code chunk number 27: crl2fig
###################################################
pairs(CRL2)


###################################################
### code chunk number 28: FMEdyna.Rnw:435-439
###################################################
Coll <- collin(SnsBact)
Coll
Coll [Coll[,"collinearity"] < 20 & Coll[ ,"N"] == 4, ]
collin(SnsBact, parset = 1:5)


###################################################
### code chunk number 29: FMEdyna.Rnw:520-525
###################################################
Dat<- data.frame(name = c("Obs1", "Obs1", "Obs2", "Obs2"),
             time = c(1, 2, 1, 2), val = c(50, 150, 1, 2),
             err = c(5, 15, 0.1, 0.2))
Mod <- data.frame(time = 0:3, Obs1 = rep(4, 4), Obs2 = 1:4)
modCost(mod = Mod, obs = Dat, y = "val")


###################################################
### code chunk number 30: FMEdyna.Rnw:530-531
###################################################
modCost(mod = Mod, obs = Dat, y = "val", err = "err")


###################################################
### code chunk number 31: FMEdyna.Rnw:540-548
###################################################
Data <- matrix (nc=2,byrow=2,data=
c(  2,  0.14,    4,  0.21,    6,  0.31,    8,  0.40,
   10,  0.69,   12,  0.97,   14,  1.42,   16,  2.0,
   18,  3.0,    20,  4.5,    22,  6.5,    24,  9.5,
   26, 13.5,    28, 20.5,    30,  29 , 35, 65, 40, 61)
)
colnames(Data) <- c("time", "Bact")
head(Data)


###################################################
### code chunk number 32: FMEdyna.Rnw:559-567
###################################################
Objective <- function(x, parset = names(x)) {
  pars[parset] <- x
  tout    <- seq(0, 50, by = 0.5)
  ## output times
  out <- solveBact(pars, tout)
  ## Model cost
  return(modCost(obs = Data, model = out))
}


###################################################
### code chunk number 33: FMEdyna.Rnw:575-577
###################################################
Coll <- collin(sF <- sensFun(func = Objective, parms = pars, varscale = 1))
Coll


###################################################
### code chunk number 34: coll
###################################################
plot(Coll, log = "y")
abline(h = 20, col = "red")


###################################################
### code chunk number 35: collfig
###################################################
plot(Coll, log = "y")
abline(h = 20, col = "red")


###################################################
### code chunk number 36: FMEdyna.Rnw:606-607
###################################################
collin(sF,parset=1:2)


###################################################
### code chunk number 37: FMEdyna.Rnw:616-619
###################################################
print(system.time(Fit <- modFit(p = c(gmax = 0.5, eff = 0.5),
                  f = Objective, lower = c(0.0, 0.0))))
summary(Fit)


###################################################
### code chunk number 38: FMEdyna.Rnw:625-632
###################################################
init <- solveBact(pars)

pars[c("gmax", "eff")] <- Fit$par
out   <- solveBact(pars)

Cost  <- modCost(obs = Data, model = out)
Cost


###################################################
### code chunk number 39: fit
###################################################
plot(out, init, xlab = "time, hour", ylab = "molC/m3", lwd = 2, 
   obs = Data, obspar = list(cex = 2, pch = 18)) 
legend ("bottomright", lwd = 2, col = 1:2, lty = 1:2, c("fitted", "original"))


###################################################
### code chunk number 40: fitfig
###################################################
plot(out, init, xlab = "time, hour", ylab = "molC/m3", lwd = 2, 
   obs = Data, obspar = list(cex = 2, pch = 18)) 
legend ("bottomright", lwd = 2, col = 1:2, lty = 1:2, c("fitted", "original"))


###################################################
### code chunk number 41: res
###################################################
plot(Cost, xlab = "time", ylab = "", main = "residuals")


###################################################
### code chunk number 42: resfig
###################################################
plot(Cost, xlab = "time", ylab = "", main = "residuals")


###################################################
### code chunk number 43: FMEdyna.Rnw:696-703
###################################################
SF<-summary(Fit)
SF
SF[]
Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/2
MCMC <- modMCMC(p = coef(Fit), f = Objective, jump = covIni,
              var0 = Var0, wvar0 = 1)


###################################################
### code chunk number 44: mcmcplot
###################################################
 plot(MCMC, Full = TRUE)


###################################################
### code chunk number 45: mcmcfig
###################################################
 plot(MCMC, Full = TRUE)


###################################################
### code chunk number 46: mcmcplot2
###################################################
pairs(MCMC)


###################################################
### code chunk number 47: mcmcfig2
###################################################
pairs(MCMC)


###################################################
### code chunk number 48: FMEdyna.Rnw:745-746
###################################################
MC <- as.mcmc(MCMC$pars)


###################################################
### code chunk number 49: cumuplot
###################################################
cumuplot(MC)


###################################################
### code chunk number 50: cumuplot
###################################################
cumuplot(MC)


###################################################
### code chunk number 51: FMEdyna.Rnw:768-770
###################################################
cov(MCMC$pars)
covIni


###################################################
### code chunk number 52: dist
###################################################
par(mfrow = c(2, 2))
Minmax <- data.frame(min = c(1, 2), max = c(2, 3))
rownames(Minmax) <- c("par1", "par2")
Mean   <- c(par1 = 1.5, par2 = 2.5)
Covar  <- matrix(nr = 2, data = c(2, 2, 2, 3))
plot(Unif(Minmax, 100), main = "Unif", xlim = c(1, 2), ylim = c(2, 3))
plot(Grid(Minmax, 100), main = "Grid", xlim = c(1, 2), ylim = c(2, 3))
plot(Latinhyper(Minmax, 5), main = "Latin hypercube", xlim = c(1, 2),
     ylim = c(2, 3))
grid()
plot(Norm(parMean = Mean, parCovar = Covar, num = 1000),
   main = "multi normal")


###################################################
### code chunk number 53: distfig
###################################################
par(mfrow = c(2, 2))
Minmax <- data.frame(min = c(1, 2), max = c(2, 3))
rownames(Minmax) <- c("par1", "par2")
Mean   <- c(par1 = 1.5, par2 = 2.5)
Covar  <- matrix(nr = 2, data = c(2, 2, 2, 3))
plot(Unif(Minmax, 100), main = "Unif", xlim = c(1, 2), ylim = c(2, 3))
plot(Grid(Minmax, 100), main = "Grid", xlim = c(1, 2), ylim = c(2, 3))
plot(Latinhyper(Minmax, 5), main = "Latin hypercube", xlim = c(1, 2),
     ylim = c(2, 3))
grid()
plot(Norm(parMean = Mean, parCovar = Covar, num = 1000),
   main = "multi normal")


