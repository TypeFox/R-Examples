### R code from vignette source 'FMEmcmc.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("FME")
options(prompt = "> ")
options(width=70)


###################################################
### code chunk number 2: FMEmcmc.Rnw:193-198
###################################################
mu  <- 10
std <- 1

Nfun <- function(p)
  -2*log(dnorm(p, mean = mu, sd = std))


###################################################
### code chunk number 3: FMEmcmc.Rnw:206-207
###################################################
MCMC <- modMCMC (f = Nfun, p = 9.5, niter = 2000, jump = 5)


###################################################
### code chunk number 4: FMEmcmc.Rnw:211-213
###################################################
MCMC <- modMCMC (f = Nfun, p = 9.5, niter = 2000, jump = 5, updatecov = 10)
summary(MCMC)


###################################################
### code chunk number 5: FMEmcmc.Rnw:219-222
###################################################
MCMC2 <- modMCMC (f = Nfun, p = 9.5, lower = 9, niter = 2000, jump = 5,
  updatecov = 10)
summary(MCMC2)


###################################################
### code chunk number 6: FMEmcmc.Rnw:231-235
###################################################
pri <- function(p) -2*log(dnorm(p, 8, 1))
MCMC3 <- modMCMC (f = Nfun, p = 9.5, niter = 2000, jump = 5,
  updatecov = 10, prior = pri)
summary(MCMC3)


###################################################
### code chunk number 7: FMEmcmc.Rnw:240-244
###################################################
summary(MCMC4 <- modMCMC(f = Nfun, p = 1, niter = 2000, jump = 5,
  updatecov = 10, prior = pri, ntrydr = 2))

MCMC4$count


###################################################
### code chunk number 8: hist1
###################################################
par(mfrow = c(2,2))
hist(MCMC$pars, xlab="x", freq = FALSE, main = "unconstrained", xlim = c(6, 14))
hist(MCMC2$pars, xlab="x", freq = FALSE, main = "x>9", xlim = c(6, 14))
hist(MCMC3$pars, xlab="x", freq = FALSE, main = "pri(x)~N(8,1)", xlim = c(6, 14))
plot(MCMC3, mfrow = NULL, main = "AM")
mtext(outer = TRUE, line = -1.5, "N(10,1)", cex = 1.25)


###################################################
### code chunk number 9: hist1
###################################################
par(mfrow = c(2,2))
hist(MCMC$pars, xlab="x", freq = FALSE, main = "unconstrained", xlim = c(6, 14))
hist(MCMC2$pars, xlab="x", freq = FALSE, main = "x>9", xlim = c(6, 14))
hist(MCMC3$pars, xlab="x", freq = FALSE, main = "pri(x)~N(8,1)", xlim = c(6, 14))
plot(MCMC3, mfrow = NULL, main = "AM")
mtext(outer = TRUE, line = -1.5, "N(10,1)", cex = 1.25)


###################################################
### code chunk number 10: FMEmcmc.Rnw:277-283
###################################################
mu  <- 1:4
std <- 1

NL <- function(p)  {
  -2*sum(log(dlnorm(p, mean = mu, sd = std)))
}


###################################################
### code chunk number 11: FMEmcmc.Rnw:289-291
###################################################
MCMCl <- modMCMC (f = NL, p = rep(1, 4), niter = 10000,
  outputlength = 1000, jump = 5)


###################################################
### code chunk number 12: logp1
###################################################
plot(MCMCl)


###################################################
### code chunk number 13: logp1
###################################################
plot(MCMCl)


###################################################
### code chunk number 14: FMEmcmc.Rnw:313-315
###################################################
MCMCl <- modMCMC (f = NL, p = rep(1, 4), niter = 5000, 
   outputlength = 1000, jump = 5, updatecov = 100, ntrydr = 2)


###################################################
### code chunk number 15: logp
###################################################
plot(MCMCl)


###################################################
### code chunk number 16: logp
###################################################
plot(MCMCl)


###################################################
### code chunk number 17: hist
###################################################
hist(MCMCl)


###################################################
### code chunk number 18: hist
###################################################
hist(MCMCl)


###################################################
### code chunk number 19: FMEmcmc.Rnw:345-347
###################################################
MCMCl$pars <- log(MCMCl$pars)
summary(MCMCl)


###################################################
### code chunk number 20: FMEmcmc.Rnw:375-378
###################################################
Banana <- function (x1, x2) {
  return(x2 - (x1^2+1))
}


###################################################
### code chunk number 21: FMEmcmc.Rnw:383-390
###################################################
pmultinorm <- function(vec, mean, Cov) {
  diff <- vec - mean
  ex   <- -0.5*t(diff) %*% solve(Cov) %*% diff
  rdet   <- sqrt(det(Cov))
  power  <- -length(diff)*0.5
  return((2.*pi)^power / rdet * exp(ex))
}


###################################################
### code chunk number 22: FMEmcmc.Rnw:394-400
###################################################
BananaSS <- function (p)
{
  P <- c(p[1], Banana(p[1], p[2]))
  Cov <- matrix(nr = 2, data = c(1, 0.9, 0.9, 1))
 -2*sum(log(pmultinorm(P, mean = 0, Cov = Cov)))
}


###################################################
### code chunk number 23: FMEmcmc.Rnw:412-415
###################################################
MCMC <- modMCMC(f = BananaSS, p = c(0, 0.5), jump = diag(nrow = 2, x = 5),
                niter = 2000)
MCMC$count


###################################################
### code chunk number 24: FMEmcmc.Rnw:422-425
###################################################
MCMC2 <- modMCMC(f = BananaSS, p = c(0, 0.5), jump = diag(nrow = 2, x = 5),
                 updatecov = 100, niter = 2000)
MCMC2$count


###################################################
### code chunk number 25: FMEmcmc.Rnw:432-435
###################################################
MCMC3 <- modMCMC(f = BananaSS, p = c(0, 0.5), jump = diag(nrow = 2, x = 5),
                 ntrydr = 2, niter = 2000)
MCMC3$count


###################################################
### code chunk number 26: FMEmcmc.Rnw:444-449
###################################################
print(system.time(
MCMC4 <- modMCMC(f = BananaSS, p = c(0, 0.5), jump = diag(nrow = 2, x = 5),
                 updatecov = 100, ntrydr = 2, niter = 2000)
))
MCMC4$count


###################################################
### code chunk number 27: banana
###################################################
par(mfrow = c(4, 2))
par(mar = c(2, 2, 4, 2))
plot(MCMC , mfrow = NULL, main = "MH")
plot(MCMC2, mfrow = NULL, main = "AM")
plot(MCMC3, mfrow = NULL, main = "DR")
plot(MCMC4, mfrow = NULL, main = "DRAM")
mtext(outer = TRUE, side = 3, line = -2, at = c(0.05, 0.95),
    c("y1", "y2"), cex = 1.25)
par(mar = c(5.1, 4.1, 4.1, 2.1))


###################################################
### code chunk number 28: banana2
###################################################
par(mfrow = c(2, 2))
xl <- c(-3, 3)
yl <- c(-1, 8)
plot(MCMC$pars,  main = "MH", xlim = xl, ylim = yl)
plot(MCMC2$pars, main = "AM", xlim = xl, ylim = yl)
plot(MCMC3$pars, main = "DR", xlim = xl, ylim = yl)
plot(MCMC4$pars, main = "DRAM", xlim = xl, ylim = yl)


###################################################
### code chunk number 29: banana2
###################################################
par(mfrow = c(2, 2))
xl <- c(-3, 3)
yl <- c(-1, 8)
plot(MCMC$pars,  main = "MH", xlim = xl, ylim = yl)
plot(MCMC2$pars, main = "AM", xlim = xl, ylim = yl)
plot(MCMC3$pars, main = "DR", xlim = xl, ylim = yl)
plot(MCMC4$pars, main = "DRAM", xlim = xl, ylim = yl)


###################################################
### code chunk number 30: FMEmcmc.Rnw:489-493
###################################################
trans <- cbind(MCMC4$pars[ ,1], Banana(MCMC4$pars[ ,1], MCMC4$pars[ ,2]))
colMeans(trans)     # was:c(0,0)
apply(trans, 2, sd) # was:1
cov(trans)          # 0.9 off-diagonal


###################################################
### code chunk number 31: FMEmcmc.Rnw:522-528
###################################################
Reaction <- function (k, times)
{
  fac <- k[1]/(k[1]+k[2])
  A   <- fac + (1-fac)*exp(-(k[1]+k[2])*times)
  return(data.frame(t=times,A=A))
}


###################################################
### code chunk number 32: FMEmcmc.Rnw:532-536
###################################################
Data     <- data.frame(
  times = c(2,     4,     6,     8,     10   ),
  A     = c(0.661, 0.668, 0.663, 0.682, 0.650))
Data


###################################################
### code chunk number 33: FMEmcmc.Rnw:544-546
###################################################
Prior <- function(p)
    return( sum(((p - c(2, 4))/200)^2 ))


###################################################
### code chunk number 34: FMEmcmc.Rnw:551-556
###################################################
residual <- function(k) return(Data$A - Reaction(k,Data$times)$A)

Fit <- modFit(p = c(k1 = 0.5, k2 = 0.5), f = residual, 
              lower = c(0, 0), upper = c(1, 1))
(sF <- summary(Fit))


###################################################
### code chunk number 35: FMEmcmc.Rnw:573-580
###################################################
mse <- sF$modVariance
Cov <- sF$cov.scaled * 2.4^2/2

print(system.time(
MCMC <- modMCMC(f = residual, p = Fit$par, jump = Cov, lower = c(0, 0),
                var0 = mse, wvar0 = 1, prior = Prior, niter = 2000)
))


###################################################
### code chunk number 36: ABMCMC
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 37: ABMCMC
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 38: FMEmcmc.Rnw:602-604
###################################################
MCMC2<- modMCMC(f = residual, p = Fit$par, jump = Cov, updatecov = 100, 
         lower = c(0, 0), var0 = mse, wvar0 = 1, prior = Prior, niter = 2000) 


###################################################
### code chunk number 39: ABMCMC2
###################################################
plot(MCMC2, Full = TRUE)


###################################################
### code chunk number 40: ABMCMC2
###################################################
plot(MCMC2, Full = TRUE)


###################################################
### code chunk number 41: ABMCMC3
###################################################
pairs(MCMC2)


###################################################
### code chunk number 42: ABMCMC3
###################################################
pairs(MCMC2)


###################################################
### code chunk number 43: FMEmcmc.Rnw:640-641
###################################################
sR <- sensRange(f=Reaction,times=seq(0,10,0.1),parInput=MCMC2$pars)


###################################################
### code chunk number 44: sr
###################################################
plot(summary(sR), xlab = "time", ylab = "Conc")
points(Data)


###################################################
### code chunk number 45: sr
###################################################
plot(summary(sR), xlab = "time", ylab = "Conc")
points(Data)


###################################################
### code chunk number 46: FMEmcmc.Rnw:677-683
###################################################
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour

Obs2<- data.frame(x=c(   20,  55,   83,  110,  138,  240,  325),   # mg COD/l
                   y=c(0.05,0.07,0.09,0.10,0.11,0.122,0.125))   # 1/hour



###################################################
### code chunk number 47: FMEmcmc.Rnw:688-689
###################################################
Model <- function(p,x) return(data.frame(x = x, y = p[1]*x/(x+p[2])))


###################################################
### code chunk number 48: FMEmcmc.Rnw:697-701
###################################################
Residuals  <- function(p) {
   cost <- modCost(model = Model(p, Obs$x), obs = Obs, x = "x")
   modCost(model = Model(p, Obs2$x), obs = Obs2, cost = cost, x = "x")
}


###################################################
### code chunk number 49: FMEmcmc.Rnw:705-708
###################################################
print(system.time(
P      <- modFit(f = Residuals, p = c(0.1, 1))
))


###################################################
### code chunk number 50: obs
###################################################
plot(Obs, xlab = "mg COD/l", ylab = "1/hour", pch = 16, cex = 1.5)
points(Obs2, pch = 18, cex = 1.5, col = "red")
lines(Model(p = P$par, x = 0:375))


###################################################
### code chunk number 51: obs
###################################################
plot(Obs, xlab = "mg COD/l", ylab = "1/hour", pch = 16, cex = 1.5)
points(Obs2, pch = 18, cex = 1.5, col = "red")
lines(Model(p = P$par, x = 0:375))


###################################################
### code chunk number 52: FMEmcmc.Rnw:732-733
###################################################
Covar   <- summary(P)$cov.scaled * 2.4^2/2


###################################################
### code chunk number 53: FMEmcmc.Rnw:748-754
###################################################
s2prior <- P$ms

print(system.time(
MCMC <- modMCMC(f = Residuals, p = P$par, jump = Covar, niter = 1000,
           var0 = s2prior, wvar0 = 0.1, lower = c(0, 0))
))


###################################################
### code chunk number 54: Monmcmc
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 55: Monmcmc
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 56: FMEmcmc.Rnw:777-783
###################################################
varprior <- P$var_ms_unweighted

print(system.time(
MCMC2 <- modMCMC(f = Residuals, p = P$par, jump = Covar, niter = 1000,
                var0 = varprior, wvar0 = 0.1, lower = c(0, 0))
))


###################################################
### code chunk number 57: Monmcmc2
###################################################
plot(MCMC2, Full = TRUE, which = NULL)


###################################################
### code chunk number 58: Monmcmc2
###################################################
plot(MCMC2, Full = TRUE, which = NULL)


###################################################
### code chunk number 59: FMEmcmc.Rnw:802-804
###################################################
summary(MCMC)
summary(MCMC2)


