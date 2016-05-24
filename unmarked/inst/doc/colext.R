### R code from vignette source 'colext.Rnw'

###################################################
### code chunk number 1: colext.Rnw:1-3
###################################################
options(width=70)
options(continue=" ")


###################################################
### code chunk number 2: colext.Rnw:377-416
###################################################
M <- 250                                # Number of sites
J <- 3                                  # num secondary sample periods
T <- 10                                 # num primary sample periods

psi <- rep(NA, T)                       # Occupancy probability
muZ <- z <- array(dim = c(M, T))        # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T))        # Detection histories

set.seed(13973)
psi[1] <- 0.4                           # Initial occupancy probability
p <- c(0.3,0.4,0.5,0.5,0.1,0.3,0.5,0.5,0.6,0.2)
phi <- runif(n=T-1, min=0.6, max=0.8)   # Survival probability (1-epsilon)
gamma <- runif(n=T-1, min=0.1, max=0.2) # Colonization probability

# Generate latent states of occurrence
# First year
z[,1] <- rbinom(M, 1, psi[1])           # Initial occupancy state
# Later years
for(i in 1:M){                          # Loop over sites
   for(k in 2:T){                        # Loop over years
      muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1]
      z[i,k] <- rbinom(1, 1, muZ[k])
   }
}

# Generate detection/non-detection data
for(i in 1:M){
   for(k in 1:T){
      prob <- z[i,k] * p[k]
      for(j in 1:J){
         y[i,j,k] <- rbinom(1, 1, prob)
      }
   }
}

# Compute annual population occupancy
for (k in 2:T){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   }


###################################################
### code chunk number 3: sim
###################################################
plot(1:T, colMeans(z), type = "b", xlab = "Year",
     ylab = "Proportion of sites occupied",
     col = "black", xlim=c(0.5, 10.5), xaxp=c(1,10,9),
     ylim = c(0,0.6), lwd = 2, lty = 1,
     frame.plot = FALSE, las = 1, pch=16)

psi.app <- colMeans(apply(y, c(1,3), max))
lines(1:T, psi.app, type = "b", col = "blue", lty=3, lwd = 2)
legend(1, 0.6, c("truth", "observed"),
       col=c("black", "blue"), lty=c(1,3), pch=c(16,1))


###################################################
### code chunk number 4: colext.Rnw:457-458
###################################################
library(unmarked)


###################################################
### code chunk number 5: colext.Rnw:467-468
###################################################
yy <- matrix(y, M, J*T)


###################################################
### code chunk number 6: colext.Rnw:476-478
###################################################
year <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
               nrow(yy), T, byrow=TRUE)


###################################################
### code chunk number 7: colext.Rnw:493-498
###################################################
simUMF <- unmarkedMultFrame(
    y = yy,
    yearlySiteCovs = list(year = year),
    numPrimary=T)
summary(simUMF)


###################################################
### code chunk number 8: colext.Rnw:577-578
###################################################
plogis(-0.813)


###################################################
### code chunk number 9: colext.Rnw:647-652 (eval = FALSE)
###################################################
## m1 <- colext(psiformula = ~1,   # First-year occupancy
##     gammaformula = ~ year-1,    # Colonization
##     epsilonformula = ~ year-1,  # Extinction
##     pformula = ~ year-1,        # Detection
##     data = simUMF)


###################################################
### code chunk number 10: colext.Rnw:654-655 (eval = FALSE)
###################################################
## m1


###################################################
### code chunk number 11: colext.Rnw:736-741 (eval = FALSE)
###################################################
## nd <- data.frame(year=c('01','02','03','04','05','06','07','08','09'))
## E.ext <- predict(m1, type='ext', newdata=nd)
## E.col <- predict(m1, type='col', newdata=nd)
## nd <- data.frame(year=c('01','02','03','04','05','06','07','08','09','10'))
## E.det <- predict(m1, type='det', newdata=nd)


###################################################
### code chunk number 12: yearlysim (eval = FALSE)
###################################################
## op <- par(mfrow=c(3,1), mai=c(0.6, 0.6, 0.1, 0.1))
## 
## with(E.ext, {   # Plot for extinction probability
##   plot(1:9, Predicted, pch=1, xaxt='n', xlab='Year',
##     ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
##     ylim=c(0,1), col=4)
##   axis(1, at=1:9, labels=nd$year[1:9])
##   arrows(1:9, lower, 1:9, upper, code=3, angle=90, length=0.03, col=4)
##   points((1:9)-0.1, 1-phi, col=1, lwd = 1, pch=16)
##   legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
##          cex=0.8)
##   })
## 
## with(E.col, {	# Plot for colonization probability
##   plot(1:9, Predicted, pch=1, xaxt='n', xlab='Year',
##     ylab=expression(paste('Colonization probability ( ', gamma, ' )')),
##     ylim=c(0,1), col=4)
##   axis(1, at=1:9, labels=nd$year[1:9])
##   arrows(1:9, lower, 1:9, upper, code=3, angle=90, length=0.03, col=4)
##   points((1:9)-0.1, gamma, col=1, lwd = 1, pch=16)
##   legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
##          cex=0.8)
##   })
## 
## with(E.det, {   # Plot for detection probability: note 10 years
##   plot(1:10, Predicted, pch=1, xaxt='n', xlab='Year',
##     ylab=expression(paste('Detection probability ( ', p, ' )')),
##     ylim=c(0,1), col=4)
##   axis(1, at=1:10, labels=nd$year)
##   arrows(1:10, lower, 1:10, upper, code=3, angle=90, length=0.03, col=4)
##   points((1:10)-0.1, p, col=1, lwd = 1, pch=16)
##   legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
##          cex=0.8)
##   })
## 
## par(op)


###################################################
### code chunk number 13: colext.Rnw:881-898
###################################################
turnover <- function(fm) {
    psi.hat <- plogis(coef(fm, type="psi"))
    if(length(psi.hat) > 1)
        stop("this function only works if psi is scalar")
    T <- getData(fm)@numPrimary
    tau.hat <- numeric(T-1)
    gamma.hat <- plogis(coef(fm, type="col"))
    phi.hat <- 1 - plogis(coef(fm, type="ext"))
    if(length(gamma.hat) != T-1 | length(phi.hat) != T-1)
        stop("this function only works if gamma and phi T-1 vectors")
    for(t in 2:T) {
        psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] +
            (1-psi.hat[t-1])*gamma.hat[t-1]
        tau.hat[t-1] <- gamma.hat[t-1]*(1-psi.hat[t-1]) / psi.hat[t-1]
        }
    return(tau.hat)
    }


###################################################
### code chunk number 14: colext.Rnw:981-995 (eval = FALSE)
###################################################
## 
## chisq <- function(fm) {
##     umf <- getData(fm)
##     y <- getY(umf)
##     sr <- fm@sitesRemoved
##     if(length(sr)>0)
##         y <- y[-sr,,drop=FALSE]
##     fv <- fitted(fm, na.rm=TRUE)
##     y[is.na(fv)] <- NA
##     sum((y-fv)^2/(fv*(1-fv)))
##     }
## 
## set.seed(344)
## pb.gof <- parboot(m0, statistic=chisq, nsim=100)


###################################################
### code chunk number 15: gof (eval = FALSE)
###################################################
## plot(pb.gof, xlab=expression(chi^2), main="", col=gray(0.95),
##      xlim=c(7300, 7700))


###################################################
### code chunk number 16: colext.Rnw:1039-1041
###################################################
data(crossbill)
colnames(crossbill)


###################################################
### code chunk number 17: colext.Rnw:1078-1081
###################################################
DATE <- as.matrix(crossbill[,32:58])
y.cross <- as.matrix(crossbill[,5:31])
y.cross[is.na(DATE) != is.na(y.cross)] <- NA


###################################################
### code chunk number 18: colext.Rnw:1093-1096
###################################################
sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE


###################################################
### code chunk number 19: colext.Rnw:1106-1112
###################################################
years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross,
    siteCovs=crossbill[,2:3], yearlySiteCovs=list(year=years),
    obsCovs=list(date=DATE),
    numPrimary=9)


###################################################
### code chunk number 20: colext.Rnw:1137-1139 (eval = FALSE)
###################################################
## # A model with constant parameters
## fm0 <- colext(~1, ~1, ~1, ~1, umf)


###################################################
### code chunk number 21: colext.Rnw:1141-1143 (eval = FALSE)
###################################################
## # Like fm0, but with year-dependent detection
## fm1 <- colext(~1, ~1, ~1, ~year, umf)


###################################################
### code chunk number 22: colext.Rnw:1145-1147 (eval = FALSE)
###################################################
## # Like fm0, but with year-dependent colonization and extinction
## fm2 <- colext(~1, ~year-1, ~year-1, ~1, umf)


###################################################
### code chunk number 23: colext.Rnw:1149-1151 (eval = FALSE)
###################################################
## # A fully time-dependent model
## fm3 <- colext(~1, ~year-1, ~year-1, ~year, umf)


###################################################
### code chunk number 24: colext.Rnw:1153-1155 (eval = FALSE)
###################################################
## # Like fm3 with forest-dependence of 1st-year occupancy
## fm4 <- colext(~forest, ~year-1, ~year-1, ~year, umf)


###################################################
### code chunk number 25: colext.Rnw:1157-1160 (eval = FALSE)
###################################################
## # Like fm4 with date- and year-dependence of detection
## fm5 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2),
##               umf, starts=c(coef(fm4), 0, 0))


###################################################
### code chunk number 26: colext.Rnw:1162-1165 (eval = FALSE)
###################################################
## # Same as fm5, but with detection in addition depending on forest cover
## fm6 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2) +
##               forest, umf)


###################################################
### code chunk number 27: cov (eval = FALSE)
###################################################
## op <- par(mfrow=c(1,2), mai=c(0.8,0.8,0.1,0.1))
## 
## nd <- data.frame(forest=seq(0, 100, length=50))
## E.psi <- predict(fm6, type="psi", newdata=nd, appendData=TRUE)
## 
## with(E.psi, {
##     plot(forest, Predicted, ylim=c(0,1), type="l",
##          xlab="Percent cover of forest",
##          ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
##     lines(forest, Predicted+1.96*SE, col=gray(0.7))
##     lines(forest, Predicted-1.96*SE, col=gray(0.7))
##     })
## 
## nd <- data.frame(date=seq(-2, 2, length=50),
##                  year=factor("2005", levels=c(unique(years))),
##                  forest=50)
## E.p <- predict(fm6, type="det", newdata=nd, appendData=TRUE)
## E.p$dateOrig <- E.p$date*sd.DATE + mean.DATE
## 
## with(E.p, {
##     plot(dateOrig, Predicted, ylim=c(0,1), type="l",
##          xlab="Julian date", ylab=expression( italic(p) ),
##          cex.lab=0.8, cex.axis=0.8)
##     lines(dateOrig, Predicted+1.96*SE, col=gray(0.7))
##     lines(dateOrig, Predicted-1.96*SE, col=gray(0.7))
##     })
## par(op)


