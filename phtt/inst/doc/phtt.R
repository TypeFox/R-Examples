### R code from vignette source 'phtt.Snw'

###################################################
### code chunk number 1: phtt.Snw:149-176
###################################################
## Install package
## install.packages("phtt")
## Load package
library("phtt")
## Load Data
data("Cigar")
N <- 46
T <- 30
## Dependent variable:
## Cigarette-Sales per Capita
l.Consumption    <- log(matrix(Cigar$sales, T, N))
## Independent variables:
## Consumer Price Index
cpi              <- matrix(Cigar$cpi, T, N)
## Real Price per Pack of Cigarettes 
l.Price          <- log(matrix(Cigar$price, T, N)/cpi)
## Real Disposable Income per Capita  
l.Income         <- log(matrix(Cigar$ndi, T, N)/cpi)

pdf("Cigar.pdf")
scl <- 1.6
par(mfrow=c(1,3), mar=c(6, 5, 5, 2.1))
matplot(l.Consumption,type="l", main="Log's of\nCigar-Consumption",ylab="",xlab="Time",cex.main=scl,cex.lab=scl,cex.axis=scl)
matplot(l.Price, type="l",      main="Log's of\nreal Prices",    ylab="",xlab="Time",cex.main=scl,cex.lab=scl, cex.axis=scl)
matplot(l.Income, type="l",     main="Log's of\nreal Income",    ylab="",xlab="Time",cex.main=scl,cex.lab=scl, cex.axis=scl)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
graphics.off()


###################################################
### code chunk number 2: phtt.Snw:291-292
###################################################
args(KSS)


###################################################
### code chunk number 3: phtt.Snw:316-324
###################################################
library("phtt")
data("Cigar")
N <- 46
T <- 30
l.Consumption   <- log(matrix(Cigar$sales, T, N))
cpi             <- matrix(Cigar$cpi,       T, N)
l.Price         <- log(matrix(Cigar$price, T, N)/cpi)
l.Income        <- log(matrix(Cigar$ndi,   T, N)/cpi)


###################################################
### code chunk number 4: phtt.Snw:329-331
###################################################
Cigar.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income) 
(Cigar.KSS.summary <- summary(Cigar.KSS))


###################################################
### code chunk number 5: phtt.Snw:335-340
###################################################
## Figure 2:
pdf("KSSM1.pdf")
scl <- 1
plot(Cigar.KSS.summary,cex.main=scl,cex.lab=scl,cex.axis=scl)
graphics.off()


###################################################
### code chunk number 6: phtt.Snw:344-345 (eval = FALSE)
###################################################
## plot(Cigar.KSS.summary)


###################################################
### code chunk number 7: phtt.Snw:437-438
###################################################
args(OptDim)


###################################################
### code chunk number 8: phtt.Snw:448-449
###################################################
OptDim(Obj = l.Consumption, criteria = "PC1")


###################################################
### code chunk number 9: phtt.Snw:454-456
###################################################
(OptDim.obj <- OptDim(Obj = l.Consumption, criteria = c("PC3",  "ER",  
                      "GR", "IPC1", "IPC2", "IPC3"), standardize = TRUE))


###################################################
### code chunk number 10: phtt.Snw:461-464
###################################################
pdf("OptDimv.pdf")
plot(OptDim.obj)
graphics.off()


###################################################
### code chunk number 11: phtt.Snw:468-469 (eval = FALSE)
###################################################
## plot(OptDim.obj) 


###################################################
### code chunk number 12: phtt.Snw:485-486 (eval = FALSE)
###################################################
## KSS(formula = l.Consumption ~ -1 + l.Price + l.Income, consult.dim = TRUE) 


###################################################
### code chunk number 13: phtt.Snw:683-684
###################################################
args(Eup)


###################################################
### code chunk number 14: phtt.Snw:699-702
###################################################
d.l.Consumption  <- diff(l.Consumption)
d.l.Price        <- diff(l.Price)
d.l.Income       <- diff(l.Income)


###################################################
### code chunk number 15: phtt.Snw:706-708
###################################################
(Cigar.Eup <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
                  dim.criterion = "PC3"))


###################################################
### code chunk number 16: phtt.Snw:724-725
###################################################
summary(Cigar.Eup)


###################################################
### code chunk number 17: phtt.Snw:731-734
###################################################
pdf("EupPlot.pdf")
plot(summary(Cigar.Eup))
graphics.off()


###################################################
### code chunk number 18: phtt.Snw:738-739 (eval = FALSE)
###################################################
## plot(summary(Cigar.Eup))


###################################################
### code chunk number 19: phtt.Snw:804-807
###################################################
Cigar2.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income,
                  additive.effects = "individual") 
Cigar2.KSS.summary <- summary(Cigar2.KSS)


###################################################
### code chunk number 20: phtt.Snw:810-813 (eval = FALSE)
###################################################
## Cigar2.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income,
##                   additive.effects = "individual") 
## (Cigar2.KSS.summary <- summary(Cigar2.KSS))


###################################################
### code chunk number 21: phtt.Snw:845-846 (eval = FALSE)
###################################################
## plot(Cigar2.KSS.summary)


###################################################
### code chunk number 22: phtt.Snw:849-855
###################################################
pdf("KSSM2.pdf")
scl <- 1.75
par(mar=c(6, 5, 5, 2.1))
plot(Cigar2.KSS.summary,cex.main=1.25,cex.lab=scl,cex.axis=scl)
par(mar=c(5.1, 4.1, 4.1, 2.1))
graphics.off()


###################################################
### code chunk number 23: phtt.Snw:902-903 (eval = FALSE)
###################################################
## checkSpecif(obj1, obj2, level = 0.05)


###################################################
### code chunk number 24: phtt.Snw:917-922 (eval = FALSE)
###################################################
## twoways.obj     <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
##                        factor.dim = 0, additive.effects = "twoways")
## not.twoways.obj <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
##                        factor.dim = 2, additive.effects = "none")  
## checkSpecif(obj1 = twoways.obj, obj2 = not.twoways.obj, level = 0.01)


###################################################
### code chunk number 25: phtt.Snw:963-966
###################################################
Eup.obj <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
               additive.effects = "twoways")
checkSpecif(Eup.obj, level = 0.01)


###################################################
### code chunk number 26: phtt.Snw:970-973
###################################################
KSS.obj <- KSS(l.Consumption ~  -1 + l.Price + l.Income, 
               additive.effects = "twoways")
checkSpecif(KSS.obj, level = 0.01)


###################################################
### code chunk number 27: phtt.Snw:1020-1031 (eval = FALSE)
###################################################
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[1]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[2]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[3]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[4]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[5]
## 
## coef(Cigar2.KSS)$Common.factors[,2]
## lambda_i1 <- coef(Cigar2.KSS)$Ind.loadings.param[,1]
## order(abs(lambda_i1),decreasing=TRUE)[1]
## coef(Cigar2.KSS)$Ind.loadings.param[7,1]
## round(range(coef(Cigar2.KSS)$Ind.loadings.param[-7,1]), digits=2)


###################################################
### code chunk number 28: phtt.Snw:1034-1041
###################################################
pdf("Factor2.pdf")
scl <- 1
par(mfrow=c(1,2))
matplot(coef(Cigar2.KSS)$Common.factors[,1]%*%t(coef(Cigar2.KSS)$Ind.loadings.param[,1]),type="l",main="Variance of time-varying indiv. effects\n in direction of the 1. common factor",xlab="Time",ylab="",cex.main=scl,cex.lab=scl,cex.axis=scl)
matplot(coef(Cigar2.KSS)$Common.factors[,2]%*%t(coef(Cigar2.KSS)$Ind.loadings.param[,2]),type="l",main="Variance of time-varying indiv. effects\n in direction of the 2. common factor",xlab="Time",ylab="",cex.main=scl,cex.lab=scl,cex.axis=scl)
par(mfrow=c(1,1))
graphics.off()


