###################################################
## R script for 'fda.usc' manuscript 
## submitted to Journal of Statistical Software
## Authors: Manuel Oviedo and Manuel Febrero-Bande
##
###################################################

###################################################
### chunk number 1: tecator dataset
###################################################
library("fda.usc")
data(tecator)
names(tecator)
absorp <- tecator$absorp.fdata
Fat20 <- ifelse(tecator$y$Fat < 20, 0, 1) * 2 + 2
absorp$names$main <- ""

dev.new( width = 150,height = 110, units = "mm")
par(mfrow = c(1 , 2))
# Figure 1 (left panel)
plot(absorp, col = Fat20)
absorp.d1 <- fdata.deriv(absorp, nderiv = 1)
# Figure 1 (right panel)
plot(absorp.d1, col = Fat20)


###################################################
### chunk number 2: Convert fdata class to fd class
###################################################
class( absorp.fd <- fdata2fd( absorp, type.basis = "fourier", nbasis = 15) )
class( absorp.fdata <- fdata(absorp.fd) )

###################################################
### chunk number 3:  phoneme data and smoothing
###################################################
data(phoneme)
learn <- phoneme$learn
l <- c(0 ,2 ^ seq(-2, 9, len = 30))
nb <- seq(7, 31, by = 2)
out0 <- min.basis(learn, lambda = l, numbasis = nb)
out1 <- min.np(learn, type.S = S.NW, par.CV = list(criteria = "GCV"))
out2 <- min.np(learn, type.S = S.LLR, par.CV = list(criteria = "GCV"))
out3 <- min.np(learn, type.S = S.KNN, h = 3:25, Ker = Ker.unif)

###################################################
### chunk number 4:  plot GCV criteria
###################################################
dev.new(width = 150, height = 110, units = "mm")
par(mfrow = c(1,2))
contour(nb, l, out0$gcv, ylab = "Lambda", xlab = "Number of basis", 
main = "GCV criteria by min.basis()")
plot(out1$h, out1$gcv, type = "l", main = "GCV criteria  by min.np() ", 
xlab = "Bandwidth (h) values",ylab = "GCV criteria", col = 3, lwd = 2)
legend(x = 3, y = 4.3, legend = c("Ker.norm-S.NW", "Ker.norm-S.LLR", "Ker.unif-S.KNN"),
box.col = "white", lwd = c(2, 2, 2), col = c(3, 4, 2),cex = 0.75)
lines(out2$h,out2$gcv, col = 4, lwd = 2)
lines(out3$h,out3$gcv, col = 2, lwd = 2)

###################################################
### chunk number 5: smoothing a fdata curve
###################################################
dev.new( width = 150, height = 100, units = "mm")
ind <- 11
nam <- expression( paste("Phoneme curve"[11]) )
plot(learn[ind, ], main = nam, lty = 2, lwd = 2, col = 8)
legend(x = 70, y = 19, legend = c("Curve","Bspline basis",
"Ker.norm-S.NW", "Ker.norm-S.LLR", "Ker.unif-S.KNN"),
lty = c(2, 1, 1, 1, 1), lwd = 2, col = c(8, 1, 3, 4, 2), box.col = "white")
lines(out0$fdata.est[ind, ], col = 1, lty = 1, lwd = 2)
lines(out1$fdata.est[ind, ], col = 3, lty = 1, lwd = 2)
lines(out2$fdata.est[ind, ], col = 4, lty = 1, lwd = 2)
lines(out3$fdata.est[ind, ], col = 2, lty = 1, lwd = 2)
lines(out0$fdata.est[ind, ], col = 1, lty = 1, lwd = 2)

###################################################
### chunk number 6:  metric and semimetric measures
###################################################
glearn <- phoneme$classlearn
mdist <- metric.lp(learn)

###################################################
### chunk number 7:  dendogram
###################################################
ind <- c(120:130,220:230)
dev.new()
a <- mdist[ind,ind]
b <- as.dist(a)
c2 <- hclust(b)
c2$labels <- phoneme$classlearn[ind]
dev.new( width = 110, height = 80, units = "mm")
plot(c2, main = "Dendogram", xlab = "Class of each leaf")

###################################################
### chunk number 8: figure depth
###################################################
data(poblenou)
nox <- poblenou$nox
working <- poblenou$nox[poblenou$df$day.festive == 0 &
 as.integer(poblenou$df$day.week) < 6]
nonworking <- poblenou$nox[poblenou$df$day.festive == 1 |
as.integer(poblenou$df$day.week) > 5]
# Centrality measures (working)

dev.new( width = 180,height = 150, units = "mm")
par( mfrow=c(2, 2) )
plot(func.mean(working), ylim = c(10, 170),
main = "Centrality measures in working days")
legend(x = 11, y = 170, cex = 1, box.col = "white", lty = 1:5,
col = c(1:5), legend = c("mean","trim.mode","trim.RP",
"median.mode","median.RP"))
lines(func.trim.mode(working, trim = 0.15), col = 2, lty = 2)
lines(func.trim.RP(working, trim = 0.15), col = 3, lty = 3)
lines(func.med.mode(working, trim = 0.15), col = 4, lty = 4)
lines(func.med.RP(working, trim = 0.15), col = 5, lty = 5)

# Centrality measures (non-working)
plot(func.mean(nonworking), ylim = c(10,170),
main = "Centrality measures in non-working days")
legend(x = 11, y = 170, cex = 1, box.col = "white",lty = 1:5,
col = c(1:5), legend = c("mean","trim.mode","trim.RP",
"median.mode","median.RP"))
lines(func.trim.mode(nonworking, trim = 0.15),col = 2, lty = 2)
lines(func.trim.RP(nonworking, trim = 0.15),col = 3, lty = 3)
lines(func.med.mode(nonworking, trim = 0.15),col = 4, lty = 4)
lines(func.med.RP(nonworking, trim = 0.15),col = 5, lty = 5)

# Measures of dispersion   (working)
plot(func.var(working),
main = "Dispersion measures in working days", ylim = c(100 ,5500))
legend(x = 11, y = 5300,cex = 1, box.col = "white", lty = 1:3, col = 1:3,
legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(working,trim = 0.15), col = 2, lty = 2)
lines(func.trimvar.RP(working,trim = 0.15), col = 3, lty = 3)

# Measures of dispersion   (non-working)
plot(func.var(nonworking),
main = "Dispersion measures in non-working days", ylim = c(100, 5500))
legend(x = 11, y = 5300, cex = 1, box.col = "white", lty = 1:3, col = 1:3,
legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(nonworking, trim = 0.15), col = 2, lty = 2)
lines(func.trimvar.RP(nonworking, trim = 0.15), col = 3, lty = 3)

###################################################
### chunk number 9:  plot funtional bootstrap
###################################################
# This takes a lot
dev.new(width = 150, height = 110, units = "mm")
par(mfrow = c(1, 2))
out.boot1 <- fdata.bootstrap(working, statistic = func.trim.RP, nb = 1000, draw = TRUE)
out.boot2 <- fdata.bootstrap(nonworking,statistic = func.trim.RP ,nb = 1000, draw = TRUE)

###################################################
### chunk number 10: functional outlier detection
###################################################
# This take a lot
nox <- poblenou$nox
working <- poblenou$nox[poblenou$df$day.festive == 0 &
as.integer(poblenou$df$day.week) < 6]
nonworking <- poblenou$nox[poblenou$df$day.festive == 1 |
as.integer(poblenou$df$day.week) > 5]

out1 <- outliers.depth.trim(working, dfunc=depth.FM,
nb=1000, smo=0.1, trim=0.06)
out2 <- outliers.depth.trim(nonworking, dfunc=depth.FM,
nb=1000, smo=0.1, trim=0.1)

###################################################
### chunk number 11: Figure outlier detection
###################################################
dev.new( width=150, height=110, units="mm")                                     
par(mfrow=c(1, 2))
plot(working, ylim=c(0,400), col="gray", lty=1,
main="NOx - Working days")
lines(working[out1[[1]]], col=2, lty=2, lwd=2)
plot(nonworking, ylim=c(0,400), col="gray", lty=1,
main="NOx - Non working days")
lines(nonworking[out2[[1]]], col=2, lty=2, lwd=2)

###################################################
### chunk number 12:  functional regression model
###################################################
ind <- 1:165
tt <- absorp[["argvals"]]
# training data
y <- tecator$y$Fat[ind]
X <- absorp[ind, ]
X.d1 <- fdata.deriv(X, nbasis = 19, nderiv = 1)
X.d2 <- fdata.deriv(X, nbasis = 19, nderiv = 2)

###################################################
### chunk number 13:   fregre.basis
###################################################
rangett <- absorp$rangeval
basis1 <- create.fourier.basis(rangeval = rangett,
nbasis=5)
res.basis1 <- fregre.basis(X.d1, y, basis.x = basis1)

###################################################
### chunk number 14: summary.fregre.fd
###################################################
dev.new( width = 160, height = 130,units = "mm")
summary(res.basis1)

###################################################
### chunk number 15: fregre.pc
###################################################
res.pc1 <- fregre.pc(X.d1, y, l = 1:6)

###################################################
### chunk number 16:  fregre.pc.cv
###################################################
res.pc2 <- fregre.pc.cv(X.d2, y, kmax = 7)
res.pc2$pc.opt
res.pc2$MSC

###################################################
### chunk number 17: fregre.pls
###################################################
fregre.pls(X.d1, y, l = 1:5)
res.pls1 <- fregre.pls.cv(X.d1, y)$fregre.pls

###################################################
### chunk number 18: influence measures
###################################################
# This take a lot
res.infl <- influence.fdata(res.basis1)
# Influence Measure for FPC regression
mat <- cbind(res.infl$DCP, res.infl$DCE, res.infl$DP)
colnames(mat)=c("CPi", "CEi", "Pi")
dev.new( width=140, height = 110, units = "mm")
pairs(mat)

###################################################
### chunk number 19: functional beta estimation
###################################################
res.boot1 <- fregre.bootstrap(res.basis1, nb = 1000,
kmax.fix = TRUE, alpha = 0.999)
res.boot2 <- fregre.bootstrap(res.pc1, nb = 1000,
kmax.fix = TRUE, alpha = 0.999)
res.boot3 <- fregre.bootstrap(res.pls1, nb = 1000,
kmax.fix = TRUE, alpha = 0.999)

dev.new( width = 170, height = 130, units = "mm")
par(mfrow=c(1, 3))
yl <- c(-200, 220)
out  <-  res.boot1$norm.boot > quantile(res.boot1$norm.boot, 0.999)
plot(res.boot1$betas.boot, 
col="grey", main = "5 Fourier basis elements", ylim = yl, lty=1, lwd=2)
lines(res.basis1$beta.est, col = 4, lwd = 1, lty = 1)
lines(res.boot1$betas.boot[out], lwd = 2, col = 2, lty = 2)
out  <-  res.boot2$norm.boot > quantile(res.boot2$norm.boot, 0.999)
plot(res.boot2$betas.boot,
col="grey", main="1st 6 FPC", ylim = yl, lty = 1, lwd = 2)
lines(res.pc1$beta.est,col = 4, lwd = 1, lty = 1)
lines(res.boot2$betas.boot[out], lwd = 2, col = 2, lty = 2)
out  <- res.boot3$norm.boot > quantile(res.boot3$norm.boot, 0.999)
plot(res.boot3$betas.boot, 
col = "grey", main = "1st 5 FPLS", ylim = yl, lty = 1, lwd = 2)
lines(res.pls1$beta.est, col = 4, lwd = 1, lty = 1)
lines(res.boot3$betas.boot[out], col = 2, lwd=2, lty = 2)

###################################################
### chunk number 20:  fregre.lm
###################################################
ind <- 1:165
dataf <- as.data.frame(tecator$y[ind, ])
newdataf <- as.data.frame(tecator$y[-ind, ])
ldata=list("df" = dataf, "X" = X, "X.d1" = X.d1, "X.d2" = X.d2)
basis.x <- list("X" = basis1)
basis2 <- create.bspline.basis(rangeval = rangett, nbasis = 5)
basis.b <- list("X" = basis2)
# functional data with basis representation
f2<- Fat ~ Water + X.d2
basis.x1 <- list("X.d1" = basis1)
basis.b1 <- list("X.d1" = basis2)
res.lm2 <- fregre.lm(f2, ldata, basis.x = basis.x1, basis.b = basis.b1)

###################################################
### chunk number 21:   fregre.np
###################################################
# Functional nonparametric regression
fregre.np(X.d1, y)

###################################################
### chunk number 22: fregre.plm
###################################################
# Functional semi-parametric regression
res.plm2 <- fregre.plm(Fat ~ Water + X.d2,ldata)

###################################################
### chunk number 23: prediction I
###################################################
newy <- matrix(tecator$y$Fa[-ind], ncol = 1)
newX.d1 <- fdata.deriv(absorp[-ind , ], nbasis = 19,nderiv = 1)
newX.d2 <- fdata.deriv(absorp[-ind, ], nbasis = 19,nderiv = 2)
res.basis2 <- fregre.basis.cv(X.d2, y, type.basis = "bspline")
pred.basis2 <- predict.fregre.fd(res.basis2, newX.d2)
res.np2 <- fregre.np.cv(X.d1, y, metric = semimetric.deriv)
pred.np2 <- predict.fregre.fd(res.np2, newX.d1)

# Functional and nonfunctional covariate
newldata <- list("df" = newdataf, "X.d1" = newX.d1, "X.d2" = newX.d2)
f1 <- Fat ~ Water + X.d1
basis.x1 <- list(X.d1 = basis1)
basis.b1 <- list(X.d1 = basis2)
res.lm1 <- fregre.lm(f1, ldata, basis.x = basis.x1,
basis.b <- basis.b1)
pred.lm1<-predict.fregre.lm(res.lm1, newldata)
res.plm1 <- fregre.plm(f1, ldata)
pred.plm1 <- predict.fregre.plm(res.plm1, newldata)

###################################################
### chunk number 24: prediction III
###################################################
pred.basis1 <- predict.fregre.fd(res.basis1, newX.d1)
pred.pc1 <- predict.fregre.fd(res.pc1, newX.d1)
pred.pls1 <- predict.fregre.fd(res.pls1, newX.d1)
res.np1 <- fregre.np.cv(X.d1, y)
pred.np1 <- predict.fregre.fd(res.np1, newX.d1)
res.pls2 <- fregre.pls.cv(X.d2, y, criteria = "CV")$fregre.pls
pred.pls2 <- predict.fregre.fd(res.pls2, newX.d2)
pred.pc2 <- predict.fregre.fd(res.pc2$fregre.pc, newX.d2)
pred.lm2 <- predict.fregre.lm(res.lm2, newldata)
pred.plm2 <- predict.fregre.plm(res.plm2, newldata)

###################################################
### chunk number 25: Table 1,  Goodness of fit
###################################################
GOF.predict <- function(model,pred,newy){
if (any(class(model)=="fregre.fd")) {
 df <- model$df
 r2 <- model$r2
 sr2 <- model$sr2
   }
else {
 if (any(class(model) == "lm")) {
 smr <- summary(model)
 df <- smr$df[1]
 r2 <- smr$r.squared
 sr2 <-   sum(model$residuals^2) / (length(model$residuals) - df)
 }
}
 MEP <- ((1 / length(newy)) * sum((newy - pred)^2)) / var(newy)
 out <- cbind(round(df, 1), round(r2, 3), round(sr2, 3), round(MEP, 4))
 colnames(out) <- c("df", "R-squared", "Sr2", "MEP")
out
}

tabl <- matrix(NA, nrow = 12, ncol = 4)
tabl[1, ] <- GOF.predict(res.basis1, pred.basis1, newy)
tabl[2, ] <- GOF.predict(res.basis2, pred.basis2, newy)
tabl[3, ] <- GOF.predict(res.pc1, pred.pc1, newy)
tabl[4, ] <- GOF.predict(res.pc2$fregre.pc, pred.pc2, newy)
tabl[5, ] <- GOF.predict(res.pls1, pred.pls1, newy)
tabl[6, ] <- GOF.predict(res.pls2, pred.pls2, newy)
tabl[7, ] <- GOF.predict(res.lm1, pred.lm1, newy)
tabl[8, ] <- GOF.predict(res.lm2, pred.lm2, newy)
tabl[9, ] <- GOF.predict(res.np1, pred.np1, newy)
tabl[10, ] <- GOF.predict(res.np2, pred.np2, newy)
tabl[11, ] <- GOF.predict(res.plm1, pred.plm1, newy)
tabl[12, ] <- GOF.predict(res.plm2, pred.plm2,newy)
colnames(tabl) <- c("df", "R-squared", "Sr2", "MEP")
tabl
###################################################
###################################################
