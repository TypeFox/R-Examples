### R code from vignette source 'FLAM_fuel.Rnw'

###################################################
### code chunk number 1: pkg-attach
###################################################
library(FDboost)


###################################################
### code chunk number 2: load-data
###################################################
data(fuelSubset)
fuel <- fuelSubset
str(fuel)

# # normalize the wavelength to 0-1
# fuel$nir.lambda0 <- (fuel$nir.lambda - min(fuel$nir.lambda)) / 
#   (max(fuel$nir.lambda) - min(fuel$nir.lambda)) 
# fuel$uvvis.lambda0 <- (fuel$uvvis.lambda - min(fuel$uvvis.lambda)) / 
#   (max(fuel$uvvis.lambda) - min(fuel$uvvis.lambda))

# compute first derivatives as first order differences
fuel$dUVVIS <- t(apply(fuel$UVVIS, 1, diff))
fuel$dNIR <- t(apply(fuel$NIR, 1, diff)) 

# get the wavelength for the derivatives
fuel$duvvis.lambda <- fuel$uvvis.lambda[-1]
fuel$dnir.lambda <- fuel$nir.lambda[-1]
# fuel$duvvis.lambda0 <- fuel$uvvis.lambda0[-1]
# fuel$dnir.lambda0 <- fuel$nir.lambda0[-1]


###################################################
### code chunk number 3: humidity-model (eval = FALSE)
###################################################
## modH2O <- FDboost(h2o ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4) 
##                     + bsignal(NIR, nir.lambda, knots=40, df=4)
##                     + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4) 
##                     + bsignal(dNIR, dnir.lambda, knots=40, df=4), 
##                     timeformula=~bols(1), data=fuel)
## 
## set.seed(212)
## cvmH2O <- suppressWarnings(cvrisk(modH2O, grid=seq(100, 5000, by=100), 
##                               folds=cv( model.weights(modH2O), 
##                               type = "bootstrap", B = 10), mc.cores=10))
## 
## par(mfrow=c(1,2))
## plot(cvmH2O)
## 
## modH2O[mstop(cvmH2O)]
## #modH2O[2400]
##   
## #### create new variable of predicted h2o
## h2o.fit <- modH2O$fitted()
##   
## plot(fuel$h2o, h2o.fit)
## abline(0,1)


###################################################
### code chunk number 4: plot-data
###################################################
# pdf("NIR_UVVIS.pdf", width=7, height=7)
jpeg("NIR_UVVIS.jpg", width=1500, height=1500)
par(mfrow=c(2,2), mar=c(4, 4, 1, 1), cex=1.5)

# generate colors depending on heat value for equidistant cuts
quants <- seq(from=min(fuel$heatan), to=max(fuel$heatan), l=11) 
cats <- cut(fuel$heatan, quants, include.lowest = TRUE)
pall <- heat.colors(12, alpha = 0.5)[1:10]
cols <- pall[cats]

## plot heatan
with(fuel, hist(heatan, breaks=quants, col=pall, 
                xlab="heat value [MJ]", main=""))

## plot heat values versus predicted humidity
with(fuel, plot(heatan~h2o, col=cols, pch=20, lwd=2,
                xlab="predicted humidity [%]", ylab="heat value [MJ]"))
with(fuel, points(heatan~h2o, col=1))

## plot the two spectra
with(fuel, matplot(uvvis.lambda, t(UVVIS), col=cols,
      lwd=1, lty=1, ylab="UV-VIS", xlab="wavelength [nm]", type="l"))
with(fuel, matplot(nir.lambda, t(NIR), col=cols,
      lwd=1, lty=1, ylab="NIR", xlab="wavelength [nm]", type="l"))
dev.off()


###################################################
### code chunk number 5: model-spec
###################################################
formula <- formula(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                   + bsignal(NIR, nir.lambda, knots=40, df=4.41))

## do a model fit:
mod <- FDboost(formula, timeformula=~bols(1), data=fuel)
mod <- mod[198]


###################################################
### code chunk number 6: cv-model-spec (eval = FALSE)
###################################################
## ## get optimal mstop and do bootstrapping for coefficient estimates
## set.seed(2703)
## val <- validateFDboost(mod, 
##                        folds=cv(model.weights(mod), type = "bootstrap", B = 50),
##                        grid = 10:500, mc.cores=10)
## 
## mopt <- val$grid[which.min(colMeans(val$oobrisk))]
## print(mopt)
## 
## ## use optimal mstop
## mod <- mod[mopt] # 198


###################################################
### code chunk number 7: model-spec-plot
###################################################
par(mfrow=c(1,2))
plot(mod, which=1, lwd=2, lty=5, rug=FALSE,
     ylab="", xlab="wavelength [nm]")

plot(mod, which=2, lwd=2, lty=5, rug=FALSE, 
     ylab="", xlab="wavelength [nm]")

# plot with bootstrapped coefficient functions
if(FALSE){
  pdf("spec_valCoef.pdf")
  par(mar=c(5, 3, 1, 1), cex.axis=1.5, cex.lab=1.5)
  plot(mod, which=1, lwd=2, col="white", main="", lty=5, rug=FALSE,
     ylab="", xlab="wavelength [nm]", ylim=range(val$coefCV[[1]]$value) )
  plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=1, add=TRUE)
  plot(mod, which=1, lwd=2, col=4, main="", lty=5, rug=FALSE, add=TRUE)

  legend("topright", legend=c("data", "BS", "mean BS", "5, 95% BS"), 
       lty=c(5,1,1,2), col=c(4,8,1,2), cex=1.5,
       lwd=c(2,1,2,2))

  plot(mod, which=2, lwd=2, col="white", main="", lty=5, rug=FALSE, 
     ylab="", xlab="wavelength [nm]", ylim=range(val$coefCV[[2]]$value) )
  plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=2, add=TRUE)
  plot(mod, which=2, lwd=2, col=4, main="", lty=5, rug=FALSE, add=TRUE)
  dev.off()
}


