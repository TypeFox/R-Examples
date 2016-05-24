#--------------------------------------------------------------------
# aquifer.R (npsp package demo)
#--------------------------------------------------------------------
#   Wolfcamp aquifer data
#
#   (c) R. Fernandez-Casal         Last revision: May 2014
#--------------------------------------------------------------------
library(npsp)
# windows(7, 5, record = TRUE)

# ?aquifer
str(aquifer)
summary(aquifer)

# Scatter plot with a color scale
with(aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer data"))


#------------------------------------------
# Linear binning

bin <- binning(aquifer[,1:2], aquifer$head, nbin = c(41,41), set.NA = TRUE) 
simage(bin, main = 'Binning averages')
points(bin$data$x, col = 'darkgray')

#------------------------------------------
# Density estimation

h.den <- diag(c(55,55)) # h.cv(as.bin.den(bin))$h 
den <- np.den(bin, h = h.den, degree = 0)

plot(den, main = 'Estimated log(density)')

# Index with grid nodes far from data
na.index <- log(den$est) < -15


#------------------------------------------
# Trend estimation

lp <- locpol(bin, h = diag(100, 2), hat.bin = TRUE)  # np.svariso.corr(lp, ...) : 'lp' must have a '$locpol$hat' component

# Delete grid nodes far from data
lp$est[na.index] <- NA

# Perspective plot with a color scale
spersp(lp, main = 'Trend estimates', zlab = 'piezometric-head levels', theta = 120)            

#------------------------------------------
# Variogram estimation

lp.resid <- residuals(lp) # lp.resid <- lp$data$y - predict(lp)
# maxlag <- 0.5*sqrt(sum(diff(apply(aquifer[,1:2], 2, range))^2))

esvar <- np.svariso(aquifer[,1:2], lp.resid, maxlag = 150, nlags = 60, h = 60)
svm <- fitsvar.sb.iso(esvar)  # dk = 2

plot(svm, main = "Nonparametric semivariogram and fitted model")

 
#------------------------------------------
# Note that the direct use of the residuals introduces a bias in the estimation 
# of the variogram. This bias is usually negative and higher at large lags 
# (e.g. Cressie, 1993, section 3.4.3).
# A correction for this bias is proposed in:
# Fernandez-Casal R. and Francisco-Fernandez M. (2013) 
# Nonparametric bias-corrected variogram estimation under non-constant trend. 
# Stoch. Environ. Res. Ris. Assess (SERRA), 1-14, doi:10.1007/s00477-013-0817-8.
# A similar algorithm (fully nonparametric) is implemented in 'np.svariso.corr'

esvar2 <- np.svariso.corr(lp, maxlag = 150, nlags = 60, h = 60, plot = TRUE)

svm2 <- fitsvar.sb.iso(esvar2)  # dk = 2

plot(svm2, main = "Nonparametric bias-corrected semivariogram and fitted models", legend = FALSE) 
plot(svm, add = TRUE, lty = 2)
legend("bottomright", legend = c("NP estimates", "fitted model", 'uncorrected'),
            lty = c(NA, 1, 2), pch = c(1, NA, NA), lwd = c(1, 2, 1))


#------------------------------------------
# Bandwidth selection and trend re-estimation?

bin2 <- binning(aquifer[,1:2], aquifer$head, nbin = c(15,15)) # to speed computations... 

h.cv(bin2, ncv = 2)   # ncv >= 2 is recommended for sparse data (when linear binning is used)

cov.dat <-  varcov(svm2, coords = aquifer[,1:2])
hcv.data(bin2, objective = "GCV", cov = cov.dat) # GCV criterion (Francisco-Fernandez and Opsomer, 2005) for spatially correlated data



#------------------------------------------
# Kriging

library(sp)
library(gstat)

spdf <- SpatialPointsDataFrame(aquifer[,1:2], data.frame(y = aquifer$head, r = lp.resid))
newdata <- SpatialPoints(coords(lp))

# Simple kriging of residuals
krig <- krige(r ~ 1, locations = spdf, newdata = newdata, model = as.vgm(svm), beta = 0)
krig.grid <- data.grid(kpred = lp$est + krig@data$var1.pred, ksd = sqrt(krig@data$var1.var), 
        grid = lp$grid)

krig2 <- krige(r ~ 1, locations = spdf, newdata = newdata, model = as.vgm(svm2), beta = 0)
krig2.grid <- data.grid(kpred = lp$est + krig2@data$var1.pred, ksd = sqrt(krig2@data$var1.var), 
        grid = lp$grid)

scale.color <- jet.colors(64)
scale.range <- c(1100, 4100)
# 1x2 plot with some room for the legend...
old.par <- par(mfrow = c(1,2), omd = c(0.01, 0.9, 0.05, 0.95), plt= c(0.08, 0.94, 0.1, 0.8))
spersp(krig.grid, main = 'Kriging predictions', col = scale.color, legend = FALSE, theta = 120)
spersp(krig2.grid, main = 'Kriging predictions \n (bias-corrected)', col = scale.color, legend = FALSE, theta = 120)
par(old.par)
splot(slim = scale.range, col = scale.color, legend.shrink = 0.6, add = TRUE)

old.par <- par(mfrow = c(1,2), omd = c(0.05, 0.85, 0.05, 0.95))
scale.range <- c(125, 200)
scale.range <- range(krig.grid$ksd, krig2.grid$ksd, finite = TRUE)
image( krig.grid, 'ksd', zlim = scale.range, main = 'Kriging sd', col = scale.color)
with(aquifer, points(lon, lat, cex = 0.75))
image( krig2.grid, 'ksd', zlim = scale.range, main = 'Kriging sd (bias-corrected)', col = scale.color)
with(aquifer, points(lon, lat, cex = 0.75))
par(old.par)
splot(slim = scale.range, col = scale.color, add = TRUE)

# NOTE: To reproduce results in SERRA paper use 'data(wolfcamp)' in package 'geoR' 
# (note also that the multiplicative Epanechnikov kernel was used in that work).
# Results obtained with 'aquifer' data set are comparable with those in Cressie (1993, section 4.1). 
