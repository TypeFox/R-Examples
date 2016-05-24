#--------------------------------------------------------------------
#   variogram.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   svar.bin        S3 class and methods
#   np.svar         S3 class and methods
#   svariso()
#   np.svariso()    S3 generic
#   as.variogram()  S3 generic
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2014
#--------------------------------------------------------------------

# Demonstration of variogram estimation tools (for geoR users)

#--------------------------------------------------------------------
# np.svariso
#--------------------------------------------------------------------

# Load packages
library(npsp)
library(geoR)

# Stationary data (data(s100) geoR)
summary(s100)
plot(s100)

# Empirical variogram
vario.geor <- variog(s100, max.dist=0.6) # geoR variog()
str(vario.geor)

# Local linear variogram 
svarnp <- np.svariso(s100$coords, s100$data, h = 0.15, maxlag = 0.6)
str(svarnp)

# Graphical comparison
oldpar <- par(mfrow=c(1,2))
plot(vario.geor, main="Empirical semivariogram")
plot(svarnp, main = "Nonparametric semivariogram")
par(oldpar) 


#--------------------------------------------------------------------
# as.variogram
#--------------------------------------------------------------------

svar.geor <- as.variogram(svarnp)

vario.ols <- variofit(svar.geor, ini = c(1, 0.5), weights = "equal")  # OLS
vario.ols
vario.wls <- variofit(svar.geor, ini = c(1, 0.5), weights = "cressie")  # WLS
vario.wls

plot(svar.geor, main = "Nonparametric estimates and fitted models")
lines(vario.ols, max.dist = 0.6)
lines(vario.wls, lty = 2, max.dist = 0.6)
legend(0.3, 0.3, legend = c("OLS","WLS"), lty = c(1, 2))


#--------------------------------------------------------------------
# svariso
#--------------------------------------------------------------------

# Empirical variogram
vario.geor <- variog(s100, max.dist=0.6)
varior.geor <- variog(s100, estimator.type = "modulus", max.dist=0.6)
str(vario.geor)

# Linearly binned variogram
svar <- svariso(s100$coords, s100$data, maxlag = 0.6, nlags = 13)
svarr <- svariso(s100$coords, s100$data, maxlag = 0.6, nlags = 13, estimator = "modulus")
str(svar)


oldpar <- par(mfrow = c(2,2))
plot(vario.geor, main = "Empirical semivariogram")
plot(varior.geor, main = "Robust semivariogram")
plot(svar, main = "Binned semivariogram")
plot(svarr, main = "Binned robust semivariogram")
par(oldpar) 


#--------------------------------------------------------------------
# locpol.svar.bin

svarnp <- locpol(svar, h = 0.1)
svarnp2 <- locpol(svar, h = 0.2)  # Equivalent to locpol(svarnp, h = 0.2) 

plot(svar, main = "Binned and nonparametric semivariograms")
plot(svarnp, add = TRUE)
plot(svarnp2, add = TRUE, lty = 2)
legend("bottomright", legend = c("h = 1","h = 2"), lty = c(1, 2))

