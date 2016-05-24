## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.

data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
matplot(wavelengths, t(gasoline$NIR), type = "l", 
        lty = 1, xlab = "Wavelength (nm)", ylab = "1/R")

data(wines, package = "ChemometricsWithRData")
colnames(wines)

table(vintages)

pairs(wines[,1:3], pch = wine.classes, col = wine.classes)

data(Prostate2000Raw, package = "ChemometricsWithRData")
plot(Prostate2000Raw$mz, Prostate2000Raw$intensity[,1],
     type = "h", xlab = "m/z", ylab = "Intensity",
     main = "Prostate data")

table(Prostate2000Raw$type)
