## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("kohonen")) {
  kohonen.present <- FALSE
  cat("Package kohonen not available - some code may not run.\nInstall it by typing 'install.packages(\"kohonen\")'")
} else {
  kohonen.present <- TRUE
}

data(wines, package = "ChemometricsWithRData")
wines.sc <- scale(wines)

set.seed(7)
wines.som <- som(wines.sc, somgrid(5, 4, "hexagonal"))
plot(wines.som, type = "codes")

for (i in c(1, 8, 11, 13)) 
  plot(wines.som, "property", property = wines.som$codes[,i],
       main = colnames(wines)[i])

plot(wines.som, type = "mapping",
     col = as.integer(vintages), pch = as.integer(vintages))
plot(wines.som, type = "dist.neighb")
plot(wines.som, "changes")
plot(wines.som, "quality")

summary(wines.som)

## Next lines take too long for a demo
##
## data(Prostate2000Raw, package = "ChemometricsWithRData")
## X <- t(Prostate2000Raw$intensity)
## types <- Prostate2000Raw$type
## prostate.som <- som(X, somgrid(7, 5, "hexagonal"))
## plot(prostate.som, "mapping", col = as.integer(types),
##      pch =  as.integer(types),
##      main = "Prostate data")
## legend("bottom", legend = levels(types),
##        col = 1:3, pch = 1:3, ncol = 3)
 

## cols <- c(2,1,3)
## units <- c(7,21,35)
## par(mfrow = c(3,1))
## for (i in 1:3)
##   plot(prostate.som$codes[units[i],], type = "l", 
##        col = cols[i], main = paste("Unit", units[i]), 
##        ylab = "codebook vector")
