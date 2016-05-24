## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("ptw")) {
  ptw.present <- FALSE
  cat("Package ptw not available - some code may not run.\nInstall it by typing 'install.packages(\"ptw\")'")
} else {
  ptw.present <- TRUE
}
if (!require("dtw")) {
  dtw.present <- FALSE
  cat("Package dtw not available - some code may not run.\nInstall it by typing 'install.packages(\"dtw\")'")
} else {
  dtw.present <- TRUE
}
if (!require("signal")) {
  signal.present <- FALSE
  cat("Package signal not available - some code may not run.\nInstall it by typing 'install.packages(\"signal\")'")
} else {
  signal.present <- TRUE
}

## smoothing
data(Prostate2000Raw, package = "ChemometricsWithRData")
## Following lines are slooooow
## prostate.array <- array(t(Prostate2000Raw$intensity),
##                         c(2, 327, 10523))
## prostate <- apply(prostate.array, c(2,3), mean)
## dim(prostate)

prostate <- rowsum(t(Prostate2000Raw$intensity), 
                   group = rep(1:327, each = 2),
                   reorder = FALSE) / 2
dim(prostate)

prostate.type <- Prostate2000Raw$type[seq(1, 654, by = 2)]
x <- cbind(prostate[1,1:500],
           Prostate2000Raw$intensity[1:500, 1:2])
matplot(Prostate2000Raw$mz[1:500], x, type = "l",
        col = c(1, "gray", "gray"), lty = c(1,2,2),
        xlab = "m/z", ylab = "response")

rmeans <- rowMeans(embed(prostate[1,1:500], 5))
plot(Prostate2000Raw$mz[1:500], prostate[1,1:500],
     type = "l", xlab = "m/z", ylab = "response",
     main = "running means", col = "gray")
lines(Prostate2000Raw$mz[3:498], rmeans, type = "l")
plot(Prostate2000Raw$mz[1:500], prostate[1,1:500],
     type = "l", xlab = "m/z", ylab = "response",
     main = "running median", col = "gray")
lines(Prostate2000Raw$mz[1:500],
      runmed(prostate[1,1:500], k = 5), type = "l")

mznew <- colMeans(matrix(Prostate2000Raw$mz[1:500], nrow = 5))
xnew <- colMeans(matrix(prostate[1, 1:500], nrow = 5))
plot(Prostate2000Raw$mz[1:500], prostate[1, 1:500], 
     type = "l", xlab = "m/z", ylab = "response",
     main = "binning", col = "gray")
lines(mznew, xnew)

## baseline removal
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
nir.diff <- t(apply(gasoline$NIR, 1, diff))
matplot(wavelengths[-1] + 1, t(nir.diff),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "1/R (1st deriv.)", lty = 1, col = 1)

if (signal.present) {
  nir.deriv <- apply(gasoline$NIR, 1, sgolayfilt, m = 1)
}

nir.msc <- msc(gasoline$NIR)

if (ptw.present) {
  data(gaschrom, package = "ptw")
  x <- gaschrom[1,]
  lsection <- 200
  xmat <- matrix(x, nrow=lsection)
  ymin <- apply(xmat, 2, min)
  plot(x, type = "l", col = "gray", ylim = c(20, 50),
       xlab = "Time", ylab = "I")
  lines(rep(ymin, each = lsection))
  
  ## Error in the book: this line is plotted in the left plot in Figure
  ## 3.6, and not as the right plot, as the text on page 20 states.
  minlocs <- seq(lsection/2 + 1, length(x), len = length(ymin))
  bsln.loess <- loess(ymin ~ minlocs)
  lines(predict(bsln.loess, 1:length(x)), lwd = 2)
  
  ## This is the right plot in Figure 3.6
  plot(x, col = "gray", type = "l", ylim = c(20, 50),
       xlab = "Time", ylab = "I")
  lines(asysm(x), lwd = 2)
  
  ## Peak alignment
  data(lcms, package = "ptw")
  plot(lcms[1,,2], type = "l", xlab = "Time", ylab = "I")
  lines(lcms[1,,3], type = "l", col = "gray")
  
  sref <- lcms[1,,2]
  ssamp <- lcms[1,,3]
  lcms.warp <- ptw(sref, ssamp, init.coef = c(0, 1, 0))
  summary(lcms.warp)
  
  plot(time, sref, type = "l", xlim = c(time[600], time[1300]),
       xlab = "Time", ylab = "I", main = "Parametric time warping")
  lines(time, ssamp + 1e6, lty = 2, col = "blue")
  lines(time, lcms.warp$warped.sample + 2e6, col = "green")
  legend("topleft", col = c("black", "blue", "green"), lty = c(1,2,1),
         legend = c("Reference", "Sample", "Warped sample"))
  
  lcms.warpRMS <- ptw(sref, ssamp, optim.crit = "RMS")
  lcms.warpRMS$warp.coef
  
  lcms.warp2 <- ptw(sref, ssamp, init = c(0, 1, 0, 0))
  lcms.warp3 <- ptw(sref, ssamp, init = c(0, 1, 0, 0, 0))
  lcms.warp4 <- ptw(sref, ssamp, init = c(0, 1, 0, 0, 0, 0))
  allwarps <- list(lcms.warp, lcms.warp2, lcms.warp3, lcms.warp4)
  wccs <- sapply(allwarps, function(x) x$crit.value)
  wccs <- round(wccs*1000) / 1000
  allwarp.funs <- sapply(allwarps, function(x) x$warp.fun)
  warpings <- allwarp.funs - 1:length(sref)
  matplot(warpings, type = "l", lty = rep(c(1,2), 2), 
          col = rep(c(1,"blue"), each = 2))
  legend("topleft", lty = rep(c(1,2), 2),
         col = rep(c(1,"blue"), each = 2),
         legend = paste("Degree", 2:5, " - WCC =", wccs))
  
  sref <- lcms[,,2]
  ssamp <- lcms[,,3]
  traces <- select.traces(sref, criterion = "var")
  lcms.warpglobal <- 
    ptw(sref, ssamp, warp.type = "global",
        selected.traces = traces$trace.nrs[1:10])
  summary(lcms.warpglobal)
  
  ## Next lines take a long time to complete...
  ## sample2.indiv.warp <-
  ##   ptw(lcms[,,3], lcms[,,2],
  ##       init.coef = c(0, 1, 0, 0, 0, 0, 0, 0, 0))
  ## sample2.global.warp <-
  ##   ptw(lcms[,,3], lcms[,,2], init.coef = c(0, 1, 0),
  ##       warp.type = "global")   ## error in book: "multiple" should be "global"
  ## sample2.final.warp <-
  ##   ptw(lcms[,,3], lcms[,,2],
  ##       init.coef = c(sample2.global.warp$warp.coef, 0, 0))
  ## plot(time, colSums(lcms[,,3]), col = "green", lwd = 2,
  ##      type = "l", main = "PTW (indiv.)", ylab = "I")
  ## lines(time, colSums(sample2.indiv.warp$warped.sample))
  ## legend("topleft", legend = c("Sample", "Reference"),
  ##        lty = 1, col = c("black", "green"), lwd = c(1,2))

  if (dtw.present) {
    ## dtw
    sref <- lcms[1,,2]
    ssamp <- lcms[1,,3]
    warpfun <- dtw(ssamp, sref, keep.internals = TRUE) # necessary for the plot
    plot(warpfun)
    abline(-20, 1, col = "gray", lty = 2)
    
    ## The code for this plot is not mentioned explicitly in the book
    contour(warpfun$costMatrix[1:500,1:500],
            x=1:500, y=1:500, drawlabels = FALSE,
            xlab="Query index",ylab="Reference index")
    lines(warpfun$index1, warpfun$index2, col="red", lwd=3)
    
    wx2 <- warp(warpfun)
    plot(time, lcms[1,,2], type = "l", xlab = "Time", ylab = "I",
         xlim = c(time[600], time[1300]), main = "Dynamic time warping")
    lines(time, lcms[1,,3] + 1e6, lty = 2, col = "blue")
    lines(time, lcms[1, wx2, 3] + 2e6, col = "green")
    legend("topleft", col = c("black", "blue", "green"), lty = c(1,2,1),
           legend = c("Reference", "Sample", "Warped sample"))
    
    sample2.dtw <- matrix(0, 100, 2000)
    for (i in 1:100) {                 
      warpfun.dtw <- dtw(lcms[i,,3], lcms[i,,2])
      new.indices <- warp(warpfun.dtw, index.reference = FALSE)
      sample2.dtw[i,] <- lcms[i,new.indices,3]
    }                                                                        
    warp.dtw.gl <- dtw(t(lcms[,,3]), t(lcms[,,2]))
    samp.aligned <- lcms[,warp(warp.dtw.gl),3]    
    
    plot(time, colSums(lcms[,,2]), col = "gray", lwd = 3,
         type = "l", main = "DTW (indiv.)", ylab = "I")
    lines(time, colSums(sample2.dtw), col = 1)
    legend("topleft", legend = c("Sample", "Reference"),
           lty = 1, col = c("black", "gray"), lwd = c(1,3))
    
    plot(time, colSums(lcms[,,2]), lwd = 3, col = "gray",
         type = "l", main = "DTW (global)", ylab = "I")
    lines(time, colSums(samp.aligned), col = 1)
    legend("topleft", legend = c("Sample", "Reference"),
           lty = 1, col = c("black", "gray"), lwd = c(1,3))
  }
}

## Peak picking
prostate.mz <- Prostate2000Raw$mz
pks10 <- pick.peaks(rmeans, 10)
plot(prostate.mz[3:498], rmeans, type = "l",
     xlab = "m/z", ylab = "Response")
abline(v = prostate.mz[pks10 + 2], col = "gray")
pks40 <- pick.peaks(rmeans, 40)
plot(prostate.mz[3:498], rmeans, type = "l",
     xlab = "m/z", ylab = "Response", main = "span = 40")
abline(v = prostate.mz[pks40 + 2], col = "gray")

range(apply(prostate[1:10,], 1, max))

range(rowSums(prostate[1:10,]))

prost10.rangesc <- sweep(prostate[1:10,], MARGIN = 1,
                         apply(prostate[1:10,], 1, max),
                         FUN = "/")
apply(prost10.rangesc, 1, max)

range(rowSums(prost10.rangesc))

prost10.lengthsc <- sweep(prostate[1:10,], MARGIN = 1,
                          apply(prostate[1:10,], 1,
                                function(x) sqrt(sum(x^2))),
                          FUN = "/")
range(apply(prost10.lengthsc, 1, max))

range(rowSums(prost10.lengthsc))

prost10.varsc <- sweep(prostate[1:10,], MARGIN = 1,
                       apply(prostate[1:10,], 1, sd),
                       FUN = "/")
range(apply(prost10.varsc, 1, max))

range(rowSums(prost10.varsc))

NIR.mc <- t(sweep(gasoline$NIR, 2, colMeans(gasoline$NIR)))

NIR.mc <- scale(gasoline$NIR, scale = FALSE)
matplot(wavelengths, t(NIR.mc),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "1/R (mean-centered)", lty = 1, col = 1)

data(wines, package = "ChemometricsWithRData")
apply(wines, 2, range)
wines.sc <- scale(wines)
wines.mc <- scale(wines, scale = FALSE)
boxplot(wines.mc ~ col(wines.mc), 
        main = "Mean-centered wine data")
boxplot(wines.sc ~ col(wines.sc), 
        main = "Autoscaled wine data")
