writeAverage <- function(filename,
                         ntimes, nwave, scans,
                         fileout = paste(filename, "Average.ivo", sep=""),
                         calibration = 1:nwave, wexplicit= FALSE) {

  s1 <- as.matrix(read.table(filename))
  cat("Read", filename, "\n")
  s19 <- list()
  cnt <- 1
  for(i in 1:scans){
    s19[[i]] <- s1[cnt:(cnt+ntimes-1), ]
    cnt <- cnt + ntimes 
  }
  sm <- s19[[1]]
  if(scans > 1) {
    for(i in 2:scans)
      sm <- sm + s19[[i]]
  }
  sm <- sm / scans
   cat("Computed average of", scans, "scans", "each having",
       nwave, "wavelengths and", ntimes, "times \n")
  times <- sm[,1]
  sm <- sm[,-1]
  
  zz <- file(fileout, "w")  
  cat("\n\n", file = zz)
  
  if(wexplicit){
    cat("Wavelength explicit\n", file = zz)
    cat(paste("Intervalnr  ",nwave,"\n", sep=""), file=zz)
    write.table(sm, file = zz, 
                col.names=calibration, row.names=times, quote=FALSE)
    fstr <- "in the 'wavelength explicit' format"
  }
  else {
    cat("Time explicit\n", file = zz)
    cat(paste("Intervalnr  ",ntimes,"\n", sep=""), file=zz)
    write.table(t(sm), file = zz, 
                row.names=calibration, col.names=times, quote=FALSE)
    fstr <- "in the 'time explicit' format"
  }
  
  close(zz)
  cat("Wrote the averaged file", fileout, fstr, "\n")
}
# either copy this file into R or use
# source("writeAverage.R")
# in R, when the file is in your working directory
#
# call the function with for example:
#
# writeAverage(filename = "rc682nm10S1Scans9.dat", ntimes = 275, nwave = 256,
#              scans = 9, fileout="xx.ivo")
#
# or just (because it's not necessary to give names of the arguments
#
# writeAverage("rc682nm10S1Scans9.dat", 275, 256, 9, "xx.ivo")
#
# Note that if 'fileout' is not given, then a file name is made automatically
#
# Note also that if you want the "ivo" file that is written to have a
# calibration in it (that is, if you know what the wavelength labels should
# be), then you can add the argument "calibration" - for example:
#
#  writeAverage(filename = "rc682nm10S1Scans9.dat", ntimes = 275, nwave = 256,
#              scans = 9, fileout="xx.ivo",
#              calibration = 425.7082 + 1.1849 * 1:256)
# 
# or just
#  writeAverage("rc682nm10S1Scans9.dat", 275, 256, 9, "xx.ivo",
#              425.7082 + 1.1849 * 1:256)
#
# to get an output file in wavelength explicit format, set wexplicit=TRUE,
# as in, for example:
#
#  writeAverage(filename = "rc682nm10S1Scans9.dat", ntimes = 275, nwave = 256,
#              scans = 9, fileout="xx.ivo",
#              calibration = 425.7082 + 1.1849 * 1:256,
#              wexplicit = TRUE)
