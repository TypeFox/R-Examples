#' Returns predictor vector for design matrix
#' @description Returns predictor vector for design matrix from 44 astronomical angular velocities.
#' @param xi Transit index
#' @param ma Max number of predictors.
#' @param ivar Dummy.
#' @param tdiff Length of input time series.
#' @return The predictor vector. Values between -1, 1.


Funcs <- function (xi, ma = 89, ivar = 1, tdiff) {

  rad <- 0.017453292519943
  xi  <- rad * xi

  omegas     <- vector()
  omegas[1]  <- 0.054809904
  omegas[2]  <- 0.115308512
  omegas[3]  <- 0.904885870
  omegas[4]  <- 1.020194382
  omegas[5]  <- 1.809771741
  omegas[6]  <- 2.040388764
  omegas[7]  <- 11.597841752
  omegas[8]  <- 11.713150263
  omegas[9]  <- 13.468112100
  omegas[10] <- 13.522922004
  omegas[11] <- 13.583420612
  omegas[12] <- 13.638230516
  omegas[13] <- 13.693040419
  omegas[14] <- 15.563310768
  omegas[15] <- 23.426300526
  omegas[16] <- 24.215877885
  omegas[17] <- 25.181262364
  omegas[18] <- 25.236072267
  omegas[19] <- 25.290882171
  omegas[20] <- 25.351380779
  omegas[21] <- 27.045844008
  omegas[22] <- 27.161152519
  omegas[23] <- 27.221651128
  omegas[24] <- 27.276461031
  omegas[25] <- 27.331270935
  omegas[26] <- 36.949222530
  omegas[27] <- 37.738799889
  omegas[28] <- 38.704184367
  omegas[29] <- 38.758994271
  omegas[30] <- 38.813804174
  omegas[31] <- 38.874302783
  omegas[32] <- 38.989611294
  omegas[33] <- 40.799383035
  omegas[34] <- 49.451950152
  omegas[35] <- 50.472144534
  omegas[36] <- 52.281916275
  omegas[37] <- 52.512533298
  omegas[38] <- 54.552922062
  omegas[39] <- 62.185294797
  omegas[40] <- 63.995066538
  omegas[41] <- 66.035455302
  omegas[42] <- 75.708216801
  omegas[43] <- 77.748605565
  omegas[44] <- 100.944289068

  afunc              <- vector()
  afunc[1]           <- 1
  rayleigh.criterion <- 0

  for(i in seq(2,ma,2)) {
    afunc[i]     <- cos(omegas[i / 2] * xi)
    afunc[i + 1] <- sin(omegas[i / 2] * xi)
  }

  hh <- 0
  for(h in 0 : 44){
    if(h < hh) next
    for(hh in ((h + 1) : 44)) {
      if(hh > 44) break
      if(h == 0) {
        rayleigh.criterion <- omegas[hh] * tdiff
      } else {
        rayleigh.criterion <- (omegas[hh] - omegas[h]) * tdiff
      }
      if((rayleigh.criterion > 360)) { break
      } else {
        afunc[2 * hh]     <- 0
        afunc[2 * hh + 1] <- 0
      }
    }
  }

  return(afunc)
}