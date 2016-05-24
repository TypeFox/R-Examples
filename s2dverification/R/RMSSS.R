RMSSS <- function(var_exp, var_obs, posloop = 1, posRMS = 2) {
  #
  #  Enlarge var_exp & var_obs & clim to 10 dim + move posloop & posRMS  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind !=posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
      stop("var_exp & var_obs must have same dimensions except along posloop")
    }
  }
  if (dimsvar[posRMS] < 2 ) {
    stop("At least 2 values required to compute RMSE")
  }
  enlvarexp <- Enlarge(var_exp, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  nexp <- dimsvar[posloop]
  nobs <- dim(var_obs)[posloop]
  posaperm <- numeric(10)
  posaperm[1] <- posloop
  posaperm[2] <- posRMS
  posaperm[3:10] <- seq(1, 10)[-c(posloop, posRMS)]
  permvarexp <- aperm(enlvarexp, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp)
  #
  #  RMSSS and its pvalue computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlRMSSS <- array(dim = c(nexp, nobs, 2, dimsaperm[3:10]))

  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      dif1 <- array(dim = dimsaperm[-1])
      dif2 <- array(dim = dimsaperm[-1])
      dif1[, , , , , , , , ] <- permvarexp[jexp, , , , , , , , 
                                , ] - permvarobs[jobs, , , , , , , , , ]
      dif2[, , , , , , , , ] <- permvarobs[jobs, , , , , , , , , ]
      rms1 <- Mean1Dim(dif1 ** 2, 1, narm = TRUE) ** 0.5
      rms2 <- Mean1Dim(dif2 ** 2, 1, narm = TRUE) ** 0.5
      rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(abs(
        rms2), na.rm = TRUE) / 1000
      enlRMSSS[jexp, jobs, 1, , , , , , , , ] <- 1 - (rms1 / rms2)
      eno1 <- Eno(dif1, 1)
      eno2 <- Eno(dif2, 1)
      F <- (eno2 * (rms2) ** 2 / (eno2 - 1)) / (eno1 * (rms1) ** 2 / (eno1 - 1))
      for (j3 in 1:dimsaperm[3]) {
        for (j4 in 1:dimsaperm[4]) {
          for (j5 in 1:dimsaperm[5]) {
            for (j6 in 1:dimsaperm[6]) {
              for (j7 in 1:dimsaperm[7]) {
                for (j8 in 1:dimsaperm[8]) {
                  for (j9 in 1:dimsaperm[9]) {
                    for (j10 in 1:dimsaperm[10]) {
                      l1 <- eno1[j3, j4, j5, j6, j7, j8, j9, j10]
                      l2 <- eno2[j3, j4, j5, j6, j7, j8, j9, j10]
                      if (is.na(l1) == FALSE & is.na(l2) == FALSE & l1 > 2 & l2 > 2) {
                        enlRMSSS[jexp, jobs, 2, j3, j4, j5, j6, j7, j8, j9, 
                                 j10] <- 1 - pf(F[j3, j4, j5, j6, j7, j8, j9, 
                                                j10], l1 - 1, l2 - 1)
                      } else {
                        enlRMSSS[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9,
                                 j10] <- NA
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  RMSSS <- array(dim = c(nexp, nobs, 2, dimsvar[-c(posloop, posRMS)]))
  RMSSS[] <- enlRMSSS
  #
  #  Output
  # ~~~~~~~~
  #
  RMSSS
}
