RMS <- function(var_exp, var_obs, posloop = 1, posRMS = 2, compROW = NULL, 
                limits = NULL) {       
  #
  #  Remove data along compROW dim if there is at least one NA between limits
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (is.null(compROW) == FALSE) {
    if (is.null(limits) == TRUE) {
      limits <- c(1, dim(var_obs)[compROW])
    }
    outrows <- (is.na(Mean1Dim(var_obs, compROW, narm = FALSE, limits))) 
    outrows <- InsertDim(outrows, compROW, dim(var_obs)[compROW])
    var_obs[which(outrows)] <- NA
  }
  
  #
  #  Enlarge var_exp & var_obs to 10 dim + move posloop & posRMS to 1st & 2nd 
  #  pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind != posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
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
  #  RMS & its confidence interval computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlrms <- array(dim = c(nexp, nobs, 3, dimsaperm[3:10]))
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      dif <- array(dim = dimsaperm[-1])
      dif[, , , , , , , , ] <- permvarexp[jexp, , , , , , , , 
                               , ] - permvarobs[jobs, , , , , , , , , ]
      enlrms[jexp, jobs, 2, , , , , , , , ] <- Mean1Dim(dif ** 2, 1, 
                                               narm = TRUE) ** 0.5
      eno <- Eno(dif, 1)
      for (j3 in 1:dimsaperm[3]){
        for (j4 in 1:dimsaperm[4]){
          for (j5 in 1:dimsaperm[5]){
            for (j6 in 1:dimsaperm[6]){
              for (j7 in 1:dimsaperm[7]){
                for (j8 in 1:dimsaperm[8]){
                  for (j9 in 1:dimsaperm[9]){
                    for (j10 in 1:dimsaperm[10]){
                      ndat <- length(sort(dif[, j3, j4, j5, j6, j7, j8, j9, 
                              j10]))
                      enlrms[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- (eno[j3, j4, j5, j6, j7, j8, j9, 
                             j10] * enlrms[jexp, jobs, 2, j3, j4, j5, j6, j7, 
                             j8, j9, j10] ** 2 / qchisq(0.975, eno[j3, j4, j5, 
                             j6, j7, j8, j9, j10] - 1)) ** 0.5
                      enlrms[jexp, jobs, 3, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- (eno[j3, j4, j5, j6, j7, j8, j9, 
                             j10] * enlrms[jexp, jobs, 2, j3, j4, j5, j6, j7, 
                             j8, j9, j10] ** 2 / qchisq(0.025, eno[j3, j4, j5, 
                             j6, j7, j8, j9, j10] - 1)) ** 0.5
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

  rms <- array(dim = c(nexp, nobs, 3, dimsvar[-c(posloop, posRMS)]))
  rms[] <- enlrms
  #  
  #  Output
  # ~~~~~~~~
  #
  rms
}
