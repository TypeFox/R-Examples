RatioRMS <- function(var_exp1, var_exp2, var_obs, posRMS = 1) {
  #
  #  Enlarge var_exps & var_obs to 10 dim + move posRMS to 1st pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp1)
  for (iind in 1:length(dimsvar)) {
    if (dim(var_exp2)[iind] != dimsvar[iind] | 
        dim(var_obs)[iind] != dimsvar[iind]) {
      stop("all input vars should have the same dimensions")
    }
  } 
  enlvarexp1 <- Enlarge(var_exp1, 10)
  enlvarexp2 <- Enlarge(var_exp2, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  posaperm <- 1:10
  posaperm[1] <- posRMS
  posaperm[posRMS] <- 1
  permvarexp1 <- aperm(enlvarexp1, posaperm)
  permvarexp2 <- aperm(enlvarexp2, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp1)
  #
  #  RMS ratio and its pvalue computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlratiorms <- array(dim = c(2, dimsaperm[2:10]))
  dif1 <- permvarexp1 - permvarobs
  dif2 <- permvarexp2 - permvarobs
  rms1 <- Mean1Dim(dif1 ** 2, 1, narm = TRUE) ** 0.5
  rms2 <- Mean1Dim(dif2 ** 2, 1, narm = TRUE) ** 0.5
  rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(
    abs(rms2), na.rm = TRUE) / 1000
  enlratiorms[1, , , , , , , , , ] <- (rms1 / rms2)
  eno1 <- Eno(dif1, 1)
  eno2 <- Eno(dif2, 1)
  F <- (eno1 * (rms1) ** 2 / (eno1 - 1)) / (eno2 * (rms2) ** 2 / (eno2 - 1))
  F[which(F < 1)] <- 1 / F[which(F < 1)]
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    l1 <- eno1[j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    l2 <- eno2[j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (is.na(l1) == FALSE & is.na(l2) == FALSE & l1 > 2 & l2 > 2) {
                      enlratiorms[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                        j10] <- (1 - pf(F[j2, j3, j4, j5, j6, j7, j8, j9, j10], 
                                 l1 - 1, l2 - 1)) * 2
                    } else {
                      enlratiorms[1, j2, j3, j4, j5, j6, j7, j8, j9, j10] <- NA
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
  
  dimsvar[posRMS] <- 2
  ratiorms <- array(dim = dimsvar)
  ratiorms[] <- aperm(enlratiorms, posaperm)
  #  
  #  Output
  # ~~~~~~~~
  #
  ratiorms
}
