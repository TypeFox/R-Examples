Corr <- function(var_exp, var_obs, posloop = 1, poscor = 2, compROW = NULL, 
                 limits = NULL, siglev = 0.95, method = 'pearson') {
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
  #  Enlarge var_exp & var_obs to 10 dim + move posloop & poscor to 1st & 2nd 
  #  pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind != posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
      stop("var_exp & var_obs must have same dimensions except along posloop")
    }
  }
  if (dimsvar[poscor] < 3 ) {
    stop("At least 3 values required to compute correlation")
  }
  if (method != "kendall" && method != "spearman" && method != "pearson") {
    stop("Wrong correlation method")
  }
  enlvarexp <- Enlarge(var_exp, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  nexp <- dimsvar[posloop]
  nobs <- dim(var_obs)[posloop]
  posaperm <- numeric(10)
  posaperm[1] <- posloop
  posaperm[2] <- poscor
  posaperm[3:10] <- seq(1, 10)[-c(posloop, poscor)]
  permvarexp <- aperm(enlvarexp, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp)
  #

  # Check the siglev arguments:
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (siglev > 1 | siglev < 0) {
    stop("siglev need to be higher than O and lower than 1")
  }
  #
  #  Loop to compute correlation for each grid point
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlCORR <- array(dim = c(nexp, nobs, 4, dimsaperm[3:10]))
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      for (j3 in 1:dimsaperm[3]) {
        for (j4 in 1:dimsaperm[4]) {
          for (j5 in 1:dimsaperm[5]) {
            for (j6 in 1:dimsaperm[6]) {
              for (j7 in 1:dimsaperm[7]) {
                for (j8 in 1:dimsaperm[8]) {
                  for (j9 in 1:dimsaperm[9]) {
                    for (j10 in 1:dimsaperm[10]) {
                      tmp1 <- permvarexp[jexp, , j3, j4, j5, j6, j7, j8, j9, 
                                         j10]
                      tmp2 <- permvarobs[jobs, , j3, j4, j5, j6, j7, j8, j9,
                                         j10]
                      if (length(sort(tmp1)) > 0 & length(sort(tmp2)) > 2) {
                        toto <- cor(tmp1, tmp2, use = "pairwise.complete.obs", method=method)
                        enlCORR[jexp, jobs, 2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- toto
                        #eno <- min(Eno(tmp2, 1), Eno(tmp1, 1))
                        if(method=="kendall" | method=="spearman"){
                          eno <- Eno(rank(tmp2), 1)                          
                        } else if(method == "pearson"){
                          eno <- Eno(tmp2, 1)                          
                        } 
                        #t <- qt(0.95, eno - 2)
                        t <- qt(siglev, eno - 2)
                        enlCORR[jexp, jobs, 4, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- sqrt((t * t) / ((t * t) + eno - 2))
                        enlCORR[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- tanh(atanh(toto) + qnorm(1-(1-siglev)/2) / sqrt(
                                #j10] <- tanh(atanh(toto) + qnorm(0.975) / sqrt(
                                        eno - 3))
                        enlCORR[jexp, jobs, 3, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- tanh(atanh(toto) + qnorm((1-siglev)/2) / sqrt(
                                #j10] <- tanh(atanh(toto) + qnorm(0.025) / sqrt(
                                        eno - 3))
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
  #
  CORR <- array(dim = c(nexp, nobs, 4, dimsvar[-c(posloop, poscor)]))
  CORR[] <- enlCORR
  #  
  #  Output
  # ~~~~~~~~
  #
  CORR
}
