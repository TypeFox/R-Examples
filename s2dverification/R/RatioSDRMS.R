RatioSDRMS <- function(var_exp, var_obs) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 7 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 4 | length(dimobs) < 4) {
    stop("At least 4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)")
  }
  for (jn in 3:max(length(dimexp), length(dimobs))) {
    if (dimexp[jn] != dimobs[jn]) {
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 7)
  var_obs <- Enlarge(var_obs, 7)
  
  #
  #  Ratio RMSE / SD and its significance level
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ensmeanexp <- Mean1Dim(var_exp, 2)
  ensmeanobs <- Mean1Dim(var_obs, 2)
  dimrms <- c(dimexp[1], dimobs[1], dimexp[4:length(dimexp)])
  dimratiormssd <- c(dimexp[1], dimobs[1], 2, dimexp[4:length(dimexp)])
  if (length(dimrms) < 6) {
    dimrms <- c(dimrms, array(1, dim = (6 - length(dimrms))))
  }
  if (length(dimratiormssd) < 7) {
    dimenlratiormssd <- c(dimratiormssd, 
                          array(1, dim = (7 - length(dimratiormssd))))
  } else {
    dimenlratiormssd <- dimratiormssd
  }
  dif <- var_exp - InsertDim(ensmeanexp, 2, dimexp[2])
  std <- apply(array(dif, dim = c(dimexp[1], dimexp[2] * dimexp[3],
               dimrms[3:6])), c(1, 3, 4, 5, 6), sd, na.rm = TRUE)
  enosd <- apply(Eno(dif, 3), c(1, 3, 4, 5, 6), sum, na.rm = TRUE)
  rms <- array(dim = dimrms)
  enlratiormssd <- array(dim = dimenlratiormssd)
  for (jexp in 1:dimexp[1]) {
    for (jobs in 1:dimobs[1]) {
      dif <- ensmeanexp[jexp, , , , , ] - ensmeanobs[jobs, , , , , ]
      rms[jexp,jobs, , , , ] <- Mean1Dim(dif ** 2, 1, narm = TRUE) ** 0.5
      enorms <- array(Eno(dif, 1), dim = dimrms[3:6])
      enlratiormssd[jexp, jobs, 1, , , , ] <- std[jexp, , , , ] / rms[jexp, 
                                              jobs, , , , ]
      for (jltime in 1:dimrms[3]) {
        for (jlev in 1:dimrms[4]) {
          for (jlat in 1:dimrms[5]) {
            for (jlon in 1:dimrms[6]) {
              l1 <- enosd[jexp, jltime, jlev, jlat, jlon]
              l2 <- enorms[jltime, jlev, jlat, jlon]
              F <- (enosd[jexp, jltime, jlev, jlat, jlon] * (std[jexp, jltime, 
                   jlev, jlat, jlon]) ** 2 / (enosd[jexp, jltime, jlev, jlat, 
                   jlon] - 1)) / (enorms[jltime, jlev, jlat, jlon] * (rms[jexp, 
                   jobs, jltime, jlev, jlat, jlon]) ** 2 / (enorms[jltime, 
                   jlev, jlat, jlon] - 1))
              if (is.na(F) == FALSE & is.na(l1) == FALSE & is.na(l2) == FALSE & l1 > 2 & l2 > 2) {
                enlratiormssd[jexp, jobs, 2, jltime, jlev, jlat, 
                              jlon] <- 1 - pf(F, l1 - 1, l2 - 1)
              } else {
                enlratiormssd[jexp, jobs, 1, jltime, jlev, jlat, jlon] <- NA
              }
            }
          }
        }
      }
    }
  }
  ratiosdrms <- array(dim = dimratiormssd)
  ratiosdrms[] <- enlratiormssd
  #  
  #  Output
  # ~~~~~~~~
  #
  ratiosdrms
}
