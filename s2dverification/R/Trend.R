Trend <- function(var, posTR = 2, interval = 1) {
  # 
  #  Enlarge the size of var to 10 and move posTR to first position 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  outdim[posTR] <- 4
  posaperm <- 1:10
  posaperm[posTR] <- 1
  posaperm[1] <- posTR
  enlvar <- aperm(enlvar, posaperm)
  dimsaperm <- outdim[posaperm]
  # 
  #  Loop on all dimensions to compute trends 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  enltrend <- array(dim = dimsaperm)
  enldetrend <- array(dim = dim(enlvar))
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    tmp <- enlvar[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (length(sort(tmp)) > 0) {
                      mon <- seq(tmp) * interval
                      lm.out <- lm(tmp ~ mon, na.action = na.omit)
                      enltrend[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- confint(lm.out)[2, 1]
                      enltrend[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[2]
                      enltrend[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- confint(lm.out)[2, 2]
                      enltrend[4, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[1]
                      enldetrend[is.na(tmp) == FALSE, j2, j3, j4, j5, j6, j7, j8, 
                        j9, j10] <- tmp[is.na(tmp) == FALSE] - lm.out$fitted.values
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
  #  Back to the original dimensions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  detrended <- array(dim = dimsvar)
  dimsvar[posTR] <- 4
  trend <- array(dim = dimsvar)
  trend[] <- aperm(enltrend, posaperm)
  detrended[] <- aperm(enldetrend, posaperm)
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(trend = trend, detrended = detrended))
}
