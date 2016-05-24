Eno <- function(obs, posdim) {
  dimsvar <- dim(obs)
  if (is.null(dimsvar)) {
    dimsvar <- length(obs)
  }
  enlobs <- Enlarge(obs, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  posaperm <- 1:10
  posaperm[posdim] <- 1
  posaperm[1] <- posdim
  enlobs <- aperm(enlobs, posaperm)
  dimsaperm <- outdim[posaperm]
  #
  #  Loop on all dimensions to compute effective number of observations 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enleno <- array(dim = c(1, dimsaperm[2:10]))
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    tmp <- enlobs[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (length(sort(tmp)) > 1 ) {
                      n <- length(sort(tmp))
                      a <- acf(tmp, lag.max = n - 1, plot = FALSE, 
                               na.action = na.pass)$acf[2:n, 1, 1]
                      s <- 0
                      for (k in 1:(n - 1)) {
                        s <- s + (((n - k) / n) * a[k])
                      }
                      enleno[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- min(n / (1 + (2 * s)), n)
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
  #dimsvar <- dimsvar[-posdim]
  if (length(dimsvar) == 1) {
    dimsvar <- 1
  } else {
    dimsvar <- dimsvar[-posdim]
  }
  effnumobs <- array(dim = dimsvar)
  effnumobs[] <- aperm(enleno, posaperm)
  #
  #  Outputs
  # ~~~~~~~~~
  #
  effnumobs
}
