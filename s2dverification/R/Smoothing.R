Smoothing <- function(var, runmeanlen = 12, numdimt = 4) {
  # 
  #  Enlarge the number of dimensions of var to 10 --> enlvar 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  u <- IniListDims(dimsvar, 10)
  enlvar <- Enlarge(var, 10)
  
  #
  #  Smoothing 
  # ~~~~~~~~~~~
  smt_enlano <- array(dim = c(dimsvar, array(1, dim = (10 - length(dimsvar)))))

  nmr1 <- floor(runmeanlen / 2)
  nltime <- dimsvar[numdimt]
  if (nmr1 > 0) { 
    for (jtime in (1+nmr1):(nltime - nmr1)) {
      # First, averaging the two extreme values
      u[[numdimt]] <- jtime - nmr1
      left <- enlvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]], 
                   u[[8]], u[[9]], u[[10]]]
      u[[numdimt]] <- jtime + nmr1
      right <- enlvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]],
                    u[[8]], u[[9]], u[[10]]]
      u[[numdimt]] <- jtime
      smt_enlano[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]], u[[8]],
               u[[9]], u[[10]]] <- (left + right) / (2 * (2 * nmr1))
      for (k in - (nmr1 - 1):(nmr1 - 1)) {
        # Second, adding the equally-weighted values around the centered
        u[[numdimt]] <- jtime + k
        mid <- enlvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]],
                      u[[8]], u[[9]], u[[10]]]
        u[[numdimt]] <- jtime
        smt_enlano[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]], 
                   u[[8]], u[[9]], u[[10]]] <- smt_enlano[u[[1]], u[[2]],
                                               u[[3]], u[[4]], u[[5]], u[[6]],
                                               u[[7]], u[[8]], u[[9]], 
                                               u[[10]]] + (mid / (2 * nmr1))
      }
    }
  }else{
    smt_enlano <- enlvar
  }
  #
  #  Reduce the number of dimensions to the original one 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  smt_ano <- array(dim = dimsvar)
  smt_ano[] <- smt_enlano
  #
  #  Outputs
  # ~~~~~~~~~
  #
  smt_ano
}
