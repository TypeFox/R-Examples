Season <- function(var, posdim = 4, monini, moninf, monsup) {
  while (monsup < moninf) {
    monsup <- monsup + 12
  }
  #
  #  Enlarge the size of var to 10 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  ntime <- dimsvar[posdim]
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  u <- IniListDims(outdim, 10)
  v <- IniListDims(outdim, 10)
  #
  #  Compute the seasonal means 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ind <- 1:ntime
  months <- ((ind - 1) + monini - 1) %% 12 + 1
  years <- ((ind - 1) + monini - 1) %/% 12

  for (jmon in moninf:monsup) {
    u[[posdim]] <- ind[which(months == ((jmon - 1) %% 12 + 1))]
    ind0 <- u[[posdim]][1]
    indf <- u[[posdim]][length(u[[posdim]])]
    if (indf > (ntime - (monsup - jmon))) {
      u[[posdim]] <- u[[posdim]][-which(u[[posdim]] == indf)]
    }
    if (ind0 < (jmon - moninf + 1)) {
      u[[posdim]] <- u[[posdim]][-which(u[[posdim]] == ind0)]
    } 
    if (jmon == moninf) { 
      nseas <- length(u[[posdim]])
      dimsvar[posdim] <- nseas
      outdim[posdim] <- nseas
      enlvarout <- array(0, dim = outdim)
    }
    v[[posdim]] <- 1:nseas
    enlvarout[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]], v[[7]], v[[8]], 
              v[[9]], v[[10]]] <- enlvarout[v[[1]], v[[2]], v[[3]], v[[4]],
                                            v[[5]], v[[6]], v[[7]], v[[8]], 
                                            v[[9]], v[[10]]] + enlvar[u[[1]], 
                                  u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], 
                                  u[[7]], u[[8]], u[[9]], u[[10]]]
  }
  varout <- array(dim = dimsvar)
  varout[] <- enlvarout
  varout <- varout / (monsup - moninf + 1)
  #
  #  Outputs
  # ~~~~~~~~~
  #
  varout
}
