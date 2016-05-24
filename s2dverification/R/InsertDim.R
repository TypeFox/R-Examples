InsertDim <- function(var, posdim, lendim) {
  #
  #  Initialisation of the output var with the required dimension length 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  outdim <- lendim
  if (posdim > 1) {
    outdim <- c(dimsvar[1:(posdim - 1)], outdim)
  }
  if (posdim <= length(dimsvar)) {
    outdim <- c(outdim, dimsvar[posdim:length(dimsvar)])
  }
  tmpvar <- array(dim = c(outdim, array(1, dim = (10 - length(outdim)))))
  #
  #  Duplicate the matrix along the required (posdim)th dimension
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  u <- IniListDims(outdim, 10)
  for (jindex in 1:lendim) {
    u[[posdim]] <- jindex
    tmpvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]], u[[8]], 
           u[[9]], u[[10]]] <- var
  }
  #
  #  Reduce the number of dimensions to the required one
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  outvar <- array(dim = outdim)
  outvar[] <- tmpvar
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar
}
