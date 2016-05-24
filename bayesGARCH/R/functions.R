## LIBRARY garchLib
## functions.bayes
## David Ardia

## formSmpl
## __input__
## MCMC : and object of the class MCMC.list or a list of matrices containing the MCMC output
## l.bi = 0 : the lenth of the burn in phase
## batch.size : size of the batching smpl
## __output__
## r : a matrix of the posterior sample
## formSmpl
## __input__
## MCMC : and object of the class MCMC.list or a list of matrices containing the MCMC output
## l.bi = 0 : the lenth of the burn in phase
## batch.size = 1 : size of the batching smpl
## __output__
## r : a matrix of the posterior sample (object of class mcmc)
"formSmpl" <- function(MCMC, l.bi = 0, batch.size = 1){
  if (missing(MCMC))
    stop ("'MCMC' is missing")
  if (!is.list(MCMC)){
    if (is.matrix(MCMC))
      MCMC <- list(MCMC)
    else
      stop ("'MCMC' is neither a list or a matrix")
  }
  if (!is.vector(l.bi) || length(l.bi) != 1)
    stop ("'l.bi' must be a scalar")
  if (l.bi < 0)
    stop ("'l.bi' must be positive")
  if (!is.vector(batch.size) || length(l.bi) != 1)
    stop ("'batch.size' must be a scalar")
  if (batch.size < 1)
    stop ("'batch.size' must be larger or equal than 1")
  r <- NULL
  n.chain <- length(MCMC)
  l.chain <- nrow(MCMC[[1]])
  if (l.bi >= l.chain)
    stop ("'l.bi >= l.chain'")
  theM <- ((l.bi + 1):l.chain)
  seqk <- seq(from = 1, to = length(theM), by = batch.size)
  theM <- theM[seqk]
  for (i in 1:n.chain)
    r <- rbind(r, MCMC[[i]][theM, ])
  cat("\nn.chain: ", n.chain, "\nl.chain: ", l.chain, "\nl.bi: ", l.bi, "\nbatch.size: ",
      batch.size, "\nsmpl size: ", nrow(r), "\n")
  dimnames(r) <- list(1:nrow(r), colnames(MCMC[[1]]))

  mcmc(r, start = 1)
}
