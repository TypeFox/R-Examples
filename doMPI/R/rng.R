#
# This function follows the outline presented in section 6 of
# the vignette for the "parallel" package written by R-Core.
#
setRngDoMPI <- function(cl, seed=NULL) {
  # save the current value of .Random.seed so it can be restored
  saveseed <- if (exists('.Random.seed', where=globalenv(), inherits=FALSE))
    get('.Random.seed', pos=globalenv(), inherits=FALSE)

  # set RNG to L'Ecuyer in order to generate .Random.seed values
  # to send to the workers, saving the previous value
  saverng <- RNGkind()
  if (saverng[1] != "L'Ecuyer-CMRG")
    RNGkind("L'Ecuyer-CMRG")

  tryCatch({
    # call set.seed if seed is not NULL
    s <- if (! is.null(seed)) {
      set.seed(seed)
      get('.Random.seed', pos=globalenv(), inherits=FALSE)
    } else {
      nextRNGStream(get('.Random.seed', pos=globalenv(), inherits=FALSE))
    }

    # send a .Random.seed value to each worker in the cluster
    for (i in seq(length=clusterSize(cl))) {
      sendToWorker(cl, i, list(seed=s))
      s <- nextRNGStream(s)
    }
  },
  finally={
    # restore the local RNG and .Random.seed
    RNGkind(saverng[1], saverng[2])
    if (is.null(saveseed))
      rm('.Random.seed', pos=globalenv())
    else
      assign('.Random.seed', saveseed, pos=globalenv())
  })

  invisible(NULL)
}
