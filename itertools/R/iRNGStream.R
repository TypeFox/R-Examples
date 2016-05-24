iRNGStream <- function(seed) {
  # Convert a single number into the appropriate vector for "L'Ecuyer-CMRG"
  if (length(seed) == 1) {
    seed <- convseed(seed)
  }

  # Error checking: this will throw an error right away if the seed is bad
  nextRNGStream(seed)

  # Define the "Next Element" function for the iterator
  nextEl <- function() (seed <<- nextRNGStream(seed))

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

iRNGSubStream <- function(seed) {
  # Convert a single number into the appropriate vector for "L'Ecuyer-CMRG"
  if (length(seed) == 1) {
    seed <- convseed(seed)
  }

  # Error checking: this will throw an error right away if the seed is bad
  nextRNGSubStream(seed)

  # Define the "Next Element" function for the iterator
  nextEl <- function() (seed <<- nextRNGSubStream(seed))

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

convseed <- function(iseed) {
  saveseed <- if (exists('.Random.seed', where=.GlobalEnv, inherits=FALSE))
    get('.Random.seed', pos=.GlobalEnv, inherits=FALSE)

  saverng <- RNGkind("L'Ecuyer-CMRG")

  tryCatch({
    set.seed(iseed)
    get('.Random.seed', pos=.GlobalEnv, inherits=FALSE)
  },
  finally={
    RNGkind(saverng[1], saverng[2])
    if (is.null(saveseed))
      rm('.Random.seed', pos=.GlobalEnv)
    else
      assign('.Random.seed', saveseed, pos=.GlobalEnv)
  })
}
