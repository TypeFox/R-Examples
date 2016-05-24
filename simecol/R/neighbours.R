## Basic Functions for Rectangular Cellular Automata
## Details: see documentation of package simecol

eightneighbours <- function(x){
  if (!is.matrix(x))    stop("x must be a matrix")
  n <- dim(x)[1]
  m <- dim(x)[2]
  y <- rep(0, length(x))
  z <- .C("eightneighbours", as.integer(n), as.integer(m),
          as.double(x), y=as.double(y), PACKAGE="simecol")$y
  dim(z) <- c(n, m)
  z
}

neighbours <- function(x, state = NULL, wdist = NULL, tol = 1e-4, bounds = 0){
  if (!is.matrix(x))    stop("x must be a matrix")
  if (!is.null(state) & !is.numeric(state)) stop("state must be numeric or NULL")
  if (!is.null(wdist) & !is.numeric(wdist)) stop("wdist must be numeric or NULL")
  if (!is.numeric(tol)) stop("tol must be numeric")

  #if (length(bounds) %% 4)
  #  warning("length of 'bounds' argument must be either one or four")
  bounds <- rep(bounds, length.out = 4)
  ## pack this into an integer bit mask
  bound <- sum(bounds * c(1L, 2L, 4L, 8L))

  n <- dim(x)[1]
  m <- dim(x)[2]

  y <- rep(0, length(x))
  ## if wdist not given do the same as eightneighbours
  if (is.null(wdist)) wdist <- matrix(c(1,1,1,1,0,1,1,1,1), nrow=3)

  ndist <- dim(wdist)[1]
  mdist <- dim(wdist)[2]
  if (mdist != ndist) stop ("wdist must be a sqare matrix")

  if ((ndist > n) || (ndist > m))
    stop("dimensions of weight matrix must be smaller dimensions of state matrix")

  ## default: all nonzero states in matrix counted
  ## we simply set all nonzero states to 1 and check against 1
  if (is.null(state)) {
    state <- 1
    x[x != 0] <- 1
  }

  if (bound == 0) {
    ## sligtly faster version
    z <- .C("neighbours", as.integer(n), as.integer(m),
            as.double(x), y = as.double(y),
            as.integer(ndist), as.double(wdist),
            as.double(state[1]), as.double(tol[1]),
            PACKAGE = "simecol")$y
  } else {
    ## more general version
    z <- .C("xneighbours", as.integer(n), as.integer(m),
            as.double(x), y = as.double(y),
            as.integer(ndist), as.double(wdist),
            as.double(state[1]), as.double(tol[1]),
            as.integer(bound),
            PACKAGE = "simecol")$y
  }
  dim(z) <- c(n, m)
  z
}

## aliases
eightneighbors <- eightneighbours
neighbors      <- neighbours
