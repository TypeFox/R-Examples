#############################################################
#
#	localdepth.simp.exact.hyperspheres function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: October, 26, 2011
#	Version: 0.2
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.simp.exact.hyperspheres <- function(x, y=NULL, tau, use=c('volume', 'diameter'), nsamp='all', nmax=1, tol=0) {
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, ncol=1)
  if (is.vector(y))
    y <- matrix(y, ncol=1)
  x <- as.matrix(x)
  y <- as.matrix(y)  
  normalize <- function(x) x/sqrt(x%*%x)
  x <- t(apply(x, 1, normalize))
  y <- t(apply(y, 1, normalize))
  nx <- nrow(x)
  ny <- nrow(y)
  nc <- ncol(x)
  if (ncol(y)!=nc)
    stop("the number of columns in 'x' and 'y' must be the same")
  nt <- choose(nx, nc)
  if (nt > .Machine$integer.max)
    nt <- .Machine$integer.max
  use <- match.arg(use)
  if (use=='volume' & nc!=3) {
    warning("The option use='volume' is available only on the sphere, we use 'diameter' instead")
    use <- "diameter"
  }
  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp) && nsamp > nt) {
      warning("Since 'nsamp' is greater than the number of simplex the 'all' method is used")
      nsamp <- 'all'
  }
  if (use=='diameter') nuse=0
  if (use=='volume') nuse=1
  
  if (is.character(nsamp) && nsamp=='all') {
    z <- .Fortran("ldsehs",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nx),
      as.integer(ny),
      as.integer(nuse),
      as.double(tol),
      depth=double(ny),
      localdepth=double(ny),
      PACKAGE = "localdepth")
  } else if (is.numeric(nsamp)){
    nt <- round(nt*nmax)
    z <- .Fortran("ldseahs",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nsamp),
      as.integer(nx),
      as.integer(ny),
      as.integer(nuse),
      as.double(tol),
      depth=double(ny),
      localdepth=double(ny),                  
      nd=double(1),
      nld=double(1),
      PACKAGE = "localdepth") 
  } else {
    stop("the argument 'nsamp' must be either 'all' or a positive number")
  }
  result <- list()
  if (is.numeric(nsamp)) {
    #### Sistemare la divisione per il numero di simplessi
    result$localdepth <- z$localdepth/z$nd
    result$depth <- z$depth/z$nd
    result$max.localdepth <- max(result$localdepth)
    result$max.depth <- max(result$depth)  
    result$num <- c(z$nd,z$nld)
  } else {
    result$localdepth <- z$localdepth/nt
    result$depth <- z$depth/nt
    result$max.localdepth <- max(result$localdepth)
    result$max.depth <- max(result$depth)
    result$num <- c(nt,nt)
  }
  result$call <- match.call()
  result$tau <- tau
  result$use <- use
  result$tol <- tol
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- nsamp
  result$method <- 'simplicial.hyperspheres'
  class(result) <- 'localdepth'
  return(result)
}
