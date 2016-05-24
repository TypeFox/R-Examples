rgp <- function(n, coord, cov.mod = "powexp", mean = 0, nugget = 0, 
                sill = 1, range = 1, smooth = 1, grid = FALSE,
                control = list()){

  dist.dim <- ncol(coord)

  if (is.null(dist.dim)){
    dist.dim <- 1
    n.site <- length(coord)
  }

  else
    n.site <- nrow(coord)
  
  if (grid && is.null(dim(coord)))
    stop("'grid' cannot be 'TRUE' if you specify univariate coordinates")

  if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel", "fbm")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel', 'fbm'")

  if (!is.null(control$method) && !(control$method %in% c("exact", "tbm", "circ")))
    stop("the argument 'method' for 'control' must be one of 'exact', 'tbm' or 'circ'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  if (cov.mod == "bessel")
    cov.mod.num <- 4
  if (cov.mod == "fbm")
    cov.mod.num <- 6

  ##Check if a regular grid is specified
  if (grid){
    reg.grid <- .isregulargrid(coord[,1], coord[,2])
    steps <- reg.grid$steps
    reg.grid <- reg.grid$reg.grid
    ngrid <- nrow(coord)
  }

  else
    reg.grid <- FALSE
    
  if (is.null(control$method)){
    ##Identify the most accurate method for simulation if not specified

    if (reg.grid)
      method <- "circ"
    
    else if (grid && (n.site^dist.dim > 500))
      method <- "tbm"
    
    else if ((n.site > 500) && (dist.dim > 1))
      method <- "tbm"
    
    else
      method <- "exact"
  }

  else
    method <- control$method

  if (is.null(control$nlines))
    nlines <- 1000

  else
    nlines <- control$nlines

  gp <- switch(method,
               "tbm" = .tbmgp(n, coord, cov.mod.num, nugget, sill, range,
                 smooth, grid, nlines = nlines),
               "exact" = .exactgp(n, coord, cov.mod.num, nugget, sill, range,
                 smooth, grid),
               "circ" = .circgp(n, ngrid, steps, dist.dim, cov.mod.num, nugget,
                 sill, range, smooth))

  return(mean + gp)
}

.exactgp <- function(n, coord, cov.mod, nugget, sill, range, smooth, grid){
  dist.dim <- ncol(coord)
  n.site <- nrow(coord)

  if (is.null(dist.dim)){
    dist.dim <- 1
    n.site <- length(coord)
  }

  if (grid)
    n.effsite <- n.site^dist.dim

  else
    n.effsite <- n.site

  gp <- .C("direct", as.integer(n), as.integer(n.site), grid, as.integer(cov.mod), 
           as.double(coord), as.integer(dist.dim), as.double(nugget), 
           as.double(sill), as.double(range), as.double(smooth), ans = double(n.effsite * n), 
           PACKAGE = "SpatialExtremes")$ans
  
  if (grid){
    if ((n == 1) && (dist.dim == 2))
      gp <- matrix(gp, n.site, n.site)
    
    else
      gp <- array(gp, c(rep(n.site, dist.dim), n))
  }

  else
    gp <- matrix(gp, ncol = n.site, nrow = n)
  
  return(gp)
}

.tbmgp <- function(n, coord, cov.mod, nugget, sill, range, smooth, grid,
                   nlines = 1000){

  n.site <- nrow(coord)
  dim <- ncol(coord)

  if (grid)
    ans <- double(n.site^dim * n)

  else
    ans <- double(n.site * n)
  
  gp <- .C("tbm", as.integer(n), as.integer(n.site), as.integer(dim),
           as.integer(cov.mod), grid, as.double(coord), as.double(nugget),
           as.double(sill), as.double(range), as.double(smooth),
           as.integer(nlines), ans = ans, PACKAGE = "SpatialExtremes")$ans

  if (grid){
    if ((n == 1) && (dim == 2))
      gp <- matrix(gp, n.site, n.site)
    
    else
      gp <- array(gp, c(rep(n.site, dim), n))
  }

  else
    gp <- matrix(gp, nrow = n, ncol = n.site)

  return(gp)
}

.circgp <- function(n, n.grid, steps, dim, cov.mod, nugget, sill, range, smooth){
  
  gp <- .C("circemb", as.integer(n), as.integer(n.grid), as.double(steps),
           as.integer(dim), as.integer(cov.mod), as.double(nugget),
           as.double(sill), as.double(range), as.double(smooth),
           ans = double(n * n.grid^2), PACKAGE = "SpatialExtremes")$ans

  if (n == 1)
    return(matrix(gp, n.grid, n.grid))

  else
    return(array(gp, c(n.grid, n.grid, n)))
}
