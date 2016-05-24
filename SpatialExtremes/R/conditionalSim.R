condrgp <- function(n, coord, data.coord, data, cov.mod = "powexp",
                    mean = 0, sill = 1, range = 1,
                    smooth = 1, grid = FALSE, control = list()){

  if (cov.mod == "caugen")
    stop("The generalized Cauchy covariance family isn't implemented")
  
  if (is.null(coord) & !is.null(data.coord))
    stop("'coord' and 'data.coord' don't match")

  if (!is.null(dim(coord)) & (any(ncol(coord) != ncol(data.coord))))
    stop("'coord' and 'data.coord' don't match")

  if (!is.null(dim(data)))
    if (nrow(data) != 1)
      stop("You can supply only one set of conditional observations")

  data <- as.numeric(data)

  new.coord <- coord

  if (grid){
    if (is.null(dim(coord)))
      stop("You cannot use 'grid = TRUE' for 1 dimensional processes")
    
    dummy <- NULL
    for (i in 1:nrow(new.coord))
      dummy <- rbind(dummy, cbind(new.coord[,1], new.coord[i,2]))

    new.coord <- dummy
  }
  
  if (is.null(dim(coord))){
    n.cond <- length(data.coord)
    n.loc <- length(new.coord)
    new.coord <- c(data.coord, new.coord)
  }

  else{
    n.cond <- nrow(data.coord)
    n.loc <- nrow(coord)
    new.coord <- rbind(data.coord, new.coord)
  }
    
  uncond <- rgp(n, new.coord, cov.mod = cov.mod, mean = mean, nugget = 0,
                sill = sill, range = range, smooth = smooth, control = control)

  weights <- kriging(data, data.coord, coord, cov.mod = cov.mod, sill = sill,
                     range = range, smooth = smooth, grid = grid,
                     only.weights = TRUE)$weights

  if (grid)
    ans <- array(NA, c(n.loc, n.loc, n))

  else
    ans <- matrix(NA, n, n.loc)
  
  res <- data - t(uncond[,1:n.cond])

  if (!grid)
    res <- t(res)
  
  ans[] <- uncond[,-(1:n.cond)] + res %*% weights

  if (grid & (n == 1))
    ans <- matrix(ans, n.loc)

  return(list(coord = coord, cond.sim = ans, data.coord = data.coord,
              data = data, cov.mod = cov.mod, grid = grid))  
}

condrmaxlin <- function(n, coord, data.coord, data, cov.mod = "gauss",
                        ..., grid = FALSE, p = 10000){

  data <- as.numeric(data)

  if (cov.mod != "gauss")
    stop("Currently only the 'discrete' Smith model is implemented")
  
  if (is.null(dim(data.coord))){
    dim <- 1
    n.cond.site <- length(data.coord)

    if (!is.null(dim(coord)))
      stop("'coord' and 'data.coord' don't match. Please check!")

    n.sim.site <- length(coord)

    if (grid){
      grid <- FALSE

      warning("You cannot use 'grid = TRUE' for one dimensional simulations.
 Setting it to 'FALSE'.")
    }
  }

  else {
    dim <- ncol(data.coord)
    n.cond.site <- nrow(data.coord)

    if (is.null(dim(coord)) || (ncol(coord) != dim))
      stop("'coord' and 'data.coord' don't match. Please check!")

    if (grid){
      dummy <- list()
      for (i in 1:dim)
        dummy <- c(dummy, list(coord[,i]))
      
      coord <- as.matrix(expand.grid(dummy))
    }
    
    n.sim.site <- nrow(coord)
  }

  if (length(data) != n.cond.site)
    stop("'data.coord' and 'data' don't match. Please check!")


  if (dim == 1)
    param <- "var"

  else if (dim == 2)
    param <- c("cov11", "cov12", "cov22")

  else
    stop("You can only use this function for one or two dimensional processes")

  if (!all(param %in% names(list(...))))
    stop(paste("You must specify ", param, sep="'"))

  param <- unlist(list(...)[param])

  if (dim == 1){
    bounds <- c(min(coord) - 4.1 * sqrt(param), max(coord) + 4.1 * sqrt(param))
    delta <- diff(bounds) / (p - 1)
    coord.grid <- seq(bounds[1], bounds[2], length = p) + delta / 2
  }

  else if (dim == 2){
    dummy <- sqrt(max(param["cov11"], param["cov22"]))
    bounds <- apply(coord, 2, range) + 4.1 * dummy * c(-1, 1)
    p <- ceiling(sqrt(p))
    delta <- (bounds[2,] - bounds[1,]) / (p - 1)
        
    coord.grid <- cbind(seq(bounds[1,1], bounds[2,1], length = p) + 0.5 * delta[1],
                        seq(bounds[1,2], bounds[2,2], length = p) + 0.5 * delta[2])
    coord.grid <- as.matrix(expand.grid(coord.grid[,1], coord.grid[,2]))
    p <- p^2
  }

  areaPixel <- prod(delta)

  ## Compute the "design matrix" for the max-linear model (only based
  ## on the conditionning points)
  dsgn.mat.cond <- .C("maxLinDsgnMat", as.double(data.coord), as.double(coord.grid),
                      as.integer(n.cond.site), as.integer(p), as.double(areaPixel),
                      as.integer(dim), as.double(param),
                      dsgnMat = double(p * n.cond.site), PACKAGE = "SpatialExtremes")$dsgnMat

  ## Get the realisation of the Zs (from the conditional distribution)
  Z <- .C("rcondMaxLin", as.double(data), as.double(dsgn.mat.cond),
          as.integer(p), as.integer(n.cond.site), as.integer(n), Z = double(p * n),
          PACKAGE = "SpatialExtremes")$Z

  ## Check if the conditional observation are honored
  sim.cond <- .C("maxLinear", as.integer(1), as.double(dsgn.mat.cond), as.double(Z),
                 as.integer(n.cond.site), as.integer(p), as.integer(FALSE),
                 sim = double(n.cond.site), PACKAGE = "SpatialExtremes")$sim

  honored <- all.equal(data, sim.cond)
  if (!isTRUE(honored))
    warning(paste("Some conditional observations aren't honored!\n The maximum absolute difference was ",
                  round(max(abs(data - sim.cond)), 2), " for site #", which.max(abs(data - sim.cond)), sep = ""))

  ## Compute the "design matrix" for the max-linear model (only for
  ## the points where we want our simulation)
  dsgn.mat.sim <- .C("maxLinDsgnMat", as.double(coord), as.double(coord.grid),
                     as.integer(n.sim.site), as.integer(p), as.double(areaPixel),
                     as.integer(dim), as.double(param),
                     dsgnMat = double(p * n.sim.site), PACKAGE = "SpatialExtremes")$dsgnMat

  ## Get the realisations at the desired locations
  ans <- .C("maxLinear", as.integer(n), as.double(dsgn.mat.sim), as.double(Z),
            as.integer(n.sim.site), as.integer(p), as.integer(grid),
            sim = double(n.sim.site * n), PACKAGE = "SpatialExtremes")$sim

  if (grid){
    if (n == 1)
      ans <- matrix(ans, sqrt(n.sim.site), sqrt(n.sim.site))

    else
      ans <- array(ans, c(sqrt(n.sim.site), sqrt(n.sim.site), n))
  }

  else
    ans <- matrix(ans, nrow = n, ncol = n.sim.site)
  
  return(ans)
}
