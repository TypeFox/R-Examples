rmaxlin <- function(n, coord, cov.mod = "gauss", dsgn.mat, grid = FALSE,
                    p = 500, ...){
  ## This function generates realisation from a max-linear model

  if (missing(dsgn.mat)){
    if (cov.mod != "gauss")
      stop("Only cov.mod = 'gauss' is currently implemented")

    if (missing(coord))
      stop("You must supply either a design matrix or 'cov.mod' and 'coord'")
    
    if (is.null(dim(coord))){
      dim <- 1
      n.site <- length(coord)

      if (grid){
        grid <- FALSE

        warning("You cannot use 'grid = TRUE' for one dimensional simulations.
 Setting it to 'FALSE'.")
      }

      param <- "var"
    }
  
    else {
      dim <- ncol(coord)
      
      if (grid){
        dummy <- list()
        for (i in 1:dim)
          dummy <- c(dummy, list(coord[,i]))
      
        coord <- as.matrix(expand.grid(dummy))
      }
    
      n.site <- nrow(coord)
      param <- c("cov11", "cov12", "cov22")
    }

    ##Build the design matrix from the discretized Smith model
    if (!all(param %in% names(list(...))))
      stop("You must specify ", param)

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

    ## Compute the "design matrix" for the max-linear model (only
    ## based on the conditionning points)
    dsgn.mat <- .C("maxLinDsgnMat", as.double(coord), as.double(coord.grid),
                   as.integer(n.site), as.integer(p), as.double(areaPixel),
                   as.integer(dim), as.double(param), dsgnMat = double(p * n.site),
                   PACKAGE = "SpatialExtremes")$dsgnMat
  }

  else {
    if (!is.matrix(dsgn.mat) || any(dsgn.mat<0))
      stop("'dsgn.mat' must be a numeric matrix with non negative entries")

    n.site <- nrow(dsgn.mat)
    p <- ncol(dsgn.mat)
    grid <- FALSE
  }

  Z <- rgev(n * p, 1, 1, 1)

  ans <- .C("maxLinear", as.integer(n), as.double(dsgn.mat), as.double(Z),
            as.integer(n.site), as.integer(p), as.integer(grid),
            sim = double(n.site * n), PACKAGE = "SpatialExtremes")$sim

  if (grid){
    if (n == 1)
      ans <- matrix(ans, sqrt(n.site), sqrt(n.site))

    else
      ans <- array(ans, c(sqrt(n.site), sqrt(n.site), n))
  }

  else
    ans <- matrix(ans, nrow = n, ncol = n.site)
  
  return(ans)
}
