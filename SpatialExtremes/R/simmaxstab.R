rmaxstab <- function(n, coord, cov.mod = "gauss", grid = FALSE,
                     control = list(), ...){
  
  if (!(cov.mod %in% c("gauss","whitmat","cauchy","powexp","bessel",
                       "iwhitmat", "icauchy", "ipowexp", "ibessel",
                       "gwhitmat", "gcauchy", "gpowexp", "gbessel",
                       "twhitmat", "tcauchy", "tpowexp", "tbessel",
                       "brown")))
    stop("'cov.mod' must be one of 'gauss', '(i/g/t)whitmat', '(i/g/t)cauchy', '(i/g/t)powexp', '(i/g/t)bessel' or 'brown'")
  
  if (!is.null(control$method) && !(control$method %in% c("direct", "tbm", "circ")))
    stop("the argument 'method' for 'control' must be one of 'direct', 'tbm' and 'circ'")
  
  if (cov.mod == "gauss")
    model <- "Smith"
  
  else if (cov.mod == "brown")
    model <- "Brown-Resnick"
  
  else if (cov.mod %in% c("whitmat","cauchy","powexp","bessel"))
    model <- "Schlather"
  
  else {
    if (substr(cov.mod, 1, 1) == "i")
      model <- "iSchlather"

    else if (substr(cov.mod, 1, 1) == "g")
      model <- "Geometric"

    else
      model <- "Extremal-t"

    cov.mod <- substr(cov.mod, 2, 10)
  }

  dist.dim <- ncol(coord)
  
  if (is.null(dist.dim))
    dist.dim <- 1

  if (dist.dim > 2)
    stop("Currently this function is only available for R or R^2")

  if ((dist.dim == 1) && grid){
    warning("You cannot use 'grid = TRUE' in dimension 1. Ignored.")
    grid <- FALSE
  }

  if (model == "Smith"){
    if ((dist.dim == 1) && (!("var" %in% names(list(...)))))
      stop("For one dimensional simulations, you must specify 'var' not 'cov11', 'cov12', 'cov22'.")
    
    if ((dist.dim == 2) && (!all(c("cov11", "cov12", "cov22") %in% names(list(...)))))
      stop("You must specify 'cov11', 'cov12', 'cov22'")

    ##Get the model parameters
    var <- list(...)$var
    cov11 <- list(...)$cov11
    cov12 <- list(...)$cov12
    cov22 <- list(...)$cov22
    cov13 <- list(...)$cov13
    cov23 <- list(...)$cov23
    cov33 <- list(...)$cov33
  }

  else if (model == "Schlather"){
    if (!all(c("nugget", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'nugget', 'range', 'smooth'")
    
    nugget <- list(...)$nugget
    range <- list(...)$range
    smooth <- list(...)$smooth
  }

  else if (model == "iSchlather"){
    if (!all(c("alpha", "nugget", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'alpha', 'nugget', 'range', 'smooth'")
    
    nugget <- list(...)$nugget
    range <- list(...)$range
    smooth <- list(...)$smooth
    alpha <- list(...)$alpha
  }

  else if (model == "Geometric") {
    if (!all(c("sigma2", "nugget", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'sigma2', 'nugget', 'range', 'smooth'")
    
    nugget <- list(...)$nugget
    range <- list(...)$range
    smooth <- list(...)$smooth
    sigma2 <- list(...)$sigma2
  }

  else if (model == "Extremal-t"){
    if (!all(c("DoF", "nugget", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'DoF', 'nugget', 'range', 'smooth'")
    
    nugget <- list(...)$nugget
    range <- list(...)$range
    smooth <- list(...)$smooth
    DoF <- list(...)$DoF
  }
  
  else if (model == "Brown-Resnick"){
    range <- list(...)$range
    smooth <- list(...)$smooth
  }

  if (dist.dim !=1){
    n.site <- nrow(coord)
    coord.range <- apply(coord, 2, range)
    center <- colMeans(coord.range)
    edge <- max(apply(coord.range, 2, diff))
  }

  else {
    n.site <- length(coord)
    coord.range <- range(coord)
    center <- mean(coord.range)
    edge <- diff(coord.range)
  }

  
  cov.mod <- switch(cov.mod, "gauss" = "gauss", "whitmat" = 1, "cauchy" = 2,
                    "powexp" = 3, "bessel" = 4)

  if (grid){
    ans <- rep(-1e10, n * n.site^dist.dim)
    reg.grid <- .isregulargrid(coord[,1], coord[,2])
    steps <- reg.grid$steps
    reg.grid <- reg.grid$reg.grid
  }

  else
    ans <- rep(-1e10, n * n.site)

  ##Identify which simulation technique is the most adapted or use the
  ##one specified by the user --- this is useless for the Smith model.
  if (is.null(control$method)){
    if (grid && reg.grid)
      method <- "circ"
    
    else if ((length(ans) / n) > 600)
      method <- "tbm"
    
    else
      method <- "direct"
  }

  else
    method <- control$method

  if (method == "tbm"){
    if (is.null(control$nlines))
      nlines <- 1000

    else
      nlines <- control$nlines
  }      
  
  if (model == "Smith")
    ans <- switch(dist.dim,
                  .C("rsmith1d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), as.double(var), ans = ans,
                     PACKAGE = "SpatialExtremes")$ans,
                  .C("rsmith2d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), grid, as.double(cov11), as.double(cov12),
                     as.double(cov22), ans = ans, PACKAGE = "SpatialExtremes")$ans)
  
  else if (model == "Schlather"){

    if (is.null(control$uBound))
      uBound <- 3.5

    else
      uBound <- control$uBound

    if (method == "direct")
      ans <- .C("rschlatherdirect", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                as.double(range), as.double(smooth), as.double(uBound), ans = ans,
                PACKAGE = "SpatialExtremes")$ans
    
    else if (method == "circ")
      ans <- .C("rschlathercirc", as.integer(n), as.integer(n.site), as.double(steps),
                as.integer(dist.dim), as.integer(cov.mod), as.double(nugget), as.double(range),
                as.double(smooth),  as.double(uBound), ans = ans, PACKAGE = "SpatialExtremes")$ans
    
    else     
      ans <- .C("rschlathertbm", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                as.double(range), as.double(smooth), as.double(uBound), as.integer(nlines),
                ans = ans, PACKAGE = "SpatialExtremes")$ans
  }

  else if (model == "Geometric"){

    if (is.null(control$uBound))
      uBound <- exp(3.5 * sqrt(sigma2) - 0.5 * sigma2)

    else
      uBound <- control$uBound

    if (method == "direct")
      ans <- .C("rgeomdirect", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                as.double(nugget), as.double(range), as.double(smooth),
                as.double(uBound), ans = ans, PACKAGE = "SpatialExtremes")$ans

    else if (method == "circ")
      ans <- .C("rgeomcirc", as.integer(n), as.integer(n.site), as.double(steps),
                as.integer(dist.dim), as.integer(cov.mod), as.double(sigma2),
                as.double(nugget), as.double(range), as.double(smooth),  as.double(uBound),
                ans = ans, PACKAGE = "SpatialExtremes")$ans
    
    else      
      ans <- .C("rgeomtbm", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                as.double(nugget), as.double(range), as.double(smooth), as.double(uBound),
                as.integer(nlines), ans = ans, PACKAGE = "SpatialExtremes")$ans
  }

  else if (model == "Extremal-t"){
    if (is.null(control$uBound))
      uBound <- 3^DoF

    else
      uBound <- control$uBound

    if (method == "direct")
      ans <- .C("rextremaltdirect", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                as.double(range), as.double(smooth), as.double(DoF), as.double(uBound),
                ans = ans, PACKAGE = "SpatialExtremes")$ans

    else if (method == "circ")
      ans <- .C("rextremaltcirc", as.integer(n), as.integer(n.site), as.double(steps),
                as.integer(dist.dim), as.integer(cov.mod), as.double(nugget), as.double(range),
                as.double(smooth), as.double(DoF), as.double(uBound), ans = ans,
                PACKAGE = "SpatialExtremes")$ans
    
    else      
      ans <- .C("rextremalttbm", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(nugget),
                as.double(range), as.double(smooth), as.double(DoF), as.double(uBound),
                as.integer(nlines), ans = ans, PACKAGE = "SpatialExtremes")$ans
  }

  else if (model == "Brown-Resnick"){
    coord <- scale(coord, scale = FALSE)

    if (is.null(control$max.sim))
      max.sim <- 1000

    else
      max.sim <- control$max.sim
    
    if (is.null(control$uBound))
      uBound <- 10

    else
      uBound <- control$uBound

    if (is.null(control$sim.type))
      sim.type <- 1

    else
      sim.type <- control$sim.type
    
    if (dist.dim == 1)
      bounds <- range(coord)
    
    else
      bounds <- apply(coord, 2, range)

    if (is.null(control$nPP))
      nPP <- 15

    else
      nPP <- control$nPP

    if (sim.type == 6){
      idx.sub.orig <- getsubregions(coord, bounds, range, smooth, dist.dim)
      n.sub.orig <- length(idx.sub.orig)
    }
    
    else
      idx.sub.orig <- n.sub.orig <- 0
    
    if (method == "direct")
      ans <- .C("rbrowndirect", as.double(coord), as.double(bounds),
                as.integer(n), as.integer(n.site), as.integer(dist.dim),
                as.integer(grid), as.double(range), as.double(smooth),
                as.double(uBound), as.integer(sim.type), as.integer(max.sim),
                as.integer(nPP), as.integer(idx.sub.orig), as.integer(n.sub.orig),
                ans = ans, PACKAGE = "SpatialExtremes")$ans
  }

  if (grid){
    if (n == 1)
      ans <- matrix(ans, n.site, n.site)

    else
      ans <- array(ans, c(n.site, n.site, n))
  }

  else
    ans <- matrix(ans, nrow = n, ncol = n.site)

  return(ans)
}

getsubregions <- function(coord, bounds, range, smooth, dist.dim){

  if (dist.dim == 1){
    coord <- matrix(coord, ncol = 1)
    bounds <- matrix(bounds, ncol = 1)
  }
  
  h.star <- 2^(1/smooth) * range
  n.windows <- 1 + floor(diff(bounds) / h.star)

  sub.bounds <- list()
  for (i in 1:dist.dim)
    sub.bounds <- c(sub.bounds,
                    list(seq(bounds[1,i], bounds[2, i], length = n.windows[i] + 1)))

  sub.centers <- list()
  for (i in 1:dist.dim)
    sub.centers <- c(sub.centers, list(0.5 * (sub.bounds[[i]][-1] +
                                              sub.bounds[[i]][-(n.windows + 1)])))

  sub.origins <- as.matrix(expand.grid(sub.centers))

  n.sub.origins <- prod(n.windows)
  idx.sub.orig <- rep(NA, n.sub.origins)
  for (i in 1:n.sub.origins){
    dummy <- sqrt(colSums((t(coord) - sub.origins[i,])^2))
    idx.sub.orig[i] <- which.min(dummy)
  }
    

  return(idx.sub.orig = idx.sub.orig)
}
  
  
