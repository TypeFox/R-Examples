distance <- function(coord, vec = FALSE){
  ##This function computes the distance between each pair of locations

  ##coord is a matrix giving the coordinates (1 row = 1 station)
  if (is.null(dim(coord))){
    n.site <- length(coord)
    dist.dim <- 1
  }

  else{
    n.site <- nrow(coord)
    dist.dim <- ncol(coord)
  }

  n.pairs <- n.site * (n.site - 1) / 2

  if (vec){
    dist <- .C("distance", as.double(coord), as.integer(dist.dim),
               as.integer(n.site), vec, dist = double(dist.dim * n.pairs),
               PACKAGE = "SpatialExtremes")$dist
    dist <- matrix(dist, ncol = dist.dim, nrow = n.pairs)
  }

  else
    dist <- .C("distance", as.double(coord), as.integer(dist.dim),
               as.integer(n.site), vec, dist = double(n.pairs),
               PACKAGE = "SpatialExtremes")$dist

  return(dist)
}

gev2frech <- function(x, loc, scale, shape, emp = FALSE){

  if (emp){
    probs <- rank(x, na.last = "keep") / (length(x) + 1)
    x <- - 1 / log(probs)
    return(x)
  }

  if (shape == 0)
    exp((x - loc)/scale)

  else
    pmax(1 + shape * (x - loc) / scale, 0)^(1/shape)
}

frech2gev <- function(x, loc, scale, shape){
  if (shape == 0)
    scale * log(pmax(x, 0)) + loc

  else
    loc + scale * (pmax(x, 0)^shape - 1) / shape
}

.qgev <- function(p, loc = 1, scale = 1, shape = 1,
                  lower.tail = TRUE){

    if ((min(p, na.rm = TRUE) <= 0) || (max(p, na.rm = TRUE) >=1))
      stop("'p' must contain probabilities in (0,1)")

    if (min(scale) < 0)
      warning("There are some invalid scale GEV parameters")

    if (length(p) != 1)
      stop("invalid p")

    if (!lower.tail)
      p <- 1 - p

    n <- length(loc)

    ans <- .C("gev", as.double(p), as.integer(n), as.double(loc),
              as.double(scale), as.double(shape), quant = double(n),
              PACKAGE = "SpatialExtremes")$quant

    return(ans)
}

.frech2gev <- function(x, loc, scale, shape)
  ##This is a specific function that do the same as frech2gev but
  ##without any checking. x MUST be a vector. Use with caution.
  loc + scale * (pmax(x, 0)^shape - 1) / shape

.isregulargrid <- function(x, y, tol.x = 1e-4, tol.y = 1e-4){
  ##This function check if the grid defined by x and y is regular
  ##i.e. the spacings along the axis is constant -- but not
  ##necessarily the same for the x and y-axis

  x.diff <- diff(x)
  y.diff <- diff(y)

  eps.x <- diff(range(x.diff))
  eps.y <- diff(range(y.diff))

  if ((eps.x <= tol.x) && (eps.y <= tol.y)){
    reg.grid <- TRUE
    steps <- c(mean(x.diff), mean(y.diff))
  }

  else{
    reg.grid <- FALSE
    steps <- NA
  }

  return(list(reg.grid = reg.grid, steps = steps))
}

.useloglink <- function(formula)
  return(substr(formula[2], 1, 3) == "log")

.getTrendSurfCoeffNames <- function(formula){
  formula.terms <- terms(formula)
  names <- attributes(formula.terms)$term.labels

  if (attributes(formula.terms)$intercept)
    names <- c("Intercept", names)

  return(names)
}

logit <- function(p, rep_one_by = 0.999, rep_zero_by = 10^-3, inv = FALSE){

    if (inv)
        return(1 / (1 + exp(-p)))
    
    p[p > rep_one_by] <- rep_one_by
    p[p<rep_zero_by] <- rep_zero_by

    return(log(p / (1-p)))
}
