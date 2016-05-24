################################################################################
## window.ramps/window.predict.ramps - Time windows for ramps MCMC output
################################################################################

.onLoad <- function(lib, pkg) {
   setAs("matrix", "data.frame", function(from) as.data.frame(from))
}

window.ramps <- function(x, iter, ...)
{
   idx <- match(iter, as.numeric(rownames(x$params)))
   idx <- sort(unique(idx[!is.na(idx)]))

   if (length(idx) == 0) stop("No matching MCMC iterations")

   x$params <- as.mcmc(x$params[idx, , drop = FALSE])
   x$z <- as.mcmc(x$z[idx, , drop = FALSE])
   x$loglik <- x$loglik[idx]
   x$evals <- x$evals[idx]
   x$control$iter <- idx

   x
}

window.predict.ramps <- function(x, iter, ...)
{
   idx <- match(iter, as.numeric(rownames(x)))
   idx <- sort(unique(idx[!is.na(idx)]))

   if (length(idx) == 0) stop("No matching MCMC iterations")

   structure(x[idx, , drop = FALSE],
        coords = attr(x, "coords"),
        class = c("predict.ramps", "matrix"))
}


################################################################################
## expand.chain - Resumed sampling for a ramps object
################################################################################

expand.chain <- function(object, n)
{
   if (class(object) != "ramps") stop("Object must be of class 'ramps'")

   nr <- nrow(object$params)
   control <- object$control

   ## Set initial values to last sample
   control$beta$init <- params2beta(object$params[nr,], control)
   control$phi$init <- params2phi(object$params[nr,], control)
   val <- prop.table(params2kappa(object$params[nr,], control))
   control$sigma2.e$init <- kappa2kappa.e(val, control)
   control$sigma2.z$init <- kappa2kappa.z(val, control)
   control$sigma2.re$init <- kappa2kappa.re(val, control)

   ## Additional iterations to sample
   control$expand <- max(control$iter)
   inc <- diff(c(0, control$iter))[nr]
   control$iter <- seq(inc, n, by = inc)

   val <- ramps.engine(object$y, object$xmat, object$kmat, object$wmat,
                       object$correlation, object$etype, object$ztype,
                       object$retype, object$weights, control)

   object$control$iter <- c(object$control$iter, control$expand + control$iter)
   object$params <- as.mcmc(rbind(object$params, val$params))
   rownames(object$params) <- object$control$iter
   object$z <- as.mcmc(rbind(object$z, val$z))
   rownames(object$z) <- object$iter
   object$loglik <- c(object$loglik, val$loglik)
   object$evals <- c(object$evals, val$evals)

   object
}


################################################################################
## genUSStateGrid/genUSStateSites - Grid point generation for Monte Carlo
##    integration
################################################################################

genUSStateGrid <- function(state, incr = NULL, resolution = NULL) {
  mymap <- map('state', state, plot=FALSE)
  lon.range <- mymap$range[1:2]
  lat.range <- mymap$range[3:4]

  ## determine incr
  if (is.null(incr)) {
    if (is.null(resolution)) stop("One of 'incr' and 'resolution' must be specified")
    if (length(resolution) != 2) stop("Resolution needs to be vector of length 2.")
    lon.incr <- diff(lon.range) / (resolution[1] - 1)
    lat.incr <- diff(lat.range) / (resolution[2] - 1)
    incr <- c(lon.incr, lat.incr)
  }

  incr <- rep(incr, length.out = 2)
  ## generate points grid covering Iowa
  lon <- seq(lon.range[1], lon.range[2], by = incr[1]) 
  lat <- seq(lat.range[1], lat.range[2], by = incr[2])
  coords <- expand.grid(lon, lat)

  ## get county ids
  counties <- map.where('county', coords[,1], coords[,2])
  ## restrict to points in state
  inState <- substr(counties, 1, nchar(state)) == state
  inState <- inState & (!is.na(inState))
  conames <- map("county", state, namesonly=TRUE, plot=FALSE) 
  blockid <- pmatch(counties, conames, duplicates.ok=TRUE)

  grid <- data.frame(lon=coords[,1], lat=coords[,2], id=blockid,
                     county=conames[blockid])
  grid <- grid[inState,]
  grid
}

genUSStateSites <- function(state, nsites) {
  mymap <- map('state', state, plot=FALSE)
  lon.range <- mymap$range[1:2]
  lat.range <- mymap$range[3:4]

  nneed <- nsites
  coords <- NULL
  while (nneed > 0) {
    lat <- runif(nsites * 2, min = lat.range[1], max = lat.range[2])
    lon <- runif(nsites * 2, min = lon.range[1], max = lon.range[2])
    statename <- map.where('state', lon, lat)
    inState <- statename == state
    newsites <- cbind(lon, lat)[inState,]
    nthis <- sum(inState)
    if (nthis > nneed) {
      coords <- cbind(coords, newsites[1:nneed,])
      break
    } else {
      coords <- cbind(coords, newsites)
      nneed <- nsites - nthis
    }
  }
  coords
}


################################################################################
## Parameter extraction
################################################################################

params2phi <- function(params, control) {
   idx <- seq(length.out = length(control$phi))
   if (is.matrix(params)) params[,idx,drop=FALSE]
   else params[idx]
}

params2kappa <- function(params, control) {
   idx <- seq(length(control$phi) + 1, length.out = length(control$sigma2.e)
              + length(control$sigma2.z) + length(control$sigma2.re))
   if (is.matrix(params)) params[,idx,drop=FALSE]
   else params[idx]
}

params2beta <- function(params, control) {
   p <- length(control$beta)
   idx <- seq(length(params) - p + 1, length.out = p)
   if (is.matrix(params)) params[,idx,drop=FALSE]
   else params[idx]
}

kappa2kappa.e <- function(kappa, control) {
   idx <- seq(length.out = length(control$sigma2.e))
   if (is.matrix(kappa)) kappa[,idx,drop=FALSE]
   else kappa[idx]
}

kappa2kappa.z <- function(kappa, control) {
   idx <- seq(length(control$sigma2.e) + 1, length.out = length(control$sigma2.z))
   if (is.matrix(kappa)) kappa[,idx,drop=FALSE]
   else kappa[idx]
}

kappa2kappa.re <- function(kappa, control) {
   idx <- seq(length(control$sigma2.e) + length(control$sigma2.z) + 1,
              length.out = length(control$sigma2.re))
   if (is.matrix(kappa)) kappa[,idx,drop=FALSE]
   else kappa[idx]
}

sigma2init <- function(control)
{
   c(control$sigma2.e$init, control$sigma2.z$init, control$sigma2.re$init)
}

sigma2tuning <- function(control)
{
   c(control$sigma2.e$tuning, control$sigma2.z$tuning, control$sigma2.re$tuning)
}

sigma2shape <- function(control)
{
   c(control$sigma2.e$shape, control$sigma2.z$shape, control$sigma2.re$shape)
}

sigma2scale <- function(control)
{
   c(control$sigma2.e$scale, control$sigma2.z$scale, control$sigma2.re$scale)
}


################################################################################
## Statistical functions
################################################################################

rdirichlet <- function(n, alpha) {
   ## taken from MCMCpack
   ## Markov Chain Monte Carlo Package (MCMCpack)
   ## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn

   l <- length(alpha)
   x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
   sm <- x %*% rep(1, l)

   x / as.vector(sm)
}

runif.simplex <- function(n) {
   ## draws uniformly from standard simplex with nvert vertices
   rdirichlet(1, rep(1, n))
}

# Multivariate Normal Random Number Generator
rmvnorm2 <- function(n, mu, sigma)
{
   mu <- as.vector(mu)
   m <- ncol(sigma)

   if (nrow(sigma) != m) {
      stop("sigma must be a square matrix")
   }
   if (length(mu) != nrow(sigma)) {
      stop("mu and sigma have non-conforming size")
   }

   matrix(rnorm(n * m), nrow=n, ncol=m) %*% chol(sigma) +
          matrix(mu, nrow=n, ncol=m, byrow=T)
}


################################################################################
## Data management function
################################################################################

inbounds <- function(x, bounds)
{
   if (length(x) != nrow(bounds)) {
      warning("Number of supplied parameter values must be ", nrow(bounds),
              " instead of ", length(x))
      return(rep(FALSE, length(x)))
   }

   type <- bounds[,3]
   lower <- ifelse(type %in% c(2,4), x >= bounds[,1], x > bounds[,1])
   upper <- ifelse(type %in% c(3,4), x <= bounds[,2], x < bounds[,2])

   val <- lower & upper
   if (!all(val)) {
      lchar <- c("(", "[", "(", "[")
      uchar <- c(")", ")", "]", "]")
      support <- apply(bounds[,c(1,2),drop=FALSE], 1, paste, collapse=", ")
      warning("Out of bounds parameter value ",
              paste(x, ifelse(val, " = ", " != "),
                    lchar[type], support, uchar[type], sep="", collapse=", "))
   }

   val
}


unique.sites <- function(x)
{
   type <- class(x)

   x <- as.matrix(x)
   coords <- unique(x)

   n <- ncol(x)
   y <- merge(cbind(x, 1:nrow(x)), cbind(coords, 1:nrow(coords)), by = 1:n)
   idx <- y[order(y[, n+1]), n+2]
   map <- Matrix(0, length(idx), nrow(coords))
   map[seq(idx) + nrow(map) * (idx - 1)] <- 1

   list(coords = as(coords, type), map = map)
}


################################################################################
## I/O functions
################################################################################

name.ext <- function(name, ext)
{
   if (is.null(ext)) name <- NULL
   else if (length(ext) > 1) name <- paste(name, ext, sep="")

   name
}

write.header <- function(header, file, n = 1)
{
    if(n > 1)  header <- paste(header, 1:n, sep="")

    if(length(header) && !is.null(file))
        cat(paste(c("iter", header), collapse="\t"), "\n", sep="", file=file)
}

write.params <- function(iter, params, file)
{
    if(length(params) && !is.null(file))
        cat(paste(c(iter, params), collapse="\t"), "\n", sep="", file=file, append=T)
}

print.iter <- function(iter)
{
   if (iter == 0) cat(iter)
   else if (iter %% 250 == 0) cat(iter, "\n")
   else if (iter %% 50 == 0) cat(iter)
   else if (iter %% 5 == 0) cat(".")
}
