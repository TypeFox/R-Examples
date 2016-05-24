

"covariog" <-  function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", uvec = "default", bins.lim = "default",
                        estimator.type = c("poisson", "not-poisson"), max.dist = NULL, pairs.min = 2)
{
  call.fc <- match.call()
  estimator.type <- match.arg(estimator.type)
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  n.data <- nrow(coords)
  if(n.data != nrow(data)) stop("Dimension of data and coords do not match ")
  n.datasets <- ncol(data)
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n.data)
  }
  if(any(units.m <= 0)) stop("units.m must be positive")
  data1 <- as.matrix(data/units.m)
  if(n.datasets == 1)
    data1 <- as.vector(data1)
  u <- as.vector(dist(as.matrix(coords)))
  if(!is.null(max.dist))
    umax <- max(u[u < max.dist])
  else umax <- max(u)
  if(all(bins.lim == "default")) {
    if(!is.null(max.dist))
      umax <- max(u[u < max.dist])
    else umax <- max(u)
    if(all(uvec == "default"))
      uvec <- seq(0, umax, l = 15)
    nvec <- length(uvec)
    d <- 0.5 * diff(uvec[seq(length=nvec)])
    bins.lim <- c(0, (uvec[seq(length=(nvec - 1))] + d), (d[nvec - 1] + uvec[nvec]))
    if(uvec[1] == 0)
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
  }
  else {
    if(!all(uvec == "default")) warning(" Both uvec and bins.lim are specified; using values in bins.lim ")
    if(!is.null(max.dist))
      if(max.dist != max(bins.lim)) stop("conflict between max.dist and max(bins.lim)")
  }
  nbins <- length(bins.lim) - 1
  if(is.null(max.dist))
    max.dist <- max(bins.lim)
  temp.list <- list(n.data = n.data, coords = coords, data = data1, nbins = nbins, bins.lim = bins.lim, max.dist = max.dist)
  bin.f <- function(data, temp.list)
    {
      result <- .C("binitprod",
                   as.integer(temp.list$n.data),
                   as.double(as.vector(temp.list$coords[, 1])),
                   as.double(as.vector(temp.list$coords[, 2])),
                   as.double(as.vector(data)),
                   as.integer(temp.list$nbins),
                   as.double(as.vector(temp.list$bins.lim)),
                   as.double(temp.list$max.dist),
                   cbin = as.integer(rep(0, temp.list$nbins)),
                   vbin = as.double(rep(0, temp.list$nbins)), PACKAGE = "geoRglm")[c("vbin", "cbin")]
      return(result)
    }
  result <- array(unlist(lapply(as.data.frame(data1), bin.f, temp.list = temp.list)), dim = c(nbins, 2, n.datasets))
  mm2 <- colMeans(as.matrix(data1))^2
  sigma2 <- log(t(t((apply(as.matrix(data1), 2, var) * (n.data - 1))/n.data - colMeans(as.matrix(data1/units.m)))/mm2) + 1)
  if(estimator.type == "poisson")
    calcresult <- array(rbind(as.vector(sigma2), log(t(t(as.matrix(result[, 1,  ]))/mm2))), dim = c(nbins + 1, n.datasets))
  else calcresult <- array(log(t(t(as.matrix(result[, 1,  ]))/mm2)), dim = c(nbins, n.datasets))
  indp <- (result[, 2, 1] >= pairs.min)
  if(estimator.type == "poisson")
    result <- list(u = c(0, uvec[indp]), v = calcresult[c(TRUE,indp),  ], n = c(n.data, result[indp, 2, 1]), bins.lim = bins.lim)
  else result <- list(u = uvec[indp], v = calcresult[indp,  ], n = result[indp, 2, 1], bins.lim = bins.lim)
  result <- c(result, list(v0 = sigma2, estimator.type = estimator.type, n.data = n.data, call = call.fc))
  class(result) <- "covariogram"
  return(result)
}


"covariog.model.env" <- 
function(geodata, coords = geodata$coords, units.m = "default", obj.covariog, model.pars, nsim = 500, prob = c(0.025, 0.975), 
     messages)
{
  call.fc <- match.call()
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  ##
  ## reading input
  ##
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, obj.covariog$n.data)
  }
  if(any(units.m <= 0)) stop("units.m must be positive")
  if(is.null(model.pars$beta) | !is.numeric(model.pars$beta)) {
    stop("argument beta must be provided in order to perform make covariogram envelopes")
  }
  else if(length(model.pars$beta) > 1) {
    stop("only one-dimensional beta is alllowed")
  }
  else beta <- model.pars$beta
  if(!is.null(model.pars$cov.model)) {
    cov.model <- model.pars$cov.model
  }
  else cov.model <- "exponential"
  if(!is.null(model.pars$kappa))
    kappa <- model.pars$kappa
  else kappa <- 0.5
  if(!is.null(model.pars$nugget))
    nugget <- model.pars$nugget
  else nugget <- 0
  cov.pars <- model.pars$cov.pars
  ##
  ## generating simulations from the model with parameters provided
  ##
  if(messages.screen) cat(paste("covariog.env: generating", nsim, "simulations(", obj.covariog$n.data,  "points). \n"))
  simula <- .pois.log.grf(obj.covariog$n.data, grid = as.matrix(coords), units.m = units.m, beta = beta, cov.model = cov.model, 
          cov.pars = cov.pars, nugget = nugget, kappa = kappa, nsim = nsim, messages = FALSE)
  ##
  ## computing empirical covariograms for the simulations
  ##
  if(messages.screen) cat(paste("covariog.env: computing the empirical covariogram for the", nsim, "simulations\n"))
  simula.result <- covariog(simula, bins.lim = obj.covariog$bins.lim, units.m = units.m, estimator.type = obj.covariog$estimator.type,
            pairs.min = min(obj.covariog$n))
  ##
  ## computing envelopes
  ##
  if(messages.screen) cat("covariog.env: computing the envelopes\n")
  limits <- apply(simula.result$v, 1, quantile, prob = prob)
  if(length(prob) == 1) {
    res.env <- list(u = obj.covariog$u, v = limits)
  }
  else {
    res.env <- list(u = obj.covariog$u, v = t(limits))
  }
  res.env$call <- call.fc
  return(res.env)
}


"lines.covariomodel" <- 
function(x, max.dist = x$max.dist, ...)
{
  if (is.null(max.dist)) stop("argument max.dist needed for this object")
  my.l <- x
  if(is.null(my.l$nugget)) my.l$nugget <- 0
  C.f <- function(x, my.l){ 
    return(ifelse(x>0, 0, my.l$nugget)+cov.spatial(x, cov.model = my.l$cov.model, kappa = my.l$kappa, cov.pars = my.l$cov.pars))
  }
  curve(C.f(x,my.l=my.l), from = 0, to = max.dist, add=TRUE, ...)
  return(invisible())
}


"plot.covariogram" <- 
function(x, max.dist = max(x$u), ylim = "default", type = "b", envelope.obj = NULL, ...)
{
        u <- x$u[x$u <= max.dist]
	if(!is.matrix(x$v)) {
		v <- x$v[x$u <= max.dist]
		ymax <- max(x$v[x$u <= max.dist])
		ymin <- min(x$v[x$u <= max.dist & x$v >  - Inf])
	}
	else {
		v <- x$v[x$u <= max.dist, 1]
		ymax <- max(x$v[x$u <= max.dist,  ])
		temp <- x$v[x$u <= max.dist,  ]
		ymin <- min(temp[temp >  - Inf])
	}
	if(!is.null(envelope.obj)) {
		ymax <- max(envelope.obj$v, ymax)
		ymin <- min(envelope.obj$v[envelope.obj$v >  - Inf], ymin)
	}
	if(is.numeric(ylim))
		plot(u, v, xlim = c(0, max.dist), ylim = ylim, xlab = "Distance", ylab = "Covariance", type = type, ...)
	else plot(u, v, xlim = c(0, max.dist), ylim = c(ymin, ymax), xlab = "Distance", ylab = "Covariance", type = type, ...)
	abline(h = 0, lty = 3)
	if(is.matrix(x$v)) {
		for(k in seq(2,ncol(x$v))) {
			v <- x$v[x$u <= max.dist, k]
			lines(u, v, type = "b")
		}
	}
	if(!is.null(envelope.obj)) {
		if(ncol(as.matrix(envelope.obj$v)) == 1)
			lines(envelope.obj$u, envelope.obj$v, lty = 4)
		else apply(envelope.obj$v, 2, lines, x = envelope.obj$u, lty = 4)
	}
	return(invisible())
}


".pois.log.grf" <- 
function(n = NULL, grid, nx = round(sqrt(n)), ny = round(sqrt(n)), nsim = 1,
        xlims = c(0, 1), ylims = c(0, 1), units.m = "default", trend = "cte", beta = stop("beta parameter(s) needed"), 
        cov.model = c("exponential", "matern", "gaussian", "spherical", "cubic", "wave", "powered.exponential", "cauchy",
        "gneiting", "gneiting.matern", "pure.nugget"), cov.pars = stop("covariance parameters (sigmasq and phi) needed"), nugget = 0,
        kappa = 0.5, method = c("cholesky", "svd", "eigen", "circular.embedding"), messages, ...)
{
#
# function for simulation of a Poisson log-Gaussian random field.
#
# the syntax is similar to grf, but beta needs to be specified.
## beta   : vector of regression parameters
## trend  : by default a constant mean
## units.m   : n-vector of observation-times for data (default is 1 for all observations) 
#        
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  ##
  call.fc <- match.call()
  method <- match.arg(method)
  cov.model <- match.arg(cov.model)
  work <- grf(n = n, grid = grid, nx = nx, ny = ny, nsim = nsim, xlims = xlims, ylims = ylims, cov.model = cov.model, cov.pars = 
              cov.pars, nugget = nugget, lambda = 1, kappa = kappa, method = method, messages =  messages.screen, ...)
  trend.data <- unclass(trend.spatial(trend = trend, geodata = work))
  n <- nrow(work$coords)
  if((length(beta) == 1) & trend == "cte") mu <- rep(beta, n) else mu <- trend.data %*% as.vector(beta)
  if(all(units.m == "default")) units.m <- rep(1, n)
  if(any(units.m <= 0)) stop("units.m must be positive")
  lambda <- matrix(rep(units.m * exp(mu), nsim), nrow = n, ncol = nsim) * exp(as.matrix(work$data))
  zsim <- matrix(rpois(n * nsim, lambda = lambda), nrow = n, ncol = nsim)
  if(nsim == 1)
    zsim <- as.vector(zsim)
  results <- list(coords = work$coords, data = zsim, cov.model = cov.model, nugget = nugget, cov.pars = cov.pars, kappa = kappa,
                  mu = mu, units.m = units.m, messages = work$messages, method = method)
  results$.Random.seed <- work$.Random.seed
  results$call <- call.fc
  class(results) <- "geodata"
  return(results)
}



