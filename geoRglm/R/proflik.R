
"proflik.glsm" <-
  function (mcmc.obj, obj.likfit.glsm, 
            phi.values, nugget.rel.values, messages, ...)
{
  ##
  ## Checking input
  ##
  call.fc <- match.call()
  temp.list <- list()
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  ##
  cov.model <- obj.likfit.glsm$cov.model
  kappa <- obj.likfit.glsm$kappa
  aniso.pars <- obj.likfit.glsm$aniso.pars
  lambda <- obj.likfit.glsm$lambda
  ##
  if(is.null(mcmc.obj$S)){
    if(is.null(mcmc.obj$mu)) stop("mcmc.obj should include either an object mu or an object S.")
    n <- temp.list$n <- nrow(mcmc.obj$mu)
    temp.list$mu <- mcmc.obj$mu
    est.boxcox <- TRUE
  }
  else{
    n <- temp.list$n <- nrow(mcmc.obj$S)
    temp.list$z <- mcmc.obj$S 
    est.boxcox <- FALSE
  }
  if(is.null(aniso.pars)) coords <- as.matrix(mcmc.obj$geodata$coords)
  else coords <- coords.aniso(coords = as.matrix(mcmc.obj$geodata$coords), aniso.pars = aniso.pars)
  temp.list$xmat <- obj.likfit.glsm$trend.matrix
  temp.list$xmat <- unclass(trend.spatial(trend=obj.likfit.glsm$trend, geodata = mcmc.obj$geodata))
  temp.list$beta.size <- dim(temp.list$xmat)[2]
  temp.list$coords <- coords
  temp.list$cov.model <- cov.model
  temp.list$kappa <- kappa
  proflik.val <- matrix(NA,length(phi.values),length(nugget.rel.values))
  if(messages.screen) cat("proflik.glsm: computing 2-D profile likelihood for the phi and relative nugget parameters \n")
  fixed.values <- list()
  ## where to find the parameter values :
  if(est.boxcox){
    fixed.values$lambda <- lambda
    ip <- list(f.tausq.rel = FALSE, f.lambda = TRUE)
  }
  else ip <- list(f.tausq.rel = FALSE)
  ##
  temp.list$log.f.sim <- mcmc.obj$log.f.sim
  temp.list$messages.screen <- getOption("verbose") ## we want a default FALSE here
  ##
  for(i in seq(along=phi.values)){
    for(j in seq(along=nugget.rel.values)){
      ini <- c(phi.values[i],nugget.rel.values[j])      
      if(est.boxcox) lik.val[i,j] <- .lik.sim.boxcox(pars=ini, fp = fixed.values, ip = ip, temp.list = temp.list)
      else proflik.val[i,j] <- .lik.sim(pars=ini, fp = fixed.values, ip = ip, temp.list = temp.list)
    }
  }
  result <- list()
  result$rangenugget.rel <- list(range=phi.values, nugget.rel = nugget.rel.values, proflik.rangenugget.rel = proflik.val, est.rangenugget.rel = obj.likfit.glsm$loglik)
  result$n.bi <- 1
  result$n.uni <-  0
  result$method.lik <- "ML"
  result$call <- call.fc
  class(result) <- "proflik" 
  return(result)
}

