
"glsm.krige" <- function(mcmc.output, locations, borders, trend.l="cte", micro.scale=NULL, dist.epsilon= 1e-10,  output)
{
  call.fc <- match.call()
  coords <- mcmc.output$geodata$coords
  cov.model <- mcmc.output$model$cov.model
  kappa <- mcmc.output$model$kappa
  beta <- mcmc.output$model$beta
  cov.pars <- mcmc.output$model$cov.pars
  nugget <- mcmc.output$model$nugget
  aniso.pars <- mcmc.output$model$aniso.pars
  trend.d <- mcmc.output$model$trend
  lambda <- mcmc.output$model$lambda
  if(is.null(micro.scale)) micro.scale <- nugget
  else if(!is.numeric(micro.scale) || length(micro.scale)>1) stop("micro.scale must be a numeric number ")
  if(missing(borders)) borders <- mcmc.output$geodata$borders
  ##
  if(missing(output)) output <- output.glm.control()
  else output <- .output.glm.check.aux(output, fct = "glsm.krige")
  sim.predict <- output$sim.predict
  messages.screen <- output$messages.screen
  if(messages.screen) cat(" glsm.krige: Prediction for a generalised linear spatial model \n")
  ##
  ## Checking for 1D prediction
  if(missing(locations))
    stop("locations need to be specified for prediction; prediction not performed")
  else {
    locations <- .geoR.check.locations(locations)
    if(is.null(trend.l))
      stop("trend.l needed for prediction")
    if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
      krige1d <- TRUE
    else krige1d <- FALSE
  }
  ##
  beta.size <- length(beta)
  if(beta.size > 1) beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ##------------------------------------------------------------
######################## ---- prediction ----- #####################
  if(mcmc.output$model$family=="binomial"){
    krige <- list(type.krige = "sk", beta = beta, trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, 
                  cov.pars = cov.pars, kappa = kappa, nugget = nugget, micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
                  aniso.pars = aniso.pars, link = mcmc.output$model$link)
    kpl.result <- .glm.krige.aux(data = mcmc.output$simulations, coords = coords, locations = locations, borders=borders, krige = krige,
                                output = list(n.predictive = ifelse(sim.predict,1,0),
                                  signal = TRUE, messages=FALSE))
  }
  else{
    if(mcmc.output$model$family!="poisson") stop("only poisson and binomial are allowed as error distribution")
    krige <- list(type.krige = "sk", beta = beta, trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, 
                  cov.pars = cov.pars, kappa = kappa, nugget = nugget, micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
                  aniso.pars = aniso.pars, lambda = lambda)
    kpl.result <- .krige.conv.extnd(data = .BC.inv(mcmc.output$simulations, lambda), coords = coords, locations = locations, borders=borders, krige = krige,
                                   output = list(n.predictive = ifelse(sim.predict,1,0), signal = TRUE, messages = FALSE))
    
  }
  ##
  kpl.result$krige.var <- rowMeans(kpl.result$krige.var) + apply(kpl.result$predict, 1, var)
  kpl.result$mcmc.error <- sqrt(asympvar(kpl.result$predict,messages=FALSE)/ncol(kpl.result$predict))
  kpl.result$predict <- rowMeans(kpl.result$predict)
  kpl.result$beta <- NULL
  kpl.result$call <- call.fc
#######################################
  attr(kpl.result, "prediction.locations") <- call.fc$locations
  if(!is.null(locations)) attr(kpl.result, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  if(!is.null(call.fc$borders)) attr(kpl.result, "borders") <- call.fc$borders
  class(kpl.result) <- "kriging"
  return(kpl.result)
}
