
".pmixed" <- 
  function(value, parms, df)
{
  if(ncol(parms$var)==1) parms$var <- as.vector(parms$var)
  if(df == Inf){
    sd.error <- sqrt(parms$var)
    if(any(sd.error < 1e-12)) sd.error[sd.error < 1e-12] <- 1e-12
    temp <- array(pnorm((value - parms$mean)/sd.error), dim = c(nrow(parms$mean), ncol(parms$mean)))
  }
  else{
    sd.error <- sqrt(parms$var*(df - 2)/df)
    if(any(sd.error < 1e-12)) sd.error[sd.error < 1e-12] <- 1e-12
    temp <- array(pt((value - parms$mean)/sd.error, df = df), dim = c(nrow(parms$mean), ncol(parms$mean)))
  }
  temp2 <- rowMeans(temp)
  return(temp2)
}

"model.glm.control" <- 
  function(trend.d = "cte", trend.l = "cte", cov.model = "matern", kappa = 0.5, aniso.pars = NULL, lambda = 0)
{
  cov.model <- match.arg(cov.model,
              choices = c("matern", "exponential", "gaussian",
                "spherical", "circular", "cubic", "wave", "power",
                "powered.exponential", "cauchy", "gneiting",
                "gneiting.matern", "pure.nugget"))
  if(cov.model == "powered.exponential" & (kappa <= 0 | kappa > 2))
    stop("model.glm.control: for power exponential correlation model the parameter kappa must be in the interval (0,2]")
  if(cov.model == "power") stop("model.glm.control: correlation function does not exist for the power variogram")
  if(!is.null(aniso.pars)){ 
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be a vector with two elements: rotation angle (in radians) and anisotropy ratio (a number > 1)")
  }
  res <- list(trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, kappa = kappa, aniso.pars = aniso.pars, lambda = lambda)
  class(res) <- "model.geoRglm"
  return(res)
}


".model.glm.check.aux" <-
  function(model, fct)
{
  if(class(model) != "model.geoRglm"){
    if(!is.list(model))
      stop(paste(fct,": the argument model only takes a list or an output of the function model.glm.control"))
    else{
      model.names <- c("trend.d", "trend.l", "cov.model", "kappa", "aniso.pars", "lambda")
      model <- .object.match.names(model,model.names)
      if(is.null(model$trend.d)) model$trend.d <- "cte"  
      if(is.null(model$trend.l)) model$trend.l <- "cte"  
      if(is.null(model$cov.model)) model$cov.model <- "matern"  
      if(is.null(model$kappa)) model$kappa <- 0.5
      if(is.null(model$lambda)){
        if(fct=="pois.krige.bayes") model$lambda <- 0
        if(fct=="binom.krige.bayes") model$lambda <- NULL
      }
      model <- model.glm.control(trend.d = model$trend.d,
                                 trend.l = model$trend.l,
                                 cov.model = model$cov.model,
                                 kappa = model$kappa,
                                 aniso.pars = model$aniso.pars,
                                 lambda = model$lambda)
    }
  }
  return(model)
}


".multgauss" <- 
  function(cov)
{
  return(crossprod(chol(cov), rnorm(n=ncol(cov))))
}

"mcmc.control" <- 
  function(S.scale, Htrunc="default", S.start, burn.in=0, thin=10, n.iter=1000*thin, phi.start="default",  phi.scale=NULL)
{
  if(missing(S.scale)) stop("S.scale parameter must to be provided for MCMC-proposal")
  if(missing(S.start) || is.null(S.start)) S.start<- "default"
  if(is.null(Htrunc)) Htrunc <- "default"
  if(is.null(burn.in)) burn.in <- 0
  if(is.null(thin)) thin <- 10
  if(is.null(n.iter)) n.iter <- 1000*thin
  res <- list(S.scale = S.scale, Htrunc = Htrunc, S.start = S.start, burn.in = burn.in, thin = thin, n.iter = n.iter,
              phi.start = phi.start, phi.scale=phi.scale)
  class(res) <- "mcmc.geoRglm"
  return(res)
}

".mcmc.check.aux" <-
  function(mcmc.input, fct)
{
  if(class(mcmc.input) != "mcmc.geoRglm"){
    if(!is.list(mcmc.input))
      stop(paste(fct,": the argument mcmc.input only takes a list or an output of the function mcmc.control"))
    else{
      mcmc.input.names <- c("S.scale", "Htrunc", "S.start", "burn.in", "thin", "n.iter", "phi.start", "phi.scale")  
      mcmc.input <- .object.match.names(mcmc.input,mcmc.input.names)
      if(fct=="pois.krige.bayes" | fct=="binom.krige.bayes"){
        if(is.null(mcmc.input$phi.start)) mcmc.input$phi.start <- "default"
      }
      mcmc.input <- mcmc.control(S.scale = mcmc.input$S.scale,Htrunc=mcmc.input$Htrunc,S.start=mcmc.input$S.start,
                                 burn.in=mcmc.input$burn.in,thin=mcmc.input$thin,n.iter=mcmc.input$n.iter,
                                 phi.start=mcmc.input$phi.start,phi.scale=mcmc.input$phi.scale)
    }
  }
  return(mcmc.input)
}


".BC.inv" <- 
  function(z,lambda)
{
  return(BCtransform(z, lambda = lambda, inverse=TRUE)$data)
}


"output.glm.control" <-
  function(sim.posterior, sim.predict, keep.mcmc.sim, quantile, threshold, inference, messages)
{
  ##
  ## Assigning default values
  ##
  if(missing(sim.posterior) || is.null(sim.posterior)) sim.posterior <- TRUE
  if(missing(sim.predict) || is.null(sim.predict)) sim.predict <- FALSE
  if(missing(keep.mcmc.sim) || is.null(keep.mcmc.sim)) keep.mcmc.sim <- TRUE
  if(missing(quantile) || is.null(quantile)) quantile.estimator <- TRUE
  else quantile.estimator <- quantile
  if(missing(threshold)) probability.estimator <- NULL 
  else probability.estimator <- threshold
  if(missing(inference) || is.null(inference)) inference <- TRUE
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  ##
  if(is.null(quantile.estimator)) quantile.estimator <- TRUE
  else{
    if(is.numeric(quantile.estimator))
      if(any(quantile.estimator < 0) | any(quantile.estimator > 1))
        stop("quantiles indicators must be numbers in the interval [0,1]\n")
    if(!inference) {
      warning("prediction not performed; quantile.estimator is set to NULL \n")
      quantile.estimator <- NULL
    }
  }
  if(!is.null(probability.estimator)){
    if(!is.numeric(probability.estimator))
      stop("threshold must be a numeric value (or vector) of cut-off value(s)\n")
    if(!inference) {
      warning("prediction not performed; probabilitites above a threshold cannot be computed\n")
      probability.estimator <- NULL
    }
  }
  res <- list(sim.posterior = sim.posterior, sim.predict = sim.predict, keep.mcmc.sim = keep.mcmc.sim,
              quantile.estimator = quantile.estimator, probability.estimator = probability.estimator,
              inference = inference, messages.screen = messages.screen)
  class(res) <- "output.geoRglm"
  return(res)
}


".output.glm.check.aux" <-
  function(output, fct)
{
  if(class(output) != "output.geoRglm"){
    if(!is.list(output))
      stop(paste(fct,": the argument output only takes a list or an output of the function output.glm.control"))
    else{
      output.names <- c("sim.posterior","sim.predict", "keep.mcmc.sim","quantile","threshold","inference","messages.screen")
      output <- .object.match.names(output,output.names)
      if(is.null(output$sim.predict)) output$sim.predict <- FALSE
      if(is.null(output$messages.screen)) output$messages.screen <- TRUE        
      output <- output.glm.control(sim.posterior = output$sim.posterior,
                                   sim.predict = output$sim.predict,
                                   keep.mcmc.sim = output$keep.mcmc.sim, quantile = output$quantile,
                                   threshold = output$threshold, inference = output$inference,
                                   messages = output$messages.screen)
    }
  }
  return(output)
}


".object.match.names" <- function(obj,names)
{
  obj.names <- names(obj)
  new.obj <- list()
  for(i in seq(along=obj)){
    n.match <- match.arg(obj.names[i], names)
    new.obj[[n.match]] <- obj[[i]]
  }
  return(new.obj)
}
