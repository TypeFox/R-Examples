"krweights" <-
  function(coords, locations, krige)
{
  if(is.vector(coords)) coords <- cbind(coords, 0)
  coords <- as.matrix(coords)
  if(is.vector(locations)) {
    if(length(locations) == 2)
      locations <- t(as.matrix(locations))
    else
      locations <- as.matrix(cbind(locations, 0))
  }
  else locations <- as.matrix(locations)
  ## setting up the model
  if(missing(krige))
    krige <- krige.control()
  else{
    ##    if(is.null(class(krige)) || class(krige) != "krige.geoR"){
    if(length(class(krige)) == 0 || class(krige) != "krige.geoR"){
      if(!is.list(krige))
        stop("krige.conv: the argument krige only takes a list or an output of the function krige.control")
      else{
        krige.names <- c("type.krige","trend.d","trend.l","obj.model",
                         "beta","cov.model", "cov.pars",
                         "kappa","nugget","micro.scale","dist.epsilon",
                         "lambda","aniso.pars")
        krige.user <- krige
        krige <- list()
        if(length(krige.user) > 0){
          for(i in 1:length(krige.user)){
            n.match <- match.arg(names(krige.user)[i], krige.names)
            krige[[n.match]] <- krige.user[[i]]
          }
        }
        if(is.null(krige$type.krige)) krige$type.krige <- "ok"  
        if(is.null(krige$trend.d)) krige$trend.d <-  "cte"
        if(is.null(krige$trend.l)) krige$trend.l <-  "cte"
        if(is.null(krige$obj.model)) krige$obj.model <-  NULL
        if(is.null(krige$beta)) krige$beta <- NULL 
        if(is.null(krige$cov.model)) krige$cov.model <- "matern"  
        if(is.null(krige$cov.pars))
          stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
        if(is.null(krige$kappa)) krige$kappa <-  0.5
        if(is.null(krige$nugget)) krige$nugget <-  0
        if(is.null(krige$micro.scale)) krige$micro.scale <- 0  
        if(is.null(krige$dist.epsilon)) krige$dist.epsilon <-  1e-10
        if(is.null(krige$aniso.pars)) krige$aniso.pars <- NULL  
        if(is.null(krige$lambda)) krige$lambda <- 1 
        krige <- krige.control(type.krige = krige$type.krige,
                               trend.d = krige$trend.d,
                               trend.l = krige$trend.l,
                               obj.model = krige$obj.model,
                               beta = krige$beta,
                               cov.model = krige$cov.model,
                               cov.pars = krige$cov.pars,
                               kappa = krige$kappa,
                               nugget = krige$nugget,
                               micro.scale = krige$micro.scale,
                               dist.epsilon = krige$dist.epsilon, 
                               aniso.pars = krige$aniso.pars,
                               lambda = krige$lambda)
        
      }
    }
  }
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  # lambda <- krige$lambda
  # beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  # micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  ##
  ## Anisotropy correction (this should be placed AFTER trend.d/trend.l
  ##
  if(!is.null(aniso.pars)) {
    if(abs(aniso.pars[2] - 1) > 0.0001){
      coords <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
      locations <- coords.aniso(coords = locations, aniso.pars = aniso.pars)
    }
  }
  ##
  ## starting kriging calculations
  ##
  cov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                        kappa = kappa, nugget = nugget, cov.pars = cov.pars)$varcov
  v0 <- cov.spatial(loccoords(coords = coords, locations = locations),
                    cov.model = cov.model, cov.pars=cov.pars,
                    kappa = kappa)
  if(krige$type == "ok"){
    cov <- cbind(cov, 1)
    cov <- rbind(cov, c(rep(1, nrow(coords)),0))
    v0 <- c(v0, 1)
  }
  wei <- drop(solve(cov, v0))
  if(krige$type == "ok"){
    if(is.vector(wei)) wei <- wei[-(nrow(coords)+1)]
    else wei <- wei[-(ncol(coords)+1),]
  }
  return(wei)
}
