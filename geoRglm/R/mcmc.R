
"glsm.mcmc" <- function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", model, mcmc.input, messages)
{
  if(missing(geodata))
    geodata <- list(coords=coords, data=data, units.m=units.m)
  call.fc <- match.call()
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  n <- length(data)
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n)
  }
  if(any(units.m <= 0)) stop("units.m must be positive")
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  dimnames(coords) <- list(NULL, NULL)
  ##
  model <- .model.glsm.mcmc.check.aux(model)
  cov.model <- model$cov.model
  kappa <- model$kappa
  beta <- model$beta
  cov.pars <- model$cov.pars
  nugget <- model$nugget
  aniso.pars <- model$aniso.pars
  trend <- model$trend
  family <- model$family
  link <- model$link
  lambda <- model$lambda
  ##
  trend.data <- unclass(trend.spatial(trend=trend, geodata = geodata))
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  if(ncol(trend.data) != length(beta)) stop("size of beta is incompatible with trend specified") 
  ##
  ## preparing for MCMC 
  ##
  if(missing(mcmc.input)) stop("glsm.mcmc: argument mcmc.input must be given")
  mcmc.input <- .mcmc.check.aux(mcmc.input, fct="glsm.mcmc")
  ##
  mean.d <- as.vector(trend.data %*% beta)
  if(!is.null(aniso.pars)) {
    invcov <- varcov.spatial(coords = coords.aniso(coords = coords, aniso.pars = aniso.pars), cov.model = cov.model, kappa = kappa, 
                             nugget = nugget, cov.pars = cov.pars, inv = TRUE, func.inv = "cholesky",
                             try.another.decomposition = FALSE)$inverse
  }
  else {
    invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget, cov.pars = cov.pars,
                             inv = TRUE, func.inv = "cholesky", try.another.decomposition = FALSE)$inverse
  }
  ##
########################----- MCMC ------#####################
  ##
  if(model$family == "binomial"){
    simulations <- .mcmc.binom.logit(data = data, units.m = units.m, meanS = mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen)
  }
  else{
    if(lambda == 0){
      simulations <- .mcmc.pois.log(data = data, units.m = units.m, meanS = mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen)
    }
    else{
      simulations <- .mcmc.pois.boxcox(data = data, units.m = units.m, meanS = mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen, lambda=lambda)
    }
  }
  kpl.result <- list(simulations=simulations$Sdata, acc.rate=simulations$acc.rate, model=model, geodata=geodata)
  kpl.result$call <- call.fc
#######################################
  class(kpl.result) <- "glsm.mcmc"
  return(kpl.result)
}


".model.glsm.mcmc.check.aux" <-
  function(model)
{
  if(class(model)=="likGLSM"){
    model$nugget <- model$nugget.rel*model$cov.pars[1]
    return(model)
  }
  else{
    if(!is.list(model)) stop("glsm.mcmc : the argument model only takes a list ")
    model.names <- c("beta", "cov.pars", "trend", "cov.model", "kappa", "aniso.pars", "nugget", "nugget.rel", "family", "link", "lambda")
    model <- .object.match.names(model,model.names)
    if(is.null(model$beta) | !is.numeric(model$beta)) stop("glsm.mcmc : need to provide beta in the model argument")
    if(is.null(model$cov.pars)) stop("glsm.mcmc : need to provide cov.pars in the model argument")
    if(is.null(model$trend)) model$trend <- "cte"   
    if(is.null(model$cov.model)) model$cov.model <- "matern"
    model$cov.model <- match.arg(model$cov.model,
                                 choices = c("matern", "exponential","gaussian",
                                   "spherical", "circular", "cubic",
                                   "wave", "power",
                                   "powered.exponential", "cauchy", "gneiting",
                                   "gneiting.matern", "pure.nugget"))
    if(model$cov.model == "power") stop("krige.glm.control: correlation function does not exist for the power variogram")
    if(is.null(model$kappa)) model$kappa <- 0.5
    if(!is.null(model$aniso.pars))
      if(length(model$aniso.pars) != 2 | !is.numeric(model$aniso.pars))
        stop("glsm.mcmc : anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(is.null(model$nugget)){
      if(!is.null(model$nugget.rel)) model$nugget <- model$nugget.rel**model$cov.pars[1]
      else model$nugget <- 0
    }
    if(is.null(model$family)) stop("glsm.mcmc : need to provide family in the model argument")
    family <- match.arg(model$family, choices = c("poisson", "binomial"))
    if(family=="poisson"){
      if(is.null(model$lambda)){
        if(is.null(model$link)){
          model$lambda <- 0
          model$link <- "log" 
        }
        if(model$link == "canonical") model$link <- "log"
        if(model$link == "boxcox") stop("glsm.mcmc : need to provide lambda in the model argument")
        if(model$link == "log") model$lambda <- 0
        if(model$link == "id") model$lambda <- 1
      }
      if(!is.null(model$lambda)){
        if(is.null(model$link)){
          if(model$lambda > 0) model$link <- "boxcox"
          if(model$lambda == 0) model$link <- "log"
        }
        if(!is.null(model$link)){
          if(model$link == "canonical"){
            if(model$lambda > 0) model$link <- "boxcox"
            if(model$lambda == 0) model$link <- "log"
          }
          if(model$link == "boxcox" & model$lambda == 0) model$link <- "log"
          if(model$link == "log" & model$lambda > 0) warning("glsm.mcmc : value of argument lambda will be ignored since it is inconsistent with log-link ")
          if(model$link == "id" & model$lambda < 1) warning("glsm.mcmc : value of argument lambda will be ignored since it is inconsistent with identity-link ")
        }
      }
    }
    if(family=="binomial"){
      if(!is.null(model$link)){
        if(model$link != "logit" & model$link != "canonical") stop("glsm.mcmc : only the canonical logit link function is implemented ")
      }
      model$link <- "logit"
      model$lambda <- NULL
    }
    return(model)
  }
}


"create.mcmc.coda" <- function(x, mcmc.input)
{
  if(exists("mcmc.input")){
    if(!is.list(mcmc.input)) stop(" mcmc.input must be given as a list ")
    if(is.null(mcmc.input$thin)) thin <- 10
    else thin <- mcmc.input$thin
    if(is.null(mcmc.input$burn.in)) st.val <- 1
    else st.val <- mcmc.input$burn.in + 1
  }
  else {
    thin <- 10
    st.val <- 1
  }
  if(class(x)=="glsm.mcmc"){
    n <- nrow(x$simulations)
    S.names <- rep(NA,n)
    for(i in 1:n) S.names[i] <- paste("S[",i,"]",sep="")
    temp <- t(x$simulations)
    colnames(temp) <- S.names
    return(coda::mcmc(data=temp, start = st.val, thin=thin))
  }
  else{
    if(class(x)=="glm.krige.bayes"){
      n <- nrow(x$posterior$simulations)
      S.names <- rep(NA,n)
      for(i in 1:n) S.names[i] <- paste("S[",i,"]",sep="")
      if(!is.null(x$prior$phi$status) && x$prior$phi$status =="fixed"){
        temp <- t(x$posterior$simulations)
        colnames(temp) <- S.names
      }
      else{
        temp <- cbind(x$posterior$phi$sample,t(x$posterior$simulations))
        colnames(temp) <- c("phi", S.names)
      }
      return(coda::mcmc(data=temp, start = st.val, thin=thin))
    }
    else{
      if(is.matrix(x)) return(coda::mcmc(data=t(x), start = st.val, thin=thin))
      if(is.vector(x)) return(coda::mcmc(data=x, start = st.val,thin=thin))
    }
  }
}

 

"asympvar" <- 
  function(timeseries, type = "mon", lag.max = 100, messages)
{
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  ##
  if(is.vector(timeseries)) n.series <- 1
  else n.series <- nrow(timeseries) 
  if(type == "mon" | type == "all" | type == "pos") {
    if(messages.screen & type == "mon")
      cat(paste("calculating the initial monotone sequence estimate \n"))
    if(messages.screen & type == "pos") 
      cat(paste("calculating the initial positive sequence estimate \n"))
    if(messages.screen & type == "all") 
      cat(paste("calculating the initial positive sequence estimate, and the initial monotone sequence estimate \n"))      
  }
  else stop("Must specify type as either: mon, pos or all")
  len.Gamma <- floor(lag.max/2)-1
  if(is.vector(timeseries)){
     if(all(is.finite(timeseries)) && var(timeseries)>0){
       asy.gamma <- acf(timeseries, type = "covariance", plot = FALSE, lag.max = lag.max)$acf
       asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:len.Gamma))])
       asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:len.Gamma))])
       asy.Gamma <- asy.gamma1 + asy.gamma2
       ##--------- initial monotone sequence estimate -----------#
       kmaxpos <- min(c(which(asy.Gamma<0)-1, len.Gamma))
       if(type == "all" | type =="mon"){
         kmax <- min(c(which(diff(asy.Gamma)>0),kmaxpos))
         monvarest <- 2*sum(asy.Gamma[seq(length=kmax)])-asy.gamma[1]
         if(kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
       }
       ##--------- initial positive sequence estimate -----------#
       if(type == "pos" | type == "all"){
         posvarest <- 2*sum(asy.Gamma[seq(length=kmaxpos)])-asy.gamma[1]
         if (kmaxpos == len.Gamma) warning("value of argument lag.max is not suffiently long")
       } 
     } else{
       if(type == "all" | type =="mon") monvarest <- NA
       if(type == "pos" | type == "all") posvarest <- NA
     }
  }
  else{
     if(type == "all" | type == "pos") posvarest <- rep(1,n.series)
     if(type == "all" | type == "mon") monvarest <- rep(1,n.series)
     for(i in seq(length=n.series)){
        if(all(is.finite(timeseries[i,])) && var(timeseries[i,])>0){     
          asy.gamma <- acf(timeseries[i,], type = "covariance", plot = FALSE, lag.max = lag.max)$acf
          asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:len.Gamma))])
          asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:len.Gamma))])
          asy.Gamma <- asy.gamma1 + asy.gamma2
          ##--------- initial monotone sequence estimate -----------#
          kmaxpos <- min(c(which(asy.Gamma<0)-1, len.Gamma))
          if(type == "all" | type =="mon"){
            kmax <- min(c(which(diff(asy.Gamma)>0),kmaxpos))
            monvarest[i] <- 2*sum(asy.Gamma[seq(length=kmax)])-asy.gamma[1]
            if(kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
          }
          ##--------- initial positive sequence estimate -----------#
          if(type == "pos" | type == "all"){
            posvarest[i] <- 2*sum(asy.Gamma[seq(length=kmaxpos)])-asy.gamma[1]
            if (kmaxpos == len.Gamma) warning("value of argument lag.max is not suffiently long")
          }
        } else{
          if(type == "all" | type =="mon") monvarest[i] <- NA
          if(type == "pos" | type == "all") posvarest[i] <- NA
        }
     }
  }
  if(type == "pos") return(posvarest)
  if(type == "all") return(list(posvarest = posvarest, monvarest = monvarest))
  if(type == "mon") return(monvarest)
}

#### consider vectorising the while loops.

