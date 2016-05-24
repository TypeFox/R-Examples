

".mcmc.aux" <- 
  function(z, data, meanS, QQ, Htrunc, S.scale, nsim, thin, QtivQ)
{
  ##
###### ------------------------ doing the mcmc-steps ----------- ############# 
  ##
  n <- length(data)
  result <- .C("mcmc1poislog",
               as.integer(n),
               z = as.double(z),
               S = as.double(rep(0, nsim * n)),
               as.double(data),
               as.double(meanS),
               as.double(as.vector(t(QQ))),
               as.double(as.vector(QtivQ)),
               as.double(rnorm(n * nsim * thin) * sqrt(S.scale)),
               as.double(runif(nsim * thin)),
               as.double(Htrunc),
               as.double(S.scale),
               as.integer(nsim),
               as.integer(thin),
               acc.rate = as.double(1), PACKAGE = "geoRglm")[c("z", "S", "acc.rate")]
  attr(result$S, "dim") <- c(n, nsim)
  return(result)
}


".mcmc.pois.log" <- 
  function(data, units.m, meanS, invcov, mcmc.input, messages.screen)
{
  ## This is the MCMC engine for the spatial Poisson log Normal model ----
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n)
      Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  QQ <- t(chol(solve(invcov + diag(data))))
  sqrtdataQ <- sqrt(data)*QQ 
  QtivQ <- diag(n)-crossprod(sqrtdataQ)
  if(any(mcmc.input$S.start=="default")) {
    z <- as.vector(solve(QQ,ifelse(data > 0, log(data), -1.96) - meanS - log(units.m)))
  }
  else{
    if(any(mcmc.input$S.start=="random")) z <- rnorm(n)
    else{
      if(is.numeric(mcmc.input$S.start)){
        if(length(mcmc.input$S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
        else z <- as.vector(solve(QQ,mcmc.input$S.start))
      }
      else  stop(" S.start must be a vector of same dimension as data ")
    }
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
## ---------------- burn-in ----------------- ######### 
  if(burn.in > 0) {
    mcmc.output <- .mcmc.aux(z, data, meanS + log(units.m), QQ, Htrunc, S.scale, 1, burn.in, QtivQ)
    if(messages.screen) cat(paste("burn-in = ", burn.in, " is finished. Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate.burn.in <- c(burn.in, mcmc.output$acc.rate)
  }
  else mcmc.output <- list(z = z)
##### ---------- sampling periode ----------- ###### 
  if(n.iter <= 1000) {
    n.temp <- round(n.iter/thin)
    n.turn <- 1
  }
  else {
    n.temp <- round(1000/thin)
    n.turn <- round(n.iter/1000)
  }
  n.sim <- n.turn * n.temp
  Sdata <- matrix(NA, n, n.sim)
  acc.rate <- matrix(NA, n.turn, 2)
  for(i in seq(length=n.turn)) {
    mcmc.output <- .mcmc.aux(mcmc.output$z, data, meanS + log(units.m), QQ, Htrunc, S.scale, n.temp, thin, QtivQ)
    Sdata[, seq((n.temp * (i - 1) + 1),(n.temp * i))] <- mcmc.output$S+meanS
    if(messages.screen) cat(paste("iter. numb.", i * n.temp * thin+burn.in, " : Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate[i,1] <-  i * n.temp * thin
    acc.rate[i,2] <- mcmc.output$acc.rate
  }
  if(messages.screen) cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  if(burn.in > 0) acc.rate <- rbind(acc.rate.burn.in,acc.rate)
  colnames(acc.rate) <- c("iter.numb", " Acc.rate")
#########
  return(list(Sdata=Sdata, acc.rate=acc.rate))
}


".mcmc.boxcox.aux" <- 
  function(z, data, units.m, meanS, QQ, Htrunc, S.scale, nsim, thin, QtivQ, lambda)
{
  ##
###### ------------------------ doing the mcmc-steps ----------- ############# 
  ##
  n <- length(data)
  result <- .C("mcmc1poisboxcox",
               as.integer(n),
               z = as.double(z),
               S = as.double(rep(0, nsim * n)),
               as.double(data),
               as.double(units.m),
               as.double(meanS),
               as.double(as.vector(t(QQ))),
               as.double(as.vector(QtivQ)),
               as.double(rnorm(n * nsim * thin) * sqrt(S.scale)),
               as.double(runif(nsim * thin)),
               as.double(Htrunc),
               as.double(S.scale),
               as.integer(nsim),
               as.integer(thin),
               as.double(lambda),
               acc.rate = as.double(1), PACKAGE = "geoRglm")[c("z", "S", "acc.rate")]
  attr(result$S, "dim") <- c(n, nsim)
  return(result)
}

".mcmc.pois.boxcox" <- 
  function(data, units.m, meanS, invcov, mcmc.input, messages.screen, lambda)
{
  ## This is the MCMC engine for the spatial Poisson - Normal model with link from the box-cox-family ----
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  fisher.l <- ifelse(data>0,data^(1-2*lambda)*units.m^(2*lambda),0)
  QQ <- t(chol(solve(invcov + diag(fisher.l)))) 
  sqrtfiQ <- sqrt(fisher.l)*QQ 
  QtivQ <- diag(n)-crossprod(sqrtfiQ)
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, (data/units.m)^lambda-1, -1.96)/lambda - meanS )       
    z <- as.vector(solve(QQ,S))
  }
  else{
    if(any(mcmc.input$S.start=="random")) z <- rnorm(n)
    else{
      if(is.numeric(mcmc.input$S.start)){
        if(length(mcmc.input$S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
        else z <- as.vector(solve(QQ,mcmc.input$S.start))
      }
      else  stop(" S.start must be a vector of same dimension as data ")
    }
  }
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n)
      Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  ## ---------------- burn-in ----------------- ######### 
  if(burn.in > 0) {
    mcmc.output <- .mcmc.boxcox.aux(z, data, units.m, meanS, QQ, Htrunc, S.scale, 1, burn.in, QtivQ, lambda)
    if(messages.screen) cat(paste("burn-in = ", burn.in, " is finished. Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate.burn.in <- c(burn.in, mcmc.output$acc.rate)
  }
  else mcmc.output <- list(z = z)
##### ---------- sampling periode ----------- ###### 
  if(n.iter <= 1000) {
    n.temp <- round(n.iter/thin)
    n.turn <- 1
  }
  else {
    n.temp <- round(1000/thin)
    n.turn <- round(n.iter/1000)
  }
  n.sim <- n.turn * n.temp
  Sdata <- matrix(NA, n, n.sim)
  acc.rate <- matrix(NA, n.turn, 2)
  for(i in seq(length=n.turn)){
    mcmc.output <- .mcmc.boxcox.aux(mcmc.output$z, data, units.m, meanS, QQ, Htrunc, S.scale, n.temp, thin, QtivQ, lambda)
    Sdata[, seq((n.temp * (i - 1) + 1),(n.temp * i))] <- mcmc.output$S+meanS
    if(messages.screen) cat(paste("iter. numb.", i * n.temp * thin+burn.in, " : Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate[i,1] <-  i * n.temp * thin
    acc.rate[i,2] <- mcmc.output$acc.rate
  }
  if(messages.screen) cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  if(burn.in > 0) acc.rate <- rbind(acc.rate.burn.in,acc.rate)
  colnames(acc.rate) <- c("iter.numb", " Acc.rate")
  remove("z")
#########
  return(list(Sdata=Sdata, acc.rate=acc.rate))
}

"krige.glm.control" <-
  function (type.krige = "sk", trend.d = "cte", trend.l = "cte", obj.model = NULL, beta, cov.model, cov.pars, kappa,
            nugget, micro.scale, dist.epsilon = 1e-10, aniso.pars, lambda)
{
  if(type.krige != "ok" & type.krige != "OK" & type.krige != "o.k." & type.krige != "O.K." & type.krige != "sk" & type.krige != "SK" & type.krige != "s.k." & type.krige != "S.K.")
    stop("pois.krige: wrong option in the argument type.krige. It should be \"sk\" or \"ok\"(if ordinary or simple kriging is to be performed)")
  if(type.krige=="OK" | type.krige=="O.K." |type.krige=="o.k.")
    type.krige <- "ok"
  if(type.krige=="SK" | type.krige=="S.K." |type.krige=="s.k.")
    type.krige <- "sk"
  ##
  if(!is.null(obj.model)){
    if(missing(beta)) beta <- obj.model$beta
    if(missing(cov.model)) cov.model <- obj.model$cov.model
    if(missing(cov.pars)) cov.pars <- obj.model$cov.pars
    if(missing(kappa)) kappa <- obj.model$kappa
    if(missing(nugget)) nugget <- obj.model$nugget
    if(missing(micro.scale)) micro.scale <- nugget
    if(missing(lambda)) lambda <- obj.model$lambda
    if(missing(aniso.pars)) aniso.pars <- obj.model$aniso.pars
  }
  else{
    if(missing(beta)) beta <- NULL
    if(missing(cov.model)) cov.model <- "matern"
    if(missing(cov.pars)) stop("covariance parameters (sigmasq and phi) should be provided")
    if(missing(kappa)) kappa <- 0.5
    if(missing(nugget)) nugget <- 0
    if(missing(micro.scale)) micro.scale <- nugget
    if(missing(lambda)) lambda <- 0
    if(missing(aniso.pars)) aniso.pars <- NULL
  }
  ##
  if(type.krige == "sk")
    if(is.null(beta) | !is.numeric(beta))
      stop(" argument beta must be provided in order to perform simple kriging")
  if(micro.scale > nugget)
    stop(" micro.scale must be in the interval [0, nugget]")
  if(!is.null(aniso.pars))
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop(" anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
  ##
  if(inherits(trend.d, "formula") | inherits(trend.l, "formula")){
    if(!inherits(trend.d, "formula") | !inherits(trend.l, "formula"))
      stop(" trend.d and trend.l must have similar specification")
  }
  else{
    if((class(trend.d)=="trend.spatial") & (class(trend.l)=="trend.spatial")){
      if(ncol(trend.d) != ncol(trend.l))
        stop("pois.krige: trend.d and trend.l do not have the same number of columns")
    }
    else{
      if(trend.d != trend.l)
        stop(" trend.l is different from trend.d")
    }
  }
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential","gaussian",
                           "spherical", "circular", "cubic",
                           "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(cov.model == "power") stop("krige.glm.control: correlation function does not exist for the power variogram")
  res <- list(type.krige = type.krige,
              trend.d = trend.d, trend.l = trend.l, 
              beta = beta,
              cov.model = cov.model, 
              cov.pars = cov.pars, kappa = kappa,
              nugget = nugget,
              micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
              aniso.pars = aniso.pars, lambda = lambda)
  class(res) <- "krige.geoRglm"
  return(res)
}

".krige.glm.check.aux" <-
  function(krige,fct)
{
  if(class(krige) != "krige.geoRglm"){
    if(!is.list(krige))
      stop(paste(fct,": the argument krige only takes a list or an output of the function krige.glm.control"))
    else{
      krige.names <-c("type.krige","trend.d","trend.l","obj.model","beta","cov.model",
                      "cov.pars","kappa","nugget","micro.scale","dist.epsilon","lambda","aniso.pars")
      krige <- .object.match.names(krige,krige.names)
      if(is.null(krige$type.krige)) krige$type.krige <- "sk"  
      if(is.null(krige$trend.d)) krige$trend.d <-  "cte"
      if(is.null(krige$trend.l)) krige$trend.l <-  "cte"
      if(is.null(krige$cov.model)) krige$cov.model <- "matern"
      if(is.null(krige$kappa)) krige$kappa <-  0.5
      if(is.null(krige$nugget)) krige$nugget <-  0
      if(is.null(krige$micro.scale)) krige$micro.scale <- krige$nugget
      if(is.null(krige$dist.epsilon)) krige$dist.epsilon <-  1e-10
      krige <- krige.glm.control(type.krige = krige$type.krige,	
                                 trend.d = krige$trend.d, trend.l = krige$trend.l,
                                 obj.model = krige$obj.model,
                                 beta = krige$beta, cov.model = krige$cov.model,
                                 cov.pars = krige$cov.pars, kappa = krige$kappa,
                                 nugget = krige$nugget, micro.scale = krige$micro.scale,
                                 dist.epsilon = krige$dist.epsilon, 
                                 aniso.pars = krige$aniso.pars)
    }
  }
  return(krige)
}



"pois.krige" <- 
function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", locations = NULL,  borders, mcmc.input, krige, output)
{
  if(missing(geodata))
    geodata <- list(coords=coords, data=data, units.m=units.m)
  if(missing(borders))
    borders <- geodata$borders
  call.fc <- match.call()
  n <- length(data)
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n)
  }
  if(any(units.m <= 0)) stop("units.m must be positive")
  if(missing(krige)) stop("must provide object krige")
  krige <- .krige.glm.check.aux(krige,fct="pois.krige")
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  trend.d <- krige$trend.d
  trend.l <- krige$trend.l
  dist.epsilon <- krige$dist.epsilon
  lambda <- krige$lambda
  if(krige$type.krige == "ok") beta.prior <- "flat"
  if(krige$type.krige == "sk") beta.prior <- "deg"
  if(missing(output)) output <- output.glm.control()
  output <- .output.glm.check.aux(output, fct="pois.krige")
  sim.predict <- output$sim.predict
  messages.screen <- output$messages.screen
  ##
  if(is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  dimnames(coords) <- list(NULL, NULL)
  ##
  if(is.null(locations)) {
    if(messages.screen) cat(paste("locations need to be specified for prediction; prediction not performed \n"))
  }
  else {
    locations <- .geoR.check.locations(locations)
    if(is.null(trend.l))
      stop("trend.l needed for prediction")
    ## Checking for 1D prediction 
    if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
      krige1d <- TRUE
    else krige1d <- FALSE
  }
  trend.data <- unclass(trend.spatial(trend=trend.d, geodata = geodata))
  beta.size <- ncol(trend.data)
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  if(beta.prior == "deg")
    if(beta.size != length(beta))
      stop("size of mean vector is incompatible with trend specified") 
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ## preparing for MCMC 
  ##
  if(missing(mcmc.input)) stop("pois.krige: argument mcmc.input must be given")
  mcmc.input <- .mcmc.check.aux(mcmc.input, fct="pois.krige")
  ##
  if(beta.prior == "deg") mean.d <-  as.vector(trend.data %*% beta)
  else mean.d <- rep(0,n)
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
  if(beta.prior == "flat") {
    ivtt <- invcov%*%trend.data
    invcov <- invcov-ivtt%*%.solve.geoR(crossprod(trend.data, ivtt),t(ivtt))
  }
  if(lambda == 0){
    intensity <- .mcmc.pois.log(data = data, units.m = units.m, meanS = mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen)
    acc.rate <- intensity$acc.rate
    intensity <- exp(intensity$Sdata)
  }
  else{
    intensity <- .mcmc.pois.boxcox(data=data, units.m=units.m, meanS=mean.d, invcov=invcov, mcmc.input=mcmc.input, messages.screen=messages.screen, lambda=lambda)
    acc.rate <- intensity$acc.rate
    intensity <- .BC.inv(intensity$Sdata, lambda)    
  }
  ##
  ##------------------------------------------------------------
######################## ---- prediction ----- #####################
  if(!is.null(locations)) {
    krige <- list(type.krige = krige$type.krige, beta = beta, trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, 
                  cov.pars = cov.pars, kappa = kappa, nugget = nugget, micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
                  aniso.pars = aniso.pars, lambda = lambda)
    kpl.result <- .krige.conv.extnd(data = intensity, coords = coords, locations = locations, borders=borders, krige = krige,
                                   output = list(n.predictive = ifelse(sim.predict,1,0), signal = TRUE, messages = FALSE))
    remove(list = c("intensity"))
    kpl.result$krige.var <- rowMeans(kpl.result$krige.var) + apply(kpl.result$predict, 1, var) 
    if(nrow(locations) > 1) kpl.result$mcmc.error <- sqrt(asympvar(kpl.result$predict)/ncol(kpl.result$predict))
    else kpl.result$mcmc.error <- sqrt(asympvar(as.vector(kpl.result$predict), messages = FALSE)/length(as.vector(kpl.result$predict)))
    kpl.result$predict <- rowMeans(kpl.result$predict)
    if(beta.prior == "flat") {
      kpl.result$beta.est <- rowMeans(kpl.result$beta)
      names(kpl.result$beta.est) <- beta.names
    }
    kpl.result$beta <- NULL
  }
  else{
    if(beta.prior == "flat") {
      ## GLS
      beta.est <- .solve.geoR(crossprod(trend.data, ivtt),t(ivtt))%*%rowMeans(log(intensity))
      kpl.result <- list(intensity=intensity, beta.est = beta.est, acc.rate=acc.rate)
    }
    else kpl.result <- list(intensity=intensity, acc.rate=acc.rate)
  }
  kpl.result$call <- call.fc
#######################################
  attr(kpl.result, "prediction.locations") <- call.fc$locations
  if(!is.null(locations)) attr(kpl.result, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  if(!is.null(call.fc$borders)) attr(kpl.result, "borders") <- call.fc$borders
  class(kpl.result) <- "kriging"
  return(kpl.result)
}
