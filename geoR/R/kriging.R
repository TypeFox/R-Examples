"krige.conv" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            locations, borders, krige, output)
{
  if(missing(geodata))
    geodata <- list(coords = coords, data = data)
  if(missing(borders))
    borders <- geodata$borders
  locations <- .check.locations(locations)
  ## Checking for 1D prediction 
  krige1d <- FALSE
  if(length(locations) > 2){
  if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
    krige1d <- TRUE
  }
  call.fc <- match.call()
  base.env <- sys.frame(sys.nframe())
  ##
  ## reading input
  ##
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
  lambda <- krige$lambda
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  ##
  ## reading output options
  ##
  if(missing(output))
    output <- output.control()
  else{
    ##    if(is.null(class(output)) || class(output) != "output.geoR"){
    if(length(class(krige)) == 0 || class(output) != "output.geoR"){
      if(!is.list(output))
        stop("krige.conv: the argument output can take only a list or an output of the function output.control")
      else{
        output.names <- c("n.posterior","n.predictive","moments","n.back.moments","simulations.predictive",
                          "mean.var","quantile","threshold","signal","messages.screen")
        output.user <- output
        output <- list()
        if(length(output.user) > 0){
          for(i in 1:length(output.user)){
            n.match <- match.arg(names(output.user)[i], output.names)
            output[[n.match]] <- output.user[[i]]
          }
        }
        if(is.null(output$n.posterior)) output$n.posterior <- 1000 
        if(is.null(output$n.predictive)) output$n.predictive <- NULL
        if(is.null(output$moments)) output$moments <- TRUE
        if(is.null(output$n.back.moments)) output$n.back.moments <- 1000 
        if(is.null(output$simulations.predictive)){
          if(is.null(output$n.predictive)) output$simulations.predictive <- NULL
          else
            output$simulations.predictive <- ifelse(output$n.predictive > 0, TRUE, FALSE)
        }
        if(is.null(output$mean.var)) output$mean.var <- NULL
        if(is.null(output$quantile)) output$quantile <- NULL
        if(is.null(output$threshold)) output$threshold <- NULL
        if(is.null(output$sim.means)) output$sim.means <- NULL
        if(is.null(output$sim.vars)) output$sim.vars <- NULL
        if(is.null(output$signal)) output$signal <- NULL
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior,
                                 n.predictive = output$n.predictive,
                                 moments = output$moments,
                                 n.back.moments = output$n.back.moments, 
                                 simulations.predictive = output$simulations.predictive,
                                 mean.var = output$mean.var,
                                 quantile = output$quantile,
                                 threshold = output$threshold,
                                 sim.means = output$sim.means,
                                 sim.vars = output$sim.vars,
                                 signal = output$signal,
                                 messages = output$messages.screen)
      }
    }
  }
  signal <- ifelse(is.null(output$signal), FALSE, output$signal)
  messages.screen <- output$messages.screen
  n.predictive <- output$n.predictive
  n.back.moments <- output$n.back.moments
  ##
  n.predictive <- ifelse(is.null(n.predictive), 0, n.predictive)
  simulations.predictive <- ifelse(is.null(output$simulations.predictive), FALSE, TRUE)
  #keep.simulations <- ifelse(is.null(output$keep.simulations), TRUE, FALSE)
  mean.estimator <- output$mean.estimator
  sim.means <- output$sim.means
  if(is.null(sim.means))
    sim.means <- ifelse(simulations.predictive, TRUE, FALSE)
  sim.vars <- output$sim.vars
  if(is.null(sim.vars)) sim.vars <- FALSE
  if(is.null(mean.estimator) & simulations.predictive)
    mean.estimator <- TRUE
  quantile.estimator <- output$quantile.estimator
  probability.estimator <- output$probability.estimator
  if(!is.null(probability.estimator)){
    if(length(probability.estimator) > 1 &
       length(probability.estimator) != nrow(locations))
      stop("krige.conv: probability.estimator must either have length 1, or have length = nrow(locations)\n")
  }
  if(simulations.predictive & n.predictive == 0) n.predictive <- 1000
  ##
  ## checking input
  ##
  if(krige$type.krige == "ok") beta.prior <- "flat"
  if(krige$type.krige == "sk") beta.prior <- "deg"
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("krige.conv: coordinates provided as a vector, assuming one spatial dimension")
  }
  coords <- as.matrix(coords)
  ##
  ## selecting locations inside the borders 
  ## and this will also used later for 
  ## values of trend.l if the case
  ##
  if(!is.null(borders)){
    nloc0 <- nrow(locations)
    ind.loc0  <- .geoR_inout(locations, borders)
#    locations <- locations.inside(locations, borders)
    if(nrow(locations) == 1)
    locations <- locations[ind.loc0,,drop=FALSE]
    else
    locations <- locations[ind.loc0,,drop=TRUE]
    if(nrow(locations) == 0)
      stop("\nkrige.conv: there are no prediction locations inside the borders")
    if(messages.screen)
      cat("krige.conv: results will be returned only for prediction locations inside the borders\n")
  }
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  ##
  ## building the trend matrix
  ##
  if(messages.screen){
    if(mode(krige$trend.d) == "numeric")
      cat("krige.conv: model with covariates matrix provided by the user")
    else
      cat(switch(as.character(krige$trend.d)[1],
                 "cte" = "krige.conv: model with constant mean",
                 "1st" = "krige.conv: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = "krige.conv: model with mean given by a 2nd order polynomial on the coordinates",
                 "krige.conv: model with mean defined by covariates provided by the user"))
    cat("\n")
  }
  if(class(krige$trend.d) == "trend.spatial")
    trend.d <- unclass(krige$trend.d)
  else
    trend.d <- unclass(trend.spatial(trend=krige$trend.d, geodata = geodata))
  if (nrow(trend.d) != nrow(coords)) 
    stop("coords and trend.d have incompatible sizes")
  beta.size <- ncol(trend.d)
  if(beta.prior == "deg")
    if(beta.size != length(beta))
      stop("size of mean vector is incompatible with trend specified") 
  if(class(krige$trend.l) == "trend.spatial")
    trend.l <- unclass(krige$trend.l)
  else
    trend.l <- unclass(trend.spatial(trend=krige$trend.l,
                                     geodata = list(coords =
    locations)))
  if(!is.null(borders))
    if(nrow(trend.l) == nloc0)
      trend.l <- trend.l[ind.loc0,,drop=FALSE]
  if (nrow(trend.l) != nrow(locations)) 
    stop("locations and trend.l have incompatible sizes")
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ## Anisotropy correction (this should be placed AFTER trend.d/trend.l
  ##
  if(!is.null(aniso.pars)) {
#    if((abs(aniso.pars[1]) > 0.001) & (abs(aniso.pars[2] - 1) > 0.001)){
    if(abs(aniso.pars[2] - 1) > 0.0001){
      if(messages.screen)
        cat("krige.conv: anisotropy correction performed\n")
      coords <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
      locations <- coords.aniso(coords = locations, aniso.pars = aniso.pars)
    }
  }
  ##
  ## Box-Cox transformation
  ##
  if(!isTRUE(all.equal(lambda, 1))) {
    if(messages.screen)
      cat("krige.conv: performing the Box-Cox data transformation\n")
    data <- BCtransform(x=data, lambda = lambda)$data
  }
  ## 
  ## setting covariance parameters
  ##
  if(is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    cpars <- c(1, phi)
  }
  else {
    stop("current version of krige.conv does not accept nested covariance models\n") 
    ##    sigmasq <- cov.pars[, 1]
    ##    phi <- cov.pars[, 2]
    ##    cpars <- cbind(1, phi)
  }
  ##  sill.partial <- micro.scale + sum(sigmasq)
  sill.partial <- sum(sigmasq)
  if(sill.partial < 1e-16){
    tausq.rel <- 0
    tausq.rel.micro <- 0
  }
  else{
    tausq.rel <- nugget/sum(sigmasq)
    tausq.rel.micro <- micro.scale/sum(sigmasq)
  }
  n <- length(data)
  ni <- nrow(trend.l)
  ##
  ## starting kriging calculations
  ##
  kc <- list()
  Vcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                         kappa = kappa, nugget = tausq.rel,
                         cov.pars = cpars)$varcov
  ivtt <- solve(Vcov,trend.d)
  ttivtt <- crossprod(ivtt,trend.d)
  if(beta.prior == "flat")
    beta.flat <- drop(solve(ttivtt,crossprod(ivtt,as.vector(data))))
  remove("ivtt")
  v0 <- loccoords(coords = coords, locations = locations)
  if(n.predictive > 0){
    ## checking if there are data points coincident with prediction locations
    loc.coincide <- apply(v0, 2, function(x, min.dist){any(x < min.dist)},
                          min.dist=krige$dist.epsilon)
    if(any(loc.coincide)) loc.coincide <- (1:ni)[loc.coincide]
    else loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data, dist.eps){return(data[x < dist.eps])}
      data.coincide <- apply(v0[, loc.coincide, drop=FALSE], 2, temp.f, data=data, dist.eps=krige$dist.epsilon)
    }
    else data.coincide <- NULL
  }
  else remove("locations")
  ## using nugget interpreted as microscale variation or measurement error
  nug.factor <- ifelse(signal, tausq.rel.micro, tausq.rel)
  ## covariances between data and prediction locations
  #v0 <- ifelse(v0 < krige$dist.epsilon, 1+nug.factor,
  #             cov.spatial(obj = v0, cov.model = cov.model, 
  #                         kappa = kappa, cov.pars = cpars))
  ind.v0 <- which(v0<krige$dist.epsilon)
  v0 <- cov.spatial(obj = v0, cov.model = cov.model, 
                    kappa = kappa, cov.pars = cpars)
  v0[ind.v0] <- 1+nug.factor
  ivv0 <- solve(Vcov,v0)
  ##tv0ivv0 <- diag(crossprod(v0,ivv0))
  tv0ivv0 <- colSums(v0 * ivv0)
  remove("Vcov", "ind.v0")
  b <- crossprod(cbind(data,trend.d),ivv0)
  if(n.predictive == 0) remove("v0", "ivv0")
#  gc(verbose = FALSE)
  tv0ivdata <- drop(b[1,])
  b <- t(trend.l) -  b[-1,, drop=FALSE]
  if(beta.prior == "deg") {
    kc$predict <- tv0ivdata + drop(crossprod(b,beta))
    kc$krige.var <- sill.partial * drop(1+nug.factor - tv0ivv0)
    kc$beta.est <- "Simple kriging performed (beta provided by user)"
  }
  if(beta.prior == "flat"){
    kc$predict <- tv0ivdata + drop(crossprod(b,beta.flat))
    #bitb <- drop(diag(crossprod(b,solve(ttivtt,b))))
    bitb <- colSums(b * solve(ttivtt,b))
    kc$krige.var <- sill.partial * drop(1+nug.factor - tv0ivv0 + bitb)
    kc$beta.est <- beta.flat
    names(kc$beta.est) <- beta.names
    if(n.predictive == 0) remove("bitb")
  }
  if(n.predictive == 0) remove("b","tv0ivv0")
  remove("tv0ivdata")
#  gc(verbose = FALSE)
  kc$distribution <- "normal"
  kc$krige.var[kc$krige.var < 1e-8] <- 0
#  if(any(round(kc$krige.var, digits=12) < 0))
  if(any(kc$krige.var < 0))
    cat("krige.conv: negative kriging variance found! Investigate why this is happening.\n")
  ##
  ## ########### Sampling from the resulting distribution ###
  ##
  if(n.predictive > 0) {
    if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)){
      warning(".Random.seed not initialised. Creating it with runif(1)")
      runif(1)
    }
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    if(messages.screen)
      cat("krige.conv: sampling from the predictive distribution (conditional simulations)\n")
    if(length(cov.pars) > 2){
      if(beta.prior == "flat"){
        ok.add.var <- drop(bitb)
        tv0ivv0 <- tv0ivv0 + ok.add.var
      }
      varcov <- (varcov.spatial(coords = locations,
                                cov.model = cov.model,
                                cov.pars = cov.pars,
                                kappa = kappa, nugget = nugget)$varcov) -
                                  tv0ivv0
      remove("tv0ivv0")
#      gc(verbose=FALSE)
## can this bit be improved in efficiency/robustness ?
      kc$simulations <-  kc$predict +
        crossprod(chol(varcov), matrix(rnorm(ni * n.predictive),
                                       ncol=n.predictive))
    }
    else{
      coincide.cond <- ((isTRUE(all.equal(nugget, 0)) | !signal) & (!is.null(loc.coincide)))
      if(coincide.cond){
        nloc <- ni - length(loc.coincide)
        ind.not.coincide <- -(loc.coincide) 
        v0 <- v0[,ind.not.coincide, drop=FALSE]
        b <- b[,ind.not.coincide, drop=FALSE]
      }
      else{
        nloc <- ni
        ind.not.coincide <- TRUE
      }
      Dval <- 1.0 + nug.factor
      if(beta.prior == "deg")
        vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
      else
        vbetai <- matrix(.solve.geoR(ttivtt),
                         ncol = beta.size, nrow = beta.size)
      df.model <- ifelse(beta.prior == "deg", n, n-beta.size)
      kc$simulations <- matrix(NA, nrow=ni, ncol=n.predictive)
      if(nloc > 0){
        ## re-write this without inverse!!!
        invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                                 kappa = kappa, nugget = tausq.rel,
                                 cov.pars = cpars, inv = TRUE,
                                 only.inv.lower.diag = TRUE)
        kc$simulations[ind.not.coincide,] <- 
          .cond.sim(env.loc = base.env, env.iter = base.env,
                    loc.coincide = loc.coincide,
                    coincide.cond = coincide.cond, 
                    tmean = kc$predict[ind.not.coincide],
                    Rinv = invcov,
                    mod = list(beta.size = beta.size, nloc = nloc,
                      Nsims = n.predictive, n = n, Dval = Dval,
                      df.model = df.model, s2 = sill.partial,
                      cov.model.number = .cor.number(cov.model),
                      phi = phi, kappa = kappa, nugget=nugget),
                    vbetai = vbetai,
                    fixed.sigmasq = TRUE)
        remove("invcov")
      }
      remove("b", "locations")
#      gc(verbose = FALSE)
      if(coincide.cond)
        kc$simulations[loc.coincide,] <- rep(data.coincide, n.predictive)
    }
    ##
    ## Backtransforming simulations
    ##
    if(!isTRUE(all.equal(lambda, 1))){
      if(messages.screen)
        cat("krige.conv: back-transforming the simulated values\n")
      if(any(kc$simulations < -1/lambda))
        warning("Truncation in the back-transformation: there are simulated values less than (- 1/lambda) in the normal scale.")
      kc$simulations <-
        BCtransform(x=kc$simulations, lambda=lambda, inverse=TRUE)$data
    }
    ##
    ## mean/quantiles/probabilities estimators from simulations
    ##
    if(!is.null(mean.estimator) | !is.null(quantile.estimator) |
       !is.null(probability.estimator) | !is.null(sim.means) |
       !is.null(sim.vars)){
      kc <- c(kc, statistics.predictive(simuls = kc$simulations,
                                        mean.var = mean.estimator,
                                        quantile = quantile.estimator,
                                        threshold = probability.estimator,
                                        sim.means = sim.means,
                                        sim.vars = sim.vars))
    }
    kc$.Random.seed <- seed
  }
  ##
  ## Backtransforming moments of the prediction distribution
  ## NOTE: this must be placed here, AFTER the simulations
  ##
  if(!isTRUE(all.equal(lambda, 1))){
    if(messages.screen){
      cat("krige.conv: back-transforming the predicted mean and variance\n")
      ##      if((abs(lambda) > 0.001) & (abs(lambda-0.5) > 0.001))
      if(!isTRUE(all.equal(lambda,0)) & !isTRUE(all.equal(lambda,0.5)))
        cat("krige.conv: back-transforming by simulating from the predictive.\n           (run the function a few times and check stability of the results.\n")
    }
    kc[c("predict", "krige.var")] <-
      backtransform.moments(lambda = lambda,
                            mean = kc$predict,
                            variance = kc$krige.var,
                            distribution = "normal",
                            n.simul = n.back.moments)[c("mean", "variance")]
  }
  ##
  message <- "krige.conv: Kriging performed using global neighbourhood"
  if(messages.screen) cat(paste(message, "\n"))
  ##
  kc$message <-  message
  kc$call <- call.fc
  ##
  ## Setting classes and attributes 
  ##
  attr(kc, "sp.dim") <- ifelse(krige1d, "1d", "2d")
  attr(kc, "prediction.locations") <- call.fc$locations
  attr(kc, "parent.env") <- parent.frame()
  if(!is.null(call.fc$coords))
    attr(kc, "data.locations") <- call.fc$coords
  else attr(kc, "data.locations") <- substitute(a$coords, list(a=substitute(geodata)))
  if(!is.null(call.fc$borders))
    attr(kc, "borders") <- call.fc$borders
  oldClass(kc) <- "kriging"
  return(kc)
}

"krige.control" <-
  function(type.krige = "ok",
           trend.d = "cte", trend.l = "cte",
#           trend.d, trend.l,
           obj.model = NULL,
           beta, cov.model, cov.pars, kappa,
           nugget, micro.scale = 0, dist.epsilon = 1e-10, 
           aniso.pars, lambda)
{
  if(type.krige != "ok" & type.krige != "OK" & type.krige != "o.k." & type.krige != "O.K." & type.krige != "sk" & type.krige != "SK" & type.krige != "s.k." & type.krige != "S.K.")
    stop("krige.conv: wrong option in the argument type.krige. It should be \"sk\" or \"ok\"(if ordinary or simple kriging is to be performed)")
  if(type.krige=="OK" | type.krige=="O.K." |type.krige=="o.k.")
    type.krige <- "ok"
  if(type.krige=="SK" | type.krige=="S.K." |type.krige=="s.k.")
    type.krige <- "sk"
  ##
  if(!is.null(obj.model)){
    if(any(class(obj.model) == "eyefit")){
      if(is.null(obj.model)) stop("eyefit object is NULL, you may have forgotten to save the fitted variogram model")
      if(length(obj.model) == 1) obj.model <- obj.model[[1]]
      else stop("select the eyefit model to be used (with [[ ]])")
    }
    if(missing(beta)) beta <- obj.model$beta
    if(missing(cov.model)) cov.model <- obj.model$cov.model
    if(missing(cov.pars)) cov.pars <- obj.model$cov.pars
    if(missing(kappa)) kappa <- obj.model$kappa
    if(missing(nugget)) nugget <- obj.model$nugget
    if(missing(lambda)) lambda <- obj.model$lambda
    if(missing(aniso.pars)) aniso.pars <- obj.model$aniso.pars
## added
    if(!is.null(obj.model$trend) && is.character(obj.model$trend))
    trend.d <- trend.l <- obj.model$trend
##    if(missing(trend.d)) trend.d <- obj.model$trend
  }
  else{
    if(missing(beta)) beta <- NULL
    if(missing(cov.model)) cov.model <- "matern"
    if(missing(cov.pars))
      stop("covariance parameters (sigmasq and phi) should be provided")
    if(missing(kappa)) kappa <- 0.5
    if(missing(nugget)) nugget <- 0
    if(missing(lambda)) lambda <- 1
    if(missing(aniso.pars)) aniso.pars <- NULL
  }
  ##
  if(type.krige == "sk")
    if(is.null(beta) | mode(beta) != "numeric")
      stop("\nkrige.conv: argument beta must be provided in order to perform simple kriging")
  cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
  if(micro.scale > nugget)
    stop("krige.control: micro.scale must be in the interval [0, nugget]")
  ##
  if(!is.null(aniso.pars))
    if(length(aniso.pars) != 2 | mode(aniso.pars) != "numeric")
      stop("krige.control: anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
  ##
  if((is.character(trend.d) | is.character(trend.l)) && (trend.d != trend.l))
        stop("krige.control: trend.l specifications is different from trend.d")
  if(inherits(trend.d, "formula") | inherits(trend.l, "formula")){
    if((inherits(trend.d, "formula") == FALSE) | (inherits(trend.l, "formula") == FALSE))
      stop("krige.control: trend.d and trend.l must have similar specification")
  }
  else{
##    if((!is.null(class(trend.d)) && class(trend.d) == "trend.spatial") &
##       (!is.null(class(trend.l)) && class(trend.l) == "trend.spatial")){
    if((length(class(trend.d)) > 0 && class(trend.d) == "trend.spatial") &
       (length(class(trend.l)) > 0 && class(trend.l) == "trend.spatial")){
      if(ncol(trend.d) != ncol(trend.l))
        stop("krige.bayes: trend.d and trend.l do not have the same number of columns")
    }
    else
      if(trend.d != trend.l)
        stop("krige.control: trend.l is different from trend.d")
  }
  ##
  res <- list(type.krige = type.krige,
              trend = trend.d, trend.d = trend.d, trend.l = trend.l, 
              beta = beta,
              cov.model = cov.model, 
              cov.pars = cov.pars, kappa = kappa,
              nugget = nugget,
              micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
              aniso.pars = aniso.pars, lambda = lambda)
  oldClass(res) <- "krige.geoR"
  return(res)
}

"image.kriging" <-
  function (x, locations, borders, values = x$predict,
            coords.data, x.leg, y.leg, ...) 
{
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations))
    locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(missing(borders)){
    if(!is.null(attr(x, "borders")))
      borders.arg <- borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    else
      borders.arg <- borders <- eval(x$call$geodata, envir= attr(x, "parent.env"))$borders
#    borders.arg <- borders <- NULL
  }
  else{
    borders.arg <- borders
    if(is.null(borders)) borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
  }
  if(missing(coords.data)) coords.data <- NULL
  else
    if(all(coords.data == TRUE))
      coords.data <-  eval(attr(x, "data.locations"), envir= attr(x, "parent.env"))
  if(missing(x.leg)) x.leg <- NULL
  if(missing(y.leg)) y.leg <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D'){
    do.call("plot.1d", c(list(x = values,
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d",
                                   data="prediction")))
  }
  else{
    ldots.image <- .ldots.set(ldots, type="image", data="prediction")
    locations <- .prepare.graph.kriging(locations=locations,
                                        borders=borders,
                                        borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                        values=values,
                                        xlim= ldots.image$xlim,
                                        ylim= ldots.image$yli)
    do.call("image", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.image))
    ##
    ## adding points at data locations
    ##
    if(!is.null(coords.data)) points(coords.data, pch=20)
    ##
    ## adding borders
    ##
    if(!is.null(borders.arg)) polygon(borders, lwd=2)
    ##
    ## adding the legend
    ##
    if(!is.null(x.leg) & !is.null(y.leg)){
###      if(is.null(ldots$col)) ldots$col <- heat.colors(12)
      do.call("legend.krige", c(list(x.leg=x.leg,
                                     y.leg=y.leg,
                                     values = locations$values[!is.na(locations$values)]),
                                     ldots))
    }
  }
  return(invisible())
}

"contour.kriging" <-
  function (x, locations, borders, values = x$predict, coords.data,
            filled=FALSE, ...) 
{
  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(missing(borders)){
    if(!is.null(attr(x, "borders")))
      borders.arg <- borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    else
      borders.arg <- borders <- eval(x$call$geodata, envir= attr(x, "parent.env"))$borders
      # borders.arg <- borders <- NULL
  }
  else{
    borders.arg <- borders
    if(is.null(borders)) borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
  }
  if(missing(coords.data)) coords.data <- NULL
  else
    if(all(coords.data == TRUE))
      coords.data <-  eval(attr(x, "data.locations"), envir= attr(x, "parent.env"))
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D'){
    do.call("plot.1d", c(list(x = values,
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d",
                                   data="prediction")))
  }
  else{
    if(filled)
      ldots.contour <- .ldots.set(ldots, type="filled.contour",
                                 data="prediction")
    else
      ldots.contour <- .ldots.set(ldots, type="contour",
                                 data="prediction")
    if(is.null(ldots.contour$asp)) ldots.contour$asp=1
    locations <- .prepare.graph.kriging(locations=locations,
                                       borders=borders,
                                       borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                       values=values,
                                       xlim= ldots.contour$xlim,
                                       ylim= ldots.contour$ylim)
    if(filled){
      if(is.null(ldots.contour$plot.axes)){
        ldots.contour$plot.axes  <- quote({
          axis(1)
          axis(2)
          if(!is.null(coords.data)) points(coords.data, pch=20)
          if(!is.null(borders)) polygon(borders, lwd=2)
        })
      }
      do.call("filled.contour", c(list(x=locations$x, y=locations$y,
                                       z=locations$values),
                                  ldots.contour))
    }
    else{
      do.call("contour", c(list(x=locations$x, y=locations$y,
                                z=locations$values), ldots.contour))
      ##
      ## adding borders
      ##
      if(!is.null(borders.arg)) polygon(borders, lwd=2)
      if(!is.null(coords.data)) points(coords.data, pch=20)
    }
  }
  return(invisible())
}

"persp.kriging" <-
  function(x, locations, borders, values = x$predict, ...)
{
  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(missing(borders)) borders <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D')
    do.call("plot.1d", c(list(x= values,
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d",
                                   data="prediction")))
  else{
    ldots.persp <- .ldots.set(ldots, type="persp", data="prediction")
    locations <- .prepare.graph.kriging(locations=locations,
                                       borders=borders,
                                       borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                       values=values,
                                       xlim= ldots.persp$xlim,
                                       ylim= ldots.persp$ylim)
    do.call("persp", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.persp))
  }
  return(invisible())
}

