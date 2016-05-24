"krige.bayes" <- 
  function(geodata, coords=geodata$coords, data=geodata$data,
           locations = "no", borders, model, prior, output)
{
  ##
  ## ======================= PART 1 ==============================
  ##                Reading and Checking Input
  ## =============================================================
  ##
  ## setting output object and environments
  ##
  if(missing(geodata))
    geodata <- list(coords = coords, data = data)
  if(missing(borders))
    borders <- geodata$borders
  call.fc <- match.call()
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)){
    warning(".Random.seed not initialised. Creating it with runif(1)")
    runif(1)
  }
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  locations <- .check.locations(locations)
  do.prediction <- ifelse(all(locations == "no"), FALSE, TRUE)
  if(do.prediction && any(!is.numeric(locations)))
		stop("krige.bayes: non numeric coordinates passed to the argument \"locations\"")
  base.env <- sys.frame(sys.nframe())
  message.prediction <- character()
  phidist<- list()
  "krige.bayes.counter" <- 
    function(.temp.ap, n.points){  
      if(n.points <= 50) cat(paste(.temp.ap, ", ", sep=""))
      if(n.points > 50 & n.points <= 500)
        if(.temp.ap %% 10 == 1) cat(paste(.temp.ap, ", ", sep=""))
      if(n.points > 500)
        if(.temp.ap %% 100 == 1) cat(paste(.temp.ap, ", ", sep=""))
      if(n.points == .temp.ap) cat("\n")
    }
  kb <- list(posterior = list(beta=list(), sigmasq=list(),
               phi=list(), tausq.rel=list()),
             predictive=list(mean = NULL, variance = NULL, distribution = NULL))
  oldClass(kb$posterior) <- c("posterior.krige.bayes", "variomodel")
  oldClass(kb$predictive) <- "predictive.krige.bayes"
  pred.env <- new.env()
  ##
  ## reading model input
  ##
  if(missing(model))
    model <- model.control()
  else{
##    if(is.null(class(model)) || class(model) != "model.geoR"){
    if(length(class(model)) == 0 || class(model) != "model.geoR"){
      if(!is.list(model))
        stop("krige.bayes: the argument model only takes a list or an output of the function model.control")
      else{
        model.names <- c("trend.d", "trend.l", "cov.model", "kappa", "aniso.pars", "lambda") 
        model.user <- model
        model <- list()
        if(length(model.user) > 0){
          for(i in 1:length(model.user)){
            n.match <- match.arg(names(model.user)[i], model.names)
            model[[n.match]] <- model.user[[i]]
          }
        }    
        if(is.null(model$trend.d)) model$trend.d <- "cte"  
        if(is.null(model$trend.l)) model$trend.l <- "cte"  
        if(is.null(model$cov.model)) model$cov.model <- "matern"  
        if(is.null(model$kappa)) model$kappa <- 0.5
        if(is.null(model$aniso.pars)) model$aniso.pars <- NULL 
        if(is.null(model$lambda)) model$lambda <- 1
        model <- model.control(trend.d = model$trend.d,
                               trend.l = model$trend.l,
                               cov.model = model$cov.model,
                               kappa = model$kappa,
                               aniso.pars = model$aniso.pars,
                               lambda = model$lambda)
      }
    }
  }
  cov.model <- model$cov.model
  cov.model.number <- .cor.number(cov.model)
  kappa <- model$kappa
  lambda <- model$lambda
  ##
  ## reading prior specification
  ##
  if(missing(prior))
    prior <- prior.control()
  else{
##    if(is.null(class(prior)) || class(prior) != "prior.geoR"){
    if(length(class(prior)) == 0 || class(prior) != "prior.geoR"){
      if(!is.list(prior))
        stop("krige.bayes: the argument prior only takes a list or an output of the function prior.control")
      else{
        prior.names <- c("beta.prior", "beta", "beta.var.std", "sigmasq.prior",
                         "sigmasq", "df.sigmasq", "phi.prior", "phi", "phi.discrete",
                         "tausq.rel.prior", "tausq.rel", "tausq.rel.discrete") 
        prior.user <- prior
        prior <- list()
        if(length(prior.user) > 0){
          for(i in 1:length(prior.user)){
            n.match <- match.arg(names(prior.user)[i], prior.names)
            prior[[n.match]] <- prior.user[[i]]
          }
        }
        ## DO NOT CHANGE ORDER OF THE NEXT 3 LINES
        if(is.null(prior$beta)) prior$beta <-  NULL
        if(is.null(prior$beta.prior))
          prior$beta.prior <-  c("flat", "normal", "fixed")
        if(is.null(prior$beta.var.std)) prior$beta.var.std <-  NULL
        ## DO NOT CHANGE ORDER OF THE NEXT 3 LINES
        if(is.null(prior$sigmasq)) prior$sigmasq <- NULL
        if(is.null(prior$sigmasq.prior))
          prior$sigmasq.prior <- c("reciprocal",  "uniform", "sc.inv.chisq",  "fixed") 
        if(is.null(prior$df.sigmasq)) prior$df.sigmasq <- NULL
        ## DO NOT CHANGE ORDER OF THE NEXT 3 LINES
        if(is.null(prior$phi)) prior$phi <- NULL
        if(is.null(prior$phi.prior))
          prior$phi.prior <- c("uniform", "exponential", "fixed", "squared.reciprocal","reciprocal")
        if(is.null(prior$phi.discrete)) prior$phi.discrete <- NULL
        ## DO NOT CHANGE ORDER OF THE NEXT 3 LINES
        if(is.null(prior$tausq.rel)) prior$tausq.rel <- 0
        if(is.null(prior$tausq.rel.prior))
          prior$tausq.rel.prior <- c("fixed", "uniform", "reciprocal")
        if(is.null(prior$tausq.rel.discrete))
          prior$tausq.rel.discrete <- NULL 
        prior <- prior.control(beta.prior = prior$beta.prior,
                               beta = prior$beta,
                               beta.var.std = prior$beta.var.std,
                               sigmasq.prior = prior$sigmasq.prior,
                               sigmasq = prior$sigmasq,
                               df.sigmasq = prior$df.sigmasq,
                               phi.prior = prior$phi.prior,
                               phi = prior$phi,
                               phi.discrete = prior$phi.discrete, 
                               tausq.rel.prior = prior$tausq.rel.prior,
                               tausq.rel = prior$tausq.rel,
                               tausq.rel.discrete = prior$tausq.rel.discrete)
      }
    }
  }
  kb$prior <- prior$priors.info
  kb$model <- model
  oldClass(kb$prior) <- "prior.geoR"
  ##
  if(prior$dep.prior){
    npr <- length(prior$sigmasq)
    nphipr <- nrow(as.matrix(prior$sigmasq))
    ntaupr <- ncol(as.matrix(prior$sigmasq))
  }
  else nphipr <- ntaupr <- npr <- 1
  beta <- prior$beta
  if(prior$beta.prior == "fixed") beta.fixed <- beta
  if(prior$beta.prior == "normal"){
    nbeta <- attr(prior$beta.var.std, "Size")
    betares <- list()
    for(j in 1:ntaupr){
      for(i in 1:nphipr){
        beta <- array(prior$beta, dim=c(nphipr, ntaupr, nbeta))[i,j,]
        beta.var.std <- array(prior$beta.var.std,
                              dim=c(nphipr, ntaupr, nbeta^2))[i,j,]
        beta.var.std <- matrix(beta.var.std, nbeta, nbeta)
        ind.pos <- (j-1)*nphipr + i
        betares[[ind.pos]] <-
          list(iv = .solve.geoR(beta.var.std),
               ivm = drop(.solve.geoR(beta.var.std, beta)),
               mivm = drop(crossprod(beta, .solve.geoR(beta.var.std, beta))))
      }
    }
  }
  if(prior$sigmasq.prior != "fixed") S2.prior <- prior$sigmasq
  else sigmasq.fixed <- S2.prior <- prior$sigmasq
  df.sigmasq.prior <- prior$df.sigmasq
  ##
  phi.discrete <- prior$phi.discrete
  exponential.par <- prior$phi  
  ##
  tausq.rel.fixed <- tausq.rel <- prior$tausq.rel
  exponential.tausq.rel.par <- prior$tausq.rel  
  if(tausq.rel.fixed > 2)
    print("WARNING: relative (NOT absolute) nugget should be specified.")
  tausq.rel.discrete <- prior$tausq.rel.discrete
  ##
  ## checking data configuration
  ##
  n <- length(data)
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("krige.bayes: vector of coordinates: assuming one spatial dimension (transect)")
  }
  coords <- as.matrix(coords)
  dists.env <- new.env()
  assign("data.dist", as.vector(dist(coords)), envir=dists.env)
  data.dist.range <- range(get("data.dist", envir=dists.env))
  data.dist.min <- data.dist.range[1]
  data.dist.max <- data.dist.range[2]
  if(round(1e12*data.dist.min) == 0)
    stop("krige.bayes: this function does not allow two data at same location")
  ##
  ## reading output options
  ##
  if(missing(output)) output <- output.control()
  else{
    if(length(class(output)) == 0 || class(output) != "output.geoR"){
      if(!is.list(output))
        stop("krige.bayes: the argument output only takes a list or an output of the function output.control")
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
  n.posterior <- output$n.posterior
  messages.screen <- output$messages.screen
  if(!do.prediction) {
    if(prior$beta.prior != "fixed" & prior$sigmasq.prior != "fixed"  &
       prior$phi.prior != "fixed" & output$messages.screen){
      cat("krige.bayes: no prediction locations provided.\n")
      cat("             Only samples of the posterior for the parameters will be returned.\n")
    }
  }
  else{
    n.predictive <- output$n.predictive
    if(is.null(n.predictive))
      n.predictive <- ifelse(prior$phi.prior == "fixed", 0, n.posterior)
    simulations.predictive <- output$simulations.predictive
    if(is.null(simulations.predictive))
      simulations.predictive <- ifelse(prior$phi.prior == "fixed", FALSE, TRUE)
    keep.simulations <- output$keep.simulations
    if(is.null(keep.simulations))
      keep.simulations <- simulations.predictive
    mean.estimator <- output$mean.estimator
    if(is.null(mean.estimator))
      mean.estimator <- ifelse(simulations.predictive, TRUE, FALSE)
    moments <- output$moments    
    if(is.null(moments) | prior$phi.prior == "fixed")
      moments <- TRUE
    n.back.moments <- output$n.back.moments
    sim.means <- output$sim.means
    if(is.null(sim.means))
      sim.means <- ifelse(simulations.predictive, TRUE, FALSE)
    sim.vars <- output$sim.vars
    if(is.null(sim.vars)) sim.vars <- FALSE
    signal <- ifelse(is.null(output$signal), TRUE, output$signal)
    quantile.estimator <- output$quantile.estimator
    probability.estimator <- output$probability.estimator
    if(simulations.predictive & n.predictive == 0) n.predictive <- 1000
  }
  ##
  ## Box-Cox transformation
  ##
  if(abs(lambda-1) > 0.001) {
    if(messages.screen)
      cat(paste("krige.bayes: Box-Cox's transformation performed for lambda =", round(lambda,digits=3), "\n"))
    data <- BCtransform(x=data, lambda = lambda)$data
  }
  ##
  ## Building trend (covariates/design) matrices:   
  ##
  dimnames(coords) <- list(NULL, NULL)
  if(nrow(coords) != length(data))
    stop("krige.bayes: number of data is different of number of data locations (coordinates)")
  if(class(model$trend.d) == "trend.spatial")
    trend.data <- unclass(model$trend.d)
  else
    trend.data <- unclass(trend.spatial(trend=model$trend.d, geodata = geodata))
  if(nrow(trend.data) != nrow(coords))
    stop("trend specification not compatible with the length of the data") 
  beta.size <- ncol(trend.data)
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  if(prior$beta.prior == "normal" |  prior$beta.prior == "fixed"){
    if(beta.size != length(beta))
      stop("size of beta incompatible with the trend model (covariates)")
  }
  if(do.prediction) {
    locations <- .check.locations(locations)    
    ##
    ## selecting locations inside the borders 
    ##
    if(!is.null(borders)){
      nloc0 <- nrow(locations)
      ind.loc0  <- .geoR_inout(locations, borders)
      ##    locations <- locations.inside(locations, borders)
      locations <- locations[ind.loc0,,drop=TRUE]
      if(nrow(locations) == 0){
        warning("\nkrige.bayes: no prediction to be performed.\n There are no prediction locations inside the borders")
        do.prediction <- FALSE
      }
      if(messages.screen)
        cat("krige.bayes: results will be returned only for prediction locations inside the borders\n")
    }
    ##
    ## Checking for 1D prediction 
    ##
    krige1d <- ifelse((length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1), TRUE, FALSE)
    ##
    ## Checking trend specification
    ##
    if(inherits(model$trend.d, "formula") | inherits(model$trend.l, "formula")){
      if((inherits(model$trend.d, "formula") == FALSE) | (inherits(model$trend.l, "formula") == FALSE))
        stop("krige.bayes: model$trend.d and model$trend.l must have similar specification\n")
    }
    else{
##      if((!is.null(class(model$trend.d)) && class(model$trend.d) == "trend.spatial") &
##         (!is.null(class(model$trend.l)) && class(model$trend.l) == "trend.spatial")){
      if((length(class(model$trend.d)) > 0 && class(model$trend.d) == "trend.spatial") &
         (length(class(model$trend.l)) > 0 && class(model$trend.l) == "trend.spatial")){
        if(ncol(model$trend.d) != ncol(model$trend.l))
          stop("krige.bayes: trend.d and trend.l do not have the same number of columns")
      }
      else
        if(model$trend.d != model$trend.l)
          stop("krige.bayes: especification of model$trend.l and model$trend.d must be compatible")
    }
    ##
    if(messages.screen){
      cat(switch(as.character(model$trend.d)[1],
                 "cte" = "krige.bayes: model with constant mean",
                 "1st" = "krige.bayes: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = "krige.bayes: model with mean given by a 2nd order polynomial on the coordinates",
                 "krige.bayes: model with mean defined by covariates provided by the user"))
      cat("\n")
    }
    ##
    dimnames(locations) <- list(NULL, NULL)
    if(class(model$trend.l) == "trend.spatial")
      assign("trend.loc", unclass(model$trend.l), envir=pred.env)
    else
      assign("trend.loc", unclass(trend.spatial(trend=model$trend.l, geodata = list(coords = locations))), envir=pred.env)
    ni <- nrow(get("trend.loc", envir=pred.env))
    if(!is.null(borders))
      if(ni == nloc0)
        assign("trend.loc",
               get("trend.loc", envir=pred.env)[ind.loc0,,drop=FALSE],
               envir=pred.env)
     if(nrow(locations) != ni)
      stop("trend.l is not compatible with number of prediction locations")
    expect.env <- new.env()
    assign("expect", 0, envir=expect.env)
    assign("expect2", 0, envir=expect.env)
  }
  ##
  ## Anisotropy correction
  ##   (warning: this must be placed here, AFTER trend matrices be defined)
  ##
  if(!is.null(model$aniso.pars)) {
#    if((abs(model$aniso.pars[1]) > 0.001) & (abs(model$aniso.pars[2] - 1) > 0.001)){
    if(abs(model$aniso.pars[2] - 1) > 0.001){
      if(messages.screen)
        cat("krige.bayes: anisotropy parameters provided and assumed to be constants\n")
      coords <- coords.aniso(coords = coords, aniso.pars = model$aniso.pars)
      if(do.prediction)
        locations <- coords.aniso(coords = locations, 
                                  aniso.pars = model$aniso.pars)
      remove(dists.env)
      dists.env <- new.env()
      assign("data.dist", as.vector(dist(coords)), envir=dists.env)
    }
  }  
  ##
  ## Distances between data and prediction locations
  ## Must be here, AFTER anisotropy be checked
  ##
  if(do.prediction){
    assign("d0", loccoords(coords = coords, locations = locations), envir=pred.env)
    ##
    ## checking coincident data and prediction locations
    ##
    loc.coincide <- apply(get("d0", envir=pred.env), 2, function(x){any(x < 1e-12)})
    if(any(loc.coincide))
      loc.coincide <- (1:ni)[loc.coincide]
    else
      loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data){return(data[x < 1e-10])}
      data.coincide <- apply(get("d0", envir=pred.env)[,loc.coincide, drop=FALSE],
                             2,temp.f, data=data)
    }
    else
      data.coincide <- NULL
    n.loc.coincide <- length(loc.coincide)
    assign("loc.coincide", loc.coincide, envir=pred.env)
    assign("data.coincide", data.coincide, envir=pred.env)
    remove(data.coincide, loc.coincide)
#    gc(verbose=FALSE)
  }
  ##
  ## Preparing prior information on beta and sigmasq
  ##
  beta.info <- list()
  sigmasq.info <- list()
  for(i in 1:npr){
    beta.info[[i]] <-
      switch(prior$beta.prior,
             fixed = list(mivm = 0, ivm = 0, iv = Inf, beta.fixed = beta.fixed, p = 0),
             flat = list(mivm = 0, ivm = 0, iv = 0, p = beta.size),
             normal = list(mivm = betares[[i]]$mivm, ivm = betares[[i]]$ivm,
               iv = betares[[i]]$iv, p = 0))
    sigmasq.info[[i]] <-
      switch(prior$sigmasq.prior,
             fixed = list(df.sigmasq = Inf, n0S0 = 0,
               sigmasq.fixed = sigmasq.fixed),
             reciprocal = list(df.sigmasq = 0, n0S0 = 0),
             uniform = list(df.sigmasq = -2, n0S0 = 0),
             sc.inv.chisq = list(df.sigmasq = df.sigmasq.prior, n0S0 = df.sigmasq.prior*S2.prior[i]))
  }
  beta.info$p <- switch(prior$beta.prior,
                        fixed = 0,
                        flat = beta.size,
                        normal = 0)
  sigmasq.info$df.sigmasq <- switch(prior$sigmasq.prior,
                                    fixed = Inf,
                                    reciprocal = 0,
                                    uniform = -2, 
                                    sc.inv.chisq = df.sigmasq.prior)
  ##
  ## ====================== PART 2 =============================
  ##                 FIXED PHI AND TAUSQ.REL
  ## ===========================================================
  ## 
  if(prior$phi.prior == "fixed"){
    phi.fixed <- prior$phi
    ##
    ## Computing parameters of the posterior for $(\beta, \sigma^2)$ 
    ## and moments of the predictive (if applies)
    ##
    bsp <- beta.sigmasq.post(n = n, beta.info = beta.info[[1]],
                             sigmasq.info = sigmasq.info[[1]], 
                             env.dists = dists.env,
                             model = list(cov.model = model$cov.model, kappa = model$kappa),
                             xmat = trend.data, y = data,
                             phi = phi.fixed, tausq.rel = tausq.rel.fixed,
                             do.prediction.moments = (do.prediction && moments),
                             do.prediction.simulations = (do.prediction && simulations.predictive),
                             env.pred = pred.env,
                             signal = (do.prediction && signal))
    ##
    ## Preparing output of the posterior distribution
    ##    
    if(prior$beta.prior == "fixed")
      kb$posterior$beta <- list(status = "fixed", fixed.value = beta.fixed)
    else{
      if(prior$sigmasq.prior == "fixed")
        kb$posterior$beta <- list(distribution = "normal")
      else
        kb$posterior$beta <- list(distribution = "t",
                                  conditional = "normal")
      kb$posterior$beta$pars <- list(mean = bsp$beta.post,
                                     var = bsp$S2.post * bsp$beta.var.std.post)
      attr(kb$posterior$beta$pars$var, "Size") <- beta.size
      class(kb$posterior$beta$pars$var) <- "betavar"
    }        
    if(prior$sigmasq.prior == "fixed")
      kb$posterior$sigmasq <- list(status="fixed", fixed.value=sigmasq.fixed)
    else
      kb$posterior$sigmasq <- list(distribution = "sc.inv.chisq",
                                   pars = list(df = bsp$df.post,
                                     S2 = bsp$S2.post))
    kb$posterior$phi<- list(status= "fixed", fixed.value = phi.fixed)
    kb$posterior$tausq.rel <-
      list(status= "fixed", fixed.value = tausq.rel.fixed)
    ##
    ## Preparing output of the predictive distribution
    ##
    kb$predictive$mean <- bsp$pred.mean
    kb$predictive$variance <- bsp$pred.var
    kb$predictive$distribution <- ifelse(prior$sigmasq.prior == "fixed",
                                         "normal", "t")
    bsp[c("pred.mean", "pred.var")] <- NULL
    ##
    ## preparing objects for simulating from the predictive
    ##
    if(do.prediction && simulations.predictive && n.predictive > 0){
      phidist$s2 <-  as.matrix(bsp$S2.post) 
      phidist$probphitausq <-  as.matrix(1)
      phidist$beta <- array(bsp$beta.post, c(1,1,beta.size)) 
      phidist$varbeta  <- array(bsp$beta.var.std.post, c(1,1,beta.size^2))
      phi.unique <- phidist$phitausq <- t(c(phi.fixed, tausq.rel.fixed))
      df.model <- bsp$df.post 
      ind.length <- 1
      inv.lower <- array(bsp$inv.lower, dim=c(1,1,(n*(n-1)/2)))
      inv.diag <- array(bsp$inv.diag, dim=c(1,1,n))
      ind.table <- n.predictive
      phi.discrete <- phi.fixed
      tausq.rel.discrete <- tausq.rel.fixed
    }
  }
  else{
    ##
    ## ====================== PART 3 =============================
    ##                 RANDOM PHI AND TAUSQ.REL
    ## ===========================================================
    ##
    if(messages.screen)
      cat("krige.bayes: computing the discrete posterior of phi/tausq.rel\n")
    ##
    ## Preparing discrete set for phi and/or tausq.rel
    ##
    if(is.null(phi.discrete)){
      phi.discrete <- seq(0, 2 * data.dist.max, l=51)[-1]
      if(messages.screen)
        cat("krige.bayes: argument `phi.discrete` not provided, using default values\n")
    } 
    if(mode(phi.discrete) != "numeric")
      stop("non-numerical value provided in phi.discrete")
    if(length(phi.discrete) == 1)
      stop("only one value provided in phi.discrete. Use prior.phi=`fixed`")
    n.phi.discrete <- length(phi.discrete)
    n.tausq.rel.discrete <- length(tausq.rel.discrete)
    phi.names <- paste("phi", phi.discrete, sep="")
    tausq.rel.names <- paste("tausqrel", tausq.rel.discrete, sep="")
    phidist$phitausq <- as.matrix(expand.grid(phi.discrete, tausq.rel.discrete))
    if(prior$phi.prior == "user" | prior$tausq.rel.prior == "user"){
      if(prior$tausq.rel.prior == "fixed")
        phidist$phitausq <-
          cbind(phidist$phitausq, prior$priors.info$phi$probs, 1)
      else{
        if(is.null(prior$joint.phi.tausq))
          phidist$phitausq <-
            cbind(phidist$phitausq,
                  as.matrix(expand.grid(prior$priors.info$phi$probs,
                                        prior$priors.info$tausq.rel$probs)))
        else
          phidist$phitausq <-
            cbind(phidist$phitausq, as.vector(prior$joint.phi.tausq.rel), 1)
      }
    }
    dimnames(phidist$phitausq) <- list(NULL, NULL)
    ##
    ##  Degrees of freedom for the posteriors
    ##
    df.model <- ifelse(sigmasq.info$df.sigmasq == Inf, Inf,
                       (n + sigmasq.info$df.sigmasq - beta.info$p))
    ##
    ## Function to compute the posterior probabilities
    ## for each parameter sets (phi, tausq.rel)
    ## 
    phi.tausq.rel.post <- function(phinug){
      par.set <- get("parset", envir=counter.env)
      if(messages.screen)
        krige.bayes.counter(.temp.ap = par.set, n.points = ntocount)
      on.exit(assign("parset", get("parset", envir=counter.env)+1, envir=counter.env))
      phi <- phinug[1]
      tausq.rel <- phinug[2]
      if(prior$beta.prior == "normal" && npr > 1) info.id <- par.set
      else info.id <- 1
      bsp <- beta.sigmasq.post(n = n, beta.info = beta.info[[info.id]],
                               sigmasq.info = sigmasq.info[[info.id]],
                               env.dists = dists.env,
                               model = list(cov.model = model$cov.model,
                                 kappa = model$kappa),
                               xmat = trend.data, y = data, phi = phi,
                               tausq.rel = tausq.rel, dets = TRUE,
                               do.prediction.moments = (do.prediction && moments),
                               do.prediction.simulations = (do.prediction && simulations.predictive),
                               env.pred = pred.env, signal = signal)
      logprobphitausq <-  (-0.5) * log(bsp$det.XiRX) -
        (bsp$log.det.to.half) - (bsp$df.post/2) * log(bsp$S2.post)
      ##print("termos")
      ##print(c(log(bsp$det.XiRX),bsp$log.det.to.half, bsp$df.post/2,
      ##        log(bsp$S2.post),logprobphitausq))
      ##
      if(prior$phi.prior == "user"){
        if(phinug[3] > 0)
          logprobphitausq <- logprobphitausq + log(phinug[3])
        else logprobphitausq <- -Inf
      }
      if(prior$phi.prior == "reciprocal"){
        if(phi > 0) logprobphitausq <- logprobphitausq - log(phi)
        else logprobphitausq <- -Inf
      }
      if(prior$phi.prior == "squared.reciprocal"){
        if(phi > 0) logprobphitausq <- logprobphitausq - 2*log(phi)
        else logprobphitausq <- -Inf
      }
      if(prior$phi.prior == "exponential"){
        logprobphitausq <- logprobphitausq - log(exponential.par) - (phi/exponential.par)
       }
      if(prior$tausq.rel.prior == "user"){
        if(phinug[4] > 0)
          logprobphitausq <- logprobphitausq + log(phinug[4])
        else logprobphitausq <- -Inf
      }
      if(prior$tausq.rel.prior == "reciprocal"){
        if(tausq.rel > 0)
          logprobphitausq <- logprobphitausq - log(tausq.rel)
        else logprobphitausq <- -Inf
      }
#      print(c(par.set,logprobphitausq))
      ##
      ## The following is a trick to go aroud a numerical problem:
      ## the value of logprobphitausq can be too small (highly negative)
      ## such that the exponential of it can be numerically
      ## equal to zero
      if(get("add.cte", envir=counter.env) && is.finite(logprobphitausq)){
        assign("cte", logprobphitausq, envir=counter.env)
        assign("add.cte", FALSE, envir=counter.env)
      }
      ##      if(get("cte", envir=counter.env) > 0)
      ##       logprobphitausq <- logprobphitausq + get("cte", envir=counter.env)
      ##    else
      logprobphitausq <- logprobphitausq - get("cte", envir=counter.env)
      ##print(logprobphitausq)
      bsp$probphitausq <- drop(exp(logprobphitausq))
      ##print(logprobphitausq)
      ##print(bsp$probphitausq)
#      print(format(c(par.set,phi,tausq.rel,logprobphitausq,bsp$probphitausq)))
      ##
      if(do.prediction && moments){
        assign("expect", (get("expect", envir=expect.env) +
                          (bsp$pred.mean * bsp$probphitausq)),
               envir= expect.env)
        assign("expect2", (get("expect2", envir=expect.env) +
                           ((bsp$pred.var + (bsp$pred.mean^2)) *
                             bsp$probphitausq)),
               envir= expect.env)
      }
      phi.ind <- which.min(abs(phi.discrete - phi))
      nug.ind <- which.min(abs(tausq.rel.discrete - tausq.rel))
      assign("pn.ind", c(phi.ind, nug.ind), envir=fn.frame)
      assign("bsp", bsp, envir=fn.frame)
      eval(expression(phidist$s2[pn.ind[1], pn.ind[2]] <- bsp$S2.post), envir=fn.frame)
      eval(expression(phidist$probphitausq[pn.ind[1], pn.ind[2]] <- bsp$probphitausq), envir=fn.frame)
      eval(expression(phidist$beta[pn.ind[1], pn.ind[2],] <- drop(bsp$beta.post)), envir=fn.frame)
      eval(expression(phidist$varbeta[pn.ind[1], pn.ind[2],] <- drop(bsp$beta.var.std.post)), envir=fn.frame)
      if(do.prediction && simulations.predictive){
        eval(expression(inv.lower[pn.ind[1], pn.ind[2],] <- bsp$inv.lower),
             envir=fn.frame)
        eval(expression(inv.diag[pn.ind[1], pn.ind[2],] <- bsp$inv.diag),
             envir=fn.frame)
      }
      return(invisible())
    }
    ##
    ## Computing the posterior probabilities and organising results
    ##
    fn.frame <- sys.frame(sys.nframe())
    phidist$s2 <- matrix(NA, n.phi.discrete, n.tausq.rel.discrete)
    dimnames(phidist$s2) <- list(phi.names, tausq.rel.names)
    phidist$probphitausq <- matrix(NA, n.phi.discrete, n.tausq.rel.discrete)
    phidist$beta <- array(NA, dim=c(n.phi.discrete, n.tausq.rel.discrete, beta.size))
    dimnames(phidist$beta) <- list(phi.names, tausq.rel.names, beta.names)
    phidist$varbeta <- array(NA, dim=c(n.phi.discrete, n.tausq.rel.discrete, beta.size^2))
    dimnames(phidist$varbeta) <- list(phi.names, tausq.rel.names, NULL)
    if(do.prediction && simulations.predictive){
      inv.lower <- array(NA, dim=c(n.phi.discrete, n.tausq.rel.discrete, (n * (n - 1))/2))
      inv.diag <- array(NA, dim=c(n.phi.discrete, n.tausq.rel.discrete, n))
      frame.inv <- sys.frame(sys.nframe())
    }
    ##
    counter.env <- new.env()
    assign("parset", 1, envir=counter.env)
    assign("add.cte", TRUE, envir=counter.env)
    assign("cte", 0, envir=counter.env)
    if(messages.screen){
      ntocount <- nrow(phidist$phitausq)
      cat(paste("krige.bayes: computing the posterior probabilities.\n             Number of parameter sets: ", ntocount,"\n"))
    }
    temp.res <- apply(phidist$phitausq, 1, phi.tausq.rel.post)
    remove("bsp")
    if(messages.screen){
      remove(counter.env)
      cat("\n")
    }
    ##
    phidist$sum.prob <- sum(phidist$probphitausq)
    phidist$probphitausq <- phidist$probphitausq/phidist$sum.prob
    ##
    ## Preparing output of the posterior distribution
    ##
    kb$posterior$beta <- list(conditional.distribution = "normal",
                              pars = list(mean=phidist$beta,
                                var = phidist$varbeta))
    attr(kb$posterior$beta$pars$var, "Size") <- beta.size
#    class(kb$posterior$beta$pars$var) <- "betavar"
    kb$posterior$sigmasq <- list(conditional.distribution = "sc.inv.chisq",
                                 pars = list(df = df.model,
                                   S2 = drop(phidist$s2))) 
    kb$posterior$phi$distribution <- drop(rowSums(phidist$probphitausq))
    names(kb$posterior$phi$distribution) <- prior$phi.discrete
    if(prior$tausq.rel.prior != "fixed"){
      kb$posterior$tausq.rel$distribution <- drop(colSums(phidist$probphitausq))
      names(kb$posterior$tausq.rel$distribution) <- tausq.rel.discrete
    }
    else{
      kb$posterior$tausq.rel <-
        list(status= "fixed", fixed.value = tausq.rel.fixed)
    }
    if(prior$phi.prior != "fixed" | prior$tausq.rel.prior != "fixed"){
      kb$posterior$joint.phi.tausq.rel <- phidist$probphitausq
      dimnames(kb$posterior$joint.phi.tausq.rel) <-
        list(phi.names, tausq.rel.names)
    }
    ##
    ##  Sampling from the posterior distribution
    ##
    if(n.posterior > 0){
      if(messages.screen)
        cat("krige.bayes: sampling from posterior distribution\n")
      ##
      ## sampling phi and/or tausq
      ##
#      print(as.vector(phidist$probphitausq))
      n.points <- length(phidist$probphitausq)
      ind <- sample((1:n.points), n.posterior, replace = TRUE,
                    prob = as.vector(phidist$probphitausq))
      phi.sam <- phidist$phitausq[ind,  ]
      ##
      ## frequencies for sampled phi/tausq
      ##
      ind.unique <- sort(unique(ind))
      ind.length <- length(ind.unique)
      ind.table <- table(ind)
      names(ind.table) <- NULL
      ##
      phi.unique <- phidist$phitausq[ind.unique,, drop=FALSE]
      if(messages.screen) {
        cat("krige.bayes: sample from the (joint) posterior of phi and tausq.rel\n")
        print(rbind(phi = phi.unique[, 1], tausq.rel = 
                    phi.unique[, 2], frequency = ind.table))
        cat("\n")
      }
      vecpars.back.order <- order(order(ind))
      ##
      ## sampling sigmasq
      ##
      sigmasq.sam <- rinvchisq(n = n.posterior, df = df.model,
                               scale = rep(as.vector(phidist$s2)[ind.unique], ind.table))
      ##
      ## sampling beta
      ##
      if(beta.size == 1) {
        vec.beta <- rep(as.vector(phidist$beta)[ind.unique],ind.table)
        vec.vbeta <- rep(as.vector(phidist$varbeta)[ind.unique], ind.table)
        beta.sam <- vec.beta + sqrt(sigmasq.sam * vec.vbeta) * rnorm(n.posterior, mean = 0, sd = 1)
      }
      else {
        ind.beta <- matrix(phidist$beta, ncol = beta.size)[ind.unique,,drop=FALSE]
        ind.beta <- ind.beta[rep(1:ind.length, ind.table),,drop=FALSE]
        ind.vbeta <- matrix(phidist$varbeta, ncol = 
                            beta.size^2)[ind.unique,,drop=FALSE]
        ind.vbeta <- ind.vbeta[rep(1:ind.length, ind.table),,drop=FALSE] * sigmasq.sam
        ##      print("2.4: try to speed up this bit!")
        temp.res <- apply(ind.vbeta, 1, rMVnorm,
                          beta.size = beta.size)
        beta.sam <- ind.beta + t(temp.res)
        remove("temp.res")
      }
      ##
      ## summaries of the posterior
      ##
      if(beta.size == 1) {
        trend.mean <- mean(beta.sam)
        trend.median <- median(beta.sam)
      }
      else {
        trend.mean <- colMeans(beta.sam)
        trend.median <- apply(beta.sam, 2, median)
      }
      S2.mean <- mean(sigmasq.sam)
      S2.median <- median(sigmasq.sam)
      phi.marg <- rowSums(phidist$probphitausq)
      .marg <- phi.marg/(sum(phi.marg))
      phi.mean <- phi.discrete %*% phi.marg
      phi.median <- median(phi.sam[, 1])
      phi.mode <- phi.discrete[which.min(abs(phi.marg - max(phi.marg)))]
      tausq.rel.marg <- colSums(phidist$probphitausq)
      tausq.rel.marg <- tausq.rel.marg/(sum(tausq.rel.marg))
      tausq.rel.mean <- tausq.rel.discrete %*% tausq.rel.marg
      tausq.rel.median <- median(phi.sam[, 2])
      tausq.rel.mode <- tausq.rel.discrete[which.min(abs(tausq.rel.marg - max(tausq.rel.marg)))]
      ##
      ## Computing the conditional mode of beta and sigmasq;
      ## conditional on the mode of (phi, tausq.rel)
      ##
      mode.ind <- which(phidist$probphitausq == max(phidist$probphitausq))
      phi.tausq.rel.mode <- phidist$phitausq[mode.ind, 1:2, drop = FALSE]
      if(nrow(phi.tausq.rel.mode) > 1){
        cat("krige.bayes: WARNING: multiple modes for phi.tausq.rel. Using the first one to compute modes of beta and sigmasq.\n")
        cat("krige.bayes: modes found are:\n")
        print(phi.tausq.rel.mode)
        phi.tausq.rel.mode <- phi.tausq.rel.mode[1,]
      }
      if(prior$beta.prior == "normal" && npr > 1)
        info.id <- mode.ind
      else info.id <- 1
      modes <- beta.sigmasq.post(n = n, beta.info = beta.info[[info.id]],
                                 sigmasq.info = sigmasq.info[[info.id]],
                                 env.dists = dists.env,
                                 model = list(cov.model = model$cov.model,
                                   kappa = model$kappa),
                                 xmat = trend.data, y = data,
                                 phi = phi.tausq.rel.mode[1],
                                 tausq.rel = phi.tausq.rel.mode[2],
                                 dets = FALSE,
                                 do.prediction.moments = FALSE,
                                 do.prediction.simulations = FALSE,
                                 env.pred = pred.env, signal = signal)
      beta.mode.cond <- modes$beta.post
      S2.mode.cond <- modes$S2.post
      rm(modes)
      ##
      ## preparing output on posterior distribution
      ##
      if(beta.size == 1)
        kb$posterior$beta$summary <-
          c(mean = trend.mean, median = trend.median, mode.cond = beta.mode.cond)
      else kb$posterior$beta$summary <-
        cbind(mean = trend.mean, median = trend.median, mode.cond = beta.mode.cond)
      kb$posterior$sigmasq$summary <-
        c(mean = S2.mean, median = S2.median, mode.cond = S2.mode.cond)
      kb$posterior$phi$summary <-
        c(mean = phi.mean, median = phi.median, mode = phi.mode)
      if(prior$tausq.rel.prior != "fixed")
        kb$posterior$tausq.rel$summary <- c(mean = tausq.rel.mean,
                                            median = tausq.rel.median,
                                            mode = tausq.rel.mode)
      kb$posterior$sample <-
        as.data.frame(cbind(drop(as.matrix(beta.sam)[vecpars.back.order,  ]),
                            sigmasq.sam[vecpars.back.order], phi.sam[,1]))
      beta.sam <- sigmasq.sam <- NULL
      names(kb$posterior$sample) <- c(beta.names, "sigmasq", "phi")
      kb$posterior$sample$tausq.rel <- phi.sam[,2]
      phi.lev <- unique(phidist$phitausq[, 1])
      kb$posterior$phi$phi.marginal <-
        data.frame(phi = phi.lev, expected = rowSums(phidist$probphitausq),
                   sampled = as.vector(table(factor(phi.sam[, 1],
                     levels = phi.lev)))/n.posterior)
      tausq.rel.lev <- unique(phidist$phitausq[, 2])
      if(prior$tausq.rel.prior != "fixed")
        kb$posterior$tausq.rel$tausq.rel.marginal <-
          data.frame(tausq.rel = tausq.rel.lev, expected = colSums(phidist$probphitausq),
                     sampled = as.vector(table(factor(phi.sam[, 2],
                       levels = tausq.rel.lev)))/n.posterior)
    }
    ##
    ##  computing results for the predictive distribution 
    ##
    if(do.prediction){
      kb$predictive$distribution <- "obtained by numerical approximation" 
      if(messages.screen)
        cat("krige.bayes: starting prediction at the provided locations\n")
      ##
      ## defining samples to be taken from the predictive
      ##
      if(n.predictive == n.posterior) {
        include.it <- FALSE
        phi.sam <- phidist$phitausq[ind,  ]
        message.prediction <- c(message.prediction, 
                                "krige.bayes: phi/tausq.rel samples for the predictive are same as for the posterior")
        if(messages.screen)
          cat(message.prediction, "\n")
      }
      else {
##phidist,  ind
## n.predictive
        include.it <- TRUE
        ind <- sample((1:(dim(phidist$phitausq)[1])), n.predictive,
                      replace = TRUE, prob = as.vector(phidist$probphitausq))
        ind.unique <- sort(unique(ind))
        ind.length <- length(ind.unique)
        ind.table <- table(ind)
        phi.unique <- phidist$phitausq[ind.unique,, drop=FALSE]
        message.prediction <- c(message.prediction, 
                                "krige.bayes: phi/tausq.rel samples for the predictive are NOT the same as for the posterior ")
        if(messages.screen) {
          cat(message.prediction, "\n")
          cat("krige.bayes: samples and their frequencies from the distribution of  phi and tau.rel when drawing from the predictive distribution\n")
          print(rbind(phi = phi.unique[, 1], tausq.rel
                      = phi.unique[, 2], frequency = ind.table))
        }
        phi.sam <- phidist$phitausq[ind,  ]
        vecpars.back.order <- order(order(ind))
      }
      if(moments){
        if(messages.screen)
          cat("krige.bayes: computing moments of the predictive distribution\n")
        kb$predictive$mean <- get("expect", envir=expect.env)/phidist$sum.prob
        remove("expect", envir=expect.env)
        kb$predictive$variance <-
          (get("expect2", envir=expect.env)/phidist$sum.prob) -
            ((kb$predictive$mean)^2)
        remove("expect2", envir=expect.env)
      }
    }
  }
  ##
  ## Backtransforming predictions
  ##
  if((do.prediction && moments) & (abs(lambda-1) > 0.001)){
    kb$predictive <-
      backtransform.moments(lambda = lambda,
                            mean = kb$predictive$mean,
                            variance = kb$predictive$variance,
                            distribution = kb$predictive$distribution,
                            n.simul = n.back.moments)
  }
  ##
  ## ======================= PART 4 ==============================
  ##                Sampling from the predictive
  ## =============================================================
  ##
  if(do.prediction && simulations.predictive){
    if(is.R()){
      if(cov.model.number > 12)
        stop("simulation in krige.bayes not implemented for the choice of correlation function")
    }
    else
      if(cov.model.number > 10)
        stop("simulation in krige.bayes not implemented for the correlation function chosen")
    krige.bayes.aux20 <- function(phinug){
      iter <- get("counter", envir=counter.env)
      if(messages.screen & prior$phi.prior != "fixed")
        krige.bayes.counter(.temp.ap = iter, n.points = ind.length)
      phinug <- as.vector(phinug)
      phi <- phinug[1]
      tausq.rel <- phinug[2]
      phi.ind <- which.min(abs(phi.discrete - phi))
      nug.ind <- which.min(abs(tausq.rel.discrete - tausq.rel))
      v0 <- cov.spatial(obj = get("d0", envir=pred.env),
                        cov.model = cov.model, kappa = kappa,
                        cov.pars = c(1, phi))
      ## care here, reusing object b
      b <- .bilinearformXAY(X = as.vector(cbind(data, trend.data)),
                           lowerA = as.vector(inv.lower[phi.ind, nug.ind,,drop=TRUE]),
                           diagA = as.vector(inv.diag[phi.ind, nug.ind,,drop=TRUE]), 
                           Y = as.vector(v0))
      tv0ivdata <- drop(b[1,])
      b <- t(get("trend.loc", envir=pred.env)) -  b[-1,, drop=FALSE]
      ##
      tmean <- tv0ivdata + drop(crossprod(b,as.vector(phidist$
                                                 beta[phi.ind, nug.ind,  ])))
      tv0ivdata <- NULL
      Nsims <- ind.table[iter]
      if (signal) Dval <- 1.0 else Dval <-  1.0 + tausq.rel
      iter.env <- sys.frame(sys.nframe())
      ## removing this because seens redundant with part (**) below
      ##      if((tausq.rel < 1e-12) & (!is.null(get("loc.coincide", envir=pred.env))))
      ##        tmean[get("loc.coincide", envir=pred.env)] <- get("data.coincide", envir=pred.env)
      coincide.cond <- (((tausq.rel < 1e-12) | !signal) & !is.null(get("loc.coincide", envir=pred.env)))
      if(coincide.cond){
        nloc <- ni - n.loc.coincide
        ind.not.coincide <- -(get("loc.coincide", envir=pred.env))
        v0 <- v0[, ind.not.coincide, drop=FALSE]
        tmean <- tmean[ind.not.coincide]
        b <- b[,ind.not.coincide, drop=FALSE]
      }
      else{
        nloc <- ni
        ind.not.coincide <- TRUE
      }
      par.set <- get("parset", envir=counter.env)
      if(prior$beta.prior == "normal" && npr > 1)
        info.id <- par.set
      else info.id <- 1
      if(any(beta.info[[info.id]]$iv == Inf))
        vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
      else
        vbetai <- matrix(drop(phidist$varbeta[phi.ind, nug.ind,  ]),
                         ncol = beta.size, nrow = beta.size)
      simul <- matrix(NA, nrow=ni, ncol=Nsims)
      if(nloc > 0)
        simul[ind.not.coincide,] <-
          .cond.sim(env.loc = base.env, env.iter = iter.env,
                   loc.coincide = get("loc.coincide", envir=pred.env),
                   coincide.cond = coincide.cond,
                   tmean = tmean,
                   Rinv = list(lower= drop(inv.lower[phi.ind, nug.ind,]),
                     diag = drop(inv.diag[phi.ind, nug.ind,])),
                   mod = list(beta.size = beta.size, nloc = nloc,
                     Nsims = Nsims, n = n, Dval = Dval,
                     df.model = df.model,
                     s2 = phidist$s2[phi.ind, nug.ind],
                     cov.model.number = cov.model.number,
                     phi = phi, kappa = kappa, nugget=tausq.rel),
                   vbetai = vbetai,
                   fixed.sigmasq = (sigmasq.info$df.sigmasq == Inf))
      ## (**) this made previous bit redundant
      if(coincide.cond)
        simul[get("loc.coincide", envir=pred.env),] <-
          rep(get("data.coincide", envir=pred.env), Nsims)
      remove("v0", "b", "tmean")
      assign("counter", (iter + 1), envir=counter.env)
      assign("parset", get("parset", envir=counter.env)+1, envir=counter.env)
      ##
      ## Back transforming (To be include in C code???)
      ##
      if(abs(lambda - 1) > 0.001){
        return(BCtransform(x=simul, lambda=lambda, inverse=TRUE)$data)
      }
      else
        return(simul)
    }
    ##
    counter.env <- new.env()
    if(messages.screen){
      cat("krige.bayes: sampling from the predictive\n")
      if(prior$phi.prior != "fixed")
        cat(paste("             Number of parameter sets: ", ind.length,"\n"))
    }
    assign("counter", 1, envir=counter.env)
    assign("parset", 1, envir=counter.env)
     kb$predictive$simulations <- 
      matrix(unlist(apply(phi.unique, 1, krige.bayes.aux20)),
             ncol = n.predictive)
    remove("inv.lower", "inv.diag", "counter.env", "pred.env")
#    gc(verbose=FALSE)
    if(messages.screen)
      if(abs(lambda-1) > 0.001) 
        cat("krige.bayes: Box-Cox data transformation performed.\n             Simulations back-transformed to the original scale\n")
    if(messages.screen)
      cat("krige.bayes: preparing summaries of the predictive distribution\n")
    ##
    ## mean/quantiles/probabilities estimators from simulations
    ##
    kb$predictive <- c(kb$predictive, 
                       statistics.predictive(simuls=kb$predictive$simulations,
                                             mean.var = mean.estimator,
                                             quantile = quantile.estimator,
                                             threshold = probability.estimator,
                                             sim.means = sim.means, 
                                             sim.vars = sim.vars))
    ##
    ## Mean and variance of each (conditionally) simulated field 
    ##
    if(sim.means && exists("vecpars.back.order"))
      kb$predictive$sim.means[vecpars.back.order]
    if(sim.vars && exists("vecpars.back.order"))
      kb$predictive$sim.vars[vecpars.back.order]
    ##
    ## keeping or not samples from predictive
    ##
    if(keep.simulations){
      if(prior$phi.prior != "fixed")
        kb$predictive$simulations <-
          kb$predictive$simulations[, vecpars.back.order]
    }
    else{
      kb$predictive$simulations <- NULL
#      gc(verbose=FALSE)
    }
    ##
    ## recording samples from  predictive if different from the posterior
    ##
    if(prior$phi.prior != "fixed"){
      if(include.it){
        phi.lev <- unique(phidist$phitausq[, 1])
        kb$predictive$phi.marginal <-
          data.frame(phi = phi.lev,
                     expected = rowSums(phidist$probphitausq),
                     sampled = as.vector(table(factor(phi.sam[, 1],
                       levels = phi.lev)))/n.predictive)
        tausq.rel.lev <- unique(phidist$phitausq[, 2])
        if(prior$tausq.rel.prior != "fixed")
          data.frame(tausq.rel = tausq.rel.lev,
                     expected = colSums(phidist$probphitausq),
                     sampled = as.vector(table(factor(phi.sam[, 2],
                       levels = tausq.rel.lev)))/n.predictive)
        else
          kb$predictive$tausq.rel.marginal <-
            paste("fixed tausq.rel with value =", tausq.rel)
        kb$predictive$tausq.rel.marginal <-
          data.frame(tausq.rel = tausq.rel.lev,
                     expected = colSums(phidist$probphitausq),
                     sampled = as.vector(table(factor(phi.sam[, 2],
                       levels = tausq.rel.lev)))/n.predictive)
      }
    }
  }
  if(!do.prediction) kb$predictive <- "no prediction locations provided"
  kb$.Random.seed <- seed
  kb$max.dist <- data.dist.max
  kb$call <- call.fc
  attr(kb, "prediction.locations") <- call.fc$locations
  attr(kb, "parent.env") <- parent.frame()
  if(!is.null(call.fc$coords))
    attr(kb, "data.locations") <- call.fc$coords
  else attr(kb, "data.locations") <- substitute(a$coords, list(a=substitute(geodata)))
  if(do.prediction) attr(kb, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  if(!is.null(call.fc$borders))
    attr(kb, "borders") <- call.fc$borders
  oldClass(kb) <- c("krige.bayes", "variomodel")
  #if(messages.screen) cat("krige.bayes: done!\n")
  return(kb)
}

".values.krige.bayes" <- 
  function (obj, values.to.plot, number.col, messages.screen)
{
  if(messages.screen)
    switch(values.to.plot,
           mean = cat("mapping the means of the predictive distribution\n"),
           variance = cat("mapping the variances of the predictive distribution\n"),
           mean.simulations = cat("mapping the means of simulations from the predictive distribution\n"),
           variance.simulations = cat("mapping the variances of simulations from the predictive distribution\n"),
           median = cat("mapping the medians of the predictive distribution\n"),
           uncertainty = cat("mapping the uncertainty of the predictive distribution\n"),
           quantiles = cat("mapping a quantile of the predictive distribution\n"),
           probabilities = cat("mapping a probability of beeing bellow threshold of the predictive distribution\n"),
           simulation = cat("mapping a simulation from the predictive distribution\n"),
           stop("wrong specification for values to plot")
           )
  switch(values.to.plot,
         mean = obj$predictive$mean,
         variance = obj$predictive$variance,
         mean.simulations = obj$predictive$mean.simulations,
         variance.simulations = obj$predictive$variance.simulations,
         median = obj$predictive$median,
         uncertainty = obj$predictive$uncertainty,
         quantiles = {
           if(!is.vector(obj$predictive$quantiles))
             if(is.null(number.col)) stop("argument number.col must be provided")
             else as.matrix(obj$predictive$quantiles)[,number.col]
           else as.matrix(obj$predictive$quantiles)[,1]
         },
         probabilities = {
           if(!is.vector(obj$predictive$probab)){
             if(is.null(number.col)) stop("argument number.col must be provided")
             else as.matrix(obj$predictive$probab)[,number.col]
           }
           else as.matrix(obj$predictive$probab)[,1]
         },
         simulation = {
           if(is.null(number.col)) stop("argument number.col must be provided")
           as.matrix(obj$predictive$simulations)[,number.col]
         },
         stop("wrong specification for values to plot")
         )
}

".prepare.graph.krige.bayes" <-
  function (obj, locations, borders, borders.obj=NULL,
            values.to.plot, number.col, xlim, ylim, messages, ...) 
{
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  if (mode(values.to.plot) != "numeric")
    values.to.plot <- .values.krige.bayes(obj=obj, values.to.plot=values.to.plot,
                                          number.col=number.col,
                                          messages.screen=messages.screen)
  locations <- locations[order(locations[, 2], locations[,1]), ]
  xx <- as.numeric(levels(as.factor(locations[, 1])))
  nx <- length(xx)
  yy <- as.numeric(levels(as.factor(locations[, 2])))
  ny <- length(yy)
  values.loc <- rep(NA, nrow(locations))
  if(length(values.loc) == length(values.to.plot)) values.loc <- values.to.plot
  if(!is.null(borders.obj)){
    borders.obj <- as.matrix(as.data.frame(borders.obj))
    inout.vec  <- .geoR_inout(locations, borders.obj)
    values.loc[inout.vec] <- values.to.plot
    rm("inout.vec")
  }
  if (!is.null(borders)){
    borders <- as.matrix(as.data.frame(borders))
    dimnames(borders) <- list(NULL, NULL)
    if(!(!is.null(borders.obj) && identical(borders,borders.obj))){
      inout.vec  <- .geoR_inout(locations, borders)
      if(length(values.loc[inout.vec]) == length(values.to.plot))
        values.loc[inout.vec] <- values.to.plot
      values.loc[!inout.vec] <- NA
      rm("inout.vec")
    }
  }
  ##
  if (missing(xlim) || is.null(xlim))
    if(is.null(borders)) xlim <- NULL
    else xlim <- range(borders[,1]) 
  if (missing(ylim) || is.null(ylim))
    if(is.null(borders)) ylim <- NULL
    else ylim <- range(borders[,2])
  coords.lims <- set.coords.lims(coords=locations, xlim=xlim, ylim=ylim)
  coords.lims[,1] <- coords.lims[,1] + c(-.025, .025) * diff(coords.lims[,1])
  coords.lims[,2] <- coords.lims[,2] + c(-.025, .025) * diff(coords.lims[,2])
  return(list(x=xx, y=yy, values = matrix(values.loc,ncol=ny), coords.lims=coords.lims))
}

"image.krige.bayes" <-
  function (x, locations, borders, 
            values.to.plot = c("mean", "variance",
              "mean.simulations", "variance.simulations",
              "quantiles", "probabilities", "simulation"),
            number.col, coords.data,
            x.leg, y.leg, messages, ...) 
{
  ldots <- list(...)
###  ldots <- match.call(expand.dots = FALSE)$...
  ##  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations))
    locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(mode(values.to.plot) != "numeric")
    values.to.plot <-
      match.arg(values.to.plot,
                choices = c("mean", "variance",
                  "mean.simulations", "variance.simulations",
                  "quantiles", "probabilities", "simulation"))
  if(missing(borders)){
    if(!is.null(attr(x, "borders"))) borders.arg <- borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    else
      borders.arg <- borders <- eval(x$call$geodata, envir= attr(x, "parent.env"))$borders
    #borders.arg <- borders <- NULL
  }
  else{
    borders.arg <- borders
    if(is.null(borders)) borders <- eval(attr(x, "borders"), envir= attr(x, "parent.env"))
  }
  if(missing(number.col)) number.col <- NULL
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
    if (mode(values.to.plot) != "numeric")
      values.to.plot <- .values.krige.bayes(obj=x, values.to.plot=values.to.plot,
                                            number.col=number.col,
                                            messages.screen=messages)
    do.call("plot.1d", c(list(x = list(coords=locations, data = values.to.plot),
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d", data="prediction")))
  }
  else{
    ldots.image <- .ldots.set(ldots, type="image", data="prediction")
    locations <- .prepare.graph.krige.bayes(obj=x,
                                           locations=locations,
                                           borders=borders,
                                           borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                           values.to.plot=values.to.plot,
                                           number.col = number.col,
                                           xlim= ldots.image$xlim,
                                           ylim= ldots.image$ylim,
                                           messages=messages)
    do.call("image", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.image))
    if(!is.null(coords.data)) points(coords.data)
    if(!is.null(borders.arg)) polygon(borders, lwd=2)
    if(is.null(ldots$col)) ldots$col <- heat.colors(12)
    if(!is.null(x.leg) & !is.null(y.leg)){
      do.call("legend.krige", c(list(x.leg=x.leg, y.leg=y.leg,
                                     values=locations$values),
                                     ldots))
    }
  }
  return(invisible())
}

"contour.krige.bayes" <-
  function (x, locations, borders, 
            values.to.plot = c("mean", "variance",
              "mean.simulations", "variance.simulations",
              "quantiles", "probabilities", "simulation"),
            filled= FALSE, number.col, coords.data,
            x.leg, y.leg, messages, ...) 
{
  ldots <- list(...)
#  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations))
    locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(mode(values.to.plot) != "numeric")
    values.to.plot <-
      match.arg(values.to.plot,
                choices = c("mean", "variance",
                  "mean.simulations", "variance.simulations",
                  "quantiles", "probabilities", "simulation"))
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
  if(missing(number.col)) number.col <- NULL
  if(missing(coords.data)) coords.data <- NULL
  else
    if(all(coords.data == TRUE))
      coords.data <-  eval(attr(x, "data.locations"), envir= attr(x, "parent.env"))
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D'){
    if (mode(values.to.plot) != "numeric")
      values.to.plot <- .values.krige.bayes(obj=x, values.to.plot=values.to.plot,
                                            number.col=number.col,
                                            messages.screen=messages)
    do.call("plot.1d", c(list(x = list(coords=locations, data = values.to.plot),
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d", data="prediction")))
  }
  else{
    if(filled)
      ldots.contour <- .ldots.set(ldots, type="filled.contour",
                                 data="prediction")
    else
      ldots.contour <- .ldots.set(ldots, type="contour",
                                 data="prediction")
    if(is.null(ldots.contour$asp)) ldots.contour$asp=1
    locations <- .prepare.graph.krige.bayes(obj=x, locations=locations,
                                           borders=borders,
                                           borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                           values.to.plot=values.to.plot,
                                           number.col = number.col,
                                           xlim = ldots.contour$xlim,
                                           ylim = ldots.contour$ylim,
                                           messages=messages)
    if(filled){
      if(is.null(ldots.contour$plot.axes)){
        ldots.contour$plot.axes <- quote({
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
                                z=locations$values),
                           ldots.contour))
      ##
      ## adding borders
      ##
      if(!is.null(borders.arg)) polygon(borders, lwd=2)
      if(!is.null(coords.data)) points(coords.data, pch=20)
    }
  }
  return(invisible())
}

"persp.krige.bayes" <-
  function (x, locations, borders, 
            values.to.plot = c("mean", "variance",
              "mean.simulations", "variance.simulations",
              "quantiles", "probabilities", "simulation"), number.col, messages, ...) 
{
  ldots <- list(...)
#  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  attach(x, pos=2, warn.conflicts=FALSE)
  on.exit(detach(2))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"), envir= attr(x, "parent.env"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2) stop("locations must be a matrix or data-frame with two columns")
  if(mode(values.to.plot) != "numeric"){
    values.to.plot <- match.arg(values.to.plot,
                                choices = c("mean", "variance",
                                  "mean.simulations",
                                  "variance.simulations",
                                  "quantiles", "probabilities",
                                  "simulation"))
  }
  if(missing(borders)) borders <- NULL
  if(missing(number.col)) number.col <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D'){
    if (mode(values.to.plot) != "numeric")
      values.to.plot <- .values.krige.bayes(obj=x, values.to.plot=values.to.plot,
                                            number.col=number.col,
                                            messages.screen=messages)
    do.call("plot.1d", c(list(x = list(coords=locations, data = values.to.plot),
                              x1vals = unique(round(locations[,1], digits=12))),
                         .ldots.set(ldots, type="plot.1d", data="prediction")))
  }
  else{
    ldots.persp <- .ldots.set(ldots, type="persp", data="prediction")
    locations <- .prepare.graph.krige.bayes(obj=x, locations=locations,
                                           borders=borders,
                                           borders.obj = eval(attr(x, "borders"), envir= attr(x, "parent.env")),
                                           values.to.plot=values.to.plot,
                                           xlim= ldots.persp$xlim,
                                           ylim= ldots.persp$ylim,
                                           number.col = number.col,
                                           messages=messages)
    do.call("persp", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.persp))
  }
  return(invisible())
}

"model.control" <-
  function(trend.d = "cte", trend.l = "cte", cov.model = "matern",
           kappa = 0.5, aniso.pars = NULL, lambda = 1)
{
  cov.model <-
    match.arg(cov.model, choices = .geoR.cov.models)
  if(cov.model == "powered.exponential" & (kappa <= 0 | kappa > 2))
    stop("model.control: for power exponential correlation model the parameter kappa must be in the interval \\(0,2\\]")
  ##  if(any(cov.model == c("exponential", "gaussian", "spherical",
  ##           "circular", "cubic", "wave", "powered.exponential",
  ##           "cauchy", "gneiting", "pure.nugget")))
  ##    kappa <- NULL
  if(!is.null(aniso.pars)) 
    if(length(aniso.pars) != 2 | mode(aniso.pars) != "numeric")
      stop("model.control: anisotropy parameters must be a vector with two elements: rotation angle (in radians) and anisotropy ratio (a number > 1)")
  res <- list(trend.d = trend.d, trend.l = trend.l,
              cov.model = cov.model,
              kappa=kappa, aniso.pars=aniso.pars, lambda=lambda)
  oldClass(res) <- "model.geoR"
  return(res)
}

"post2prior" <- function(obj)
{
  if(length(class(obj)) == 0)
    stop("post.prior: argument must be an object of the class `krige.bayes` or `posterior.krige.bayes`")
  if(any(class(obj) == "krige.bayes")) obj <- obj$posterior
  if(all(class(obj) != "posterior.krige.bayes"))
    stop("post.prior: argument must be an object of the class `krige.bayes` or `posterior.krige.bayes`")
  ##
  ## beta
  ##
  if(!is.null(obj$beta$status) &&
     obj$beta$status == "fixed"){
    beta.prior <- "fixed"
    beta <- obj$beta$fixed.value
    beta.var.std <- NULL
  }
  else{
    beta.prior <- obj$beta$conditional.distribution
    beta <- obj$beta$pars$mean
    beta.var.std <- obj$beta$pars$var
  }
  ##
  ## sigmasq
  ##
  if(!is.null(obj$sigmasq$status) &&
     obj$sigmasq$status == "fixed"){
    sigamsq.prior <- "fixed"
    sigmasq <- obj$sigmasq$fixed.value
    df.sigmasq <- NULL
  }
  else{
    sigmasq.prior <- obj$sigmasq$conditional.distribution
    sigmasq <- obj$sigmasq$pars$S2
    df.sigmasq <-  obj$sigmasq$pars$df
  }
  ##
  ## phi
  ##
  if(!is.null(obj$phi$status) &&
     obj$phi$status == "fixed"){
    phi.prior <- "fixed"
    phi <- obj$phi$fixed.value
    phi.discrete <- NULL
  }
  else{
    phi.prior <- obj$phi$phi.marginal[,"expected"]
    phi <- NULL
    phi.discrete <- obj$phi$phi.marginal[,"phi"]
  }
  ##
  ## tausq.rel
  ##
  if(!is.null(obj$tausq.rel$status) &&
     obj$tausq.rel$status == "fixed"){
    tausq.rel.prior <- "fixed"
    tausq.rel <- obj$tausq.rel$fixed.value
    tausq.rel.discrete <- NULL 
  }
  else{
    tausq.rel.prior <- obj$tausq.rel$tausq.rel.marginal[,"expected"]
    tausq.rel <- 0
    tausq.rel.discrete <- obj$tausq.rel$tausq.rel.marginal[,"tausq.rel"]
  }
  ##
  res <- prior.control(beta.prior = beta.prior, beta = beta,
                       beta.var.std = beta.var.std,
                       sigmasq.prior = sigmasq.prior, 
                       sigmasq = sigmasq,  df.sigmasq = df.sigmasq,
                       phi.prior = phi.prior, 
                       phi = phi, phi.discrete = phi.discrete, 
                       tausq.rel.prior = tausq.rel.prior,
                       tausq.rel = tausq.rel,
                       tausq.rel.discrete = tausq.rel.discrete)
  res$joint.phi.tausq.rel <- obj$joint.phi.tausq.rel
  res$dep.prior <- TRUE
  return(res)
}

"prior.control" <-
  function(beta.prior = c("flat", "normal", "fixed"),
           beta = NULL, beta.var.std = NULL,
           sigmasq.prior = c("reciprocal",  "uniform", "sc.inv.chisq",  "fixed"),
           sigmasq = NULL,  df.sigmasq = NULL,
           phi.prior = c("uniform", "exponential", "fixed", "squared.reciprocal","reciprocal"),
           phi = NULL, phi.discrete = NULL, 
           tausq.rel.prior = c("fixed", "uniform", "reciprocal"),
           tausq.rel = 0,
           tausq.rel.discrete = NULL)
{
  ##
  ## 1. Checking parameters for the priors
  ##
  ##
  ## beta
  ##
  beta.prior <- match.arg(beta.prior)
  if(beta.prior == "fixed" & is.null(beta))
    stop("prior.control: argument beta must be provided with fixed prior for this parameter")
  if(beta.prior == "normal"){
    if(is.null(beta) | is.null(beta.var.std))
      stop("prior.control: arguments `beta` and `beta.var.std` must be provided with normal prior for the parameter beta")
  }
  ##
  ## sigmasq
  ##
  sigmasq.prior <- match.arg(sigmasq.prior)
  if(sigmasq.prior == "fixed" & is.null(sigmasq))
    stop("prior.control: argument `sigmasq' must be provided with fixed prior for the parameter sigmasq")
  if(sigmasq.prior == "sc.inv.chisq")
    if(is.null(sigmasq) | is.null(df.sigmasq))
      stop("prior.control: arguments `sigmasq` and `df.sigmasq' must be provided for this choice of prior distribution")
  if(!is.null(sigmasq))
    if(any(sigmasq < 0))
      stop("prior.control: negative values not allowed for `sigmasq'")
  ##
  ## phi
  ##
  if(!is.null(phi) && length(phi) > 1)
    stop("prior.control: length of phi must be one. Use phi.prior and phi.discrete to specify the prior for phi or enter a single fixed value for phi")
  if(mode(phi.prior) == "numeric"){
    phi.prior.probs <- phi.prior
    phi.prior <- "user"
    if(is.null(phi.discrete))
      stop("prior.control: argument phi.discrete with support points for phi must be provided\n")
    if(length(phi.prior.probs) != length(phi.discrete))
      stop("prior.control: user provided phi.prior and phi.discrete have incompatible dimensions\n")
    if(round(sum(phi.prior.probs), digits=8) != 1)
      stop("prior.control: prior probabilities provided for phi do not add up to 1")
  }
  else
    phi.prior <- match.arg(phi.prior,
                           choices = c("uniform", "exponential", "fixed",
                             "squared.reciprocal","reciprocal"))
  if(phi.prior == "fixed"){
    if(is.null(phi)){
      stop("prior.control: argument `phi` must be provided with fixed prior for this parameter")
    }
    phi.discrete <- phi
  }
  else{
    if(phi.prior == "exponential" && (is.null(phi) | (length(phi) > 1)))
      stop("prior.control: argument `phi` must be provided when using the exponential prior for the parameter phi")
    ## instead of commented below the probability at zero is set to zero
    ##    if(any(phi.prior == c("reciprocal", "squared.reciprocal")) &
    ##       any(phi.discrete == 0)){
    ##      warning("degenerated prior at phi = 0. Excluding value phi.discrete[1] = 0")
    ##      phi.discrete <- phi.discrete[phi.discrete > 1e-12]
    ##    }
#    if(is.null(phi.discrete))
#      stop("prior.control: argument phi.discrete with support points for phi must be provided\n")
#    else{
    if(!is.null(phi.discrete)){
      discrete.diff <- diff(phi.discrete)
      if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
        stop("prior.control: the current implementation requires equally spaced values in the argument `phi.discrete`\n")
    } 
  }
  if(any(phi.discrete < 0)) stop("prior.control: negative values not allowed for parameter phi")
  ##
  ## tausq
  ##
  if(length(tausq.rel) > 1)
    stop("prior.control: length of tausq.rel must be one. Use tausq.rel.prior and tausq.rel.discrete to specify the prior for tausq.rel or enter a single fixed value for tausq.rel")
  if(mode(tausq.rel.prior) == "numeric"){
    tausq.rel.prior.probs <- tausq.rel.prior
    tausq.rel.prior <- "user"
    if(is.null(tausq.rel.discrete))
      stop("prior.control: argument tausq.rel.discrete with support points for tausq.rel must be provided\n")
    if(length(tausq.rel.prior.probs) != length(tausq.rel.discrete))
      stop("prior.control: user provided tausq.rel.prior and tausq.rel.discrete have incompatible dimensions\n")
    if(round(sum(tausq.rel.prior.probs), digits=8) != 1)
      stop("prior.control: prior probabilities for tausq.rel provided do not add up to 1")
  }
  else
    tausq.rel.prior <- match.arg(tausq.rel.prior, choices = c("fixed", "uniform", "reciprocal"))
  if(tausq.rel.prior == "fixed"){
    if(is.null(tausq.rel) | mode(tausq.rel) != "numeric")
      stop("prior.control: argument `tausq.rel` must be provided with fixed prior for the parameter tausq.rel")
    tausq.rel.discrete <- tausq.rel
  }
  else{
    if(is.null(tausq.rel.discrete))
      stop("prior.control: argument `tausq.rel.discrete` must be provided with chosen prior for the parameter tausq.rel")  
    discrete.diff <- diff(tausq.rel.discrete)
    if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
      stop("prior.control: the current implementation requires equally spaced values in the argument `tausq.rel.discrete`\n")
  }
  if(any(tausq.rel.discrete < 0))
    stop("prior.control: negative values not allowed for parameter tausq.rel")
  ##
  ## Further checks on dimensions
  ##
  if(phi.prior != "fixed"){
    if(mode(phi.discrete) == "numeric"){
      if(is.null(tausq.rel.discrete)) nsets <- length(phi.discrete)
      else nsets <- length(phi.discrete) * length(tausq.rel.discrete)
    }
    else nsets <- 0
    if(sigmasq.prior == "sc.inv.chisq"){
      if(length(sigmasq) == nsets) dep.prior <- TRUE
      else dep.prior <- FALSE
    }
    else dep.prior <- FALSE
    if(beta.prior == "normal"){
      if(dep.prior){
        if(((length(beta)/nsets)^2) != (length(beta.var.std)/nsets))
          stop("prior.control: beta and beta.var.std have incompatible dimensions")
      }
      else{
        if((length(beta))^2 != length(beta.var.std))
          stop("prior.control: beta and beta.var.std have incompatible dimensions")
#        if(exists("trySilent")){
        if(inherits(try(.solve.geoR(beta.var.std), silent=TRUE), "try-error"))
          stop("prior.control: singular matrix in beta.var.std")
        if(inherits(try(chol(beta.var.std), silent=TRUE), "try-error"))
          stop("prior.control: no Cholesky decomposition for beta.var.std")
#        }
#        else{
#          error.now <- options()$show.error.messages
#          if (is.null(error.now) | error.now) 
#            on.exit(options(show.error.messages = TRUE))
#          options(show.error.messages = FALSE)
#          if(inherits(try(.solve.geoR(beta.var.std)), "try-error"))
#            stop("prior.control: singular matrix in beta.var.std")
#          if(inherits(try(chol(beta.var.std)), "try-error"))
#            stop("prior.control: no Cholesky decomposition for beta.var.std")
#        }
        if(any(beta.var.std != t(beta.var.std)))
          stop("prior.control: non symmetric matrix in beta.var.std")
      }
    }
  }
  else dep.prior <- FALSE
  if(!dep.prior & beta.prior == "normal"){
    attr(beta.var.std, "Size") <- length(beta)
  }
  ##
  ip <- list(beta=list(), sigmasq=list(), phi=list(), tausq.rel=list())
  ##
  if(beta.prior == "fixed"){
    ip$beta$status <- "fixed"
    ip$beta$fixed.value <- beta 
  }
  else{
    ip$beta$status <- "random"
    ip$beta$dist <- beta.prior
    if(beta.prior == "flat")
      ip$beta$pars <- c(0, +Inf)
    if(beta.prior == "normal"){
      if(length(beta) == 1)
        ip$beta$pars <- c(mean=beta, var.std=beta.var.std)
      else
        ip$beta$pars <- list(mean=beta, var.std=beta.var.std)
    }
  }
  ##
  if(sigmasq.prior == "fixed"){
    ip$sigmasq$status <- "fixed"
    ip$sigmasq$fixed.value <- sigmasq 
  }
  else{
    ip$sigmasq$status <- "random"
    ip$sigmasq$dist <-  sigmasq.prior
    if(sigmasq.prior == "reciprocal")
      ip$sigmasq$pars <- c(df=0, var=+Inf)
    if(sigmasq.prior == "uniform")
      ip$sigmasq$pars <- c(df=-2, var=+Inf)
    if(sigmasq.prior == "sc.inv.chisq")
      ip$sigmasq$pars <- c(df=df.sigmasq, var=sigmasq)
  }
  ##
  if(phi.prior == "fixed"){
    ip$phi$status <- "fixed"
    ip$phi$fixed.value <- phi
  }
  else{
    ip$phi$status <- "random"
    ip$phi$dist <- phi.prior
    if(is.null(phi.discrete))
      ip$phi$probs <- NULL
    else{
      pd <- as.vector(phi.discrete)
      names(pd) <- NULL
      ip$phi$probs <-
        switch(phi.prior,
               uniform = rep(1, length(pd)),
               exponential = dexp(pd, rate=1/phi),
               squared.reciprocal = ifelse((pd > 0), 1/(pd^2),0),
               reciprocal = ifelse((pd > 0), 1/pd, 0),
               user = phi.prior.probs)
      names(ip$phi$probs) <- phi.discrete
    }
    if(phi.prior == "exponential")
      ip$phi$pars <- c(ip$phi$pars, exp.par=phi)
#    else
    ip$phi$probs <- ip$phi$probs/sum(ip$phi$probs)
  }
  ##
  if(tausq.rel.prior == "fixed"){
    ip$tausq.rel$status <- "fixed"
    ip$tausq.rel$fixed.value <- tausq.rel 
  }
  else{
    ip$tausq.rel$status <- "random"
    ip$tausq.rel$dist <- tausq.rel.prior

    if(is.null(tausq.rel.discrete))
      ip$tausq.rel$probs <- NULL
    else{
      td <- as.vector(tausq.rel.discrete)
      names(td) <- NULL
      ip$tausq.rel$probs <-
        switch(tausq.rel.prior,
               uniform = rep(1, length(td)),
               reciprocal = ifelse((td > 0), 1/td, 0),
               user = tausq.rel.prior.probs)
      names(ip$tausq.rel$probs) <- tausq.rel.discrete
    }
    ip$tausq.rel$probs <- ip$tausq.rel$probs/sum(ip$tausq.rel$probs)
  }
  ## checking valid options for random/fixed parameters
  if(ip$phi$status == "random")
    if(ip$beta$status == "fixed" | ip$sigmasq$status == "fixed")
      stop("random phi with fixed sigmasq and/or beta not implemented")
  ##if(ip$sigmasq$status == "random")
  ##  if(ip$beta$status == "fixed")
  ##    stop("random sigmasq with fixed beta not implemented")
  ##
  res <- list(beta.prior = beta.prior, beta = beta,
              beta.var.std = beta.var.std,
              sigmasq.prior = sigmasq.prior,
              sigmasq = sigmasq, df.sigmasq = df.sigmasq,
              phi.prior = phi.prior, phi = phi,
              phi.discrete = phi.discrete,  
              tausq.rel.prior = tausq.rel.prior,
              tausq.rel = tausq.rel,
              tausq.rel.discrete = tausq.rel.discrete, 
              priors.info = ip, dep.prior = dep.prior)
  oldClass(res) <- "prior.geoR"
  return(res)
}

"output.control" <-
  function(n.posterior, n.predictive, moments, n.back.moments, 
           simulations.predictive, mean.var, quantile,
           threshold, sim.means, sim.vars, signal, messages)
{
  ##
  ## Assigning default values
  ##
  if(missing(n.posterior)) n.posterior <- 1000
  if(missing(n.predictive)) n.predictive <- NULL
  if(missing(moments)) moments <- TRUE
  if(missing(n.back.moments)) n.back.moments <- 1000
  if(missing(simulations.predictive)){
    if(is.null(n.predictive)) simulations.predictive <- NULL
    else
      simulations.predictive <- ifelse(n.predictive > 0, TRUE, FALSE)
  }
  if(missing(mean.var)) mean.estimator <- NULL
  else mean.estimator <- mean.var
  if(missing(quantile))  quantile.estimator <- NULL
  else quantile.estimator <- quantile
  if(missing(threshold)) probability.estimator <- NULL 
  else probability.estimator <- threshold
  if(missing(sim.means)) sim.means <- NULL 
  if(missing(sim.vars)) sim.vars <- NULL 
  if(missing(signal)) signal <- NULL
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  ##
  ##
  ##
  if(!is.null(quantile.estimator) | !is.null(probability.estimator) | !is.null(mean.estimator) | !is.null(sim.means) | !is.null(sim.vars)){
    if(is.null(simulations.predictive)) keep.simulations <- FALSE
    else  keep.simulations <- ifelse(simulations.predictive, TRUE, FALSE)
    simulations.predictive <- TRUE
  }
  else keep.simulations <- NULL
  ##
  if(!is.null(quantile.estimator)){
    if(mode(quantile.estimator) == "numeric")
      if(any(quantile.estimator < 0) | any(quantile.estimator > 1))
        stop("output.control: quantiles indicators must be numbers in the interval [0,1]\n")
    if(all(quantile.estimator == TRUE))
      quantile.estimator <- c(0.025, 0.5, 0.975)
  }
  if(!is.null(probability.estimator)){
    if(mode(probability.estimator) != "numeric")
      stop("output.control: probability.estimator must be a numeric value (or vector) of cut-off value(s)\n")
  }
  res <- list(n.posterior = n.posterior, n.predictive = n.predictive,
              moments = moments, n.back.moments = n.back.moments,
              simulations.predictive = simulations.predictive,
              keep.simulations = keep.simulations,
              mean.estimator = mean.estimator,
              quantile.estimator = quantile.estimator,
              probability.estimator = probability.estimator,
              sim.means = sim.means, sim.vars = sim.vars,
              signal = signal, messages.screen = messages.screen)
  oldClass(res) <- "output.geoR"
  return(res)
}

"beta.sigmasq.post" <-
  function(n, beta.info, sigmasq.info, env.dists,
           model, xmat, y, phi, tausq.rel, do.prediction.moments,
           do.prediction.simulations,
           dets = FALSE, env.pred = NULL, signal)
{
  ##-----------------------------------------------------------------
  ##
  ## dists.env is an environment containing the objetc "data.dist" 
  ## with the distances between pairs os points (output of dists())
  ##
  ## sigmasq.info should contain:
  ##        df.sigmasq: df in prior for sigmasq
  ##                  df.sigmasq = 0 : reciprocal prior for sigmasq
  ##                  df.sigmasq = Inf : fixed sigmasq
  ##        n0S0 : sum of squares in prior for sigmasq
  ## beta.info should contain:
  ##        mivm : computed from the prior of beta: m\prime V^{-1} m
  ##        ivm  : computed from the prior of beta: V^{-1} m
  ##        iv   : computed from the prior of beta: V^{-1}
  ##                  iv = 0 : flat prior for beta
  ##                  iv = Inf : fixed beta
  ##        p    : degrees of freedom correction
  ##                  p = 0  : beta fixed or w/ normal prior
  ##                  p = beta.size : flat prior for beta  
  ## might contain
  ##        beta.fixed
  ##        sigmasq.fixed
  ##-----------------------------------------------------------------
  ##
  R <- varcov.spatial(dists.lowertri = get("data.dist", envir=env.dists),
                      cov.model = model$cov.model,
                      kappa = model$kappa, nugget = tausq.rel,
                      cov.pars = c(1, phi), det = dets)
  iRy.x <- solve(R$varcov,cbind(y, xmat))
  yiRy <- crossprod(y,iRy.x[,1])
  xiRy.x <- crossprod(xmat,iRy.x)
  ##
  ## 1. Computing parameters of posterior for beta
  ##
  if(any(beta.info$iv == Inf)){
    beta.post <- beta.info$beta.fixed
    beta.var.std.post <- 0
    inv.beta.var.std.post <- Inf
  }
  else{
    inv.beta.var.std.post <- drop(beta.info$iv + xiRy.x[,-1])
    beta.var.std.post <- .solve.geoR(inv.beta.var.std.post)
    beta.post <- drop(beta.var.std.post %*% (beta.info$ivm + xiRy.x[,1]))
  }
  ##
  ## 2. Computing parameters of posterior for sigmasq
  ##
  if(sigmasq.info$df.sigmasq == Inf){
    S2.post <- sigmasq.info$sigmasq.fixed
    df.post <- Inf
  }
  else{
    df.post <- n + sigmasq.info$df.sigmasq - beta.info$p
    ##
    if(any(beta.info$iv == Inf)){
      S2.post <- sigmasq.info$n0S0 + yiRy -
        2*crossprod(beta.post, xiRy.x[,1]) +
        crossprod(beta.post, (xiRy.x[,-1] %*% beta.post))
    }
    else
      S2.post <- sigmasq.info$n0S0 + beta.info$mivm + yiRy -
        crossprod(beta.post, (inv.beta.var.std.post %*% beta.post))
    S2.post <- drop(S2.post/df.post)
  }
  ##
  res <- list(beta.post = beta.post,
              beta.var.std.post = beta.var.std.post,
              df.post = df.post, S2.post = S2.post)
  if(dets){
    res$log.det.to.half <- R$log.det.to.half
    res$det.XiRX <- det(xiRy.x[,-1, drop=FALSE])
  }
  ##
  if(do.prediction.moments){
    env.r0 <- new.env()
    assign("r0",cov.spatial(obj = get("d0", envir=env.pred),
                            cov.model = model$cov.model,
                            kappa = model$kappa, cov.pars = c(1, phi)),
           envir=env.r0)
    ## care here, reusing b
    b <- crossprod(get("r0", envir=env.r0),iRy.x)
    riRy <- b[,1, drop=FALSE]
    b <- get("trend.loc", envir=env.pred) -  b[,-1, drop=FALSE]
    ##
    res$pred.mean <- drop(riRy + b %*% beta.post)
    if((tausq.rel < 1e-12) & (!is.null(get("loc.coincide", envir=env.pred))))
      res$pred.mean[get("loc.coincide", envir=env.pred)] <- get("data.coincide", envir=env.pred)
    ##
#    R.riRr.bVb <- 1 - .diagquadraticformXAX(X = get("r0", envir=env.r0),
#                                           lowerA = iR$lower.inverse,
#                                           diagA = iR$diag.inverse)
    R.riRr.bVb <- 1 - colSums(get("r0", envir=env.r0) *
                              solve(R$varcov,get("r0", envir=env.r0)))
    remove("env.r0")
     if(all(beta.info$iv != Inf))
       R.riRr.bVb <- R.riRr.bVb +
         .diagquadraticformXAX(X = t(b),
                              lowerA=beta.var.std.post[lower.tri(beta.var.std.post)],
                              diagA = diag(beta.var.std.post))
    ##
    nug.factor <- ifelse(signal, 0, tausq.rel)
    res$pred.var <- S2.post * (nug.factor + R.riRr.bVb)
    if(((tausq.rel < 1e-12) | signal) & !is.null(get("loc.coincide", envir=env.pred)))
      res$pred.var[get("loc.coincide", envir=env.pred)] <- 0
    res$pred.var[res$pred.var < 1e-16] <- 0
    if(sigmasq.info$df.sigmasq != Inf)
      res$pred.var <- (df.post/(df.post-2)) * res$pred.var
  }
  ##
  if(do.prediction.simulations){
    ## check how to do without inverse!!!
    iR <- varcov.spatial(dists.lowertri = get("data.dist", envir=env.dists),
                         cov.model = model$cov.model,
                         kappa = model$kappa, nugget = tausq.rel,
                         cov.pars = c(1, phi), inv = TRUE,
                         only.inv.lower.diag = TRUE, det = dets)
    res$inv.diag <- iR$diag.inverse
    res$inv.lower <- iR$lower.inverse
  }
  return(res)
}

"sample.prior" <-
  function(n, kb.obj = NULL, prior = prior.control())
{
  call.fn <- match.call()
  ##
  if(!is.null(kb.obj))
    prior <- eval(kb.obj$call$prior)
  ##
  ## Checking for improper priors
  ##
  if(prior$beta.prior == "flat")
    stop("sampling is not possible: improper prior for beta")
  if(any(prior$sigmasq.prior == c("reciprocal", "uniform")))
    stop("sampling is not possible: improper prior for sigmasq")
  ##
  ## preparing output object
  ##
  beta.size <- length(prior$beta)
  if(beta.size == 1)
    beta.name <- "beta"
  else
    beta.name <- paste("beta", (0:(beta.size-1)), sep="")
  simul <- as.data.frame(matrix(0, nrow=n, ncol = beta.size+3))
  names(simul) <- c(beta.name, "sigmasq","phi","tausq.rel")
  ##
  ## Sampling phi and tausq.rel
  ##
  if(prior$phi.prior == "fixed" & prior$tausq.rel.prior == "fixed"){
    simul$phi <- rep(prior$phi, n)
    simul$tausq <- rep(prior$tausq.rel, n)
  }
  else{
    ##
    ## Building the discrete prior distribution
    ##
    phi.discrete <- prior$phi.discrete
    tausq.rel.discrete <- prior$tausq.rel.discrete
    both.discrete <- expand.grid(phi.discrete, tausq.rel.discrete)
    "prob.discrete" <- function(phi.d, tausq.rel.d, prior){
      if(all(prior$phi.prior == "user"))
        pd <- prior$priors.info$phi$probs
      else pd <- phi.d
      if(all(prior$tausq.rel.prior == "user"))
        td <- prior$priors.info$tausq.rel$probs
      else td <- tausq.rel.d
      probs <- switch(prior$phi.prior,
                      user = outer(pd, td, function(x,y){x}),
#                      uniform = outer(pd, td, function(x,y){x-x+1}),
                      uniform = matrix(1, nrow=length(pd), ncol=length(td)),
                      reciprocal = outer(pd, td, function(x,y){ifelse(x>0, 1/x, 0.0)}),
                      squared.reciprocal = outer(pd, td, function(x,y){ifelse(x>0, 1/(x^2), 0.0)}),
                      exponential = outer(pd, td, function(x,y){(1/prior$exponential.par) * exp(x^(1/prior$exponential.par))}),
#                      fixed = outer(pd, td, function(x,y){x-x+1}))
                      fixed = matrix(1, nrow=length(pd), ncol=length(td))
                      )
      if(prior$tausq.rel.prior == "user")
        probs <- t(t(probs) * td)
      if(prior$tausq.rel.prior == "reciprocal"){
        probs <- t(probs) * 1/td
        probs[td == 0] <- 0
        probs <- t(probs)
      }
      return(probs/sum(probs))
    }
    both.discrete$probs <- as.vector(prob.discrete(phi.d = phi.discrete,
#                                                   tausq.discrete = tausq.discrete,
                                                   tausq.rel.d = tausq.rel.discrete,
                                                   prior = prior))
    n.points <- nrow(both.discrete)
    ind <- sample((1:n.points), n, replace = TRUE,
                  prob = both.discrete$probs)
    simul$phi <- both.discrete[ind, 1]
    simul$tausq.rel <- both.discrete[ind, 2]
  }
  ##
  if(prior$sigmasq.prior == "fixed")
    simul$sigmasq <- rep(prior$sigmasq, n)
  else
    simul$sigmasq <- rinvchisq(n, df = prior$df.sigmasq,
                               scale = prior$sigmasq)
  ##
  if(prior$beta.prior == "fixed")
    simul[,1:beta.size] <- matrix(rep(prior$beta, rep(n, beta.size)), ncol=beta.size)
  else{
    if(beta.size == 1){
      simul$beta <- rnorm(n, mean = prior$beta,
                          sd=sqrt(simul$sigmasq * prior$beta.var.std))
    }
    else{
      "sample.beta" <- function(sigmasq, beta, beta.var.std){
        cov.values <- sigmasq * beta.var.std
        cov.svd <- svd(cov.values)
        cov.decomp <- cov.svd$u %*% (t(cov.svd$u) * sqrt(cov.svd$d))
        zsim <- beta + drop(cov.decomp %*% rnorm(length(beta)))
        return(zsim)
      }
      simul[,1:beta.size] <-
        t(sapply(simul$sigmasq, sample.beta, beta = prior$beta,
                 beta.var.std = prior$beta.var.std))
    }
  }
  attr(simul, "Call") <- call.fn
  return(simul) 
}

"sample.posterior" <-
  function(n, kb.obj)
{
  call.fn <- match.call()
  ##
  if(length(class(kb.obj)) == 0)
    stop("kb.obj must be an object with an output of krige.bayes")
  if(any(class(kb.obj) == "krige.bayes")) post <- kb.obj$posterior
  if(any(class(kb.obj) == "posterior.krige.bayes")) post <- kb.obj
  if(all(class(kb.obj) != "krige.bayes") &
     all(class(kb.obj) != "posterior.krige.bayes"))
    stop("kb.obj must be an object with an output of krige.bayes")
  ##
  ## preparing data frame to store the output 
  ##
  if(length(dim(post$beta$pars$mean)) == 2) beta.size <- 1
  else beta.size <- dim(post$beta$pars$mean)[3]
  if(beta.size == 1)
    beta.name <- "beta"
  else
    beta.name <- paste("beta", (0:(beta.size-1)), sep="")
  simul <- as.data.frame(matrix(0, nrow=n, ncol = beta.size+3))
  names(simul) <- c(beta.name, "sigmasq","phi","tausq.rel")
  ##
  ## Sampling phi and tausq.rel
  ##
  if(post$phi$status == "fixed" & post$tausq.rel$status == "fixed"){
    simul$phi <- rep(post$phi$fixed.value, n)
    simul$tausq.rel <- rep(post$tausq.rel$fixed.value, n)
    ind <- 1
    phi.discrete <- post$phi$fixed.value
  }
  else{
    ##
    ## sampling phi and tausq.rel
    ##
    n.points <- length(post$joint.phi.tausq.rel)
    phi.discrete <- post$phi$phi.marginal[,"phi"]
    tausq.rel.discrete <- post$tausq.rel$tausq.rel.marginal[,"tausq.rel"]
    phi.tau.grid <- expand.grid(phi.discrete, tausq.rel.discrete)
    ind <- sample((1:n.points), n, replace = TRUE,
                  prob = as.vector(post$joint.phi.tausq.rel))
    simul$phi <- phi.tau.grid[ind, 1]
    simul$tausq.rel <- phi.tau.grid[ind, 2]
  }
  ##
  ## sampling sigmasq
  ##
  if(post$sigmasq$status == "fixed")
    simul$sigmasq <- rep(post$sigmasq$fixed.value, n)
  else
    simul$sigmasq <- rinvchisq(n, df = post$sigmasq$pars$df,
                               scale = post$sigmasq$pars$S2[ind])
  ##
  ## sampling beta
  ##
  if(post$beta$status == "fixed")
    simul[,1:beta.size] <-
      matrix(rep(post$beta$fixed.value, rep(n, beta.size)), ncol=beta.size)
  else{
    beta.size <- length(post$beta)
    if(beta.size == 1)
      simul$beta <- rnorm(n, mean = post$beta$pars$mean[ind],
                          sd=sqrt(post$beta$pars$var[ind]))
    else{
      if(post$phi$status == "fixed" & post$tausq.rel$status == "fixed"){
        simul[,1:beta.size] <-
          MASS::mvrnorm(n=n, mu = post$beta$pars$mean,
                  Sigma = matrix(post$beta$pars$var, ncol=beta.size))
      }
      else{
        "simula.betavec" <- function(i, nphi){
          nc <- ceiling(i/nphi)
          nr <- i %% nphi
          if(nr == 0) nr <- nphi
          beta.sim <-
            MASS::mvrnorm(n=n, mu = post$beta$pars$mean[nr,nc,],
                    Sigma = matrix(post$beta$pars$var[nr,nc,],
                      ncol=beta.size))
          return(beta.sim)
        }
        simul[,1:beta.size] <- t(sapply(ind, simula.betavec,
                                 nphi = length(phi.discrete)))
      }
    }
  }
  names(simul) <- c(beta.name, c("sigmasq", "phi", "tausq.rel"))
  attr(simul, "Call") <- call.fn
  return(simul) 
}

"statistics.predictive" <-
  function(simuls, mean.var = TRUE, quantile, threshold, sim.means, sim.vars)
{
  results <- list()
  if(missing(quantile) || is.null(quantile)) quantile.estimator <- NULL
  else quantile.estimator <- quantile
  if(missing(threshold) || is.null(threshold)) probability.estimator <- NULL
  else probability.estimator <- threshold
  ##
  if(!is.null(mean.var) && mean.var){
    results$mean.simulations <- drop(rowMeans(simuls))
    results$variance.simulations <- drop(apply(simuls, 1, var))
  }
  if(!is.null(quantile.estimator)) {
    results$quantiles.simulations <-
      drop(apply(simuls, 1, quantile, probs = quantile.estimator))
    if(length(quantile.estimator) > 1) {
      results$quantiles.simulations <-
        as.data.frame(t(results$quantiles))
      names(results$quantiles.simulations) <-
        paste("q", 100 * quantile.estimator, sep = "")
    }
  }
  if(!is.null(probability.estimator)) {
    nsims <- ncol(simuls)
    "prob.cutoff" <- function(x, thres, nsims){
      return(sapply(thres, FUN = function(cut){sum(x <= cut)/nsims}))
    }
    results$probabilities.simulations <-
      drop(apply(simuls, 1, prob.cutoff,
                 thres = probability.estimator, nsims = nsims))
    if(length(threshold) > 1){
      results$probabilities.simulations <-
        as.data.frame(t(results$probabilities.simulations))
      names(results$probabilities.simulations) <-
        paste("threshold", probability.estimator, sep = "")
    }
  }
  if(!missing(sim.means)){
    if(!is.null(sim.means) && sim.means){
      results$sim.means <- drop(colMeans(simuls))
      results$variance.simulations <- drop(apply(simuls, 1, var))
    }
  }
  if(!missing(sim.vars)){
    if(!is.null(sim.vars) && sim.vars){
      results$sim.vars <- drop(apply(simuls, 2, var))
    }
  }
  return(results)
}

"rMVnorm" <-
  function(cov.values, beta.size)
{
  ##
  ## This function produces a sample from  a multivariate normal distribution 
  ## mean is 0 and cov.values is a vector of length beta.size^2
  ##
  cov.values <- matrix(cov.values, ncol = beta.size)
  cov.svd <- svd(cov.values)
  cov.decomp <- cov.svd$u %*% (t(cov.svd$u) * sqrt(cov.svd$d))
  zsim <- as.vector(cov.decomp %*% rnorm(beta.size))
  return(zsim)
}

"plot.krige.bayes" <-
  function(x, phi.dist = TRUE, tausq.rel.dist = TRUE, add = FALSE,
           type = c("bars", "h", "l", "b", "o", "p"), thin, ...)
{
  if(length(class(krige.bayes)) > 0 && all(class(x) != "krige.bayes"))
    stop("object x must be of the class `krige.bayes`")
  if(missing(thin)) thin <- c(1,1)
  if(length(thin) == 1) thin <- rep(thin, 2)
  ##
  type <- match.arg(type)
  ldots <- list(...)
  if(is.null(ldots$col)){
    if(type == "bars") col <- 0:1
    else col <- "black"
  }
  else col <- ldots$col
  if(type != "bars"){
    if(is.null(ldots$lty)) lty <- 1
    else lty <- ldots$lty
    if(is.null(ldots$lwd)) lwd <- 1:2
    else lwd <- ldots$lwd
  }
  ##
  if(phi.dist){
    if(x$prior$phi$status == "fixed")
      cat("parameter `phi` is fixed\n")
    else{
      phi.vals <- x$posterior$phi$phi.marginal[,"phi"]
      phi.off <- 0.1 * diff(phi.vals[1:2])
      ## aqui
      phi.table <- rbind(x$prior$phi$probs, x$posterior$phi$dist)
      colnames(phi.table) <- phi.vals
      ## thining
      nphi <- length(phi.vals)
      ind <- seq(1, nphi, by = thin[1])
      phi.vals <- phi.vals[ind]
      phi.table <- phi.table[,ind]
      if(is.null(ldots$ylim)) phi.ylim <- c(0, 1.1*max(phi.table))
      else phi.ylim <- ldots$ylim
      if(type == "bars")
        barplot(phi.table, legend.text=c("prior", "posterior"),
                beside=TRUE, col=col, ylim = phi.ylim, 
                xlab = expression(phi), ylab = "density")
      else{
        if(type=="h")
          phi.vals <- cbind(phi.vals - phi.off, 
                            phi.vals + phi.off)
        matplot(phi.vals, t(phi.table), type = type,
                lwd = lwd, lty = lty, ylim = phi.ylim, 
                col = col, xlab = expression(phi),
                ylab = "density", add = add)
      }
    }
  }
  if(tausq.rel.dist){
    if(x$prior$tausq.rel$status == "fixed")
      cat("parameter `tausq.rel` is fixed\n")
    else{
      tausq.rel.vals <- x$posterior$tausq.rel$tausq.rel.marginal[,"tausq.rel"]
      tausq.rel.off <- 0.1 * diff(tausq.rel.vals[1:2])
      tausq.rel.table <- rbind(x$prior$tausq.rel$probs,
                               x$posterior$tausq.rel$dist)
      colnames(tausq.rel.table) <-  tausq.rel.vals
      ## thining
      ntausq.rel <- length(tausq.rel.vals)
      ind <- seq(1, ntausq.rel, by = thin[2])
      tausq.rel.vals <- tausq.rel.vals[ind]
      tausq.rel.table <- tausq.rel.table[,ind]
      if(is.null(ldots$ylim)) tau.ylim <- c(0, 1.1*max(tausq.rel.table))
      else tau.ylim <- ldots$ylim
      if(type == "bars")
        barplot(tausq.rel.table, legend.text=c("prior", "posterior"),
                beside=TRUE, col=col, ylim = tau.ylim, 
                xlab = expression(tau[rel]^2), ylab = "density")
      else{
        if(type=="h")
          tausq.rel.vals <- cbind(tausq.rel.vals - tausq.rel.off,
                                  tausq.rel.vals + tausq.rel.off)
        matplot(tausq.rel.vals, t(tausq.rel.table), type = type, lwd = lwd, lty = lty,
                col = col, ylim = tau.ylim, 
                xlab = expression(tau[rel]^2), ylab = "density", add = add)
      }
    }
  }
  return(invisible())
}

##"lines.posterior.krige.bayes" <-
##  function(x, parameter = c("beta", "sigmasq", "phi", "tausq.rel"), ...)#
##{
##  if(parameter == "beta"){
##    attach(x$posterior$beta, pos=1)
##    
##  return(invisible())#
##}

  
"print.betavar" <-
  function(x, ...)
{
  size <- attr(x, "Size")
  x <- matrix(x, size, size)
  betavar <- matrix(NA, size, size)
  betavar[row(betavar) >= col(betavar)] <- x[row(betavar) >= col(betavar)]
  if(size > 1){
    labels <- paste("beta", 0:(size-1), sep="")
    dimnames(betavar) <- list(labels, labels)
  }
  print(betavar, na="")
  return(invisible(x))
}


"hist.krige.bayes" <-
  function(x, pars, density.est = TRUE,
           histogram = TRUE, ...)
{
  Ldots <- list(...)
  if(missing(pars) | (!missing(pars) && pars == -1)){
    ppars <- names(x$posterior$sample)
    np <- length(ppars) 
    if(x$prior$tausq.rel$status != "random") ppars <- ppars[-np]
    if(x$prior$phi$status != "random") ppars <- ppars[-(np-1)]
    if(x$prior$sigmasq$status != "random") ppars <- ppars[-(np-2)]
    if((!missing(pars) && pars == -1)) ppars <- ppars[-1]
    pars <- ppars
  }
  res <- list(histogram = list(), density.estimation = list())
  for(ipar in pars){
    if(substr(ipar, 1,4) == "beta")
      if(nchar(ipar) == 4) xl <- expression(beta)
      else xl <- substitute(beta[N], list(N=as.numeric(strsplit(ipar, "beta")[[1]][2])))
    if(ipar == "sigmasq") xl <- expression(sigma^2)
    if(ipar == "phi") xl <- expression(phi)
    if(ipar == "tausq.rel") xl <- expression(tau[rel]^2)
    y <- as.vector(x$posterior$sample[[ipar]])
    ymax <- 0
    if(histogram){
      res$histogram[[ipar]] <- plH <- hist(y, plot = FALSE, ...)
      ymax <- max(plH$dens)
    }
    if(density.est){
      if(is.null(Ldots$width)){
       if(requireNamespace("MASS", quietly=TRUE))
          plD <- density(y,width=MASS::bandwidth.nrd(y), ...)
        else plD <- density(y, ...)
      }
      else
        plD <- lines(density(y, width = Ldots$width))
  res$density.estimation[[ipar]] <- plD 
      ymax <- max(c(ymax, plD$y))
    }
    if(histogram){
      plot(plH, ylim =c(0, ymax), freq = FALSE,
           xlab=xl, main= Ldots$main, ...)
      if(density.est) lines(plD)
    }
    else
      if(density.est) plot(plD, xlab=xl, ...)
  }
  return(invisible(res))
}

"lines.variomodel.krige.bayes" <- 
  function(x, summary.posterior, max.dist, uvec,
           posterior = c("variogram", "parameters"),  ...)
{
  my.l <- list()
  ##
  ## Setting the maximum distance to compute the variogram
  ##
  if(missing(max.dist)){
    my.l$max.dist <- x$max.dist
    if (is.null(my.l$max.dist) | mode(my.l$max.dist) != "numeric") 
      stop("a numerical value must be provided to the argument max.dist")
  }
  else my.l$max.dist <- max.dist
  ##
  ## picking the variogram model
  ##
 # if(is.null(x$call$cov.model))
 #   my.l$cov.model <- "exponential"
 # else {
    my.l$cov.model <- x$model$cov.model
    my.l$kappa <- x$model$kappa
 #   if(any(x$model$cov.model == c("matern", "powered.exponential",
 #            "cauchy", "gencauchy", "gneiting.matern")))
 #     my.l$kappa <- x$call$kappa
 #   else my.l$kappa <- NULL
 # }
  ##
  posterior <- match.arg(posterior)
  if(is.function(summary.posterior)) spost <- post.fc <- summary.posterior
  else spost <- match.arg(summary.posterior, choices = c("mode", "median", "mean"))
  ##
  if(!is.null(x$posterior$sample) & posterior == "variogram"){
    if(!is.function(spost))
      stop("summary.posterior must be a function when posterior = `variogram`")
    if(missing(uvec)) my.l$uvec <- seq(0, my.l$max.dist, l=51)
    calc.vario <- function(x, info = my.l){
      return((x[1] * (1 + x[3])) -
             cov.spatial(info$uvec, cov.model = my.l$cov.model, kappa = my.l$kappa, cov.pars = x[1:2]))
    }
    post.vario <- apply(x$posterior$sample[c("sigmasq","phi","tausq.rel")], 1, calc.vario)
    gamma.post <- drop(apply(post.vario, 1, post.fc))
    if(is.vector(gamma.post))
      lines(my.l$uvec, gamma.post, ...)
    else
      matplot(my.l$uvec, t(gamma.post), add = TRUE, ...)
  }
  else{
    if(is.function(spost))
      stop("summary.posterior must be one of `mean`, `median` or `mode` when posterior = `parameters`")
    if(spost == "mode")
      spost1 <- "mode.cond"
    else spost1 <- spost
    my.l$cov.pars <- c(x$posterior$sigmasq$summary[spost1],
                       x$posterior$phi$summary[spost])
    names(my.l$cov.pars) <- NULL
    if(mode(x$posterior$tausq.rel$summary) == "numeric")
      nugget <- x$posterior$tausq.rel$summary[spost] * my.l$cov.pars[1]
    else nugget <- 0
    names(nugget) <- NULL
    my.l$sill.total <- nugget + my.l$cov.pars[1]
    gamma.f <- function(x, my.l)
      {
        return(my.l$sill.total -
               cov.spatial(x, cov.model = my.l$cov.model, kappa = my.l$kappa,
                           cov.pars = my.l$cov.pars))
      }
    curve(gamma.f(x,my.l=my.l), from = 0, to = my.l$max.dist, add=TRUE, ...)
  }
  return(invisible())
}

"print.krige.bayes" <-
  function(x, ...)
{
  print.default(x, ...)
}

"print.posterior.krige.bayes" <-
  function(x, ...)
{
  print.default(x, ...)
}



