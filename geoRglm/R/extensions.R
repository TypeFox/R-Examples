".krige.bayes.extnd" <- 
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
    geodata <- list(coords=coords, data=data)
  call.fc <- match.call()
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  do.prediction <- ifelse(all(locations == "no"), FALSE, TRUE)
  base.env <- sys.frame(sys.nframe())
  message.prediction <- character()
  ##
  ## reading model input
  ##
  if(missing(model))
    model <- model.control()
  else{
    if(class(model) != "model.geoR"){
      if(!is.list(model))
        stop(".krige.bayes.extnd: the argument model only takes a list or an output of the function model.control")
      else{
        model.names <- c("trend.d", "trend.l", "cov.model", "kappa", "aniso.pars", "lambda") 
        model <- .object.match.names(model,model.names)
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
  ##
  ## reading prior input
  ##
  if(missing(prior))
    prior <- prior.control()
  else{
    if(class(prior) != "prior.geoR"){
      if(!is.list(prior))
        stop(".krige.bayes.extnd: the argument prior only takes a list or an output of the function prior.control")
      else{
        prior.names <- c("beta.prior", "beta", "beta.var.std", "sigmasq.prior",
                         "sigmasq", "df.sigmasq", "phi.prior", "phi", "phi.discrete",
                         "tausq.rel.prior", "tausq.rel", "tausq.rel.discrete") 
        prior <- .object.match.names(prior,prior.names)
        ## DO NOT CHANGE ORDER OF THE NEXT 3 LINES
        if(is.null(prior$beta)) prior$beta <-  NULL
        if(is.null(prior$beta.prior)) prior$beta.prior <-  c("flat", "normal", "fixed")
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
          prior$tausq.rel.prior <- c("fixed", "uniform")
        if(is.null(prior$tausq.rel.discrete)) prior$tausq.rel.discrete <- NULL 
        prior <- prior.control(beta.prior = prior$beta.prior,
                               beta = prior$beta, beta.var.std = prior$beta.var.std,
                               sigmasq.prior = prior$sigmasq.prior,
                               sigmasq = prior$sigmasq,  df.sigmasq = prior$df.sigmasq,
                               phi.prior = prior$phi.prior,
                               phi = prior$phi, phi.discrete = prior$phi.discrete, 
                               tausq.rel.prior = prior$tausq.rel.prior,
                               tausq.rel = prior$tausq.rel,
                               tausq.rel.discrete = prior$tausq.rel.discrete)
      } 
    }
  }
  ##
  kb <- list(posterior = list(beta=list(), sigmasq=list(), phi=list(), tausq.rel=list()),
             predictive=list(mean = NULL, variance = NULL, distribution = NULL),
             prior = prior$priors.info, model = model) 
  ##class(kb$posterior) <- "krige.bayes.posterior"
  ##class(kb$predictive) <- "krige.bayes.predictive"
  ##class(kb$prior) <- "krige.bayes.prior"
  pred.env <- new.env()
  ##
  beta <- prior$beta
  if(prior$beta.prior == "fixed") beta.fixed <- beta
  if(prior$beta.prior == "normal"){
    beta.var.std <- prior$beta.var.std
    inv.beta.var.std <- .solve.geoR(beta.var.std)
    betares <- list(iv = inv.beta.var.std, ivm = drop(.solve.geoR(beta.var.std, beta)),
                    mivm = drop(crossprod(beta, .solve.geoR(beta.var.std, beta))))
  }
  if(prior$sigmasq.prior == "fixed") sigmasq.fixed <- prior$sigmasq
  else S2.prior <- prior$sigmasq
  df.sigmasq <- prior$df.sigmasq
  ##
  if(prior$phi.prior != "fixed") stop(".krige.bayes.extnd: only phi fixed is allowed.")
  ##
  tausq.rel.fixed <- tausq.rel <- prior$tausq.rel
  if(prior$tausq.rel.prior != "fixed") stop(".krige.bayes.extnd: only tausq fixed is allowed.")
  ##
  ## checking data configuration
  ##
  data <- as.matrix(data)
  n.datasets <- ncol(data)
  n <- nrow(data)
  if(n.datasets == 1) stop(".krige.bayes.extnd: this function is for multiple datasets. Use krige.bayes instead.")
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning(".krige.bayes.extnd: vector of coordinates: assuming one spatial dimension (transect)")
  }
  coords <- as.matrix(coords)
  dists.env <- new.env()
  assign("data.dist", as.vector(dist(coords)), envir=dists.env)
  data.dist.range <- range(get("data.dist", envir=dists.env))
  data.dist.min <- data.dist.range[1]
  data.dist.max <- data.dist.range[2]
  if(1e12*data.dist.min < 0.5)
    stop(".krige.bayes.extnd: this function does not allow two data at same location")
  ##
  ## reading output options
  ##
  if(missing(output))
    output <- output.control()
  else{
    if(class(output) != "output.geoR"){
      if(!is.list(output))
        stop(".krige.bayes.extnd: the argument output only takes a list or an output of the function output.control")
      else{
        output.names <- c("n.posterior","n.predictive","moments","n.back.moments","simulations.predictive",
                          "mean.var","quantile","threshold","signal","messages.screen")
        output <- .object.match.names(output,output.names)
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
        if(is.null(output$signal)) output$signal <- NULL
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior,
                                 n.predictive = output$n.predictive,
                                 moments = output$moments,
                                 n.back.moments = output$n.back.moments, 
                                 simulations.predictive = output$simulations.predictive,
                                 mean.var = output$mean.var, quantile = output$quantile,
                                 threshold = output$threshold, signal = output$signal,
                                 messages = output$messages.screen)
      }
    }
  }
  n.posterior <- output$n.posterior
  messages.screen <- output$messages.screen
  if(do.prediction){
    n.predictive <- as.integer(output$n.predictive)
    if(is.null(n.predictive)) n.predictive <- as.integer(0)
    simulations.predictive <- output$simulations.predictive
    if(is.null(simulations.predictive)) simulations.predictive <- FALSE
    if(!is.null(output$signal) && output$signal) stop(".krige.bayes.extnd: prediction of the signal is not implemented")
    if(!is.null(output$probability.estimator)) stop(".krige.bayes.extnd: probability.estimator not implemented\n")
    if(!is.null(output$quantile.estimator)) stop(".krige.bayes.extnd: quantile.estimator not implemented\n")
    if(simulations.predictive && n.predictive == 0) n.predictive <- as.integer(1)
  }
  ##
  ## Box-Cox transformation
  ##
  if(abs(model$lambda-1)>0.0001) stop(".krige.bayes.extnd: Box-Cox transformation not allowed \n")    
  ##
  ## Building trend (covariates/design) matrices:   
  ##
  dimnames(coords) <- list(NULL, NULL)
  if(nrow(coords) != n) stop(".krige.bayes.extnd: number of data is different of number of data locations (coordinates)")
  trend.data <- unclass(trend.spatial(trend=model$trend.d, geodata = geodata))
  beta.size <- ncol(trend.data)
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  if((prior$beta.prior == "normal") && (beta.size != length(beta)) )
    stop(".krige.bayes.extnd: size of beta incompatible with the trend model (covariates)")
  ##
  if(do.prediction) {
    locations <- .geoR.check.locations(locations)
    if(!is.null(borders)){
      nloc0 <- nrow(locations)
      ind.loc0  <- .geoR_inout(locations, borders)
      locations <- locations[ind.loc0,,drop=TRUE]
      if(nrow(locations) == 0){
        warning("\n .krige.bayes.extnd: no prediction to be performed.\n             There are no prediction locations inside the borders")
        do.prediction <- FALSE
      }
      if(messages.screen)
        cat(".krige.bayes.extnd: results will be returned only for prediction locations inside the borders\n")
    }
    ##
    ## Checking trend specification
    ##
    if(inherits(model$trend.d, "formula") | inherits(model$trend.l, "formula")){
      if((inherits(model$trend.d, "formula") == FALSE) | (inherits(model$trend.l, "formula") == FALSE))
        stop(".krige.bayes.extnd: model$trend.d and model$trend.l must have similar specification\n")
    }
    else{
      if((class(model$trend.d) == "trend.spatial") & (class(model$trend.l) == "trend.spatial")){
        if(ncol(model$trend.d) != ncol(model$trend.l))
          stop(".krige.bayes.extnd: trend.d and trend.l do not have the same number of columns")
      }
      else
        if(model$trend.d != model$trend.l)
          stop(".krige.bayes.extnd: especification of model$trend.l and model$trend.d must be similar")
    }
    ##
    if(messages.screen){
      cat(switch(model$trend.d,
                 "cte" = ".krige.bayes.extnd: model with mean being constant",
                 "1st" = ".krige.bayes.extnd: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = ".krige.bayes.extnd: model with mean given by a 2nd order polynomial on the coordinates",
                 ".krige.bayes.extnd: model with mean defined by covariates provided by the user"))
      cat("\n")
    }
    ##
    dimnames(locations) <- list(NULL, NULL)
    if(class(model$trend.l) == "trend.spatial")
      assign("trend.loc", unclass(model$trend.l), envir=pred.env)
    else
      assign("trend.loc", unclass(trend.spatial(trend=model$trend.l, geodata = list(coords = locations))), envir=pred.env)
    ni <- nrow(get("trend.loc", envir=pred.env))
    if(!is.null(borders) && ni == nloc0)
      assign("trend.loc", get("trend.loc", envir=pred.env)[ind.loc0,,drop=TRUE], envir=pred.env)
    if(nrow(locations) != ni)
      stop("trend.l is not compatible with number of prediction locations") 
  }
  ##
  ## Anisotropy correction
  ##   (warning: this must be placed here, AFTER trend matrices be defined)
  ##
  if(!is.null(model$aniso.pars)) {
    if((abs(model$aniso.pars[1]) > 0.001) & (abs(model$aniso.pars[2] - 1) > 0.001)){
      if(messages.screen) cat(".krige.bayes.extnd: anisotropy parameters provided and assumed to be constants\n")
      coords <- coords.aniso(coords = coords, aniso.pars = model$aniso.pars)
      if(do.prediction) locations <- coords.aniso(coords = locations, aniso.pars = model$aniso.pars)
      remove("dists.env")
      dists.env <- new.env()
      assign("data.dist", as.vector(dist(coords)), envir=dists.env)
    }
  }
  ##
  ## Distances between data and prediction locations
  ## Must be here AFTER anisotropy be checked
  ##
  if(do.prediction){
    assign("d0", loccoords(coords = coords, locations = locations), envir=pred.env)
    ##
    ## checking coincident data and prediction locations
    ##
    loc.coincide <- apply(get("d0", envir=pred.env), 2, function(x){any(x < 1e-10)})
    if(any(loc.coincide))
      loc.coincide <- which(loc.coincide)
    else
      loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data){return(data[x < 1e-10,])}
      data.coincide <- t(apply(get("d0", envir=pred.env)[,loc.coincide, drop=FALSE],
                               2,temp.f, data=data))
    }
    else data.coincide <- NULL
    n.loc.coincide <- length(loc.coincide)    
  }
  ##
  ## Preparing prior information on beta and sigmasq
  ##
  beta.info <-
    switch(prior$beta.prior,
           fixed = list(mivm = 0, ivm = 0, iv = Inf, beta.fixed = beta.fixed, p = 0),
           flat = list(mivm = 0, ivm = 0, iv = 0, p = beta.size),
           normal = list(mivm = betares$mivm, ivm = betares$ivm, iv = betares$iv, p= 0))
  sigmasq.info <-
    switch(prior$sigmasq.prior,
           fixed = list(df.sigmasq = Inf, n0S0 = 0, sigmasq.fixed = sigmasq.fixed),
           reciprocal = list(df.sigmasq = 0, n0S0 = 0),
           uniform = list(df.sigmasq = -2, n0S0 = 0),
           sc.inv.chisq = list(df.sigmasq = df.sigmasq, n0S0 = df.sigmasq*S2.prior))
  ##
  ## ====================== PART 2 =============================
  ##                 FIXED PHI AND TAUSQ.REL
  ## ===========================================================
  ##
  phi.fixed <- prior$phi
  ##
  ## Computing parameters of the posterior for $\(\beta, \sigma^2)$ 
  ## and variables to be used for prediction (if applies)
  ##
  ind.datasets <- seq(length=n.datasets)
  R <- varcov.spatial(dists.lowertri = get("data.dist", envir=dists.env), cov.model = model$cov.model,
                      kappa = model$kappa, nugget = tausq.rel.fixed, cov.pars = c(1, phi.fixed))
  iRy.x <- solve(R$varcov, cbind(data, trend.data))
  yiRy <-  colSums(data * iRy.x[,ind.datasets])
  xiRy.x <- crossprod(trend.data, iRy.x)
  ##
  if(!is.matrix(xiRy.x)) xiRy.x <- matrix(xiRy.x, 1, n.datasets+beta.size)
  ##
  ## 1. Computing parameters of posterior for beta
  ##
  if(any(beta.info$iv == Inf)){ 
    beta.post <- beta.info$beta.fixed
    beta.var.std.post <- matrix(0, ncol = beta.size, nrow = beta.size)
    inv.beta.var.std.post <- Inf
  }
  else{
    inv.beta.var.std.post <- as.matrix(beta.info$iv + xiRy.x[,-ind.datasets, drop = FALSE])    
    beta.var.std.post <- .solve.geoR(inv.beta.var.std.post)
    beta.post <- beta.var.std.post %*% (beta.info$ivm + xiRy.x[,ind.datasets])
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
      S2.post <- sigmasq.info$n0S0 + yiRy - 2*crossprod(beta.post, xiRy.x[,ind.datasets, drop = FALSE]) + as.vector(crossprod(beta.post,xiRy.x[,-ind.datasets, drop = FALSE])%*%beta.post)
    }
    else{
      S2.post <- sigmasq.info$n0S0 + beta.info$mivm + yiRy - .diagquadraticformXAX(beta.post, inv.beta.var.std.post[lower.tri(inv.beta.var.std.post)], diag(inv.beta.var.std.post))     
    }
    S2.post <- drop(S2.post/df.post)
  }
  remove("yiRy","xiRy.x")
  ##
  ## Preparing output of the posterior distribution
  ##
  if(prior$beta.prior == "fixed")
    kb$posterior$beta <- list(status = "fixed", fixed.value = beta.fixed)
  else {
    if(prior$sigmasq.prior == "fixed")
      kb$posterior$beta <- list(distribution = "normal")
    else
      kb$posterior$beta <- list(distribution = "t", conditional = "normal")
    kb$posterior$beta$pars <- list(mean = beta.post,
                                   var = beta.var.std.post%o%S2.post)
  }
  if(prior$sigmasq.prior == "fixed")
    kb$posterior$sigmasq <- list(status="fixed", fixed.value=sigmasq.fixed)
  else
    kb$posterior$sigmasq <- list(distribution = "sc.inv.chisq", pars = list(df = df.post, S2 = S2.post))
  kb$posterior$phi<- list(status= "fixed", fixed.value = phi.fixed)
  kb$posterior$tausq.rel <- list(status= "fixed", fixed.value = tausq.rel.fixed)
  ##
  ## Preparing output of the predictive distribution
  ##
  if(do.prediction){
    v0 <- cov.spatial(obj = get("d0", envir=pred.env),
                      cov.model = cov.model, kappa = kappa,
                      cov.pars = c(1, phi.fixed))
    ## care here, reusing object b
    b <- crossprod(iRy.x,v0)
    remove("iRy.x")
    tv0ivdata <- t(b[ind.datasets, , drop=FALSE])
    b <- t(get("trend.loc", envir=pred.env)) - b[-ind.datasets, , drop=FALSE]
    ##
    if(any(beta.info$iv == Inf))
      kb$predictive$mean <- tv0ivdata + as.vector(crossprod(b, beta.post))
    else kb$predictive$mean <- tv0ivdata + crossprod(b, beta.post)
    if((tausq.rel.fixed < 1e-12) & (!is.null(loc.coincide))){
      kb$predictive$mean[loc.coincide,] <- data.coincide
    }
    R.riRr.bVb <- 1 - colSums(v0 * solve(R$varcov,v0))
    if(all(beta.info$iv != Inf))
      R.riRr.bVb <- R.riRr.bVb +
        .diagquadraticformXAX(X = t(b),
                              lowerA=beta.var.std.post[lower.tri(beta.var.std.post)],
                              diagA = diag(beta.var.std.post))
    ##
    kb$predictive$variance <- as.vector(tausq.rel.fixed + R.riRr.bVb)%*%t(S2.post)
    if(((tausq.rel.fixed < 1e-12) ) & !is.null(loc.coincide))
      kb$predictive$variance[loc.coincide] <- 0
    kb$predictive$variance[kb$predictive$variance < 1e-12] <- 0
    if(sigmasq.info$df.sigmasq != Inf)
      kb$predictive$variance <- (df.post/(df.post-2))*kb$predictive$variance
    kb$predictive$distribution <- ifelse(prior$sigmasq.prior == "fixed", "normal", "t")
    remove("R.riRr.bVb")
  }
  ##
  ## ======================= PART 4 ==============================
  ##                Sampling from the predictive
  ## =============================================================
  ##
  if(do.prediction && simulations.predictive){
    if(is.R()){
      if(cov.model.number > 12)
        stop("simulation in .krige.bayes.extnd not implemented for the choice of correlation function")
    }
    else
      if(cov.model.number > 10)
        stop("simulation in .krige.bayes.extnd not implemented for the choice of correlation function")
    if(messages.screen){
      cat(".krige.bayes.extnd: sampling from the predictive\n")
    }
    tmean <- kb$predictive$mean
    tv0ivdata <- NULL        ### se efter om den kan fjernes foer
    Dval <-  1.0 + tausq.rel
    coincide.cond <- any(loc.coincide)
    if(coincide.cond){
      nloc <- ni - n.loc.coincide
      ind.not.coincide <- (-loc.coincide)
      v0 <- v0[, ind.not.coincide, drop=FALSE]
      tmean <- tmean[ind.not.coincide, , drop=FALSE]
      b <- b[,ind.not.coincide, drop=FALSE]
    }
    else{
      nloc <- ni
      ind.not.coincide <- TRUE
    }
    if(n.predictive > 1){
      warning("n.predictive > 1 is not implemented, n.predictive = 1")
      n.predictive <- 1
    }
    kb$predictive$simulations <- matrix(NA, nrow=ni, ncol=n.datasets)
    if(nloc>0){
      ## next step will be do the following without using the inverse
      iR <- varcov.spatial(dists.lowertri = get("data.dist",
                             envir=dists.env),
                           cov.model = model$cov.model,
                           kappa = model$kappa, nugget = tausq.rel.fixed,
                           cov.pars = c(1, phi.fixed), inv = TRUE,
                           only.inv.lower.diag = TRUE)
      kb$predictive$simulations[ind.not.coincide,] <-
        geoR::.cond.sim(env.loc = base.env, env.iter = base.env,
                  loc.coincide = loc.coincide,
                  coincide.cond = coincide.cond, tmean = tmean,
                  Rinv = list(lower=iR$lower.inverse, diag=iR$diag.inverse),
                  mod = list(beta.size = beta.size, nloc = nloc,
                    Nsims = n.datasets, n = n,
                    Dval = Dval, df.model = df.post,
                    s2 = S2.post,
                    cov.model.number = cov.model.number,
                    phi = phi.fixed, kappa = kappa),
                  vbetai = beta.var.std.post,
                  fixed.sigmasq = (sigmasq.info$df.sigmasq == Inf))
    }
    if(coincide.cond)
      kb$predictive$simulations[loc.coincide,] <- data.coincide
  }
  if(!do.prediction) kb$predictive <- "no prediction locations provided"
  kb$.Random.seed <- seed
  kb$max.dist <- data.dist.max
  kb$call <- call.fc
  attr(kb, "prediction.locations") <- call.fc$locations
  return(kb)
}


".krige.conv.extnd" <- 
function(geodata, coords = geodata$coords, data = geodata$data,
         locations, borders, krige, output)
{
##############################################################################
  ##     An extended version of krige.conv (geoR) allowing for multivariate data 
  ##     (useful for performing kriging on simulations: output from MCMC)
##############################################################################
  ##     data   : n*m-matrix of values of observations (m datasets)
  ##     n.predictive  : only 0 and 1 is allowed
##############################################################################
  base.env <- sys.frame(sys.nframe())
  if(missing(geodata))
    geodata <- list(coords=coords, data=data)
  cl <- match.call()
  data <- as.matrix(data)
  n.datasets <- ncol(data)
  n <- nrow(data)
  ##
  ## reading input
  ##
  if(missing(krige))
    krige <- krige.control()
  else{
    if(class(krige) != "krige.geoR"){
      if(!is.list(krige))
        stop(".krige.conv.extnd: the argument krige only takes a list or an output of the function krige.control")
      else{
        krige.names <-c("type.krige","trend.d","trend.l","obj.model","beta","cov.model",
"cov.pars","kappa","nugget","micro.scale","dist.epsilon","lambda","aniso.pars")
        krige <- .object.match.names(krige,krige.names)
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
                               trend.d = krige$trend.d, trend.l = krige$trend.l,
                               obj.model = krige$obj.model,
                               beta = krige$beta, cov.model = krige$cov.model,
                               cov.pars = krige$cov.pars, kappa = krige$kappa,
                               nugget = krige$nugget, micro.scale = krige$micro.scale,
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
  ##
  ## reading output options
  ##
  if(missing(output))
    output <- output.control()
  else{
    if(class(output) != "output.geoR"){
      if(!is.list(output))
        stop(".krige.conv.extnd: the argument output only takes a list or an output of the function output.control")
      else{
        output.names <- c("n.posterior","n.predictive","moments","n.back.moments","simulations.predictive",
                          "mean.var","quantile","threshold","signal","messages.screen")
        output <- .object.match.names(output,output.names)
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
        if(is.null(output$signal)) output$signal <- NULL
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior,
                                 n.predictive = output$n.predictive,
                                 moments = output$moments,
                                 n.back.moments = output$n.back.moments, 
                                 simulations.predictive = output$simulations.predictive,
                                 mean.var = output$mean.var, quantile = output$quantile,
                                 threshold = output$threshold, signal = output$signal,
                                 messages = output$messages.screen)
      }
    }
  }
  ##
  signal <- ifelse(is.null(output$signal), FALSE, output$signal)
  messages.screen <- output$messages.screen
  n.predictive <- output$n.predictive
  n.back.moments <- output$n.back.moments
  ##
  ## checking input
  ##
  if(krige$type.krige == "ok") beta.prior <- "flat"
  else beta.prior <- "deg"
  if(n.predictive > 1) {
    warning("n.predictive redefined: n.predictive=1")
    n.predictive <- 1
  }
  if(n.datasets == 1)
    stop("Only one dataset; please use krige.conv instead")
  ##
  if(is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  locations <- .geoR.check.locations(locations)
  ## selecting locations inside the borders 
  ## and also values of trend.l if the case
  ##
  if(!is.null(borders)){
    nloc0 <- nrow(locations)
    ind.loc0  <- .geoR_inout(locations, borders)
    locations <- locations[ind.loc0,,drop=FALSE]
    if(nrow(locations) == 0)
      stop("\n .krige.conv.extnd : there are no prediction locations inside the borders")
    if(messages.screen)
      cat("krige.conv: results will be returned only for prediction locations inside the borders\n")
  }
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  ##
  if(messages.screen){
    if(is.numeric(krige$trend.d))
      cat(".krige.conv.extnd: model with covariates matrix provided by the user")
    else
      cat(switch(as.character(krige$trend.d)[1],
                 "cte" = ".krige.conv.extnd: model with mean being constant",
                 "1st" = ".krige.conv.extnd: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = ".krige.conv.extnd: model with mean given by a 2nd order polynomial on the coordinates",
                 ".krige.conv.extnd: model with mean defined by covariates provided by the user"))
    cat("\n")
  }
  trend.data <- unclass(trend.spatial(trend=krige$trend.d, geodata = geodata))
  beta.size <- ncol(trend.data)
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  trend.l <- unclass(trend.spatial(trend=krige$trend.l, geodata = list(coords = locations)))
  if(!is.null(borders) && nrow(trend.l) == nloc0) trend.l <- trend.l[ind.loc0,,drop=FALSE]
  if (nrow(trend.l) != nrow(locations)) stop("locations and trend.l have incompatible sizes")
  ni <- nrow(trend.l)
  ##
  ## Anisotropy correction (should be placed AFTER trend.d/trend.l
  ##
  if(!is.null(krige$aniso.pars)){
    if(messages.screen) cat(".krige.conv.extnd: anisotropy correction performed\n")
    coords <- coords.aniso(coords = coords, aniso.pars = krige$aniso.pars)
    locations <- coords.aniso(coords = locations, aniso.pars = krige$aniso.pars)
  }
  ##
  ## Box-Cox transformation
  ##
  if(lambda != 1) {
    if(messages.screen) cat(".krige.conv.extnd: Data transformation (Box-Cox) performed.\n")
    if(lambda == 0)
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  ##
  ## setting covariance parameters
  ##
  tausq <- nugget
  if(is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
  }
  else stop("covariance parameter should be given as a vector")
  cpars <- c(1, phi)
  tausq.rel <- tausq/sigmasq
  ## using nugget interpreted as microscale variation or measurement error
  nug.factor <- ifelse(signal, krige$micro.scale/sigmasq, tausq.rel)
  ##
  ## starting kriging calculations
  ##
  kc.result <- list()
  ##
  Vcov <- varcov.spatial(coords = coords, cov.model = cov.model,
                         kappa = kappa, nugget = tausq.rel,
                         cov.pars = cpars)$varcov
  ivtt <- solve(Vcov,trend.data)
  ttivtt <- crossprod(ivtt,trend.data)
  if(beta.prior == "flat") beta.flat <- solve(ttivtt,crossprod(ivtt, data))
  remove("ivtt")
  ## PJ: From a numerical perspective there quite a few ugly bits
  ##here. See R-help e-mails from Bates and Ripley about
  ##LS/GLS. However, changing is not simple.
  temp <- NULL
  v0 <- loccoords(coords = coords, locations = locations)
  if(n.predictive > 0){
    ## checking if there are data points coincident with prediction locations
    loc.coincide <- apply(v0, 2, function(x, min.dist){any(x < min.dist)}, min.dist=krige$dist.epsilon)
    if(any(loc.coincide)) loc.coincide <- (1:ni)[loc.coincide]
    else loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data, dist.eps){return(data[x < dist.eps,, drop=FALSE])}
      data.coincide <- apply(v0[, loc.coincide, drop=FALSE], 2, temp.f, data=data, dist.eps=krige$dist.epsilon)
    }
    else data.coincide <- NULL
  }
  ind.v0 <- which(v0<krige$dist.epsilon)
  v0 <- cov.spatial(obj = v0, cov.model = cov.model, kappa = kappa, cov.pars = cpars)
  v0[ind.v0] <- 1+nug.factor

  
  ivv0 <- solve(Vcov,v0)
  tv0ivv0 <- colSums(v0 * ivv0)
  b <- t(trend.l - crossprod(ivv0, trend.data))
  tv0ivdata <- crossprod(ivv0,data)
  if(n.predictive == 0) {
    remove(list = "v0")
    gc(verbose = FALSE)
  }
  if(beta.prior == "deg") {
    kc.result$predict <- array(tv0ivdata + rep(drop(crossprod(b,beta)), n.datasets), dim = c(ni, n.datasets))
    if(n.predictive == 0) remove(list = c("b", "tv0ivdata"))
    else remove("tv0ivdata")
    kc.result$krige.var <- sigmasq*drop(1+nug.factor - tv0ivv0)
    beta.est <- paste(".krige.conv.extnd: Simple kriging (beta provided by user)\n")
  }
  if(beta.prior == "flat") {
    kc.result$predict <- array(tv0ivdata, dim = c(ni, n.datasets)) + crossprod(b,beta.flat)
    remove(list = c("tv0ivdata"))
    bitb <- colSums(b * solve(ttivtt,b))
    kc.result$krige.var <- sigmasq*drop(1+nug.factor - tv0ivv0 + bitb)
    kc.result$beta.est <- beta.flat
    remove("beta.flat")
  }
  kc.result$krige.var[kc.result$krige.var < 1e-12] <- 0
  if(any(kc.result$krige.var < 0)) warning(".krige.conv.extnd: negative kriging variance found! Investigate why this is happening.\n")
  out.message <- ".krige.conv.extnd: Kriging performed using global neighbourhood"
  if(messages.screen) cat(paste(out.message, "\n"))
############## Sampling from the resulting distribution #####################
  if(n.predictive > 0) {
    if(messages.screen) cat(".krige.conv.extnd: sampling from the predictive distribution (conditional simulations)\n")    
    Dval <- 1. + nug.factor
    if(beta.prior == "deg") vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
    else vbetai <- matrix(.solve.geoR(ttivtt), ncol = beta.size, nrow = beta.size)
    coincide.cond <- (((round(1e12 * nugget) == 0) | !signal) & (!is.null(loc.coincide)))
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
    kc.result$simulations <- matrix(0, nrow = ni, ncol = n.datasets)
    if(nloc>0){
      ## re-write this without using inverse!!!
      invcov <- varcov.spatial(coords = coords, cov.model = cov.model,
                               kappa = kappa, nugget = tausq.rel,
                               cov.pars = cpars, 
                               inv = TRUE, only.inv.lower.diag = TRUE)
      kc.result$simulations[ind.not.coincide,  ] <-
        geoR::.cond.sim(env.loc = base.env, env.iter = base.env,
                  loc.coincide = loc.coincide,
                  coincide.cond = coincide.cond,
                  tmean = kc.result$predict[ind.not.coincide, , drop = FALSE],
                  Rinv = invcov,
                  mod = list(beta.size = beta.size, nloc = nloc,
                    Nsims = n.datasets, n = n, Dval = Dval,
                    df.model = NULL, s2 = sigmasq,
                    cov.model.number = .cor.number(cov.model),
                    phi = phi, kappa = kappa),
                  vbetai = vbetai, fixed.sigmasq = TRUE)
    }
    if(coincide.cond)
      kc.result$simulations[loc.coincide,  ] <- kc.result$predict[loc.coincide,, drop = FALSE]         
    remove(list = c("v0", "invcov", "b"))
  }
######################	Back-transforming predictions ############
  if(lambda != 1) {
    if(lambda == 0){
      predict.transf <- kc.result$predict
      if(messages.screen) cat(".krige.conv.extnd: back-transforming the predictions using formula for EXP() \n")
      kc.result$predict <- exp(predict.transf + 0.5 * kc.result$krige.var)
      kc.result$krige.var <- (exp(2 * predict.transf + kc.result$krige.var))*expm1(kc.result$krige.var)
      remove("predict.transf")
    }
    if(lambda > 0){
      ## using second order taylor-expansion + facts for N(0,1) [third moment = 0 ; fourth moment = 12].
      if(messages.screen) cat(".krige.conv.extnd: back-transforming predictions using 2. order Taylor expansion for g^{-1}() \n")
      ivBC <- .BC.inv(kc.result$predict,lambda)
      kc.result$predict <- ivBC + 0.5 * ((1-lambda)*ivBC^(1-2*lambda))*kc.result$krige.var
      kc.result$krige.var <- (ivBC^(1-lambda))^2*kc.result$krige.var + (11/4)*((1-lambda)*ivBC^(1-2*lambda))^2*kc.result$krige.var^2
      remove("ivBC")
    }
    if(lambda < 0){
      if(messages.screen) cat(".krige.conv.extnd: resulting distribution has no mean for lambda < 0 - back transformation not performed\n")
    }
    if(n.predictive > 0) {
      kc.result$simulations <- .BC.inv(kc.result$simulations,lambda)
    }
  }
  else{
    temp <- kc.result$krige.var
    kc.result$krige.var <- matrix(0, nrow = ni, ncol = n.datasets)
    kc.result$krige.var[] <- temp
  }
  kc.result <- c(kc.result, list(message = out.message, call = cl))
#######################################
  return(kc.result)
}
