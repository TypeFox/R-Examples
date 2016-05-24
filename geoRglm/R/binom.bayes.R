

".mcmc.bayes.binom.logit" <- 
  function(data, units.m, trend, mcmc.input, messages.screen, cov.model, kappa, tausq.rel, coords, ss.sigma, df, phi.prior, phi.discrete)
{
  ##
  ## This is the MCMC engine for the Bayesian analysis of a spatial binomial logit Normal model
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, log(data), -1) - ifelse(units.m-data > 0, log(units.m-data), -1) )
  }
  else{
    if(!any(mcmc.input$S.start=="random")){
      if(is.numeric(mcmc.input$S.start)){
        if(length(mcmc.input$S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
        S <- as.vector(mcmc.input$S.start)
      }
      else  stop(" S.start must be a vector of same dimension as data ")
    }
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1 && pnorm((phi.discrete[nmphi] - phi.discrete[1])/(nmphi - 1), sd = sqrt(phi.scale)) > 0.975)
      warning("Consider making the grid in phi.discrete more dense. The algorithm may have problems moving.")
  }
  messages.C <- ifelse(messages.screen,1,0) ## for C-code
  ##
  ##                                                                      
  ## ---------- sampling ----------- ###### 
  cov.model.number <- .cor.number(cov.model)
  beta.size <- if(is.vector(trend)) 1 else ncol(trend)
  n.sim <- floor(n.iter/thin)
  ## remember this rather odd coding for telling that S.start is from the prior !!!
  if(any(mcmc.input$S.start=="random")) Sdata <- as.double(as.vector(c(rep(0, n.sim*n - 1),1)))
  else Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  result <-  .C("mcmcrun4binom",
                as.integer(n),
                as.double(data),
                as.double(units.m),
                as.double(as.vector(t(trend))),
                as.integer(beta.size),
                as.integer(cov.model.number),
                as.double(kappa),
                as.double(tausq.rel),
                as.double(coords[,1]),
		as.double(coords[,2]),
                as.double(S.scale),
                as.double(phi.scale),
                as.integer(n.iter),
                as.integer(thin),
                as.integer(burn.in),
                as.integer(messages.C),
                as.double(ss.sigma),
                as.integer(df),
                as.double(phi.prior),
                as.double(phi.discrete),
                as.integer(nmphi),
                Sdata = Sdata,
                phi.sample = as.double(rep(phi, n.sim)),
		acc.rate = rep(0,floor(n.iter/1000)+1), 
		acc.rate.phi = rep(0,floor(n.iter/1000)+1), PACKAGE = "geoRglm")[c("Sdata", "phi.sample","acc.rate","acc.rate.phi" )]
  attr(result$Sdata, "dim") <- c(n, n.sim)
  if(nmphi>1) result$acc.rate <- cbind(burn.in + seq(0,floor(n.iter/1000))*1000, result$acc.rate,result$acc.rate.phi)
  else result$acc.rate <- cbind(burn.in + seq(0,floor(n.iter/1000))*1000, result$acc.rate)
  result$acc.rate.phi <- NULL
  if(burn.in==0) result$acc.rate <- result$acc.rate[-1,,drop=FALSE]
  if(nmphi>1) colnames(result$acc.rate) <- c("iter.numb", "Acc.rate", "Acc.rate.phi")
  else colnames(result$acc.rate) <- c("iter.numb", "Acc.rate")
  if(messages.screen) cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n")) 
  return(result)
}

".mcmc.bayes.conj.binom.logit" <- 
  function(data, units.m, meanS, ttvbetatt, mcmc.input, messages.screen, cov.model, kappa, tausq.rel, coords, ss.sigma, df, phi.prior, phi.discrete)
{
  ##
  ## This is the MCMC engine for the Bayesian analysis (with normal prior for beta) of a spatial binomial logit Normal model
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, log(data), -1) - ifelse(units.m-data > 0, log(units.m-data), -1) ) - meanS
  }
  else{
    if(!any(mcmc.input$S.start=="random")){
      if(is.numeric(mcmc.input$S.start)){
        if(length(mcmc.input$S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
        S <- as.vector(mcmc.input$S.start)
      }
      else  stop(" S.start must be a vector of same dimension as data ")
    }
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1 && pnorm((phi.discrete[nmphi] - phi.discrete[1])/(nmphi - 1), sd = sqrt(phi.scale)) > 0.975)
      warning("Consider making the grid in phi.discrete more dense. The algorithm may have problems moving.")
  }
  messages.C <- ifelse(messages.screen,1,0) ## for C-code
  ##
  ## ---------- sampling ----------- ###### 
  cov.model.number <- .cor.number(cov.model)
  n.sim <- floor(n.iter/thin)
  ## remeber this rather odd coding for telling that S.start is from the prior !!!
  if(any(mcmc.input$S.start=="random")) Sdata <- as.double(as.vector(c(rep(0, n.sim*n - 1),1)))
  else Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  result <-  .C("mcmcrun5binom",
                as.integer(n),
                as.double(data),
                as.double(units.m),
                as.double(as.vector(meanS)),
                as.double(as.vector(ttvbetatt)),
                as.integer(cov.model.number),
                as.double(kappa),
                as.double(tausq.rel),
                as.double(coords[,1]),
		as.double(coords[,2]),
                as.double(S.scale),
                as.double(phi.scale),
                as.integer(n.iter),
                as.integer(thin),
                as.integer(burn.in),
                as.integer(messages.C),
                as.double(ss.sigma),
                as.integer(df),
                as.double(phi.prior),
                as.double(phi.discrete),
                as.integer(nmphi),
                Sdata = Sdata,
                phi.sample= as.double(rep(phi, n.sim)),
		acc.rate = rep(0,floor(n.iter/1000)+1), 
		acc.rate.phi = rep(0,floor(n.iter/1000)+1), PACKAGE = "geoRglm")[c("Sdata", "phi.sample","acc.rate","acc.rate.phi" )]
  attr(result$Sdata, "dim") <- c(n, n.sim)
  if(nmphi>1) result$acc.rate <- cbind(burn.in + seq(0,floor(n.iter/1000))*1000, result$acc.rate,result$acc.rate.phi)
  else result$acc.rate <- cbind(burn.in + seq(0,floor(n.iter/1000))*1000, result$acc.rate)
  result$acc.rate.phi <- NULL
  if(burn.in==0) result$acc.rate <- result$acc.rate[-1,,drop=FALSE]
  if(nmphi>1) colnames(result$acc.rate) <- c("iter.numb", "Acc.rate", "Acc.rate.phi")
  else colnames(result$acc.rate) <- c("iter.numb", "Acc.rate")
  if(messages.screen) cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  return(result)
}


"binom.krige.bayes" <- 
  function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", locations = "no", borders, model, prior, mcmc.input, output){
###########
  if(missing(geodata))
    geodata <- list(coords=coords, data=data, units.m=units.m)
  if(missing(borders))
    borders <- geodata$borders
  call.fc <- match.call()
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  do.prediction <- ifelse(all(locations == "no"), FALSE, TRUE)
  ##
  ## Checking data configuration
  ##
  if(is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  n <- length(data)
  if(nrow(coords) != n) stop("number of data is different from number of data locations (coordinates)")
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n)
  }
  if(any(units.m <= 0)) stop("units.m must be positive")
  ##
  ## reading input
  ##
  if(missing(model)) model <- model.glm.control()
  else model <- .model.glm.check.aux(model, fct = "binom.krige.bayes")
  cov.model <- model$cov.model
  kappa <- model$kappa
  tausq.rel <- prior$tausq.rel
  ## reading prior input
  ##
  if(missing(prior))  stop("binom.krige.bayes: argument prior must be given ")
  else prior <- .prior.glm.check.aux(prior, fct = "binom.krige.bayes")
  beta.prior <- prior$beta.prior
  beta <- prior$beta
  beta.var <- prior$beta.var.std
  sigmasq.prior <- prior$sigmasq.prior
  if(sigmasq.prior == "fixed") sigmasq <- prior$sigmasq
  else{
    df.sigmasq <- prior$df.sigmasq
    S2.prior <- prior$sigmasq
  }
  phi.prior <- prior$phi.prior 
  phi <- prior$phi
  if(phi.prior != "fixed") phi.discrete <- prior$phi.discrete
  else phi.discrete <- phi
  ##
  ## reading output options
  ##
  if(missing(output)) output <- output.glm.control()
  else output <- .output.glm.check.aux(output, fct = "binom.krige.bayes")
  quantile.estimator <- output$quantile.estimator
  probability.estimator <- output$probability.estimator
  inference <- output$inference
  messages.screen <- output$messages.screen
  ## check == here
  data.dist <- as.vector(dist(coords))
  if(1000000000000. * min(data.dist) == 0) stop("Two coords are identical; not allowed.")
  ##
  trend.d <- model$trend.d
  if(messages.screen) {
    cat(switch(as.character(trend.d)[1],
                 "cte" = "binom.krige.bayes: model with mean being constant",
                 "1st" = "binom.krige.bayes: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = "binom.krige.bayes: model with mean given by a 2nd order polynomial on the coordinates",
                 "binom.krige.bayes: model with mean defined by covariates provided by the user"))
    cat("\n")
  }
  trend.data <- unclass(trend.spatial(trend=trend.d, geodata = geodata))
  dimnames(coords) <- list(NULL, NULL)
  dimnames(trend.data) <- list(NULL, NULL)
  beta.size <- ncol(trend.data)
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  if(beta.prior == "normal" |  beta.prior == "fixed"){
    if(beta.size != length(beta))
      stop("binom.krige.bayes: size of beta incompatible with the trend model (covariates)")
  }
  if(sigmasq.prior == "uniform"){
    if(max(units.m) == 1) warning("This choice of sigmasq.prior gives an improper posterior !!!!!!! \n")
    if(sum(ifelse(units.m>1,1,0)) < beta.size + 3) warning("This choice of sigmasq.prior may give an improper posterior !!!!!!! \n")
  }
  aniso.pars <- model$aniso.par
  if(!is.null(aniso.pars)) coords.transf <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
  else coords.transf <- coords
  ##
  # checking prediction locations
  ##
  if((inference) & (do.prediction)){
    locations <- .geoR.check.locations(locations)
    ## Checking the consistency between coords, locations, and trends
    trend.l <- model$trend.l 
    ## Checking for 1D prediction 
    if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
      krige1d <- TRUE
    else krige1d <- FALSE
    ##
    if(is.null(trend.l)) stop("trend.l needed for prediction")
    if(inherits(trend.d, "formula") | inherits(trend.l, "formula")){
      if((!inherits(trend.d, "formula")) | (!inherits(trend.l, "formula")))
        stop("trend.d and trend.l must have similar specification\n")
    }
    else{
      if((class(trend.d)=="trend.spatial") & (class(trend.l)=="trend.spatial")){
        if(ncol(trend.d) != ncol(trend.l))
          stop("trend.d and trend.l do not have the same number of columns")
      }
      else if(trend.d != trend.l) stop("trend.l is different from trend.d")
    }
    ## 
    if(nrow(unclass(trend.spatial(trend=model$trend.l, geodata = list(coords = locations)))) != nrow(locations))
      stop("binom.krige.bayes: number of points to be estimated is different of the number of trend locations")
    if(!is.null(borders)){
      ind.loc0  <- .geoR_inout(locations, borders)
      if(!any(ind.loc0)){
        warning("\n binom.krige.bayes: no prediction to be performed.\n             There are no prediction locations inside the borders")
        do.prediction <- FALSE
      }
    }
    kb.results <- list(posterior = list(), predictive = list())
  }
  else {
    if(do.prediction & messages.screen) cat(paste("need to specify inference=TRUE to make predictions \n"))
    kb.results <- list(posterior = list(), predictive = paste("prediction not performed"))
    do.prediction <- FALSE
  }
  ##
  ## ##### preparing for MCMC -------------------------------------------------------
  ##
  if(missing(mcmc.input)) stop("binom.krige.bayes: argument mcmc.input must be given")
  mcmc.input <- .mcmc.check.aux(mcmc.input, fct="binom.krige.bayes")
  ##
  if(beta.prior == "fixed" | beta.prior == "normal") mean.d <- as.vector(trend.data %*% beta)
  else mean.d <- rep(0, n)
  if(sigmasq.prior != "fixed"){
    if(beta.prior == "flat") df.model <- n - beta.size + df.sigmasq
    else df.model <- n + df.sigmasq
  }
  else df.model <- Inf
  if(beta.prior == "normal"){
    if(beta.size > 1) ttvbetatt <- trend.data%*%beta.var%*%t(trend.data)
    else ttvbetatt <- crossprod(t(trend.data))*beta.var
  } else ttvbetatt <- matrix(0, beta.size, beta.size)
  if(sigmasq.prior == "fixed"){     ### implies that phi is fixed !
    invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = tausq.rel*sigmasq,
                               cov.pars = c(sigmasq,phi), inv = TRUE, func.inv = "cholesky",
                               try.another.decomposition = FALSE)$inverse
    if(beta.prior != "fixed"){
      ivtt <- invcov%*%trend.data
      if(beta.prior == "normal") invcov <- invcov-ivtt%*%.solve.geoR(crossprod(trend.data, ivtt) + solve(beta.var), t(ivtt))
      else invcov <- invcov-ivtt%*%.solve.geoR(crossprod(trend.data, ivtt), t(ivtt))
    }
  }
  if((phi.prior == "fixed") && (sigmasq.prior != "fixed")){
    phi.prior.prob <- 1
    phi.discrete <- phi
  }
  else phi.prior.prob <-  prior$priors.info$phi$probs
  ##
############----------PART 2 ------------##############################
############-----------MCMC -------------##############################
  ##
  if(sigmasq.prior == "fixed"){ 
    log.odds <- .mcmc.binom.logit(data = data, units.m = units.m, meanS = mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen)
  }
  else {
    kb.results$posterior$phi <- list()
    ## take care re-using log.odds !
    if(beta.prior == "flat"){
      log.odds <- .mcmc.bayes.binom.logit(data=data, units.m=units.m, trend=trend.data, mcmc.input=mcmc.input, messages.screen=messages.screen, cov.model=cov.model, 
                                         kappa=kappa, tausq.rel = tausq.rel, coords=coords.transf, 
                                         ss.sigma = df.sigmasq*S2.prior, df = df.model, phi.prior = phi.prior.prob,
                                         phi.discrete = phi.discrete)
    }
    else{     
      log.odds <- .mcmc.bayes.conj.binom.logit(data=data, units.m=units.m, meanS = mean.d, ttvbetatt = ttvbetatt, mcmc.input=mcmc.input, messages.screen=messages.screen,
                                              cov.model=cov.model, kappa=kappa, tausq.rel = tausq.rel,
                                              coords=coords.transf, ss.sigma = df.sigmasq*S2.prior, df = df.model,
                                              phi.prior = phi.prior.prob, phi.discrete = phi.discrete)
    }
    kb.results$posterior$phi$sample <- log.odds$phi.sample
  }
  kb.results$posterior$acc.rate  <- log.odds$acc.rate
  log.odds <- log.odds$Sdata
  ##
##############-------------PART 3----------######################
##############------------prediction-------######################
  ##
  n.sim <- ncol(log.odds)
  if(inference) {
    if(phi.prior=="fixed") phi.posterior <- list(phi.prior=phi.prior, phi=phi)
    else  phi.posterior <- list(phi.prior=phi.prior, phi.discrete=phi.discrete, sample=kb.results$posterior$phi$sample)
    predict.temp <- .pred.aux(S=log.odds, coords=coords, locations=locations, borders=borders, model=model, prior=prior, output=output, phi.posterior=phi.posterior, link="logit")
    temp.post <- predict.temp$temp.post
    if(do.prediction) {
      temp.pred <- predict.temp$temp.pred
      kb.results$predictive$simulations <- predict.temp$pred.simulations
    }
    ##
    if(do.prediction) {
      ##
      if(!is.null(borders)){
        nloc0 <- nrow(locations)
        ind.loc0  <- .geoR_inout(locations, borders)
        locations <- locations[ind.loc0,]
      }
      ni <- nrow(locations)  
      d0mat <- loccoords(coords, locations)
      loc.coincide <- (colSums(d0mat < 1e-10) == 1)
      if((is.logical(quantile.estimator) && (quantile.estimator)) || (is.numeric(quantile.estimator))){
        predi.q <- .pred.quan.aux(temp.pred, loc.coincide, df.model, ni, quantile.estimator)
        kb.results$predictive$median <- plogis(predi.q$median)
        kb.results$predictive$uncertainty <- (plogis(predi.q$upper) - plogis(predi.q$lower))/4      
        if(is.data.frame(predi.q$quantiles)){
          names.q <- names(predi.q$quantiles)
          kb.results$predictive$quantiles <- as.data.frame(plogis(as.matrix(predi.q$quantiles)))
          names(kb.results$predictive$quantiles) <- names.q
        }
        else kb.results$predictive$quantiles <- plogis(predi.q$quantiles)
      }
      ##
      ## ------ probability estimators
      ##
      if(!is.null(probability.estimator)) {
        logit.probab <- ifelse(probability.estimator < 1, log(probability.estimator) - log(1-probability.estimator), 1e+17)
        logit.probab <- ifelse(probability.estimator > 0, logit.probab, 1e-17)
        len.p <- length(probability.estimator)
        if(len.p== 1){
          kb.results$predictive$probability <- round(.pmixed(logit.probab, temp.pred, df.model), digits = 3)
        }
        else{
          kb.results$predictive$probability <- matrix(NA,ni,len.p)
          for(ii in seq(length=len.p)){
            kb.results$predictive$probability[,ii] <- round(.pmixed(logit.probab[ii], temp.pred, df.model), digits = 3)
          }
        }
      }
      remove("temp.pred")
      ## 
      if(messages.screen) cat("binom.krige.bayes: Prediction performed \n")
    }
    else {
      kb.results$predictive <- "no locations to perform prediction were provided"
      if(messages.screen) cat("Only Bayesian estimation of model parameters \n")
    }
    ##
    ##----- calculating posterior summaries ----------------##
    ##
    if(beta.prior == "fixed") kb.results$posterior$beta <- paste("provided by user: ", beta)
    else {
      kb.results$posterior$beta <- list()
      kb.results$posterior$beta$mean <- rowMeans(temp.post$beta.mean)
      names(kb.results$posterior$beta$mean) <- beta.names
      kb.results$posterior$beta$var <- rowMeans(temp.post$beta.var, dims = 2) + var(t(temp.post$beta.mean))
      dimnames(kb.results$posterior$beta$var) <- list(beta.names,beta.names)
    }
    if(sigmasq.prior == "fixed") kb.results$posterior$sigmasq <- paste("provided by user: ", sigmasq) 
    else{
      kb.results$posterior$sigmasq <- list()
      kb.results$posterior$sigmasq$mean <- mean(temp.post$S2)*df.model/(df.model-2)
      kb.results$posterior$sigmasq$var <- (mean(temp.post$S2)*2/(df.model-4) + var(temp.post$S2))*df.model^2/(df.model-2)^2
    }
    if(phi.prior == "fixed") kb.results$posterior$phi <- paste("provided by user: ", phi) 
    else{
      kb.results$posterior$phi$mean <- mean(kb.results$posterior$phi$sample)
      kb.results$posterior$phi$var <- var(kb.results$posterior$phi$sample)
    }
    ##
    ## Simulations from the posterior of parameters.
    ##
    if(output$sim.posterior){
      if(beta.size == 1) {
        if(sigmasq.prior == "fixed") {
          if(beta.prior != "fixed")
            kb.results$posterior$beta$sample <- rnorm(n.sim) * as.vector(sqrt(temp.post$beta.var)) + as.vector(temp.post$beta.mean)
        }
        else{
          kb.results$posterior$sigmasq$sample <- rinvchisq(n.sim, df.model, temp.post$S2)
          if(beta.prior != "fixed"){
            cond.beta.sd <- sqrt((as.vector(temp.post$beta.var) * kb.results$posterior$sigmasq$sample)/temp.post$S2)
            kb.results$posterior$beta$sample <- rnorm(n.sim) * cond.beta.sd + as.vector(temp.post$beta.mean)
          }
        }
      }
      else {
        if(sigmasq.prior == "fixed"){
          if(beta.prior != "fixed")
            kb.results$posterior$beta$sample <- array(apply(temp.post$beta.var,3,.multgauss),dim=c(beta.size, n.sim))+temp.post$beta.mean
        }
        else {
          kb.results$posterior$sigmasq$sample <- rinvchisq(n.sim, df.model, temp.post$S2)
          if(beta.prior != "fixed"){
            cond.beta.var <- temp.post$beta.var *rep(kb.results$posterior$sigmasq$sample/temp.post$S2,rep(beta.size^2,n.sim))
            kb.results$posterior$beta$sample <- array(apply(cond.beta.var,3,.multgauss),dim=c(beta.size, n.sim)) + temp.post$beta.mean
          }
        }
      }
    }
    remove("temp.post")
  }
  if(output$keep.mcmc.sim) kb.results$posterior$simulations <- plogis(log.odds)  
  model$lambda <- NULL
  kb.results$model <- model
  kb.results$prior <- prior$priors.info
  kb.results$mcmc.input <- mcmc.input
  kb.results$.Random.seed <- seed
  kb.results$call <- call.fc
  attr(kb.results, "prediction.locations") <- call.fc$locations
  if(do.prediction) attr(kb.results, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  if(!is.null(call.fc$borders)) attr(kb.results, "borders") <- call.fc$borders
  class(kb.results) <- "glm.krige.bayes"
  return(kb.results)
}

