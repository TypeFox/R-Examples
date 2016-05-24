
"prepare.likfit.glsm" <- function(mcmc.output, use.intensity = FALSE)
{
### this is for the glsm.mcmc() function
  ##
  if(class(mcmc.output) != "glsm.mcmc") stop("mcmc.output must be an object of class ``glsm.mcmc'' ")
  if(sum(mcmc.output$acc.rate[,2])==0) stop("mcmc.output must not have all accept-rates equal to zero ")
  n.dat <- nrow(mcmc.output$simulations)
  n.sim <- ncol(mcmc.output$simulations)
  if(use.intensity){
    if(mcmc.output$model$family != "poisson") stop("use.intensity = TRUE is only allowed for the Poisson error distribution")
    if(any(mcmc.output$geodata$data == 0)) stop("use.intensity = TRUE is only allowed when all data are positive ")
  }
  lambda <- mcmc.output$model$lambda
  S <- mcmc.output$simulations
  trend.data <- unclass(trend.spatial(trend = mcmc.output$model$trend, geodata=mcmc.output$geodata))
  beta.size <- ifelse(is.matrix(trend.data), ncol(trend.data), 1)
  cov.model <- mcmc.output$model$cov.model
  kappa <- mcmc.output$model$kappa
  beta <- mcmc.output$model$beta
  sigmasq <- mcmc.output$model$cov.pars[1]
  phi <- mcmc.output$model$cov.pars[2]
  nugget.rel <- mcmc.output$model$nugget/sigmasq
  if(is.null(mcmc.output$model$aniso.pars)) coords <- mcmc.output$geodata$coords
  else coords <- coords.aniso(coords = mcmc.output$geodata$coords, aniso.pars = mcmc.output$model$aniso.pars)
  invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget.rel, 
                           cov.pars = c(1,phi), det = TRUE, only.inv.lower.diag = TRUE)
  SivS <- .diagquadraticformXAX(S, invcov$lower.inverse, invcov$diag.inverse) 
  if(beta.size == 1){
    DivD <- as.vector(.bilinearformXAY(trend.data,invcov$lower.inverse,invcov$diag.inverse,trend.data)) 
    SivD <- as.vector(.bilinearformXAY(S, invcov$lower.inverse, invcov$diag.inverse, trend.data))
    log.f.sim <- -invcov$log.det.to.half - 0.5*n.dat*log(sigmasq) - 0.5*(SivS-2*SivD*beta+DivD*beta^2)/sigmasq
  }
  else{
    DivD <- .bilinearformXAY(trend.data,invcov$lower.inverse,invcov$diag.inverse,trend.data)
    SivD <- .bilinearformXAY(S, invcov$lower.inverse, invcov$diag.inverse, trend.data)
    log.f.sim <-  -invcov$log.det.to.half - 0.5*n.dat*log(sigmasq) - 0.5*as.vector(SivS-2*as.vector(SivD%*%beta)+t(beta)%*%DivD%*%beta)/sigmasq
  }
  if(use.intensity){
    if(lambda==0){
      logJ <- -apply(S,2,sum)
      mu <- exp(S)
    }
    else {
      logJ <- apply(log(S*lambda+1),2,sum)*(lambda-1)/lambda
      mu <- (S*lambda+1)^(1/lambda)
    }
    return(list(family=mcmc.output$model$family, link=mcmc.output$model$link, mu = mu, geodata = mcmc.output$geodata, trend = mcmc.output$model$trend,
                aniso.pars = mcmc.output$model$aniso.pars, lambda = lambda, log.f.sim = log.f.sim + logJ))
  }
  else{
    return(list(family=mcmc.output$model$family, link=mcmc.output$model$link, S = S, geodata = mcmc.output$geodata, trend = mcmc.output$model$trend,
                aniso.pars = mcmc.output$model$aniso.pars, lambda = lambda, log.f.sim = log.f.sim))
  }
}

".func.val" <- function(SivS, SivD, DivD, beta, sigmasq, log.f)
{
  beta.size <- length(beta)
  if(beta.size == 1) ff <- exp(-0.5*(SivS-2*SivD*beta+DivD*beta^2)/sigmasq-log.f)
  else ff <- exp(-0.5*as.vector(SivS-2*as.vector(SivD%*%beta)+t(beta)%*%DivD%*%beta)/sigmasq-log.f)
  if(any(!is.finite(ff))){
    print(summary(ff))
    stop("Some function values are not finite")
  }
  return(ff)
}

".NewtonRhapson.step" <-
  function(SivS, SivD, DivD, SivDi2, ff, n.dat, beta, sigmasq, steplen)
{
  beta.size <- length(beta)
  n.sim <- length(SivS)
  n.dat <- length(SivS)
  if(any(!is.finite(ff))){
    print(summary(ff))
    stop("Some function values are not finite")
  }
  meanff <-  mean(ff)
  ## removed 1/(sigmasq^(n.dat/2) everywhere
  if(beta.size == 1){
    F1 <- mean((SivD-DivD*beta)*ff)/sigmasq
    SivDbeta <- SivD*beta
    betaDivDbeta <- DivD*beta^2
    F11 <- (mean(SivDi2*ff)-(beta*DivD)^2*meanff)/(sigmasq^2)-2*beta*DivD*F1/sigmasq-meanff*DivD/sigmasq
    F12 <- -(n.dat/2+1)*F1/sigmasq+0.5*mean((SivD-DivD*beta)*(SivS-2*SivDbeta+betaDivDbeta)*ff)/(sigmasq^3)
  }
  else{
    F1 <- (colMeans(SivD*ff)-(DivD%*%beta)*meanff)/sigmasq
    SivDbeta <- as.vector(SivD%*%beta)
    betaDivDbeta <- as.vector(t(beta)%*%DivD%*%beta)
    F11 <- colMeans(SivDi2*ff)/(sigmasq^2)-(DivD%*%beta)%*%t(DivD%*%beta)*meanff/(sigmasq^2)-((DivD%*%beta)%*%t(F1)+F1%*%t(DivD%*%beta))/sigmasq-meanff*DivD/sigmasq
    F12 <- -(n.dat/2+1)*F1/sigmasq+0.5*(colMeans(SivD*(SivS-2*SivDbeta+betaDivDbeta)*ff)-(DivD%*%beta)*mean((SivS-2*SivDbeta+betaDivDbeta)*ff))/(sigmasq^3)
  }
  F2 <- -(n.dat/2)*meanff/sigmasq + 0.5*mean((SivS-2*SivDbeta+betaDivDbeta)*ff)/(sigmasq^2)
  F22 <- (n.dat/2)*(n.dat/2+1)*meanff/(sigmasq^2)-0.5*(n.dat+2)*mean((SivS-2*SivDbeta+betaDivDbeta)*ff)/(sigmasq^3)+0.25*mean((SivS-2*SivDbeta+betaDivDbeta)^2*ff)/(sigmasq^4)
  Delta2 <- rbind(cbind(F11,F12),c(t(F12),F22))
  if(det(Delta2) != 0){
    stepvec <- solve(Delta2,c(F1,F2))
    betanew <- beta - steplen*stepvec[seq(length=beta.size)]
    sigmasqnew <- sigmasq - steplen*stepvec[beta.size+1]
    flat.message <- FALSE
  }
  else{
    flat.message <- TRUE
    betanew <- beta 
    sigmasqnew <- sigmasq 
  } 
  return(list(betanew=betanew, sigmasqnew=sigmasqnew, flat = flat.message))
}

".maxim.aux1" <-
  function(S, invcov, trend, log.f.sim, messages.screen=FALSE)
{
  n.sim <- ncol(S)
  n.dat <- nrow(S)  
  if(is.matrix(trend)) beta.size <- ncol(trend)
  else beta.size <- 1
  SivS <- .diagquadraticformXAX(S,invcov$lower.inverse,invcov$diag.inverse)
  if(beta.size == 1){
    SivD <- as.vector(.bilinearformXAY(S,invcov$lower.inverse,invcov$diag.inverse,trend))
    DivD <- as.vector(.bilinearformXAY(trend,invcov$lower.inverse,invcov$diag.inverse,trend))
    beta.hat <- SivD/DivD
    sigmasq.hat <- .diagquadraticformXAX(S-trend%*%t(beta.hat),invcov$lower.inverse,invcov$diag.inverse)/n.dat
    SivDi2 <- SivD^2
    corr1 <- mean(-0.5*(SivS-2*SivD*mean(beta.hat)+DivD*mean(beta.hat)^2)/mean(sigmasq.hat)-log.f.sim)
  }
  else{
    SivD <- .bilinearformXAY(S,invcov$lower.inverse,invcov$diag.inverse,trend)
    DivD <- .bilinearformXAY(trend,invcov$lower.inverse,invcov$diag.inverse,trend)
    beta.hat <- t(solve(DivD,t(SivD))) ####### might be improved (GLS)
    sigmasq.hat <- .diagquadraticformXAX(S-trend%*%t(beta.hat),invcov$lower.inverse,invcov$diag.inverse)/n.dat
    "cp" <- function(x){return(x%*%t(x))}
    SivDi2 <- array(t(apply(SivD,1,cp)),dim=c(n.sim,beta.size,beta.size))
    corr1 <- mean(-0.5*(SivS-2*as.vector(SivD%*%colMeans(beta.hat))+t(colMeans(beta.hat))%*%DivD%*%colMeans(beta.hat))/mean(sigmasq.hat)-log.f.sim)
  }
  log.f.c <- log.f.sim + corr1
  log.hh <- (-Inf)
  for(ll in seq(length=ncol(S))){
    if(beta.size == 1) ff <- .func.val(SivS, SivD, DivD, beta.hat[ll], sigmasq.hat[ll], log.f.c)
    else ff <- .func.val(SivS, SivD, DivD, beta.hat[ll,], sigmasq.hat[ll], log.f.c)
    log.hhnew <- log(mean(ff))-(n.dat/2)*log(sigmasq.hat[ll])
    if(log.hhnew>log.hh){
      if(beta.size == 1) beta <-beta.hat[ll]
      else beta <- beta.hat[ll,]
      sigmasq <- sigmasq.hat[ll]
      log.hh <- log.hhnew 
    }
  }
  ## new correction (again for numerical purposes).
  if(beta.size == 1) corr2 <- max(-0.5*(SivS-2*SivD*beta+DivD*beta^2)/sigmasq-log.f.sim)
  else corr2 <- max(-0.5*(SivS-2*as.vector(SivD%*%beta)+t(beta)%*%DivD%*%beta)/sigmasq-log.f.sim)
  log.f.c <- log.f.sim + corr2
  ff <- .func.val(SivS, SivD, DivD, beta, sigmasq, log.f.c)
  log.hh <- log(mean(ff))-(n.dat/2)*log(sigmasq)
  ##
  test <- 1
  test2 <- 1
  steplen <- 1
  while(test>0.0000000000001 | test2 > 0 ){
    New <- .NewtonRhapson.step(SivS, SivD, DivD, SivDi2, ff, n.dat, beta, sigmasq, steplen)
    ffnew <- .func.val(SivS, SivD, DivD, New$beta, New$sigmasq, log.f.c)
    log.hhnew <- log(mean(ffnew))-(n.dat/2)*log(New$sigmasq)
    if(New$flat & messages.screen){
      cat(paste("Problems when optimising w.r.t. beta and sigmasq: likelihood is very flat \n"))
      cat(paste("likelihood value at this stage is = ",log.hhnew+corr2,"\n"))
    }
    test <- sum((beta-New$beta)^2)+(sigmasq-New$sigmasq)^2
    test2 <- log.hhnew-log.hh
    if(test2>0){
      log.hh <- log.hhnew   
      beta <- New$beta 
      sigmasq <- New$sigmasq
      ff <- ffnew  
    }
    else{
      steplen <- steplen/2 
    }
  }
  return(list(beta = beta, sigmasq = sigmasq, logh = log.hh+corr2)) 
}


".lik.sim" <- function(pars, fp, ip, temp.list)
{ 
  ## Obligatory parameter:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  messages.screen <- ifelse(is.null(temp.list$messages.screen), TRUE,temp.list$messages.screen)
  if(messages.screen) cat(paste("phi = ",phi, "tausq.rel = ",tausq.rel,"\n"))
  ##
  ## Computing likelihood
  ##
  iv <- varcov.spatial(dists.lowertri = as.vector(dist(temp.list$coords)), cov.model = temp.list$cov.model, kappa = temp.list$kappa,
                       nugget = tausq.rel, cov.pars = c(1, phi), only.inv.lower.diag = TRUE, det = TRUE)
  negloglik <- (iv$log.det.to.half - .maxim.aux1(S=temp.list$z, invcov=iv, trend = temp.list$xmat, log.f.sim = temp.list$log.f.sim, messages.screen=messages.screen)$logh)
  if(messages.screen) cat(paste("log-likelihood = ",-negloglik,"\n"))
  return(negloglik)
}


".lik.sim.boxcox" <-
  function(pars, fp, ip, temp.list)
{ 
### Function for finding m.l.e. for a given phi based on samples from mu=g^{-1}(S) ###############
### This function is only valid when all observations are positive.
  ##
  ## Obligatory parameters:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  if(ip$f.lambda) lambda <- fp$lambda
  else lambda <- pars[length(pars)]
  ##
  mu <- temp.list$mu
  ##
  messages.screen <- ifelse(is.null(temp.list$messages.screen), TRUE,temp.list$messages.screen)
  if(messages.screen) cat(paste("phi = ",phi, "tausq.rel = ",tausq.rel, "lambda= ", lambda,"\n"))
  ##
  ## computing the determinant of the transformation
  ##
  log.J.lambda <- apply(log(mu),2,sum)*(lambda-1)
  if(lambda ==0) mu <- log(mu)
  else mu <- (mu^lambda-1)/lambda
  ##
  ## Computing likelihood
  ##
  iv <- varcov.spatial(dists.lowertri = as.vector(dist(temp.list$coords)), cov.model = temp.list$cov.model, kappa = temp.list$kappa,
                       nugget = tausq.rel, cov.pars = c(1, phi), only.inv.lower.diag = TRUE, det = TRUE)
  negloglik <- (iv$log.det.to.half - .maxim.aux1(S=mu, invcov=iv, trend = temp.list$xmat, log.f.sim = temp.list$log.f.sim-log.J.lambda, messages.screen=messages.screen)$logh)
  if(messages.screen) cat(paste("log-likelihood = ",-negloglik,"\n"))
  return(negloglik)
}


"likfit.glsm" <-
  function (mcmc.obj, trend = mcmc.obj$trend,
            cov.model = "matern", 
            kappa = 0.5, ini.phi, fix.nugget.rel = FALSE, nugget.rel = 0, aniso.pars = NULL, 
            fix.lambda = TRUE, lambda = NULL, limits = pars.limits(), messages, ...)
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
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential","gaussian",
                           "spherical", "circular", "cubic",
                           "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(!is.null(kappa)){
    if(cov.model == "matern" & kappa == 0.5) cov.model <- "exponential"
  }
  ##
  if(is.null(mcmc.obj$S)){
    if(is.null(mcmc.obj$mu)) stop("mcmc.obj should include either an object mu or an object S.")
    if(is.null(lambda)) stop("must specify lambda")
    n <- temp.list$n <- nrow(mcmc.obj$mu)
    temp.list$mu <- mcmc.obj$mu
    if(mcmc.obj$family == "poisson") another.boxcox <- TRUE
    else stop("estimation of lambda is only possible for Poisson error distribution and boxcox link function")
  }
  else{
    if(!is.null(lambda)) warning("cannot use argument lambda with the given objects in mcmc.obj")
    if(!fix.lambda){
      warning("cannot estimate lambda from the given objects in mcmc.obj")
      fix.lambda <- TRUE
    }
    n <- temp.list$n <- nrow(mcmc.obj$S)
    temp.list$z <- mcmc.obj$S 
    another.boxcox <- FALSE
  }
  coords <- mcmc.obj$geodata$coords
  if (n != nrow(coords)) stop("Number of locations does not match with number of data")
  if(!is.null(aniso.pars)){
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be a vector with two elements: rotation angle (in radians) and anisotropy ratio (a number > 1)")
    coords <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
  }
  else{
    if(!is.null(mcmc.obj$aniso.pars)){
      coords <- coords.aniso(coords = coords, aniso.pars = mcmc.obj$aniso.pars)
      aniso.pars <- mcmc.obj$aniso.pars
    }
  }
  temp.list$xmat <- unclass(trend.spatial(trend = trend, geodata=mcmc.obj$geodata))
  beta.size <- temp.list$beta.size <- dim(temp.list$xmat)[2]
  ##
  temp.list$coords <- coords
  temp.list$cov.model <- cov.model
  temp.list$kappa <- kappa
  if(cov.model=="pure.nugget"){
    ini.phi <- 1
    if(!fix.nugget.rel) nugget.rel <- 0
    fix.nugget.rel <- TRUE
  }
  else if(missing("ini.phi")) stop("likfit.glsm : must specify ini.phi ")
  ini <- ini.phi
  lower.optim <- c(limits$phi["lower"])
  upper.optim <- c(limits$phi["upper"])
  fixed.values <- list()
  ##
  if(fix.nugget.rel) {
    fixed.values$tausq.rel <- nugget.rel
  }
  else {
    ini <- c(ini, nugget.rel)
    lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
    upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
  }
  if(another.boxcox){
    if(fix.lambda) {
      fixed.values$lambda <- lambda
    }
    else {
      ini <- c(ini, lambda)
      lower.optim <- c(lower.optim, limits$lambda["lower"])
      upper.optim <- c(upper.optim, limits$lambda["upper"])
    }
    ip <- list(f.tausq.rel = fix.nugget.rel, f.lambda = fix.lambda)
  }
  else ip <- list(f.tausq.rel = fix.nugget.rel)
  names(ini) <- NULL
  temp.list$log.f.sim <- mcmc.obj$log.f.sim
  temp.list$messages.screen <- messages.screen
  ##
  npars <- beta.size + 1 + ifelse(cov.model == "pure.nugget",0,1) + sum(!unlist(ip))
  if(cov.model != "pure.nugget" | !fix.lambda){
    if(messages.screen){
      cat("--------------------------------------------------------------------\n")
      cat("likfit.glsm: likelihood maximisation using the function optim.\n") 
    }
    if(another.boxcox){
      lik.optim <- optim(par = ini, fn = .lik.sim.boxcox, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                         fp = fixed.values, ip = ip, temp.list = temp.list, ...)
    }
    else{
      lik.optim <- optim(par = ini, fn = .lik.sim, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                         fp = fixed.values, ip = ip, temp.list = temp.list, ...)
    }
    ##
    if(messages.screen) 
      cat("likfit.glsm: end of numerical maximisation.\n")
    par.est <- lik.optim$par
    phi <- par.est[1]
    ##
    ## Values of the maximised likelihood
    ##
    loglik.max <-  - lik.optim$value
    ##
    ## Assigning values for estimated parameters
    ##
    if(!fix.nugget.rel){
      nugget.rel <- par.est[2]
    }
    if(!fix.lambda){
      lambda <- par.est[length(par.est)]
    }
  }
  if(cov.model == "pure.nugget") phi <- NA
  ##
  gc(verbose = FALSE)  
  ##
  ## Computing estimated beta and sigmasq
  if((is.na(phi) | phi < 1e-12))
    siv <- list(diag.inverse = rep(1/(1+nugget.rel), n), lower.inverse = rep(0,n*(n-1)/2))
  else{
    siv <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget.rel, cov.pars = c(1, phi),
                          only.inv.lower.diag = TRUE)
  }
  if(another.boxcox){
    if(lambda == 0){
      result <- .maxim.aux1(S = log(mcmc.obj$mu), invcov = siv, trend = temp.list$xmat, log.f.sim = temp.list$log.f.sim - apply(log(mcmc.obj$mu),2,sum)*(lambda-1), messages.screen=messages.screen)
    }
    else{
      result <- .maxim.aux1(S = (mcmc.obj$mu^lambda-1)/lambda, invcov = siv, trend = as.vector(temp.list$xmat),
                           log.f.sim = temp.list$log.f.sim - apply(log(mcmc.obj$mu),2,sum)*(lambda-1), messages.screen=messages.screen)
    }
    if(cov.model == "pure.nugget"){
      loglik.max <- result$logh
      lik.optim <- NULL
    }
    results <- list(family=mcmc.obj$family, link="boxcox", cov.model = cov.model, beta = result$beta, cov.pars=c(result$sigmasq, phi), nugget.rel = nugget.rel, 
                    kappa = kappa, aniso.pars = aniso.pars, lambda = lambda, trend = trend, npars=npars, loglik = loglik.max,
                    info.minimisation.function = lik.optim, call = call.fc)
  }
  else{
    result <- .maxim.aux1(S = mcmc.obj$S, invcov = siv, trend = temp.list$xmat,log.f.sim = temp.list$log.f.sim, messages.screen=messages.screen)
    if(cov.model == "pure.nugget"){
      loglik.max <- result$logh
      lik.optim <- NULL
    }
    results <- list(family=mcmc.obj$family, link=mcmc.obj$link, cov.model = cov.model, beta = result$beta, cov.pars = c(result$sigmasq, phi), nugget.rel = nugget.rel,
                    kappa = kappa, aniso.pars = aniso.pars, lambda = mcmc.obj$lambda, trend = trend, npars=npars, loglik = loglik.max,
                    info.minimisation.function = lik.optim, call = call.fc)
  }
  ##
  if(beta.size == 1) beta.name <- "beta"
  else beta.name <- paste("beta", 0:(beta.size-1), sep="")
  if(mcmc.obj$family == "poisson"){
    par.su <- data.frame(status=rep(-9, beta.size + 4))
    par.su$status <- c(rep("estimated", beta.size+2), ifelse(c(fix.nugget.rel,fix.lambda),"fixed", "estimated"))
    if(cov.model == "pure.nugget"){
      par.su$status[beta.size+2] <- ""
    }
    par.su$values <- round(c(results$beta, results$cov.pars, results$nugget.rel, results$lambda), digits=4)
    row.names(par.su) <- c(beta.name, "sigmasq", "phi", "tausq.rel", "lambda")
  }
  else{  ### binomial distribution so far only includes the logit-link. Therefore no transformation parameter
    par.su <- data.frame(status=rep(-9, beta.size + 3))
    par.su$status <- c(rep("estimated", beta.size+2), ifelse(fix.nugget.rel,"fixed", "estimated"))
    if(cov.model == "pure.nugget"){
      par.su$status[beta.size+2] <- ""
    }
    par.su$values <- round(c(results$beta, results$cov.pars, results$nugget.rel), digits=4)
    row.names(par.su) <- c(beta.name, "sigmasq", "phi", "tausq.rel")
  }
  results$parameters.summary <- par.su
  class(results) <- "likGLSM"
  return(results)
}

"print.likGLSM" <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  est.pars <- as.vector(x$parameters.summary[x$parameters.summary[,1] == "estimated",2])
  names.est.pars <- dimnames(x$parameters.summary[x$parameters.summary[,1] == "estimated",])[[1]]
  names(est.pars) <- names.est.pars
  cat("likfit.glsm: estimated model parameters:\n")
  print.default(format(est.pars, digits=digits), ...)
  cat("\n likfit.glsm : maximised log-likelihood = ")
  cat(format(x$loglik, digits=digits))
  cat("\n")
  return(invisible())
}

"summary.likGLSM" <-
  function(object, ...)
{
  names.pars <- dimnames(object$parameters.summary)[[1]]
  summ.lik <- list()
  summ.lik$method.lik <- "maximum likelihood"
  summ.lik$family <- object$family
  summ.lik$link <- object$link
  summ.lik$mean.component <- object$beta
  names(summ.lik$mean.component) <- names.pars[seq(along=object$beta)]
  summ.lik$cov.model <- object$cov.model
  summ.lik$kappa <- object$kappa
  summ.lik$aniso.pars <- object$aniso.pars
  summ.lik$spatial.component <- object$parameters.summary[c("sigmasq", "phi"),]
  summ.lik$nugget.component <- object$parameters.summary[c("tausq.rel"),, drop=FALSE]
  summ.lik$transformation  <- object$parameters.summary[c("lambda"),, drop=FALSE]
  summ.lik$likelihood <- list(log.L = object$loglik, n.params = as.integer(object$npars))
  summ.lik$estimated.pars <- dimnames(object$parameters.summary[object$parameters.summary[,1] == "estimated",])[[1]]
  summ.lik$call <- object$call
  class(summ.lik) <- "summary.likGLSM"
  return(summ.lik)
}

"print.summary.likGLSM" <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(length(class(x)) == 0 || all(class(x) != "summary.likGLSM"))
    stop("object is not of the class \"summary.likGLSM\"")
  cat("Summary of the maximum likelihood parameter estimation\n")
  cat("-----------------------------------\n")
  cat(paste("Family = ", x$family, ", Link = ", x$link, "\n" ))
  cat("\n")
  cat("Parameters of the mean component (trend):")
  cat("\n")
  print.default(format(x$mean.component, digits=digits), ...)
  cat("\n")
  ##
  cat("Parameters of the spatial component:")
  cat("\n")
  cat(paste("   correlation function:", x$cov.model))
  if(x$cov.model == "matern" | x$cov.model == "powered.exponential" |
     x$cov.model == "cauchy" | x$cov.model == "gneiting.matern"){
    cat(paste("\n          kappa = ", x$kappa))
    if(x$cov.model == "matern" & (round(x$kappa, digits=digits)  == 0.5)) cat(" (exponential)")
  }
  cat("\n")
  if(!is.null(x$aniso.pars)){
    cat(paste("\n (fixed) anisotropy parameters (angle, ratio) = (", x$aniso.pars, "( \n"))
  }
  cat(paste("\n      (estimated) variance parameter sigmasq (partial sill) = ", format(x$spatial.component[1,2], digits=digits)))
  if(x$cov.model != "pure.nugget") cat(paste("\n      (estimated) cor. fct. parameter phi (range parameter)  = ", format(x$spatial.component[2,2], digits=digits)))
  cat("\n")
  if(x$nugget.component[,1] == "estimated")
    cat(paste("\n (estimated) relative nugget = ", format(x$nugget.component[,2], digits=digits)))
  else
    cat(paste("\n (fixed) relative nugget =", x$nugget.component[,2]))
  cat("\n")
  cat("\n")
  if(x$family == "poisson" && x$link == "boxcox"){
    cat("\n")
    cat("Transformation parameter:")
    cat("\n")
    lambda <- x$transformation[,2]
    if(x$transformation[,1] == "estimated")
      cat(paste("      (estimated) Box-Cox parameter =", format(lambda, digits=digits)))
    else{
      cat(paste("      (fixed) Box-Cox parameter =", lambda))
      if(abs(lambda - 1) <  0.0001) cat(" (no transformation)")
      if(abs(lambda) < 0.0001) cat(" (log-transformation)")
    }
  }
  cat("\n")
  cat("\n")
  cat("Maximised Likelihood:")
  cat("\n")
  print(format(x$likelihood, digits=digits))
  cat("\n")
  cat("Call:")
  cat("\n")
  print(x$call)
  cat("\n")
  invisible(x)
}

