acdFit <- function(durations = NULL, model = "ACD", dist = "exponential", 
                   order = NULL, startPara = NULL,  dailyRestart = 0, optimFnc = "optim", 
                   method = "Nelder-Mead", output = TRUE, bootstrapErrors = FALSE, 
                   forceErrExpec = TRUE, fixedParamPos = NULL, bp = NULL, control = list()){  
  
  
  vecProvided <- FALSE
  if(is.data.frame(durations)){   
    
    if("adjDur" %in% colnames(durations)){
      dur <- durations$adjDur
    } else if("durations" %in% colnames(durations)){
      warning("no 'adjDur' column for diurnally adjusted durations found - used unadjusted durations instead")
      dur <- durations$durations
    } else stop("neither a 'durations' or a 'adjDur' column was found in the data.frame 'durations'")
    
    if("time" %in% colnames(durations)){
      if(!("POSIXlt" %in% class(durations$time))) durations$time <- as.POSIXlt(durations$time)
      time = durations$time
    } else time <- NULL      
    
  } else if(is.vector(durations)){
    
    dur = durations
    vecProvided <- TRUE
    time <- NULL
    
  } else stop("'durations' must be a data.frame or a vector")
  
  N <- length(dur)
  mean <- mean(dur)
  currentTime <- Sys.time()
  
  #provides the possibility of entering truncated and/or case mismatched arguments:
  model <- match.arg(toupper(model), c("ACD", "LACD1", "LACD2", "AMACD", "ABACD", "BACD", "SNIACD", "LSNIACD"))
  dist <- match.arg(tolower(dist), c("exponential", "weibull", "burr", "gengamma", "genf", "qweibull", "mixqwe", "mixqww", "mixinvgauss"))
  
  distCode <- .getDistCode(dist)
  #checks startPara and order input:
  if(length(startPara) != 0){
    if(length(order) == 0) order <- .setOrder(model)
    .checkOrderAndPara(order, startPara, distCode, model)
    paraTemp <- .seperateStartPara(startPara, model, distCode, order)
    distStartPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
  }else{
    if(length(order) != 0){
      .checkOrder(order, model)
    } else{
      order <- .setOrder(model)
    }
    startPara <- .setStartPara(model, distCode, mean, order)
    distStartPara <- startPara$distStartPara
    startPara <- startPara$startPara
  }
  
  if(model %in% c("SNIACD", "LSNIACD") && length(bp) == 0) bp <- .setBP(order[3])
  
  #checks the control list arguments:
  con <- list(newDay = 0,
              maxit = 4000,
              trace = 0,
              B = 999,
              BootRoundTosec = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  #controls that time was provided if dailyRestart == 1, and computes the vector of indices for first observation of a news day
  if(dailyRestart != 0 && length(time) != 0){
    if(length(con$newDay) == 0 || con$newDay == 0) con$newDay <- .getNewDay(time) #computes the vector of indices for first observation of a news day
  } else if(dailyRestart != 0 && length(time) == 0 && con$newDay == 0){
    warning("can only use daily restart of the conditional mean if the clocktimes of transactions or the vector newDay (as a control parameter) are provided! Estimation done assuming no daily restart.")
  }
  
  mean <- mean(dur)    
  
  #makes it possible to have fixed parameters when optimizing the log likelihood:
  fixedParam <- NULL
  if(length(fixedParamPos) != 0){
    fixedParam <- startPara[fixedParamPos]
    startPara <- startPara[!fixedParamPos]
  }  
  
  if(optimFnc == "optim"){    
    fit <- stats::optim(startPara, .getLLcall,
                 method = method, hessian=TRUE,
                 dur = dur, model = model, order=order, distCode = distCode, newDay = con$newDay,
                 mean=mean, returnMu = FALSE, breakPoints = bp, forceErrExpec = forceErrExpec, 
                 fixedParam = fixedParam, fixedParamPos = fixedParamPos, control = list(maxit = con$maxit, trace = con$trace))
    
    #if some parameters were set to be fixed, the fixed and estimated parameters are recombined:
    if(length(fixedParamPos) != 0) parTemp <- .returnfixedPara(fit$par, fixedParam, fixedParamPos)
    else parTemp <- fit$par
    
    mu <- .getLLcall(param = parTemp, dur = dur, model = model, order = order, mean = mean, distCode = distCode,  newDay = con$newDay, returnMu = TRUE, breakPoints = bp, forceErrExpec = forceErrExpec)  
  } else if(optimFnc == "nlminb"){
    fit <- stats::nlminb(start = startPara, objective = .getLLcall, dur = dur, 
                  model = model, order=order, distCode = distCode, newDay = con$newDay, 
                  mean=mean, returnMu = FALSE, breakPoints = bp, forceErrExpec = forceErrExpec,
                  fixedParam = fixedParam, fixedParamPos = fixedParamPos, 
                  control = list(trace = con$trace, iter.max = con$maxit), lower = -Inf, upper = Inf)
    hessianTemp <- stats::optimHess(fit$par, .getLLcall, 
                             dur = dur, model = model, order=order, distCode = distCode, newDay = con$newDay, 
                             mean=mean, returnMu = FALSE, breakPoints = bp)
    fit <- list(par = fit$par, hessian = hessianTemp, value = fit$objective, convergence = fit$convergence, counts = fit$evaluations[2])
    mu <- .getLLcall(param = fit$par, dur = dur, model = model, order=order, mean = mean, distCode = distCode,  newDay = con$newDay, returnMu = TRUE, breakPoints = bp, forceErrExpec = forceErrExpec)  
  } else if(optimFnc == "solnp"){
    options(warn = -1)
    utils::capture.output(
      fit <- Rsolnp::solnp(pars=startPara, fun = .getLLcall, 
                   #ineqfun = ineq, ineqUB = .999, ineqLB = 0, LB = LB, UB = UB, 
                   dur = dur, model = model, order = order, mean = mean, distCode = distCode, returnMu = FALSE, forceErrExpec = forceErrExpec,  
                   breakPoints = bp, newDay = con$newDay, control = list(outer.iter = con$maxit, trace = con$trace))
    )
    options(warn = 0)
    fit <- list(par = fit$pars, hessian = fit$hessian, 
                value = fit$values[length(fit$values)], convergence = fit$convergence, counts = fit$nfuneval)
    mu <- .getLLcall(param = fit$par, dur = dur, model = model, order=order, mean = mean, distCode = distCode,  newDay = con$newDay, returnMu = TRUE, breakPoints = bp, forceErrExpec = forceErrExpec)
  } else if(optimFnc == "optimx"){
    fit <- optimx::optimx(startPara, .getLLcall, 
                          method = method, 
                          dur = dur, model = model, order=order, distCode = distCode, newDay = con$newDay, mean=mean, forceErrExpec = forceErrExpec, 
                          returnMu = FALSE, breakPoints = bp, itnmax = con$maxit, control = list(trace = con$trace, kkt = FALSE))
    hessianTemp <- stats::optimHess(fit[1, 1:length(startPara)], .getLLcall, 
                             dur = dur, model = model, order=order, distCode = distCode, newDay = con$newDay, 
                             mean = mean, returnMu = FALSE, breakPoints = bp)
    rownames(hessianTemp) <- NULL; colnames(hessianTemp) <- NULL
    fit <- list(par = as.numeric(fit[1, 1:length(startPara)]), hessian = hessianTemp, 
                value = fit$value, convergence = fit$convcode, counts = fit$fevals)    
    mu <- .getLLcall(param = fit$par, dur = dur, model = model, order=order, mean = mean, distCode = distCode,  newDay = con$newDay, returnMu = TRUE, breakPoints = bp, forceErrExpec = forceErrExpec)      
  }
  
  ##if(model == "ACD"){
  if(bootstrapErrors){
    if(!con$BootRoundTosec){ #the simulation in the bootstraps wont be rounded to seconds
      
      bootPar <- matrix(nrow = con$B, ncol = length(fit$par))
      i <- 1
      percDone = 5
      failed = 0
      bootConverged <- rep(-99, con$B)
      
      bootStartTime <- Sys.time()
      while(i <= con$B){
        bootDur <- sim_ACD(N, model = model, param = fit$par, order = order, startX = mean, startMu = mean, errors = mu$resi, dist = dist, roundToSec = FALSE)
        bootTemp <- tryCatch(stats::optim(par = fit$par, fn = .getLLcall, dur = bootDur, model = model, order = order, mean = mean(bootDur), distCode = distCode, returnMu = FALSE,  hessian = F, control = list(maxit = con$maxit, trace = con$trace, trace = con$trace)), error = function(e) {NULL})
        bootParTemp <- bootTemp$par 
        
        if(length(bootParTemp) != 0 && all(abs(bootParTemp[-1]) < 1.5) && bootTemp$convergence == 0){
          bootPar[i, ] <- bootParTemp
          bootConverged[i] <- bootTemp$convergence
          i <- i + 1          
        } else failed <- failed + 1
        
        if ((i / con$B) >= .05 && percDone == 5) cat("Estimated time for bootstrap simulation: ", round(difftime(Sys.time(), bootStartTime, units = "secs")*20), " sec \n\nbootstrap % done: \n")
        if ((i / con$B) >= (percDone / 100)) {cat(percDone,"% "); percDone = percDone + 5}
      }      
      cat("\ntime for bootstrap simulation: ", round(difftime(Sys.time(), bootStartTime, units = "secs")), " sec")
      cat("\n\n")
      cat(failed, "of the ", con$B, " bootstrap estimations failed and were resimulated\n")
      ##cat("of the successful bootstraps ", sum(bootConverged == 0)," converged using the numerical maximazation function")
      ##cat("\n")
      bootErr <- sqrt(diag(stats::cov(bootPar)))
      bootCorr <- stats::cor(bootPar)
      bootMean <- apply(bootPar, 2, mean)
    } else{
      mu <- .getLLcall(param = fit$par, dur = dur, model = model, order=order, mean = mean, distCode = distCode, newDay = con$newDay, returnMu = TRUE) 
      bootPar <- matrix(nrow = con$B, ncol = length(fit$par))
      i <- 1
      percDone = 5
      cat("bootstrap % done: ")
      bootDurTemp <- sim_ACD((con$B*(N+50))*1.3, param = fit$par, model = model, order = order, startX = mean, startMu = mean, errors = mu$resi, roundToSec = FALSE)
      bootDur <- bootDurTemp[bootDurTemp!=0]
      if(bootDur<N) 
        while(i <= con$B){
          bootParTemp <- tryCatch(stats::optim(par = fit$par, fn = .getLLcall, model = model, x=bootDur[((i-(1+failed))*(N+50)+1):((i+failed)*(N+50))], order=order,mean=mean(bootDur), dist=distCode, returnMu=FALSE,  newDay = con$newDay, hessian = TRUE, control = list(maxit = con$maxit, trace = con$trace))$par, error = function(e) {NULL})
          if(length(bootParTemp) != 0){
            bootPar[i, ] <- bootParTemp
            i <- i + 1          
          } else failed <- failed + 1
          if ((i / con$B) >= (percDone / 100)) {cat(percDone,"% "); percDone = percDone + 5}
        }      
      cat("\n")
      bootErr <- sqrt(diag(stats::cov(bootPar)))
      bootCorr <- stats::cor(bootPar)
      bootMean <- apply(bootPar, 2, mean)
    }
  }
  
  #if QML and not preset (fixed) paramters, then robust errors will be computed:
  if(model == "ACD" && dist == "exponential" && length(fixedParamPos) == 0){
    QLMscore <- .Call("getScoreACDExp",
                      as.double(dur),
                      as.double(mu$mu),
                      as.double(fit$par),                     
                      as.integer(order),
                      as.integer(0), PACKAGE = "ACDm")
    
    sandwich <- solve(as.matrix(as.data.frame(QLMscore[3]))) %*% as.matrix(as.data.frame(QLMscore[4])) %*% solve(as.matrix(as.data.frame(QLMscore[3])))
    robustSE <- sqrt(diag(sandwich)) 
    robustCorr <- solve(diag(robustSE)) %*% sandwich %*% solve(diag(robustSE))
  }
  else{
    robustSE <- NULL
    robustCorr <- NULL
  }
  
  if(bootstrapErrors) namedParameters <- .getCoef(para = fit$par , model = model, dist = dist, hessian = fit$hessian, 
                                                 order = order, bootError = bootErr, bootCorr = bootCorr, bootMean = bootMean, 
                                                 robustCorr = robustCorr, robustSE = robustSE, fixedParam = fixedParam, fixedParamPos = fixedParamPos)
  else namedParameters <- .getCoef(para = fit$par , model = model, dist = dist, hessian = fit$hessian, order = order, 
                                  robustCorr = robustCorr, robustSE = robustSE, fixedParam = fixedParam, fixedParamPos = fixedParamPos)
  
  
  
  N <- length(dur)
  Npar <- length(startPara)
  LogLikelihood <- -fit$value
  AIC <- 2 * (Npar - LogLikelihood)
  BIC <- -2 * LogLikelihood + Npar * log(N)
  MSE <- mean((dur-mu$mu)^2)
  GoodnessOfFit <- data.frame("value" = c(LogLikelihood, AIC, BIC, MSE))
  rownames(GoodnessOfFit) <- c("LogLikelihood", "AIC", "BIC", "MSE")
  #computes the value of the unfree distribution parameter as a function of the others (if mean was forced to be 1):
  if(forceErrExpec == 1) forcedDistPara <- .returnFixedMeanPara(distCode, namedParameters$DPar)
  else{ #if forceErrExpec was false (0) this parameter was fixed at 1 in the estimation
    forcedDistPara <- 1
    names(forcedDistPara) <- names(.returnFixedMeanPara(distCode, namedParameters$DPar)) #only to get the name of the parameter
  } 
  
  if(!vecProvided) returnValue <- list(call = match.call(),
                                       durations = durations)
  else returnValue <- list(call = match.call(),
                           durations = data.frame(durations = dur))
  
  returnValue <- append(returnValue, list(
                                          muHats = mu$mu, 
                                          residuals = mu$resi,
                                          model = model,
                                          order = order,
                                          distribution = dist, 
                                          distCode = distCode,
                                          startPara = startPara,
                                          mPara = namedParameters$MPar, 
                                          dPara = namedParameters$DPar, 
                                          Npar = Npar[1],
                                          goodnessOfFit = GoodnessOfFit,
                                          parameterInference = namedParameters$Inference,
                                          forcedDistPara = forcedDistPara,
                                          comments = namedParameters$comment,
                                          hessian = fit$hessian,
                                          N = N,
                                          evals = fit$counts[1],
                                          convergence = fit$convergence,
                                          estimationTime = difftime(Sys.time(), currentTime, units = "secs"),
                                          description = paste("Estimated at", currentTime, "by user", Sys.info()[["user"]]),
                                          newDayVector = con$newDay))
  
    
  #if bootstrapp errors: adds the bootstrapp inference
  if(bootstrapErrors) returnValue <- append(returnValue, list(bootstrapEstimates = bootPar, bootConverged = bootConverged, bootErr = bootErr, bootMean = bootMean, bootCorr = namedParameters$bootCorr, bootPar = bootPar))
  #if QML (ACD and exponetial): adds the robust correlation
  if(model == "ACD" && dist == "exponential" && length(fixedParamPos) == 0) returnValue <- append(returnValue, list(robustCorr = namedParameters$robustCorr))
  #if SNIACD: adds the break points:
  if(model %in% c("SNIACD", "LSNIACD")) returnValue <- append(returnValue, list(SNIACDbp = bp))
  
  
  class(returnValue) <-  c("acdFit", class(returnValue))
  
  if(output) print(returnValue)
  
  acdFit <- returnValue
}
