sim_ACD <- function(N = 1000, 
                    model = "ACD",
                    dist = "exponential",
                    param = NULL, 
                    order = NULL,
                    Nburn = 50,
                    startX = c(1),
                    startMu = c(1),
                    errors = NULL,
                    sampleErrors = TRUE,
                    roundToSec = FALSE,
                    rm0 = FALSE,
                    diurnalFactor = FALSE,
                    splineObj = NULL,
                    open = NULL,
                    close = NULL){
  
  #provides the possibility of entering truncated and/or case mismatched arguments:
  model <- match.arg(toupper(model), c("ACD", "LACD1", "LACD2", "AMACD", "ABACD"))
  dist <- match.arg(tolower(dist), c("exponential", "weibull", "burr", "gengamma", "genf"))
  
  distCode <- .getDistCode(dist)
  #checks param and order input:
  if(length(param) != 0){
    if(length(order) == 0) order <- .setOrder(model)
    .checkOrderAndPara(order, param, distCode, model)
    paraTemp <- .seperateStartPara(param, model, distCode, order)
    distPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
  }else{
    if(length(order) != 0){
      .checkOrder(order, model)
    } else{
      order <- .setOrder(model)
    }
    paraTemp <- .setStartPara(model, distCode, 1, order)
    distPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
    param <- c(distPara, startPara)
  }
  
  if(length(errors) == 0){
    if(dist == "exponential") e <- stats::rexp(N + Nburn)
    else if(dist == "weibull"){
      e <- stats::rweibull(N + Nburn, shape = distPara, scale = 1/(gamma(1+1/distPara)))
    }
    else if(dist == "burr"){
      kappa <- distPara[1]
      sig2 <- distPara[2]
      muPara <- burrExpectation(theta = 1, kappa = kappa, sig2 = sig2)^kappa
      e <- rburr(N + Nburn, theta = 1, kappa = 1.2, sig2 = .3)
    } else if(dist == "gengamma"){
      kappa <- distPara[1]
      gammaPara <- distPara[2]
      
      e <- rgengamma(N + Nburn, gamma = gammaPara, kappa = kappa, forceExpectation = T)
    } else if(dist == "genf"){
      stop("Simulations are not available for the generelized F distribution")
    }
  } else{
    if(sampleErrors) e <- sample(errors, size = N + Nburn, replace = TRUE)
    else{
      if(length(errors) != N + Nburn) stop("the 'errors' vector needs to be of length N + Nburn if sampleErrors = FALSE")
      e <- errors
    }
  } 
    
  maxpq = max(order)
  if(maxpq > min(length(startX), length(startMu))){
    startX <- rep(startX, length.out = maxpq)
    startMu <- rep(startMu, length.out = maxpq)
  } 
  
  if(diurnalFactor){
    if(length(splineObj) == 0){ 
      splineObj <- ACDm::defaultSplineObj
      open = "10:00:00"
      close = "18:25:00"
    }
    
    knots <- c(splineObj[[1]]$knots, splineObj[[2]]$knots, splineObj[[3]]$knots, splineObj[[4]]$knots, splineObj[[5]]$knots)*60
    konst <- c(splineObj[[1]][[2]][,1], splineObj[[2]][[2]][,1], splineObj[[3]][[2]][,1], splineObj[[4]][[2]][,1], splineObj[[5]][[2]][,1])
    lin <- c(splineObj[[1]][[2]][,2], splineObj[[2]][[2]][,2], splineObj[[3]][[2]][,2], splineObj[[4]][[2]][,2], splineObj[[5]][[2]][,2])/60
    sq <- c(splineObj[[1]][[2]][,3], splineObj[[2]][[2]][,3], splineObj[[3]][[2]][,3], splineObj[[4]][[2]][,3], splineObj[[5]][[2]][,3])/60^2
    qub <- c(splineObj[[1]][[2]][,4], splineObj[[2]][[2]][,4], splineObj[[3]][[2]][,4], splineObj[[4]][[2]][,4], splineObj[[5]][[2]][,4])/60^3
    splineNewDay <- cumsum(c(0, length(splineObj[[1]]$knots), length(splineObj[[2]]$knots), length(splineObj[[3]]$knots) , length(splineObj[[4]]$knots)))
    
    opensek <- as.POSIXlt(strptime(open, "%H:%M:%S"))
    opensek <- opensek$h * 3600 + opensek$min * 60 + opensek$sec
    closesek <- as.POSIXlt(strptime(close, "%H:%M:%S"))
    closesek <- closesek$h * 3600 + closesek$min * 60 + closesek$sec
    
    temp<-.Call("sim_ACDSpline",
                as.integer(N),
                param[1:(1+order[1]+order[2])],
                order,
                startX,
                startMu,
                e,
                as.integer(Nburn),
                opensek,
                closesek,
                knots,
                konst,
                lin,
                sq,
                qub,
                splineNewDay, PACKAGE = "ACDm")   
    
    #time <- strptime("2014-01-06 00:00:00", "%Y-%m-%d %H:%M:%S") + ((temp[[1]] %/% 432000) * 604800 + (temp[[1]] %% 432000)) * 86400 
    
    if(roundToSec){
      df <- data.frame(time = strptime("2014-01-06 00:00:00", "%Y-%m-%d %H:%M:%S") + ((temp[[1]] %/% 5) * 7 + (temp[[1]] %% 5)) * 60 * 60 * 24 + ceiling(temp[[2]]))
      utils::capture.output(dur <- computeDurations(transactions = df, open = open, close = close, rm0dur = F, type = "transactions"))      
    } else { #doesnt yet work 
      df <- data.frame(time = strptime("2014-01-06 00:00:00", "%Y-%m-%d %H:%M:%S") + ((temp[[1]] %/% 5) * 7 + (temp[[1]] %% 5)) * 60 * 60 * 24 + temp[[2]])
      utils::capture.output(dur <- computeDurations(transactions = df, open = open, close = close, rm0dur = F, type = "transactions"))
    }    
    return(df)
  } else if(!diurnalFactor){
    
    cFunction <- switch(model,
                        ACD = "sim_ACDCALL",
                        LACD1 = "sim_LACD1",
                        LACD2 = "sim_LACD2",
                        AMACD = "sim_AMACD",
                        ABACD = "sim_ABACD")
    
    if(!roundToSec){
      return(.Call(cFunction,
                   as.integer(N),
                   startPara,
                   order,
                   startX,
                   startMu,
                   e,
                   as.integer(Nburn), PACKAGE = "ACDm"))
    } else{
      durTemp <- ceiling(cumsum(.Call(cFunction,
                                      as.integer(N),
                                      startPara,
                                      order,
                                      startX,
                                      startMu,
                                      e,
                                      as.integer(Nburn), PACKAGE = "ACDm")))
      if(!rm0)  return(c(durTemp[1], durTemp[-1]-durTemp[-N]))
      else{
        durTemp <- c(durTemp[1], durTemp[-1]-durTemp[-N])
        return(durTemp[durTemp != 0])
      }
    }
  }
}