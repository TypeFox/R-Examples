"fitModel" <-
  function (data, modspec = list(), datasetind = vector(), modeldiffs = list(), 
            opt = opt(),lprogress=FALSE) 
  {
    currModel <- getModel(data, modspec, modeldiffs, datasetind, opt)
    globalEnvir = as.environment(1)
    
    tr <- getTheta(currModel)
    theta <- tr$theta
    currModel <- tr$mod  
    
    currTheta <- getThetaCl(theta, currModel)
    
    iter <- opt@iter
    
    assign(".currTheta", currTheta, envir = globalEnvir)
    assign(".currModel", currModel, envir = globalEnvir)
    
    if(opt@algorithm == "nls") {
      if (lprogress) 
        nlsprogress<-capture.output(currModel@fit@nlsres$onls <- nls(~rescomp(theta=t,d=d,currModel=currModel),
                                                                     data=list(d=vector(),currModel=currModel), 
                                                                     control =
                                                                       nls.control(maxiter = iter,
                                                                                   minFactor = opt@minFactor,
                                                                                   warnOnly = TRUE,
                                                                                   printEval = FALSE),
                                                                     start = list(t = theta),
                                                                     algorithm = opt@nlsalgorithm,
                                                                     trace = TRUE))
      else  
        currModel@fit@nlsres$onls <- nls(~rescomp(theta=t,d=d,currModel=currModel),
                                         data=list(d=vector(),currModel=currModel), 
                                         control =
                                           nls.control(maxiter = iter,
                                                       minFactor = opt@minFactor,
                                                       warnOnly = TRUE,
                                                       printEval = TRUE),
                                         start = list(t = theta),
                                         algorithm = opt@nlsalgorithm,
                                         trace = TRUE)
    }
    else if(opt@algorithm == "nls.lm") {
      if(length(opt@parscale)!=0) {
        if(opt@parscale[1] == "abs")
          parscale <- 1/abs(theta)
        else
          parscale <- opt@parscale
        currModel@fit@nlsres$onls <- nls.lm(par=theta, fn = rescomp, 
                                            currModel=currModel,
                                            control=list(nprint=1,
                                                         maxiter = iter,maxfev=opt@maxfev,
                                                         diag = parscale))
      }
      else
        currModel@fit@nlsres$onls <- nls.lm(par=theta, fn = rescomp, 
                                            currModel=currModel,
                                            control=list(nprint=1,
                                                         maxiter = iter,maxfev=opt@maxfev))
    }
    else if(opt@algorithm == "optim") {
      currModel@fit@nlsres$onls <- optim(par=theta, fn = rescomp,
                                         currModel=currModel,
                                         method = opt@optimmethod,
                                         control=list(trace=100,
                                                      maxit=iter),
                                         hessian=TRUE)
      
    }
    currModel@finished <- TRUE
    if(opt@algorithm == "nls")
      resFinal <- rescomp(theta=currModel@fit@nlsres$onls$m$getPars(),
                          currModel=currModel)
    else 
      resFinal <- rescomp(theta=currModel@fit@nlsres$onls$par, currModel=currModel)
    
    currModel <- resFinal$currModel
    currTheta <- resFinal$currTheta
    
    assign(".currTheta", currTheta, envir = globalEnvir)
    assign(".currModel", currModel, envir = globalEnvir)
    
    if (opt@plot) 
      plotter(currModel@modellist[[1]], currModel, currTheta, opt)
    
    currModel@optlist[[1]] <- opt
    
    if (lprogress) 
      ret <- list(currModel = currModel, currTheta = currTheta,nlsprogress=nlsprogress)
    else
      ret <- list(currModel = currModel, currTheta = currTheta)
    
    return(ret)
  }
