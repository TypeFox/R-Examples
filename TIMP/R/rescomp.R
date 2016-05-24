"rescomp" <-
  function (theta=vector(), d=vector(), currModel=currModel, currTheta=vector()) 
  {
    if(length(currTheta) == 0) 
      currTheta <- getThetaCl(theta, currModel)
    groups <- currModel@groups 
    m <- currModel@modellist
    resid <- clpindepX <-list()
    nexp <- length(m) 
    for(i in 1:nexp) {
      clpindepX[[i]] <- if(!m[[i]]@clpdep || m[[i]]@getX) 
        getClpindepX(model = m[[i]], theta =
                       currTheta[[i]], multimodel = currModel,
                     returnX = FALSE, rawtheta= theta, dind=0)
      else matrix() 
    }
    for(i in 1:length(groups)) {
      resid[[i]] <- residPart(model = m[[1]], 
                              group = groups[[i]], multimodel = currModel, 
                              thetalist = currTheta, clpindepX = clpindepX,
                              finished = currModel@finished,
                              returnX = FALSE, rawtheta = theta) 
      if(currModel@finished){
        currModel <- fillResult(group = groups[[i]],
                                multimodel = currModel, thetalist = currTheta,
                                clpindepX = clpindepX, rlist = resid[[i]],
                                rawtheta = theta)
      }
    }
    if(currModel@finished) {
      currModel@fit@nlsres$onls$nclp <- currModel@nclp
      if(currModel@optlist[[1]]@sumnls) {
        if(class(currModel@fit@nlsres$onls) == "nls")
          class(currModel@fit@nlsres$onls) <- "timp.nls"
        else if(class(currModel@fit@nlsres$onls) == "nls.lm")
          class(currModel@fit@nlsres$onls) <- "timp.nls.lm"
        else
          class(currModel@fit@nlsres$onls) <- "timp.optim"
        currModel@fit@nlsres$sumonls <- summary(currModel@fit@nlsres$onls,
                                                currModel=currModel,
                                                currTheta=currTheta)
      }
      if(currModel@stderrclp) {
        for(i in 1:length(groups)) {
          currModel <- getStdErrClp(group = groups[[i]],
                                    multimodel = currModel, thetalist = currTheta,
                                    clpindepX = clpindepX, rlist = resid[[i]],
                                    rawtheta = theta)
        }
      }
    }
    ## if using a trilinear type model, we have cp=AE; so separate A out. 
    if(currModel@finished && currModel@trilinear){ 
      trires <- triResolve(currModel, currTheta)
      currModel <- trires$currModel
      currTheta <- trires$currTheta
    }
    if(currModel@finished && m[[1]]@mod_type == "kin") {
      if (m[[1]]@fullk) { 
        for(i in 1:nexp) {
          nocolsums <- length(m[[1]]@lightregimespec) > 0 # lightdiff (see compModel.R)
          eig <- fullKF(currTheta[[i]]@kinpar, currTheta[[i]]@kinscal, m[[1]]@kmat, currTheta[[i]]@jvec, m[[1]]@fixedkmat, m[[1]]@kinscalspecial,
                        m[[1]]@kinscalspecialspec, nocolsums)
          currTheta[[i]]@eigenvaluesK <- eig$values
        }
      }
    }
    if(currModel@finished) {
      return(list(currModel=currModel,currTheta=currTheta))
    }
    if(currModel@algorithm == "optim") ## minimize this sum 
      retval <- sum(unlist(resid))  
    else
      retval <- unlist(resid) ## nls and nls.lm want the residuals,
    ## to minimize the sum of their squares 
    retval
  }

