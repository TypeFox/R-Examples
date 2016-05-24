
#  Fit the penalized occupancy models of Hutchinson et al (2015).

computeMPLElambda = function(formula, data, knownOcc = numeric(0), starts, method = "BFGS", engine = c("C", "R")) {

  designMats <- getDesign(data, formula)
  X <- designMats$X; V <- designMats$V; y <- designMats$y
  removed <- designMats$removed.sites
  y <- truncateToBinary(y)

  ## convert knownOcc to logical so we can correctly to handle NAs.
  knownOccLog <- rep(FALSE, numSites(data))
  knownOccLog[knownOcc] <- TRUE
  if(length(removed)>0)
     knownOccLog <- knownOccLog[-removed]

  nDP <- ncol(V)
  nOP <- ncol(X)
  nP <- nDP + nOP
  if(!missing(starts) && length(starts) != nP)
      stop(paste("The number of starting values should be", nP))
  if(missing(starts)) starts <- rep(0, nP)

  LRparams = glm.fit(x=X,y=apply(y,1,max),family=binomial(),intercept=F,start=starts[1:nOP])
  naiveOcc = mean(LRparams$fitted.values)
  occuOutMLE = occu(formula,data,knownOcc, starts,
                 method = "BFGS", engine = c("C", "R"), se = TRUE)
  meanDet = mean((1+exp(-occuOutMLE[2]@estimates%*%t(V)))^-1)
  MPLElambda = sqrt(sum(diag(occuOutMLE[2]@covMat)))*(1-(1-meanDet)^(dim(y)[2]))*(1-naiveOcc) # what if there are different numbers of visits to different sites?
  return(MPLElambda)      
}

occuPEN_CV <- function(formula, data, knownOcc = numeric(0), starts,
                 method = "BFGS", engine = c("C", "R"), 
		 lambdaVec = c(0,2^seq(-4,4)),
		 pen.type = c("Bayes","Ridge"),
		 k = 5,
		 foldAssignments = NA,
		 ...) 
{
  if(!is(data, "unmarkedFrameOccu"))
        stop("Data is not an unmarkedFrameOccu object.")

  pen.type = pen.type[1]
  if (pen.type=="MPLE") stop("MPLE does not require cross-validation.")
  if (!(pen.type=="Bayes" | pen.type=="Ridge")) 
      stop("pen.type not recognized.  Choose Bayes or Ridge.")

  if (length(lambdaVec)==1) stop("Must provide more than one lambda for cross-validation.")
  
  engine <- match.arg(engine, c("C", "R"))
  designMats <- getDesign(data, formula)
  X <- designMats$X; V <- designMats$V; y <- designMats$y
  y <- truncateToBinary(y)
  J <- ncol(y)
  M <- nrow(y)
  
  if (!(length(foldAssignments)==1 & is.na(foldAssignments)[1])) { # user-supplied foldAssignments
    if (!(k==length(unique(foldAssignments)))) stop("Value of k does not match number of folds indicated in foldAssignments.")
  } else { # create foldAssignments 
  # attempt to include sites with and without observations in each fold
    foldAssignments = c(1:M)
    idxsWithObs = which(rowSums(y)>0)
    idxsWoObs = which(rowSums(y)==0)
    if (length(idxsWithObs)>0 & length(idxsWoObs)>0) {
      foldAssignments[idxsWithObs] = sample(rep(1:k,ceiling(length(idxsWithObs)/k))[1:length(idxsWithObs)])
      foldAssignments[idxsWoObs] = sample(rep(1:k,ceiling(length(idxsWoObs)/k))[1:length(idxsWoObs)])
    } else if (k<=M) {
      foldAssignments = sample(rep(1:k,ceiling(M/k)))[1:M]
    } else {
      stop("k>M. More folds than sites creates folds. Specify a smaller k.")
    }     
  }
  #print(foldAssignments)
  foldNames = unique(foldAssignments)

  if(identical(engine, "C")) {
    nll <- function(params) {
          beta.psi <- params[1:nOP]
          beta.p <- params[(nOP+1):nP]
          .Call("nll_occu",
                  yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec,
                  X.offset, V.offset, 
                  PACKAGE = "unmarked")
      }
  } else {
    nll <- function(params) { # penalize this function
        psi <- plogis(X %*% params[1 : nOP] + X.offset)
        psi[knownOccLog] <- 1
        pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1 # so that NA's don't modify likelihood
        cpmat <- matrix(cp, M, J, byrow = TRUE) #
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        -sum(loglik)
    }
  } # end if (engine)
 
  lambdaScores = lambdaVec*0 # score by held-out likelihood

  for (f in 1:k) {				  
    fold = foldNames[f]
    occuTrain = data[which(foldAssignments!=fold),] # train on NOT this fold
    occuTest = data[which(foldAssignments==fold),] # test on this fold
    
    designMats <- getDesign(occuTest, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    removed <- designMats$removed.sites
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if(is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    y <- truncateToBinary(y)
    J <- ncol(y)
    M <- nrow(y)

    ## convert knownOcc to logical so we can correctly to handle NAs.
    knownOccLog <- rep(FALSE, numSites(data))
    knownOccLog[knownOcc] <- TRUE
    if(length(removed)>0)
        knownOccLog <- knownOccLog[-removed]

    occParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nOP <- ncol(X)
    nP <- nDP + nOP

    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))
    if(missing(starts)) starts <- rep(0, nP)

    yvec <- as.numeric(t(y))
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

    # For each lambda, get parameters on the training set, and use them
    #  to compute the likelihood on the held-out test fold.
    for (la in 1:length(lambdaVec)) { 
      occuOut = occuPEN(formula, occuTrain, starts, lambda=lambdaVec[la],pen.type=pen.type) 
      ests = c(as.numeric(occuOut[1]@estimates),as.numeric(occuOut[2]@estimates))
      lambdaScores[la] = lambdaScores[la] + nll(ests)
    } # la
  } # f

  bestLambda = lambdaVec[which.min(lambdaScores)]
  #print(lambdaScores)

  occuOut = occuPEN(formula, data, starts=starts, lambda=bestLambda, pen.type=pen.type) 
 
  umfit <- new("unmarkedFitOccuPEN_CV", fitType = "occu", call = match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = occuOut@estimates, AIC = occuOut@AIC, 
		 opt = occuOut@opt,
                 negLogLike = occuOut@negLogLike,
                 nllFun = occuOut@nllFun, knownOcc = knownOccLog, 
		 pen.type = pen.type, lambdaVec = lambdaVec,
		 k = k, foldAssignments = foldAssignments,
		 lambdaScores = lambdaScores, chosenLambda = bestLambda)

  return(umfit)

} # fn: occuPEN_CV

occuPEN <- function(formula, data, knownOcc = numeric(0), starts,
                 method = "BFGS", engine = c("C", "R"),
#		 se = TRUE, 
		 lambda = 0, 
		 pen.type = c("Bayes","Ridge","MPLE"),
		 ...) 
 {
    if(!is(data, "unmarkedFrameOccu"))
        stop("Data is not an unmarkedFrameOccu object.")

    pen.type = pen.type[1]
    if (!(pen.type=="Bayes" | pen.type=="Ridge" | pen.type=="MPLE")) 
        stop("pen.type not recognized.  Choose Bayes, Ridge, or MPLE.")

    engine <- match.arg(engine, c("C", "R"))

    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y

    if (ncol(X)==1 & pen.type=="MPLE") stop("MPLE requires occupancy covariates.")

    if (ncol(X)==1 & ncol(V)==1 & pen.type=="Ridge") stop("Ridge requires covariates.")

    removed <- designMats$removed.sites
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if(is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    y <- truncateToBinary(y)
    J <- ncol(y)
    M <- nrow(y)

    ## convert knownOcc to logical so we can correctly to handle NAs.
    knownOccLog <- rep(FALSE, numSites(data))
    knownOccLog[knownOcc] <- TRUE
    if(length(removed)>0)
        knownOccLog <- knownOccLog[-removed]

    occParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nOP <- ncol(X)

    nP <- nDP + nOP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))
    if(missing(starts)) starts <- rep(0, nP)

    yvec <- as.numeric(t(y))
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

    ## need to add offsets !!!!!!!!!!!!!!
    ## and fix bug causing crash when NAs are in V

    ## compute logistic regression MPLE targets and lambda:
    if (pen.type=="MPLE") {
      LRparams = glm.fit(x=X,y=apply(y,1,max),family=binomial(),intercept=F,start=starts[1:nOP])
      MPLElambda = computeMPLElambda(formula, data, knownOcc = numeric(0), starts, method = "BFGS", engine = c("C", "R"))
      if (MPLElambda != lambda) warning("Supplied lambda does not match the computed value. Proceeding with the supplied lambda.")
    }
    
    if(identical(engine, "C")) {
        nll <- function(params) {
            beta.psi <- params[1:nOP]
            beta.p <- params[(nOP+1):nP]
  	    if (pen.type=="Bayes") {
	      penalty = sum(params^2)*lambda*0.5
  	    } else if (pen.type=="Ridge") {
	      penalty = 0
	      if (nOP>1) { penalty = penalty + sum((params[2:nOP])^2) }
	      if (nDP>1) { penalty = penalty + sum((params[(nOP+2):nP])^2) }
	      penalty = penalty*lambda*0.5
	    } else if (pen.type=="MPLE") {
	      penalty = abs(params[1:nOP]-LRparams$coefficients)
	      penalty = sum(penalty)*lambda
	    } else {
	      stop("pen.type not found")
	    }
	
            .Call("nll_occuPEN",
                  yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec,
                  X.offset, V.offset, penalty,
                  PACKAGE = "unmarked")
        }
    } else {
      nll <- function(params) { # penalize this function
          psi <- plogis(X %*% params[1 : nOP] + X.offset)
          psi[knownOccLog] <- 1
          pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
          cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
          cp[navec] <- 1 # so that NA's don't modify likelihood
          cpmat <- matrix(cp, M, J, byrow = TRUE) #
          loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
          #-sum(loglik)
 
  	  if (pen.type=="Bayes") {
	    penalty = sum(params^2)*lambda*0.5
  	  } else if (pen.type=="Ridge") {
	    penalty = 0
	    if (nOP>1) { penalty = penalty + sum((params[2:nOP])^2) }
	    if (nDP>1) { penalty = penalty + sum((params[(nOP+2):nP])^2) }
	    penalty = penalty*lambda*0.5
	  } else if (pen.type=="MPLE") {
	    penalty = abs(params[1:nOP]-LRparams$coefficients)
	    penalty = sum(penalty)*lambda
	  } else {
	    stop("pen.type not found")
	  }
	  penLL = sum(loglik) - penalty
	  return(-penLL)
        }
    } # end if (engine)

    fm <- optim(starts, nll, method = method, hessian = FALSE, ...)
    opt <- fm
    covMat <- matrix(NA, nP, nP)

    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
    names(ests) <- c(occParms, detParms)

    state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                              estimates = ests[1:nOP],
                              covMat = as.matrix(covMat[1:nOP,1:nOP]),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad")

    det <- unmarkedEstimate(name = "Detection", short.name = "p",
                            estimates = ests[(nOP + 1) : nP],
                            covMat = as.matrix(covMat[(nOP + 1) : nP,
                                                     (nOP + 1) : nP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    estimateList <- unmarkedEstimateList(list(state=state, det=det))

    umfit <- new("unmarkedFitOccuPEN", fitType = "occu", call = match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = opt,
                 negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog, 
		 pen.type = pen.type, lambda = c(lambda))

    return(umfit)
}
