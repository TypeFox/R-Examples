geem <- function(formula, id, waves=NULL, data = parent.frame(), family = gaussian, corstr = "independence", Mv = 1, weights = NULL, corr.mat = NULL, init.beta=NULL, init.alpha=NULL, init.phi = 1, scale.fix=FALSE, nodummy=FALSE, sandwich=TRUE, maxit=20, tol=0.00001){
  call <- match.call()
  
  famret <- getfam(family)
  
  if(inherits(famret, "family")){
    LinkFun <- famret$linkfun
    InvLink <- famret$linkinv
    VarFun <- famret$variance
    InvLinkDeriv <- famret$mu.eta
  }else{
    LinkFun <- famret$LinkFun
    VarFun <- famret$VarFun
    InvLink <- famret$InvLink
    InvLinkDeriv <- famret$InvLinkDeriv
  }  
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  ### First, get all the relevant elements from the arguments
  dat <- model.frame(formula, data, na.action=na.pass)
  nn <- dim(dat)[1]
  
  if(typeof(data) == "environment"){
    id <- id
    weights <- weights
    if(is.null(call$weights)) weights <- rep(1, nn)
    waves <- waves
  }
  else{
    if(length(call$id) == 1){
      subj.col <- which(colnames(data) == call$id)  
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      }else{
        id <- eval(call$id, envir=parent.frame())
      }
    }else if(is.null(call$id)){
      id <- 1:nn
    }
    
    if(length(call$weights) == 1){
      weights.col <- which(colnames(data) == call$weights)  
      if(length(weights.col) > 0){
        weights <- data[,weights.col]
      }else{
        weights <- eval(call$weights, envir=parent.frame())
      }
    }else if(is.null(call$weights)){
      weights <- rep.int(1,nn)
    }
    
    if(length(call$waves) == 1){
      waves.col <- which(colnames(data) == call$waves)  
      if(length(waves.col) > 0){
        waves <- data[,waves.col]
      }else{
        waves <- eval(call$waves, envir=parent.frame())
      }
    }else if(is.null(call$waves)){
      waves <- NULL
    }    
  }
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves
  
  if(!is.numeric(dat$waves) & !is.null(dat$waves)) stop("waves must be either an integer vector or NULL")
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  # Note that we need to assign weight 0 to rows with NAs 
  # in order to preserve the correlation structure
  na.inds <- NULL
  
  if(any(is.na(dat))){
    na.inds <- which(is.na(dat), arr.ind=T)
  }
  
  #SORT THE DATA ACCORDING TO WAVES
  if(!is.null(waves)){
    dat <- dat[order(id, waves),]
  }
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(corstr, cor.vec)
  
  if(is.na(cor.match)){stop("Unsupported correlation structure")}  
  
  if(!is.null(dat$waves)){
    wavespl <- split(dat$waves, dat$id)
    idspl <- split(dat$id, dat$id)
    
    maxwave <- rep(0, length(wavespl))
    incomp <- rep(0, length(wavespl))
    
    for(i in 1:length(wavespl)){
      maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 1
      if(maxwave[i] != length(wavespl[[i]])){
        incomp[i] <- 1
      }
    }
  
    #If there are gaps and correlation isn't exchangeable or independent
    #then we'll add some dummy rows
    if(!is.element(cor.match, c(1,3)) & (sum(incomp) > 0) & !nodummy){
      dat <- dummyrows(formula, dat, incomp, maxwave, wavespl, idspl)
      id <- dat$id
      waves <- dat$waves
      weights <- dat$weights
    }
  }

  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      if(is.factor(dat[,i])){
        dat[na.inds[,1], i] <- levels(dat[,i])[1]  
      }else{
        dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
      }
    }
  }
  
  
  includedvec <- weights>0
  
  
  inclsplit <- split(includedvec, id)
  
  dropid <- NULL
  allobs <- T
  if(any(!includedvec)){
    allobs <- F
    for(i in 1:length(unique(id))){
      if(all(!inclsplit[[i]])){
        dropid <- c(dropid, i)
      }    
    }
  }
  
  
  if(length(dropid)>0){
    dropind <- which(is.element(id, dropid))
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- id[-dropind]
  }
  nn <- dim(dat)[1]
  K <- length(unique(id))

  
  modterms <- terms(formula)
  
  X <- model.matrix(formula,dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)			
  
  p <- dim(X)[2]
  
  
  
  ### if no offset is given, then set to zero
  if(is.null(offset)){
    off <- rep(0, nn)
  }else{
    off <- offset
  }
  
  # Is there an intercept column?
  interceptcol <- apply(X==1, 2, all)
  
  ## Basic check to see if link and variance functions make any kind of sense
  linkOfMean <- LinkFun(mean(Y))
  if( any(is.infinite(linkOfMean) | is.nan(linkOfMean)) ){
    stop("Infinite or NaN in the link of the mean of responses.  Make sure link function makes sense for these data.")
  }
  if( any(is.infinite( VarFun(mean(Y))) | is.nan( VarFun(mean(Y)))) ){
    stop("Infinite or NaN in the variance of the mean of responses.  Make sure variance function makes sense for these data.")
  }
  
  if(is.null(init.beta)){
    if(any(interceptcol)){
      #if there is an intercept and no initial beta, then use link of mean of response
      init.beta <- rep(0, dim(X)[2])
      init.beta[which(interceptcol)] <- linkOfMean
    }else{
      stop("Must supply an initial beta if not using an intercept.")
    }
  }
  
  
  # Number of included observations for each cluster
  includedlen <- rep(0, K)
  len <- rep(0,K)
  uniqueid <- unique(id)
  
  tmpwgt <- as.numeric(includedvec)
  idspl <-ifelse(tmpwgt==0, NA, id)
  includedlen <- as.numeric(summary(split(Y, idspl, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
  # Get vector of cluster sizes... remember this len variable
  #len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  

  # Set the initial alpha value
  if(is.null(init.alpha)){
    alpha.new <- 0.2
    if(cor.match==4){
      # If corstr = "m-dep"
      alpha.new <- 0.2^(1:Mv)
    }else if(cor.match==5){
      # If corstr = "unstructured"
      alpha.new <- rep(0.2, sum(1:(max(len)-1)))
    }else if(cor.match==7){
      # If corstr = "userdefined"
      alpha.new <- rep(0.2, max(unique(as.vector(corr.mat))))
    }
  }else{
    alpha.new <- init.alpha	
  }
  #if no initial overdispersion parameter, start at 1
  if(is.null(init.phi)){
    phi <- 1
  }else{
    phi <- init.phi
  }
  
  beta <- init.beta
  

  
  #Set up matrix storage
  StdErr <- Diagonal(nn)
  dInvLinkdEta <- Diagonal(nn)
  Resid <- Diagonal(nn)
  
  
  # Initialize for each correlation structure
  if(cor.match == 1){
    # INDEPENDENCE
    R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
    BlockDiag <- getBlockDiag(len)$BDiag
  }else if(cor.match == 2){
    # AR-1
    tmp <- buildAlphaInvAR(len)
    # These are the vectors needed to update the inverse correlation
    a1<- tmp$a1
    a2 <- tmp$a2
    a3 <- tmp$a3
    a4 <- tmp$a4
    # row.vec and col.vec for the big block diagonal of correlation inverses
    # both are vectors of indices that facilitate in updating R.alpha.inv
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    BlockDiag <- getBlockDiag(len)$BDiag
    
  }else if(cor.match == 3){
    # EXCHANGEABLE
    # Build a block diagonal correlation matrix for updating and sandwich calculation
    # this matrix is block diagonal with all ones.  Each block is of dimension cluster size.
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    
    # Create a vector of length number of observations with associated cluster size for each observation
    n.vec <- vector("numeric", nn)
    index <- c(cumsum(len) - len, nn)
    for(i in 1:K){
      n.vec[(index[i]+1) : index[i+1]] <-  rep(len[i], len[i])
    }
  }else if(cor.match == 4){
    # M-DEPENDENT, check that M is not too large
    if(Mv >= max(len)){
      stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
    }
    
    # Build block diagonal similar to in exchangeable case, also get row indices and column
    # indices for fast matrix updating later.		
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    R.alpha.inv <- NULL
  }else if(cor.match == 5){
    # UNSTRUCTURED
    if( max(len^2 - len)/2 > length(len)){
      stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
    }
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec	
  }else if(cor.match == 6){
    # FIXED
    # check if matrix meets some basic conditions
    corr.mat <- checkFixedMat(corr.mat, len)
    
    R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
    BlockDiag <- getBlockDiag(len)$BDiag
  }else if(cor.match == 7){
    # USERDEFINED
    corr.mat <- checkUserMat(corr.mat, len)
    
    # get the structure of the correlation matrix in a way that
    # I can use later on.
    tmp1 <- getUserStructure(corr.mat)
    corr.list <- tmp1$corr.list
    user.row <- tmp1$row.vec
    user.col <- tmp1$col.vec
    struct.vec <- tmp1$struct.vec
    
    # the same block diagonal trick.
    tmp2 <- getBlockDiag(len)
    BlockDiag <- tmp2$BDiag
    row.vec <- tmp2$row.vec
    col.vec <- tmp2$col.vec
    
  }else if(cor.match == 0){
    stop("Ambiguous Correlation Structure Specification")
  }else{
    stop("Unsupported Correlation Structure")	
  }
  
  stop <- F
  converged <- F	
  count <- 0
  beta.old <- beta
  unstable <- F
  phi.old <- phi
  
  
  # Main fisher scoring loop
  while(!stop){		
    count <- count+1
    
    eta <- as.vector(X %*% beta) + off
    
    mu <- InvLink(eta)
    
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    if(!scale.fix){
      phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, includedlen, sqrtW)
    }
    phi.new <- phi
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(cor.match == 2){
      # AR-1
      alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs, sqrtW, BlockDiag)
      R.alpha.inv <- getAlphaInvAR(alpha.new, a1,a2,a3,a4, row.vec, col.vec)/phi
    }else if(cor.match == 3){
      #EXCHANGEABLE		  
      alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen, sqrtW)
      R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
    }else if(cor.match == 4){
      # M-DEPENDENT
      if(Mv==1){
        alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs, sqrtW, BlockDiag)
      }else{
        alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, Mv, included, includedlen, allobs, sqrtW)
        if(sum(len>Mv) <= p){
          unstable <- T
        }
      }
      if(any(alpha.new >= 1)){
        stop <- T
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, col.vec)/phi		
    }else if(cor.match == 5){
      # UNSTRUCTURED
      alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, id, len, StdErr, Resid,  p, BlockDiag, included, includedlen, allobs, sqrtW)
      # This has happened to me (greater than 1 correlation estimate)
      if(any(alpha.new >= 1)){
        stop <- T
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, row.vec, col.vec)/phi
    }else if(cor.match ==6){
      # FIXED CORRELATION, DON'T NEED TO RECOMPUTE
      R.alpha.inv <- R.alpha.inv*phi.old/phi
    }else if(cor.match == 7){
      # USER SPECIFIED
      alpha.new <- updateAlphaUser(Y, mu, phi, id, len, StdErr, Resid, p, BlockDiag, user.row, user.col, corr.list, included, includedlen, allobs, sqrtW)
      R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
    }else if(cor.match == 1){
      # INDEPENDENT
      R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
      alpha.new <- "independent"
    }
    
    
    beta.list <- updateBeta(Y, X, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, sqrtW)	
    beta <- beta.list$beta
    phi.old <- phi
    if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){converged <- T; stop <- T}
    if(count >= maxit){stop <- T}		
    beta.old <- beta		
  }
  biggest <- which.max(len)[1]
  index <- cumsum(len[biggest])
  if(K == 1){
    biggest.R.alpha.inv <- R.alpha.inv
    if(cor.match == 6) {
      biggest.R.alpha <- corr.mat
    }else{
      biggest.R.alpha <- solve(R.alpha.inv)
    }
  }else{
    biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
    if(cor.match == 6){
      biggest.R.alpha <- corr.mat[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
    }else{
      biggest.R.alpha <- solve(biggest.R.alpha.inv)
    }
  }

  
  eta <- as.vector(X %*% beta) + off
  if(sandwich){
    sandvar.list <- getSandwich(Y, X, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, beta.list$hess, StdErr, dInvLinkdEta, BlockDiag, sqrtW)
  }else{
    sandvar.list <- list()
    sandvar.list$sandvar <- "no sandwich"
  }
  
  if(!converged){warning("Did not converge")}
  if(unstable){warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")}
  
  
  # Create object of class geem with information about the fit
  dat <- model.frame(formula, data, na.action=na.pass)
  X <- model.matrix(formula, dat)
  
  if(is.character(alpha.new)){alpha.new <- 0}
  results <- list()
  results$beta <- as.vector(beta)
  results$phi <- phi
  results$alpha <- alpha.new
  if(cor.match == 6){
    results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat,1)!=0)])
  }
  results$coefnames <- colnames(X)
  results$niter <- count
  results$converged <- converged
  results$naiv.var <- solve(beta.list$hess)*phi.new  ## call model-based
  results$var <- sandvar.list$sandvar
  results$call <- call
  results$corr <- cor.vec[cor.match]
  results$clusz <- len
  results$FunList <- famret
  results$X <- X
  results$offset <- off
  results$eta <- eta
  results$dropped <- dropid
  results$weights <- weights
  results$y <- Y
  results$biggest.R.alpha <- biggest.R.alpha
  class(results) <- "geem"
  return(results)
}






### Simple moment estimator of dispersion parameter
updatePhi <- function(YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW){
  nn <- sum(includedlen)
  resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu))
  phi <- (1/(sum(included)-p))*crossprod(resid, resid) 
  return(as.numeric(phi))	
}

### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
updateBeta = function(YY, XX, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, sqrtW){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(XX%*%beta.new) + off
    
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)	
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*%XX, R.alpha.inv %*% sqrtW %*% StdErr %*%dInvLinkdEta %*% XX)
    esteq <- crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%XX , R.alpha.inv %*% sqrtW %*% StdErr %*% (YY - mu))
    
    update <- solve(hess, esteq)
    if(max(abs(update)/beta.new) < 100*tol){break}
    
    beta.new <- beta.new + as.vector(update)
  }
  return(list(beta = beta.new, hess = hess))
}

### Calculate the sandiwch estimator as usual.
getSandwich = function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag, sqrtW){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)			
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% sqrtW %*% StdErr %*% BlockDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% XX)
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}

