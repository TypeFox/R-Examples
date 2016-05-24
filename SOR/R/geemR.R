geemR <- function(formula, id,data = parent.frame(), family = gaussian, corstr = "independence", Mv = 1, weights = NULL, corr.mat = NULL, init.beta=NULL, init.alpha=NULL, init.phi = NULL, scale.fix=FALSE, sandwich=TRUE, maxit=20, tol=0.00001){
  call <- match.call()
  
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))  
  }
  sandwich = F
  if(is.function(family)){
    family <- family()
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family) && !is.null(family$family)){
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family)){
    if(length(match(names(family), c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv"))) == 4){
      LinkFun <- family$LinkFun
      InvLink <- family$InvLink
      VarFun <- family$VarFun
      InvLinkDeriv <- family$InvLinkDeriv
    }else{
      LinkFun <- family[[1]]
      VarFun <- family[[2]]
      InvLink <- family[[3]]
      InvLinkDeriv <- family[[4]]
    }
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  ### First, get all the relevant elements from the arguments
  if(typeof(data) == "environment"){id = id
                                    data <- model.frame(formula, data=data)
  }
  else{
    subj.col <- which(colnames(data) == call$id)  
    id <- data[,subj.col]
  }
  nn <- length(id)
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  na.inds <- NULL
  
  if(any(is.na(data))){
    na.inds <- which(is.na(data), arr.ind=T)
  }
  
  if(is.null(weights)){
    weights <- rep.int(1, nn)
  }
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      data[na.inds[,1], i] <- median(data[,i], na.rm=T)
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
    data <- data[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- data[,subj.col]
  }
  nn <- dim(data)[1]
  K <- length(unique(id))
  
  modterms <- terms(formula)
  modframe <- model.frame(formula, data)
  
  X <- model.matrix(formula,modframe)
  Y <- model.response(modframe)
  offset <- model.offset(modframe)  		
  
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
  tmpwgt[tmpwgt==0] <- NA
  includedlen <- as.numeric(summary(split(Y, id*tmpwgt, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
  # Get vector of cluster sizes... remember this len variable
  #len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(corstr, cor.vec)
  
  if(is.na(cor.match)){stop("Unsupported correlation structure")}  
  
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
    
    diag(StdErr) <- sqrt(1/VarFun(eta))
    
    if(!scale.fix){
      phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, includedlen)
    }
    phi.new <- phi
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(cor.match == 2){
      # AR-1
      alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
      R.alpha.inv <- updateAlphaInvAR(alpha.new, a1,a2,a3,a4, row.vec, col.vec)/phi
    }else if(cor.match == 3){
      #EXCHANGEABLE		  
      alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen)
      R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
    }else if(cor.match == 4){
      # M-DEPENDENT
      if(Mv==1){
        alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
      }else{
        alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, Mv, included, includedlen, allobs)
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
      alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, id, len, StdErr, Resid,  p, BlockDiag, included, includedlen, allobs)
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
      alpha.new <- updateAlphaUser(Y, mu, phi, id, len, StdErr, Resid, p, BlockDiag, user.row, user.col, corr.list, included, includedlen, allobs)
      R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
    }else if(cor.match == 1){
      # INDEPENDENT
      R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
      alpha.new <- "independent"
    }
    
    
    beta.list <- updateBetaR(Y, X, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, sqrtW)	
    beta <- beta.list$beta
    phi.old <- phi
    if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){converged <- T; stop <- T}
    if(count >= maxit){stop <- T}		
    beta.old <- beta		
  }
  biggest <- which.max(len)[1]
  index <- cumsum(len[biggest])
  biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
  eta <- as.vector(X %*% beta) + off
  if(sandwich){
    sandvar.list <- getSandwichR(Y, X, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, beta.list$hess, StdErr, dInvLinkdEta, BlockDiag, sqrtW)
  }else{
    sandvar.list <- list()
    sandvar.list$sandvar <- "no sandwich"
  }
  
  if(!converged){warning("Did not converge")}
  if(unstable){warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")}
  
  
  # Create object of class geem with information about the fit
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
  results$FunList <- FunList
  results$R.a.inv <- R.alpha.inv
  results$X <- X
  results$offset <- off
  results$eta <- eta
  results$dropped <- dropid
  class(results) <- "geem"
  return(results)
}

### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
updateBetaR = function(YY, XX, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, sqrtW){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(XX%*%beta.new) + off
    
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)  
    diag(StdErr) <- sqrt(1/VarFun(eta))
    
    hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*%XX, R.alpha.inv %*% sqrtW %*% StdErr %*%dInvLinkdEta %*% XX)
    esteq <- crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%XX , R.alpha.inv %*% sqrtW %*% StdErr %*% (YY - mu))
    
    update <- solve(hess, esteq)
    if(max(abs(update)/beta.new) < 100*tol){break}
    
    beta.new <- beta.new + as.vector(update)
  }
  return(list(beta = beta.new, hess = hess))
}

### Calculate the sandiwch estimator as usual.
getSandwichR = function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag, sqrtW){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)			
  diag(StdErr) <- sqrt(1/VarFun(eta))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% sqrtW %*% StdErr %*% BlockDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% XX)
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}


geem <- function(formula, id,data = parent.frame(), family = gaussian, corstr = "independence", Mv = 1, weights = NULL, corr.mat = NULL, init.beta=NULL, init.alpha=NULL, init.phi = NULL, scale.fix=FALSE, sandwich=TRUE, maxit=20, tol=0.00001){
  call <- match.call()
  
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))  
  }
  
  if(is.function(family)){
    family <- family()
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family) && !is.null(family$family)){
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family)){
    if(length(match(names(family), c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv"))) == 4){
      LinkFun <- family$LinkFun
      InvLink <- family$InvLink
      VarFun <- family$VarFun
      InvLinkDeriv <- family$InvLinkDeriv
    }else{
      LinkFun <- family[[1]]
      VarFun <- family[[2]]
      InvLink <- family[[3]]
      InvLinkDeriv <- family[[4]]
    }
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  ### First, get all the relevant elements from the arguments
  if(typeof(data) == "environment"){id = id
    data <- model.frame(formula, data=data)
  }
  else{
    subj.col <- which(colnames(data) == call$id)  
    id <- data[,subj.col]
  }
  nn <- length(id)
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  na.inds <- NULL
  
  if(any(is.na(data))){
    na.inds <- which(is.na(data), arr.ind=T)
  }
  
  if(is.null(weights)){
    weights <- rep.int(1, nn)
  }
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      data[na.inds[,1], i] <- median(data[,i], na.rm=T)
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
    data <- data[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- data[,subj.col]
  }
  nn <- dim(data)[1]
  K <- length(unique(id))
  
  modterms <- terms(formula)
  modframe <- model.frame(formula, data)
  
  X <- model.matrix(formula,modframe)
  Y <- model.response(modframe)
  offset <- model.offset(modframe)			
  
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
  tmpwgt[tmpwgt==0] <- NA
  includedlen <- as.numeric(summary(split(Y, id*tmpwgt, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
  # Get vector of cluster sizes... remember this len variable
  #len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(corstr, cor.vec)
  
  if(is.na(cor.match)){stop("Unsupported correlation structure")}  
  
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
      phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, includedlen)
    }
    phi.new <- phi
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(cor.match == 2){
      # AR-1
      alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
      R.alpha.inv <- updateAlphaInvAR(alpha.new, a1,a2,a3,a4, row.vec, col.vec)/phi
    }else if(cor.match == 3){
      #EXCHANGEABLE		  
      alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen)
      R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
    }else if(cor.match == 4){
      # M-DEPENDENT
      if(Mv==1){
        alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
      }else{
        alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, Mv, included, includedlen, allobs)
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
      alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, id, len, StdErr, Resid,  p, BlockDiag, included, includedlen, allobs)
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
      alpha.new <- updateAlphaUser(Y, mu, phi, id, len, StdErr, Resid, p, BlockDiag, user.row, user.col, corr.list, included, includedlen, allobs)
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
  biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
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
  results$FunList <- FunList
  results$R.a.inv <- R.alpha.inv
  results$X <- X
  results$offset <- off
  results$eta <- eta
  results$dropped <- dropid
  class(results) <- "geem"
  return(results)
}


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


### Get the necessary structures for the inverse matrix of AR-1
### returns the row and column indices of the AR-1 inverse
### a1, a2, a3, and a4 are used to compute the entries in the matrix
buildAlphaInvAR <- function(len){
  
  nn <- sum(len)
  K <- length(len)
  a1 <- a2 <- a3 <- a4 <- vector("numeric", nn)
  index <- c(cumsum(len) - len, nn)
  
  for (i in 1:K)  {
    if(len[i] > 1)  {
      a1[(index[i]+1) : index[i+1]] <- c(rep(-1,times = len[i]-1),0)
      a2[(index[i]+1) : index[i+1]] <- c(0,rep(1,times=len[i]-2),0)
      a3[(index[i]+1) : index[i+1]] <- c(1,rep(0,times=len[i]-2),1)
      a4[(index[i]+1) : index[i+1]] <- c(rep(0,times=len[i]))
    }
    else if (len[i] == 1)  {
      a1[(index[i]+1) : index[i+1]] <- 0
      a2[(index[i]+1) : index[i+1]] <- 0
      a3[(index[i]+1) : index[i+1]] <- 0
      a4[(index[i]+1) : index[i+1]] <- 1
    }
  }
  
  a1 <- a1[1:(nn-1)]
  subdiag.col <- 1:(nn-1)
  subdiag.row <- 2:nn
  
  row.vec <- c(subdiag.row, (1:nn), subdiag.col)
  col.vec <- c(subdiag.col, (1:nn), subdiag.row)
  return(list(row.vec = row.vec, col.vec= col.vec, a1 = a1, a2=a2, a3=a3, a4=a4))
}

### Get a vector with the entries of the AR-1 inverse matrix
updateAlphaInvAR <- function(alpha.new, a1,a2,a3,a4, row.vec, col.vec){
  corr.vec <- c(alpha.new*a1/(1-alpha.new^2) , ( (1+alpha.new^2)*a2 + a3)/(1-alpha.new^2) + a4, alpha.new*a1/(1-alpha.new^2))	
  return(as(sparseMatrix(i= row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}

### Returns the full inverse matrix of the correlation for EXCHANGEABLE structure
getAlphaInvEX <- function(alpha.new, diag.vec, BlockDiag){
  return(as(BlockDiag %*% Diagonal(x = (-alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)))) + Diagonal( x = ((1+(diag.vec-2)*alpha.new)/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)) + alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)))), "symmetricMatrix"))
}

### Calculate the parameter for the EXCHANGEABLE correlation structure
updateAlphaEX <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen){
  
  Resid <- StdErr %*% included %*% Diagonal(x = YY - mu)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  denom <-  phi*(crossprod(includedlen, pmax(includedlen-1, 0))/2 - p)
  alpha <- (sum(BlockDiag) - phi*(sum(includedlen)-p))/2
  alpha.new <- alpha/denom
  return(alpha.new)
}

### Calculate the parameters for the M-DEPENDENT correlation structure
updateAlphaMDEP <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, m, included, includedlen, allobs){
  
  Resid <- StdErr %*% included %*% Diagonal(x = YY - mu)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", m)
  for(i in 1:m){	
    if(sum(includedlen>i) > p){
      bandmat <- drop0(band(BlockDiag, i,i))
      
      if(allobs){alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i))-p))
      }else{alpha.new[i] <- sum( bandmat)/(phi*(length(bandmat@i)-p))}
      
    }else{
      # If we don't have many observations for a certain parameter, don't divide by p
      # ensures we don't have NaN errors.
      bandmat <- drop0(band(BlockDiag, i,i))
      
      if(allobs){alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i))))
      }else{alpha.new[i] <- sum( bandmat)/(phi*length(bandmat@i))}
      
    }
  }
  return(alpha.new)
}

### Calculate the parameter for the AR-1 correlation, also used for 1-DEPENDENT
updateAlphaAR <- function(YY, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs){
  K <- length(len)
  oneobs <- which(len == 1)
  
  resid <- diag(StdErr %*% included %*% Diagonal(x = YY - mu))
  
  len2 = len
  includedvec2 <- includedvec
  if(length(oneobs) > 0){
    index <- c(0, (cumsum(len) -len)[2:K], sum(len))
    len2 <- len[-oneobs]
    resid <- resid[-index[oneobs]]
    includedvec2 <- includedvec[-index[oneobs]]
  }
  
  
  nn <- length(resid)
  lastobs <- cumsum(len2)
  
  shiftresid1 <- resid[1:nn-1]
  shiftresid2 <- resid[2:nn]
  if(!allobs){
    shiftresid1 <- shiftresid1[-lastobs]
    shiftresid2 <- shiftresid2[-lastobs]
    s1incvec2 <- includedvec2[1:nn-1]
    s2incvec2 <- includedvec2[2:nn]
    s1incvec2 <- s1incvec2[-lastobs]
    s2incvec2 <- s2incvec2[-lastobs]
    
    alphasum <- crossprod(shiftresid1, shiftresid2)
    
    denom <- (as.vector(crossprod(s1incvec2, s2incvec2)) - p)*phi
  }else{
    alphasum <- crossprod(shiftresid1[-(cumsum(len2))], shiftresid2[-(cumsum(len2))])
    denom <- (sum(len2-1) - p)*phi
  }
  
  alpha <- alphasum/denom
  return(as.numeric(alpha))
}

### Simple moment estimator of dispersion parameter
updatePhi <- function(YY, mu, VarFun, p, StdErr, included, includedlen){
  nn <- sum(includedlen)
  resid <- diag(StdErr %*% included %*% Diagonal(x = YY - mu))
  phi <- (1/(nn-p))*crossprod(resid, resid) 
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

### Get the inverse of M-DEPENDENT correlation matrix.
getAlphaInvMDEP <- function(alpha.new, len, row.vec, col.vec){
  K <- length(len)
  N <- sum(len)
  m <- length(alpha.new)
  
  # First get all of the unique block sizes.
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sum(len^2))
  mat.inverses <- list()
  index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  
  for(i in 1:length(mat.sizes)){
    # Now create and invert each matrix
    if(mat.sizes[i] == 1){
      mat.inverses[[i]] <- 1
    }else{		
      mtmp <- min(m, mat.sizes[i]-1)
      a1 = list()
      a1[[1]] <- rep(1, mat.sizes[i])
      for(j in 1:mtmp){
        a1[[j+1]] <- rep(alpha.new[j], mat.sizes[i]-j)
      }
      
      tmp <- bandSparse(mat.sizes[i], k=c(0:mtmp),diagonals=a1, symmetric=T )
      mat.inverses[[i]] <- as.vector(solve(tmp))
    }
  }
  # Put all inverted matrices in a vector in the right order
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}

### Get a vector of elements of the inverse correlation matrix for UNSTRUCTURED
### Inversion strategy follows the same basic idea as M-DEPENDENT
getAlphaInvUnstruc <- function(alpha.new, len, row.vec, col.vec){
  K <- length(len)
  unstr.row <- NULL
  unstr.col <- NULL
  ml <- max(len)
  sl2 <- sum(len^2)
  for(i in 2:ml){
    unstr.row <- c(unstr.row, 1:(i-1))
    unstr.col <- c(unstr.col, rep(i, each=i-1))
  }
  unstr.row <- c(unstr.row, 1:ml)
  unstr.col <- c(unstr.col, 1:ml)
  xvec <- c(alpha.new, rep(1, ml))
  # Get the biggest matrix implied by the cluster sizes
  biggestMat <- forceSymmetric(sparseMatrix(i=unstr.row, j=unstr.col, x=xvec))
  
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sl2)
  mat.inverses <- list()
  index <- vector("numeric", K+1)
  index[1] <- 0
  index[2:K] <-  (cumsum(len^2) -len^2)[2:K]
  index[K+1] <- sl2
  
  for(i in 1:length(mat.sizes)){
    tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))		
  }
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
  
}

### Calculate alpha values for UNSTRUCTURED correlation
updateAlphaUnstruc <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen, allobs){
  
  Resid <- StdErr %*% included %*% Diagonal(x = YY - mu)
  
  ml <- max(len)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", sum(1:(ml-1)))
  lalph <- length(alpha.new)
  
  row.vec <- NULL
  col.vec <- NULL
  for(i in 2:ml){
    row.vec <- c(row.vec, 1:(i-1))
    col.vec <- c(col.vec, rep(i, each=i-1))
  }
  index <- cumsum(len)-len
  if(sum(includedlen == max(len)) <= p){stop("Number of clusters of largest size is less than p.")}
  for(i in 1:lalph){
    # Get all of the indices of the matrix corresponding to the correlation
    # we want to estimate.
    newrow <- index[which(len>=col.vec[i])] + row.vec[i]
    newcol <- index[which(len>=col.vec[i])] + col.vec[i]
    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- (phi*(length(newrow)-p))
    }else{denom <- (phi*(sum(bdtmp!=0)-p))}
    alpha.new[i] <- sum(bdtmp)/denom
  }
  
  return(alpha.new)
}

### Invert the FIXED correlation structure.  Again,
### uses same basic technique as M-DEPENDENT
getAlphaInvFixed <- function(mat, len){
  K <- length(len)
  mat.sizes <- sort(unique(len))
  mat.inverses <- list()
  sl2 <- sum(len^2)	
  corr.vec <- vector("numeric", sl2)
  index <- vector("numeric", K+1)
  index[1] <- 0
  index[2:K] <-  (cumsum(len^2) -len^2)[2:K]
  index[K+1] <- sl2
  for(i in 1:length(mat.sizes)){
    tmp <- mat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))
  }
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  
  return(as(getBlockDiag(len, corr.vec)$BDiag, "symmetricMatrix"))	
}

### Get a block diagonal matrix. Each block has dimension corresponding to
### each cluster size.  By default, each block is just a matrix filled with ones.
getBlockDiag <- function(len, xvec=NULL){
  K <- length(len)
  if(is.null(xvec)){
    xvec <- rep.int(1, sum(len^2))
  }
  
  row.vec <- col.vec <- vector("numeric", sum(len^2))
  add.vec <- cumsum(len) - len
  index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  for(i in 1:K){
    row.vec[(index[i] + 1):(index[i+1])] <- rep.int( (1:len[i]) + add.vec[i], len[i])
    col.vec[(index[i] + 1):(index[i+1])] <- rep( (1:len[i]) + add.vec[i], each=len[i])
  }	
  BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
  return(list(BDiag = as(BlockDiag, "symmetricMatrix"), row.vec =row.vec, col.vec=col.vec))
}

### Get the structure of the USERDEFINED correlation matrix implied by the
### corr.mat argument to geem.
getUserStructure <- function(corr.mat){
  ml <- dim(corr.mat)[1]
  
  row.vec <- NULL
  col.vec <- NULL
  for(i in 2:ml){
    row.vec <- c(row.vec, 1:(i-1))
    col.vec <- c(col.vec, rep(i, each=i-1))
  }
  
  struct.vec <- corr.mat[cbind(row.vec, col.vec)]
  
  corr.list <- vector("list", max(struct.vec))
  for(i in 1:max(struct.vec)){
    corr.list[[i]] <- which(struct.vec == i)
  }
  return(list(corr.list = corr.list, row.vec = row.vec, col.vec = col.vec, struct.vec = struct.vec))
}

### Update the alpha (possibly) vector for the USERDEFINED correlation matrix.	
updateAlphaUser <- function(YY, mu, phi, id, len, StdErr, Resid, p, BlockDiag, row.vec, col.vec, corr.list, included, includedlen, allobs){
  Resid <- StdErr %*% included %*% Diagonal(x = YY - mu)
  
  ml <- max(len)
  
  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", length(corr.list))
  
  index <- cumsum(len)-len
  
  for(i in 1:length(alpha.new)){
    newrow <- NULL
    newcol <- NULL
    for(j in 1:length(corr.list[[i]])){
      newrow <- c(newrow, index[which(len >= col.vec[corr.list[[i]]][j])] + row.vec[corr.list[[i]][j]])
      newcol <- c(newcol, index[which(len >= col.vec[corr.list[[i]]][j])] + col.vec[corr.list[[i]][j]])
    }
    
    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- phi*(length(newrow) - p)
    }else{denom <- phi*(sum(bdtmp!=0)-p)}
    alpha.new[i] <- sum(bdtmp)/denom
  }
  return(alpha.new)	
}

### Get the inverse correlation matrix for USERDEFINED.
getAlphaInvUser <- function(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec){
  K <- length(len)
  ml <- max(len)
  sl2 <- sum(len^2)
  
  # Indices for the correlation matrix for the subject
  # with the most observations.
  user.row <- c(user.row, 1:ml)
  user.col <- c(user.col, 1:ml)
  # The entries of the biggest matrix
  xvec <- rep.int(0, length(struct.vec))
  for(i in 1:length(alpha.new)){
    xvec[which(struct.vec == i)] <- alpha.new[i]
  }
  
  xvec <- c(xvec, rep(1, ml))	
  
  biggestMat <- forceSymmetric(sparseMatrix(i=user.row, j=user.col, x=xvec))
  
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sl2)
  mat.inverses <- list()
  
  for(i in 1:length(mat.sizes)){
    tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))		
  }
  
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}

### Check some conditions on the USERDEFINED correlation structure supplied.
checkUserMat <- function(corr.mat, len){
  if(is.null(corr.mat)){
    stop("corr.mat must be specified if using user defined correlation structure")
  }
  if(dim(corr.mat)[1] < max(len)){
    stop("corr.mat needs to be at least as long as the maximum cluster size.")
  }
  test.vec <- as.vector(corr.mat)
  if(any(abs(test.vec-round(test.vec)) > .Machine$double.eps )){
    stop("entries in corr.mat must be integers.")
  }
  max.val <- max(test.vec)
  min.val <- min(test.vec)
  if(!all(sort(unique(test.vec)) == min.val:max.val)){
    stop("entries in corr.mat must be consecutive integers starting at 1.")
  }
  return(corr.mat[1:max(len), 1:max(len)])	
}

### Check some conditions on the FIXED correlation structure.
checkFixedMat <- function(corr.mat, len){
  if(is.null(corr.mat)){
    stop("corr.mat must be specified if using fixed correlation structure")
  }
  if(dim(corr.mat)[1] < max(len)){
    stop("Dimensions of corr.mat must be at least as large as largest cluster")
  }
  if(!isSymmetric(corr.mat)){
    stop("corr.mat must be symmetric")
  }
  if(determinant(corr.mat, logarithm=T)$modulus == -Inf){
    stop("supplied correlation matrix is not invertible.")
  }	
  return(corr.mat[1:max(len), 1:max(len)])	
}

### summary function for geem object.
summary.geem <- function(object, ...)  {
  Coefs <- matrix(0,nrow=length(object$beta),ncol=5)
  Coefs[,1] <- c(object$beta)
  naive <- is.character(object$var)
  Coefs[,2] <- sqrt(diag(object$naiv.var))
  if(naive){Coefs[,3] <- rep(0, length(object$beta))}else{Coefs[,3] <- sqrt(diag(object$var))}
  if(naive){Coefs[,4] <- Coefs[,1]/Coefs[,2]}else{Coefs[,4] <- Coefs[,1]/Coefs[,3]}
  Coefs[,5] <- round(2*pnorm(abs(Coefs[,4]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Naive SE","Robust SE", "wald", "p")
  
  summ <- list(beta = Coefs[,1], se.model = Coefs[,2], se.robust = Coefs[,3], wald.test = Coefs[,4], p = Coefs[,5],
               alpha = object$alpha, corr = object$corr, phi = object$phi, niter = object$niter, clusz = object$clusz, coefnames = object$coefnames)
  
  class(summ) <- 'summary.geem'
  return(summ)
  #rownames(Coefs) <- c(object$coefnames)
  #print("Call: ", object$call, "\n")
  #print(Coefs)
  #cat("\n Est. Correlation: ", object$alpha, "\n")
  #cat(" Correlation Structure: ", object$corr, "\n")
  #cat(" Est. Scale Parameter: ", object$phi, "\n")
  #cat("\n Number of GEE iterations:", object$niter, "\n")
  #cat(" Number of Clusters: ", length(object$clusz), "   Maximum Cluster Size: ", max(object$clusz), "\n")
}

### print function for summary.geem object
print.summary.geem <- function(x, ...){
  Coefs <- matrix(0,nrow=length(x$coefnames),ncol=5)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Model SE","Robust SE", "wald", "p")
  Coefs[,1] <- x$beta
  Coefs[,2] <- x$se.model
  Coefs[,3] <- x$se.robust
  Coefs[,4] <- x$wald.test
  Coefs[,5] <- x$p
  
  #print("Call: ", object$call, "\n")
  print(Coefs)
  cat("\n Est. Correlation: ", x$alpha, "\n")
  cat(" Correlation Structure: ", x$corr, "\n")
  cat(" Est. Scale Parameter: ", x$phi, "\n")
  cat("\n Number of GEE iterations:", x$niter, "\n")
  cat(" Number of Clusters: ", length(x$clusz), "   Maximum Cluster Size: ", max(x$clusz), "\n")
}

### print function for geem object
print.geem <- function(x, ...){
  coefdf <- data.frame(x$beta)
  rownames(coefdf) <- x$coefnames
  colnames(coefdf) <- ""
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(coefdf))
  cat("\n Scale Parameter: ", x$phi, "\n")
  cat("\n Correlation Model: ", x$corr)
  cat("\n Estimated Correlation Parameters: ", x$alpha, "\n")
  cat("\n Number of clusters: ", length(x$clusz), "  Maximum cluster size: ", max(x$clusz), "\n")
}

### fitted function for geem object
fitted.geem <- function(object, ...){
  InvLink <- object$FunList$InvLink
  return(InvLink(object$eta))
}

predict.geem <- function(object, newdata = NULL,...){
  coefs <- object$beta
  if(is.null(newdata)){
    return(as.vector(object$X %*% object$beta))
  }else{
    if(dim(newdata)[2] != length(coefs)){warning("New observations must have the same number of rows as coefficients in the model")}
    return(as.vector(newdata %*% object$beta))
  }
}