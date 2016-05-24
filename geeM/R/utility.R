getfam <- function(family){
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))  
  }
  
  if(is.function(family)){
    family <- family()
    return(family)
  }else if(inherits(family, "family")){
    return(family)
  }else if(is.list(family)){
    if(length(match(names(family), c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv"))) == 4){
      famname <- "custom"
      LinkFun <- family$LinkFun
      InvLink <- family$InvLink
      VarFun <- family$VarFun
      InvLinkDeriv <- family$InvLinkDeriv
    }else{
      famname <- "custom"
      LinkFun <- family[[1]]
      VarFun <- family[[2]]
      InvLink <- family[[3]]
      InvLinkDeriv <- family[[4]]
    }
    
    
    FunList <- list("family"= famname, "LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv) 
    return(FunList)
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
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
  if(K == 1){
    index <- c(0, sum(len^2))
  }else{
    index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  }
  
  for(i in 1:K){
    row.vec[(index[i] + 1):(index[i+1])] <- rep.int( (1:len[i]) + add.vec[i], len[i])
    col.vec[(index[i] + 1):(index[i+1])] <- rep( (1:len[i]) + add.vec[i], each=len[i])
  }	
  BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
  
  if(!is.null(xvec)){
    testsymm <- abs(sum(skewpart(BlockDiag)))
    if(testsymm != 0) {
      warning("Correlation matrix is not computed to be exactly symmetric. Taking only the symmetric part.")
    }
  }
  return(list(BDiag = symmpart(BlockDiag), row.vec =row.vec, col.vec=col.vec))
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





### fitted function for geem object
fitted.geem <- function(object, ...){
  #InvLink <- object$FunList$InvLink
  InvLink <- object$FunList[[if(object$FunList$family == "custom") "InvLink" else "linkinv"]]
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

coef.geem <- function(object, ...){
  coefs <- object$beta
  names(coefs) <- object$coefnames
  return(coefs)
}

family.geem <- function(object,...){
  
  return(object$FunList)
}

dummyrows <- function(formula, dat, incomp, maxwave, wavespl, idspl){
  missing <- missid <- misswave <- rep(0, sum(maxwave))
  index <- 1
  for(i in 1:length(wavespl)){
    wa <- wavespl[[i]]
    index <- index+1
    for(j in 2:length(wa)){
      wdiff <- wa[j] - wa[j-1] -1
      if(wdiff > 0){
        missing[index:(index+wdiff-1)] <- (wa[j-1]+1):(wa[j]-1)
        missid[index:(index+wdiff-1)] <- idspl[[i]][1]
      }
      index <- index+wdiff+1
    }
  }
  dat2 <- as.data.frame(matrix(nrow=sum(maxwave), ncol=dim(dat)[2]))
  colnames(dat2) <- colnames(dat)
  dat2[missing==0,] <- dat
  dat2$id[missing > 0] <- missid[missing>0]
  
  dat2$weights[missing > 0] <- 0
  dat2$waves[missing > 0] <- missing[missing > 0]
  
  NAcols <- which(!is.element(names(dat2), c("id", "waves", "weights")))
  for(i in NAcols){
    dat2[missing>0, NAcols[i]] <- median(dat2[,NAcols[i]], na.rm=TRUE)
  }
  
  retdat <- model.frame(formula, dat2, na.action=na.pass)
  retdat$id <- dat2$id
  retdat$weights <- dat2$weights
  retdat$waves <- dat2$waves
  return(retdat)
}

