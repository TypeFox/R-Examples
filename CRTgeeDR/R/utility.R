#Exactly the function 'getfam' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
getfam <- function(family){
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
    return(FunList)
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
}

#Exactly the function 'getBlockDiag' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
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

#Exactly the function 'checkFixedMat' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
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




#' Fit CRTgeeDR object.
#' 
#' Fit CRTgeeDR object to a dataset
#' 
#' @param object CRTgeeDR object
#' @param ... ignored
#' @aliases fitted.CRTgeeDR fitted
#' @method fitted CRTgeeDR 
#' @export
### fitted function for CRTgeeDR object
fitted.CRTgeeDR <- function(object, ...){
  InvLink <- object$FunList$InvLink
  return(InvLink(object$eta))
}

#' Predict CRTgeeDR object.
#' 
#' Predict CRTgeeDR object to a dataset
#' 
#' @param object CRTgeeDR object
#' @param newdata dataframe, new dataset to which the CRTgeeDRneed to be used for prediction
#' @param ... ignored
#' @method predict CRTgeeDR 
#' @aliases predict.CRTgeeDR predict
#' @export
predict.CRTgeeDR <- function(object, newdata = NULL,...){
  coefs <- object$beta
  if(is.null(newdata)){
    return(as.vector(object$X %*% object$beta))
  }else{
    if(dim(newdata)[2] != length(coefs)){warning("New observations must have the same number of rows as coefficients in the model")}
    return(as.vector(newdata %*% object$beta))
  }
}


cleandata<-function(dat,type,nameY,cc,formula,print=TRUE){
  na.inds <- NULL
  if(any(is.na(dat))){
    na.inds <- which(is.na(dat), arr.ind=T)
  }
  
  if(sum(c(which(colnames(dat)%in%all.vars(formula[-2])))%in%na.inds[,2])!=0){
    if(type=="marginal model"){
      if(print)warning(paste("It exists missing data for covariates in the ",type," model \n -> Weights for these observation are set to 0",sep=""))
    }else{
      if(print)warning(paste("It exists missing data for covariates in the ",type," model \n -> Single imputation had been used for computation",sep=""))
    }
  }
  
  if(!is.null(na.inds)){
    if(cc){
      #### If missing set to median value (No problem for outcome but carefull with covariates)
      dat$weights[unique(na.inds[,1])] <- 0
      for(i in unique(na.inds[,2])){
        if(is.factor(dat[,i])){
          dat[na.inds[,1], i] <- levels(dat[,i])[1]  
        }else{
          dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
        }
      }
    }else{
      for(i in unique(na.inds[,2][!na.inds[,2] %in% which(colnames(dat)==nameY)])){
        if(is.factor(dat[,i])){
          dat[na.inds[,1], i] <- levels(dat[,i])[1]  
        }else{
          dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
        }
      }
    }
  }
  return(dat)
}

#Exactly the function 'cormax.ind' included in the package 'geesmv',
#authored by Ming Wang
#under the GPL-2 license.
cormax.ind <-
  function(n){
    matrix<-matrix(0,nrow=n,ncol=n)
    diag(matrix)<-rep(1,n)
    return(matrix)
  }


#Exactly the function 'cormax.ar1' included in the package 'geesmv',
#authored by Ming Wang
#under the GPL-2 license.
cormax.ar1 <-
  function(n, alpha){
    n.max<- max(n)
    cor.max<- diag(1,n.max)
    lowertri<- rep(0,0)
    for(j in (n.max-1):1){
      lowertri<- c(lowertri,1:j)
    }
    cor.max[lower.tri(cor.max)]<- alpha^lowertri
    cor.max[upper.tri(cor.max)]<- alpha^lowertri[length(lowertri):1]
    return(cor.max)
  }


#Exactly the function 'cormax.exch' included in the package 'geesmv',
#authored by Ming Wang
#under the GPL-2 license.
cormax.exch <-
  function(n, alpha){
    n.max<- n
    cor.max<- diag(1,n.max)
    cor.max[lower.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2)
    cor.max[upper.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2)
    return(cor.max)
  }

#Exactly the function 'cluster.size' included in the package 'geesmv',
#authored by Ming Wang
#under the GPL-2 license.
cluster.size <-
  function(id){
    clid<- unique(id)
    m<- length(unique(id))
    n<- rep(0,m)
    autotime<- rep(0,0)
    for(i in 1:m){
      n[i]<- length(which(id==clid[i]))
      autotime<- c(autotime,1:n[i])
    }
    id<- rep(1:m,n)
    return(list(m=m,n=n,id=id,autotime=autotime))
  }

#Exactly the function 'vec' included in the package 'fBasics',
#authored by Rmetrics Core Team
#under the GPL-2 license.
vec <-
  function(x)
  {
    t(t(as.vector(x)))
  }


