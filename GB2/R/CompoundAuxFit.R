pkl.cavgb2 <- function(z,lambda){
  L <- dim(lambda)[2]+1
  pkl <- exp(z%*%lambda)
  ck <- apply(pkl,1,sum)+1
  pkl <- pkl/ck
  pkl <- cbind(pkl,1-rowSums(pkl))
  colnames(pkl) <- paste("comp", 1:L, sep="")
return(pkl)
# pkl: matrix I x L
}

lambda0.cavgb2 <- function(pl0, z, w=rep(1, dim(z)[1])){
# vl0 - initial values of vl (do not depend on k!!!)
  vl0 <- vofp.cgb2(pl0)
  sw <- sum(w)
  I <- dim(z)[2]
  wz <- (w*z) 
  swz <- apply(wz, 2, sum)
# means - the vector of mean values of the aux vars
# sum(w_k*z_{ik})/sum(w_k)
  means <- swz/sw
  lambda0 <- (1/(I*means))%*%t(vl0)   
return(lambda0) 
}

logl.cavgb2 <- function(fac, z, lambda, w=rep(1, dim(fac)[1])){
  sw <- sum(w)
  pkL <- pkl.cavgb2(z,lambda)
  mixt <- rowSums(pkL*fac)
  logcompaux <- log(mixt)
  logL <- sum(w*logcompaux)/sw
return(logL)
}

scores.cavgb2 <- function(fac, z, lambda, w=rep(1, dim(fac)[1])){
# fac    : (N x L) - matrix of Gamma factors
# z      : (N x I) - matrix of auxiliary variables
# w      : (N x 1) - vector of weights
# lambda : (I x L-1)- matrix of parameters
  sw <- sum(w)
  evl <- exp(z%*%lambda)
  ck <- apply(evl,1,sum)+1
  pkl <- evl/ck #dimension N x L-1
  pkL <- cbind(pkl,1-rowSums(pkl)) #dimension N x L
  L <-dim(pkL)[2]
  denom <- rowSums(pkL*fac)
  num <- fac[,-L]
  midt <- num/as.vector(denom) - 1
  dlogL <- t(pkl*midt)%*%(z*w)/sw
return(dlogL)
}


ml.cavgb2 <- function (fac, z, lambda0, w = rep(1, dim(fac)[1]), maxiter = 100, fnscale=length(w)) 
{
    dl <- dim(lambda0)
    lambda0vec <- as.vector(lambda0)
    fn <- function(lambdavec, fac, z, w) {
        lambda <- matrix(lambdavec, nrow = dl[1], ncol = dl[2])
        return(-logl.cavgb2(fac, z, lambda, w))
    }
    gr <- function(lambdavec, fac, z, w) {
        lambda <- matrix(lambdavec, nrow = dl[1], ncol = dl[2])
        sc <- -scores.cavgb2(fac, z, lambda, w)
        sctr <- t(sc)
        asvecsc <- as.vector(sctr)
        return(asvecsc)
    }
    opt <- optim(lambda0vec, fn, gr, fac, z, w, method = "BFGS", 
        control = list(maxit = maxiter, fnscale = fnscale), 
        hessian = FALSE)
    lambdafit <- opt$par
    lambdafitm <- matrix(lambdafit, nrow = dl[1], ncol = dl[2])
    colnames(lambdafitm) <- paste("comp", 1:dl[2], sep = "")
    rownames(lambdafitm) <- colnames(z)
    return(list(lambdafitm, opt))
}