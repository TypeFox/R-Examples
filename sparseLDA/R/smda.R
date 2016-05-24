smda <- function (x, ...) UseMethod("smda")


smda.default <- function(x, y, Z = NULL, Rj = NULL, lambda=1e-6, stop, maxIte=50, Q=R-1, trace=FALSE, tol=1e-4, ...){
  ##
  ## smda performs Sparse Mixture Disciminant Analysis
  ## Solving: argmin{|(Y*theta-X*b)|_2^2 + t*|beta|_1 + lambda*|beta|_2^2}
  ##
  ## INPUT:
  ## x      : matrix of n observations down the rows and p variable columns. The
  ##          columns are assumed normalized
  ## Z      : matrix initializing the probabilities representing the groups
  ## Rj     : K length vector containing the number of subclasses in each of
  ##          the K classes
  ## lambda : the weight on the L2-norm for elastic net regression. Default: 1e-6
  ## stop   : nonzero STOP will perform
  ##          elastic net regression with early stopping. If STOP is negative, its 
  ##          absolute value corresponds to the desired number of variables. If STOP
  ##          is positive, it corresponds to an upper bound on the L1-norm of the
  ##          b coefficients. There is a one to one correspondence between stop
  ##          and t.
  ## maxIte : Maximum number of iterations. Default: 50.
  ## trace  : trace = FALSE turns printing of RSS off and trace = TRUE turns it on.
  ## tol    : Tolerance for the stopping criterion (change in RSS). Default: 1e-4
  ##
  ## OUTPUT:
  ## $beta   : The regression parameters
  ## $theta  : Optimal scores
  ## $Z      : Updated subclass probabilities
  ## $rss    : Residual Sum of Squares at each itearation
  ##
  ## Author: Line H. Clemmensen, IMM, DTU, lhc@imm.dtu.dk
  ## Based on the elastic net algorithm by Hui Zou and Trevor Hastie
  ##
  
  
  ## this is stright from nnet:::formula
  class.ind <- function(cl) {    
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  orth.Q=function(dp,Qj,theta){
    Qjp=Qj*as.vector(dp)
    qjtheta=t(Qjp)%*%theta
    theta=theta-Qj%*%qjtheta
    thetan=sqrt(apply(dp*theta^2,2,sum))
    scale(theta,FALSE,thetan)
  }
  rtheta=function(K,dp){
    jj=rnorm(K);
    jj/sqrt(sum(jj^2)*dp)
  }
  
  if(is.factor(y))
  {
    classes <- levels(y)
    factorY <- y
    y <- class.ind(y)
  } else {
    if(is.null(colnames(y))) colnames(y) <- paste("class", 1:ncol(y), sep = "")
    classes <- colnames(y)
    factorY <- factor(colnames(y)[apply(y, 1, which.max)])
  }
  
  if(is.null(Rj)) Rj <- rep(3, length(classes))
  if(length(Rj) == 1) Rj <- rep(Rj, length(classes))
  
  classKey <- rep(classes, times = Rj)
  
  subClasses <- classKey
  for(i in seq(along = classes))
  {
    tmp <- subClasses[subClasses == classes[i]]
    subClasses[subClasses == classes[i]] <- paste(tmp, seq(along = tmp), sep = ".")
  }
  
  if(!is.matrix(x)) x <- as.matrix(x)
  predNames <- colnames(x)
  
  N <- dim(x)[1]
  p <- dim(x)[2]
  K <- length(Rj) ## number of classes
  
  ## make Z from y
  if(is.null(Z))
  {
    tmp <- mda.start(x, factorY, subclasses = Rj,  start.method = "kmeans")
    Z <- matrix(0, nrow = nrow(x), ncol = sum(Rj))
    for(i in seq(along = tmp))
    {
      colIndex <- which(classKey == names(tmp)[i])
      rowIndex <- which(factorY == names(tmp)[i])
      Z[rowIndex, colIndex] <- tmp[[i]]
    }
    rm(tmp)
  }
  Isubcl <- apply(Z,1,which.max)
  colnames(Z) <- subClasses
  factorSubY <- factor(colnames(Z)[apply(Z, 1, which.max)])
  
  R <- dim(Z)[2] ## total number of subclasses
  if(Q>(R-1))stop("at most R-1 variates allowed")
  RSSold <- 1e8 
  RSS <- 1e6
  item <- 0
  Zhat <- matrix(0,N,Q)
  dpi <- apply(Z,2,sum)/N
  zdp <- scale(Z, FALSE, dpi) 
  theta <- matrix(0,R,Q)
  Ztheta <- Z%*%theta  ## N x Q
  rss <- rep(0,maxIte)
  b <- matrix(0,p,Q)
  if (length(stop)< Q){
    stop <- rep(stop[1],1,Q)
  }
  if (stop[1]<0) sparse <- "varnum" else sparse <- "penalty" 
  
  while (abs(RSSold-RSS)/RSS > tol & item < maxIte){ 
    RSSold <- RSS
    item <- item + 1
    Qj=matrix(1,R,1)
    ## 1. Estimate beta and theta
    for(j in 1:Q){
      RSS <- 1e6
      RSSold <- Inf
      ite <- 0
      thetaj=rtheta(R,dpi)
      thetaj=orth.Q(dpi,Qj,thetaj)
      while (abs(RSSold-RSS)/RSS > tol & ite < maxIte){ 
        RSSold <- RSS
        ite <- ite + 1
        ## 1. Estimate beta:    
        Yc <- Z%*%thetaj 
        beta <- solvebeta(x, Yc, paras=c(lambda, abs(stop[j])),sparse=sparse) # elasticnet to estimate beta
        yhatj=x%*%beta
        thetaj=orth.Q(dpi,Qj,drop(t(zdp)%*%yhatj))
        RSS=sum((yhatj-Yc)^2)+lambda*sum(beta^2)
        if (trace){ 
          cat('ite: ', ite, ' ridge cost: ', RSS, ' |b|_1: ', sum(abs(beta)),'\n')
        }
      }
      rss[j]=RSS
      Qj=cbind(Qj,thetaj)
      theta[,j]=thetaj
      b[,j]=beta
    }
    Zhat <- x%*%b
    Ztheta <- Z%*%theta
    RSS <- sum((Ztheta-Zhat)^2) + lambda*sum(b^2)
    rss[ite] <- RSS
    if (trace){
      cat('ite: ', item, ' ridge cost: ', RSS, ' l1-norm: ', sum(abs(b)), '\n')
    }
    
    ## 2. update parameter estimates - Expectation Maximization:
    Sigma <- matrix(0,Q,Q)
    mu <- matrix(0,Q*R,K)
    dim(mu) <- c(Q,R,K)
    for (i in 1:K){
      IK <- (sum(Rj[1:i-1])+1):(sum(Rj[1:i-1])+Rj[i])
      for (j in 1:Rj[i]){
        Ik <- which(Isubcl==IK[j])
        Ik.length <- length(Ik)
        sumZ <- sum(Z[Ik,IK[j]])
        mu[,IK[j],i] = apply(((Z[Ik,IK[j]])%*%matrix(1,1,Q))*Zhat[Ik,,drop = FALSE],2,sum)/sumZ
        Sigma = Sigma + t(Zhat[Ik,,drop = FALSE]-matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i]))*(Z[Ik,IK[j]]%*%matrix(1,1,Q)))%*%(Zhat[Ik,,drop = FALSE]-
                                                                                                                                     matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i])))/sumZ
      }
    }
    if (kappa(Sigma)>1e8){
      Sigma = Sigma + 1e-3*diag(rep(1,Q))
    }
    Sigma_inv <- solve(Sigma)
    
    for (i in 1:K){
      IK <- (sum(Rj[1:i-1])+1):(sum(Rj[1:i-1])+Rj[i])
      Dmahal_K <- matrix(0,N,Rj[i])
      for (j in 1:Rj[i]){
        Dmahal_K[,j] <- diag((Zhat-matrix(1,N,1)%*%t(matrix(mu[,IK[j],i])))%*%Sigma_inv%*%t(Zhat-
                                                                                              matrix(1,N,1)%*%t(matrix(mu[,IK[j],i]))))
      }
      sum_K <- apply(matrix(1,N,1)%*%dpi[IK]*exp(-Dmahal_K/2),1,sum)
      for (j in 1:Rj[i]){
        Z[,IK[j]] <- dpi[IK[j]]*exp(-Dmahal_K[,j]/2)/(sum_K+1e-3)
      }
      dpi[IK] <- sum(Z[,IK])
      dpi[IK] <- dpi[IK]/sum(dpi[IK])
    }
    zdp <- scale(Z, FALSE, dpi)
    # end update of parameters Expectation-Maximization.
    Ztheta <- Z%*%theta
    RSS <- sum((Ztheta-Zhat)^2)+lambda*sum(b^2)
    if (trace){  
      cat('EM update, ridge cost: ', RSS, ' l1-norm: ', sum(abs(b)), '\n')
    }
    if (ite==maxIte){warning('Forced exit. Maximum number of iterations reached in smda step.')}
  }
  
  notZero <- apply(b, 1, function(x) any(x != 0))
  b <- b[notZero,,drop = FALSE]
  origP <- ncol(x)
  x <- x[, notZero, drop = FALSE]
  varNames <- colnames(x)
  
  sl <- x %*% b
  colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
  lobj<-lda(sl, factorSubY, ...)
  
  structure(
    list(call = match.call(),
         beta = b,
         theta = theta,
         Z = Z,
         Rj = Rj,
         K = K,
         mu = mu,
         Sigma = Sigma,
         Sigma_inv = Sigma_inv,
         dpi = dpi,
         varNames = varNames,
         varIndex = which(notZero),
         origP = origP,
         rss = rss[1:ite],
         fit = lobj,
         classes = classes,
         subClasses = subClasses,
         lambda = lambda,
         stop = stop),
    class = "smda")
}


predict.smda <- function(object, newdata = NULL, ...)
{
  if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
  if(!is.null(object$varNames))
  {
    newdata <- newdata[, object$varNames, drop = FALSE]
  } else {
    if(ncol(newdata) != object$origP) stop("dimensions of training and testing X different")
    newdata <- newdata[, object$varIndex, drop = FALSE]
  }
  x <- newdata %*% object$beta
  #subPred <- predict(object$fit, newdata = x, ...)
  #calculate subclass probabilities
  Rz <- c(0,object$Rj)
  Zt <- matrix(0,dim(x)[1],sum(object$Rj))
  for (i in 1:object$K){
    IK <- (sum(object$Rj[1:i-1])+1):(sum(object$Rj[1:i-1])+object$Rj[i])
    Dmahal_K <- matrix(0,dim(x)[1],object$Rj[i])
    for (j in 1:object$Rj[i]){
      Dmahal_K[,j] <- diag((x-matrix(1,dim(x)[1],1)%*%t(matrix(object$mu[,IK[j],i])))%*%object$Sigma_inv%*%t(x-
                                                                                                               matrix(1,dim(x)[1],1)%*%t(matrix(object$mu[,IK[j],i]))))
    }
    sum_K <- apply(matrix(1,dim(x)[1],1)%*%object$dpi[IK]*exp(-Dmahal_K/2),1,sum)
    for (j in 1:object$Rj[i]){
      Zt[,IK[j]] <- object$dpi[IK[j]]*exp(-Dmahal_K[,j]/2)/(sum_K+1e-3)
    }
  }
  pr <- matrix(0,dim(x)[1],object$K)
  for (i in 1:object$K){
    pr[,i] <- apply(Zt[,(sum(Rz[1:i])+1):(sum(Rz[1:i])+Rz[i+1])],1,sum)
  }
  class <- factor(object$classes[apply(pr,1,which.max)], levels = object$classes)
  colnames(pr) <- object$classes
  colnames(Zt) <- object$subClasses
  list(class=class,
       classprob = pr/apply(pr,1,sum),
       subprob = Zt/apply(Zt,1,sum))
}


print.smda <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  
  classInfo <- paste(paste(x$classes, " (", x$Rj, ")", sep = ""), collapse = ", ")
  
  if(all(x$stop < 0))
  {
    stopVal <- paste(-x$stop[1], "variables")
  } else {
    stopVal <- paste(
      paste(format(x$stop, digits = digits),
            collapse = ", "),
      "L1 bounds")
  }
  
  cat("lambda =", format(x$lambda, digits = digits),
      "\nstop =", stopVal,
      "\nsubclasses =", classInfo,
      "\n\n")
  
  top <- if(!is.null(x$varNames)) x$varNames else paste("Predictor", x$varIndex, sep = "")
  varOrder <- if(is.matrix(x$beta)) order(apply(abs(x$beta), 1, sum)) else order(abs(x$beta))
  top <- top[varOrder]
  top <- top[1:min(5, length(top))]
  top <- paste(top, collapse = ", ")
  
  if(nrow(x$beta) > 5)
  {
    cat("Top 5 predictors (out of ",
        length(x$varIndex),
        "):\n\t",
        top,
        sep = "")
  } else {
    cat("Predictors:\t",
        top,
        "\n",
        sep = "")
  }
  cat("\n")
  invisible(x)
}
