sda <- function (x, ...) UseMethod("sda")


sda.default <- function(x, y, lambda=1e-6, stop=-p, maxIte=100, Q=K-1, trace=FALSE, tol=1e-6, ...){
    ##
    ## sda performs Sparse Linear Disciminant Analysis
    ## Solving: argmin{|(y*theta-x*b)|_2^2 + t*|beta|_1 + lambda*|beta|_2^2}
    ##
    ## InPUT:
    ## x      : matrix of n observations down the rows and p variable columns. The
    ##          columns are assumed normalized
    ## y      : matrix initializing the dummy variables representing the groups or a factor
    ## lambda : the weight on the L2-norm for elastic net regression. Default: 1e-6.
    ## stop   : If STOP is negative, its 
    ##          absolute value corresponds to the desired number of variables. If STOP
    ##          is positive, it corresponds to an upper bound on the L1-norm of the
    ##          b coefficients. There is a one to one correspondence between stop
    ##          and t.
    ## maxIte : Maximum number of iterations. Default: 100.
    ## Q      : number of discriminative directions. default Q=K
    ## trace  : trace = FALSE turns printing of RSS off and trace = TRUE turns it on.
    ## tol    : Tolerance for the stopping criterion (change in RSS). Default is 1e-6.
    ##
    ## OUTPUT:
    ## $beta    : The sparse loadings
    ## $theta   : Optimal scores
    ## $rss     : Residual Sum of Squares at each itearation
    ##
    ## Author: Line H. Clemmensen, IMM, DTU, lhc@imm.dtu.dk
    ## Modified by Trevor Hastie to become sequential
    ## Based on the elastic net algorithm by Hui Zou and Trevor Hastie
    ##

    ## this is stright from nnet:::formula
    class.ind <- function(cl) {
      Ik=diag(length(levels(cl)))
      x=Ik[as.numeric(cl),]
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
        classes <- colnames(y)
        factorY <- factor(colnames(y)[apply(y, 1, which.max)])
      }
    if(!is.matrix(x)) x <- as.matrix(x)
    x=scale(x,TRUE,FALSE)##This centering is essential for the trivial solution to disappear
    predNames <- colnames(x)
    
    n <- dim(x)[1]
    p <- dim(x)[2]
    K <- dim(y)[2]
    ones=rep(1,K)
    if(Q>(K-1))stop("at most K-1 variates allowed")
    
    
    dpi=as.vector(table(factorY)/n)
    #Dpi <- diag(dpi) ## diagonal matrix of class priors
    #Dpi_inv <- diag(1/sqrt(diag(Dpi)))
                                        # set-up stop criterion for elasticnet
    if (length(stop)< (Q)){
      stop <- rep(stop[1],1,Q)
    }
    if (stop[1]<0) sparse <- "varnum" else sparse <- "penalty" 
    Yhat <- matrix(0,n,Q)
    b <- matrix(0,p,Q)
    rss <- rep(0,Q)
    theta=matrix(0,K,Q)
    Qj=matrix(ones,K,1)
    ydp=scale(y,FALSE,dpi)
    for(j in 1:Q){
      RSS <- 1e6
      RSSold <- Inf
      ite <- 0
      thetaj=rtheta(K,dpi)
      thetaj=orth.Q(dpi,Qj,thetaj)
      while (abs(RSSold-RSS)/RSS > tol & ite < maxIte){ 
        RSSold <- RSS
        ite <- ite + 1
        ## 1. Estimate beta:    
        Yc <- y%*%thetaj 
        beta<- solvebeta(x, Yc, paras=c(lambda, abs(stop[j])),sparse=sparse) # elasticnet to estimate beta
        yhatj=x%*%beta
        thetaj=orth.Q(dpi,Qj,drop(t(ydp)%*%yhatj))
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
    
    
    if (trace){
      cat('final update, total ridge cost: ', sum(rss), ' |b|_1: ', sum(abs(b)),'\n')
    }

                                        # calcualte classes for LDA
    if (all(b==0)){ # all values in b and sl are zero
      warning("no non-zero elements - try other regularization parameters")
      b <- matrix(0,p,1)
      sl <- matrix(0,n,1)
      notZero <- matrix(TRUE,p,1) # fake NotZero for alignment purpose
      colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
      prior <- t(as.matrix(apply(y,2,sum)/n))
      colnames(prior) <- classes
      varNames <- colnames(x)
      origP <- ncol(x)
      lobj<-list(
                 prior=prior,
                 means = NULL,
                 x = x,
                 covw = matrix(0,1,1)
                 )
    }
    else{
                                        # remove predictors which are not included (do not have non-zero parameter estimates)
      notZero <- apply(b, 1, function(x) any(x != 0))
      b <- as.matrix(b[notZero,])
      origP <- ncol(x)
      x <- x[, notZero, drop = FALSE]
      varNames <- colnames(x)

### remove directions with only zero elements (this can be caused by a too high weight on L1-penalty)
      if (is.vector(b)){b<-t(as.matrix(b))}
      notZeroC <- apply(b,2,function(x) any(x!=0))
      b <- as.matrix(b[,notZeroC])
      
      
      sl <- x %*% b
      colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
      lobj<-lda(sl, factorY, ...)
    }
    

    structure(
              list(call = match.call(),
                   beta = b,
                   theta = theta,
                   varNames = varNames,
                   varIndex = which(notZero),
                   origP = origP,
                   rss = rss,
                   fit = lobj,
                   classes = classes,
                   lambda = lambda,
                   stop = stop),
              class = "sda")
  }
  
  predict.sda <- function(object, newdata = NULL, ...)
  {
    if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
    if(!is.null(object$varNames))
      {
        newdata <- newdata[, object$varNames, drop = FALSE]
      } else {
        if(ncol(newdata) != object$origP) stop("dimensions of training and testing X different")
        newdata <- newdata[, object$varIndex, drop = FALSE]
        
      }
    xnew <- newdata %*% object$beta
    pred <- predict(object$fit,xnew)
    pred
    
  }
      
   print.sda <- function(x, digits = max(3, getOption("digits") - 3), ...)
    {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
 
    classInfo <- paste(x$classes, collapse = ", ")

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
        "\nclasses =", classInfo,
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
    invisible(x)
  }
      
      



