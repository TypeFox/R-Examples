Youden3Grp.boxcox.fit <-
function (object1, xmat1, lambda, lambda2 = NULL, add.to.data1 = 0, object2, xmat2,add.to.data2 = 0, object3, xmat3,add.to.data3 = 0,...)
{
  ######This function estimate lambda by ML method which is to be used for Box-Cox transformation
  ###object1, 2, 3: data vectors x, y, z for the 3 groups
  ###xmat1, 2, 3: default a column of 1s if no other covariates are provided
  ###lambda: all 3 groups x, y, z use the same transformation [x^(lambda)-1]/lambda etc
  ##lambda2: if limit to a interval
  ##add.to.data1, 2, 3: add to x, y, z so that they are not negative,default to be zero
  ###..., can add other arguments to optimize function, eg. specify "interval=c(min(c(data1,data2,data3)), max(c(data1,data2,data3)))"
  
  ###NOTE: if lambda2 not given, search within a given interval (at Line 707) of (-5,5), will this be too narrow in some situations??
  # The method used is a combination of golden section search and successive parabolic interpolation, and was designed for use with continuous functions.
  #Convergence is never much slower than that for a Fibonacci search.  If 'f' has a continuous second derivative which is positive at the minimum (which is not at 'lower' or
  #'upper'), then convergence is superlinear, and usually of the order of about 1.324.

  
    call.fc <- match.call()

    data1 <- object1 + add.to.data1
    data2 <- object2 + add.to.data2
    data3 <- object3 + add.to.data3
    
    if (is.null(lambda2) && any(data1 <= 0)&& any(data2 <= 0)&& any(data3 <= 0)) 
        stop("Transformation requires positive data")

    data1 <- as.vector(data1)
    data2 <- as.vector(data2)
    data3 <- as.vector(data3)
    
    n1 <- length(data1)
    n2 <- length(data2)
    n3 <- length(data3)
    
    if (missing(xmat1)) 
        xmat1 <- rep(1, n1)
    if (missing(xmat2)) 
        xmat2 <- rep(1, n2)
    if (missing(xmat3)) 
        xmat3 <- rep(1, n3)
    
    xmat1 <- as.matrix(xmat1)
    xmat2 <- as.matrix(xmat2)
    xmat3 <- as.matrix(xmat3)
    
    if (any(xmat1[, 1] != 1)) 
        xmat1 <- cbind(1, xmat1)
    if (any(xmat2[, 1] != 1)) 
        xmat2 <- cbind(1, xmat2)
    if (any(xmat3[, 1] != 1)) 
        xmat3 <- cbind(1, xmat3)
    
    xmat1 <- xmat1[!is.na(data1), , drop = FALSE]
    xmat2 <- xmat2[!is.na(data2), , drop = FALSE]
    xmat3 <- xmat3[!is.na(data3), , drop = FALSE]
    
    data1 <- data1[!is.na(data1)]
    data2 <- data2[!is.na(data2)]
    data3 <- data3[!is.na(data3)]
    
    n1 <- length(data1)
    n2 <- length(data2)
    n3 <- length(data3)
    
    beta1.size <- ncol(xmat1)
    beta2.size <- ncol(xmat2)
    beta3.size <- ncol(xmat3)
    
    if (nrow(xmat1) != length(data1)) 
        stop("xmat1 and data1 have incompatible lengths")
    if (nrow(xmat2) != length(data2)) 
        stop("xmat3 and data3 have incompatible lengths")
    if (nrow(xmat3) != length(data3)) 
        stop("xmat3 and data3 have incompatible lengths")
    

    
    if (all(data1 > 0)) 
        absmin1 <- 0
    else absmin1 <- abs(min(data1)) + 1e-05 * diff(range(data1))
    if (all(data2 > 0)) 
        absmin2 <- 0
    else absmin2 <- abs(min(data2)) + 1e-05 * diff(range(data2))
    
    if (all(data3 > 0)) 
        absmin3 <- 0
    else absmin3 <- abs(min(data3)) + 1e-05 * diff(range(data3))

    lik.method <- "ML"##we  only use ML method
    
    if (!is.null(lambda2)) {
        if (missing(lambda)) 
            lambda.ini <- seq(-2, 2, by = 0.2)
        else lambda.ini <- lambda
        lambda2.ini <- 0
        if (isTRUE(lambda2)) 
            lambda2.ini <- min(absmin1,absmin2,absmin3)
        else if (mode(lambda2) == "numeric") 
            lambda2.ini <- lambda2
        lambdas.ini <- as.matrix(expand.grid(lambda.ini, lambda2.ini))
        if (length(as.matrix(lambdas.ini)) > 2) {
            lamlik <- apply(lambdas.ini, 1, negloglik.boxcox.3Grp, data1 = data1 + absmin1, xmat1 = xmat1, data2 = data2 + absmin2, xmat2 = xmat2,data3 = data3 + absmin3, xmat3 = xmat3,lik.method = lik.method)
            lambdas.ini <- lambdas.ini[which(lamlik == min(lamlik)), ]
          }
        lambdas.ini <- unname(drop(lambdas.ini))
        lik.lambda <- optim(par = lambdas.ini, fn = negloglik.boxcox.3Grp, 
            method = "L-BFGS-B", lower = c(-Inf, min(absmin1,absmin2,absmin3)), data1 = data1 + absmin1, xmat1 = xmat1,
                            data2 = data2 + absmin2, xmat2 = xmat2,data3 = data3 + absmin3, xmat3 = xmat3,lik.method = lik.method)

      }
    else {
        ##changed (-5,5) to (-20,20)
        lik.lambda <- optimize(negloglik.boxcox.3Grp, interval = c(-20, 20), data1 = data1 + absmin1, xmat1 = xmat1,data2 = data2 + absmin2, xmat2 = xmat2,data3 = data3 + absmin3, xmat3 = xmat3,lik.method = lik.method,tol=0.01)##08/06/2009,change tol from default (around 0.0001) to 0.1 for speed
        ######02/25/2010 change the argument interval=c(min(c(data1,data2,data3)), max(c(data1,data2,data3))) instead of (-20,20) in the optimize function 
        #lik.lambda <- optimize(negloglik.boxcox.3Grp, interval =c(min(c(data1,data2,data3)), max(c(data1,data2,data3))) , data1 = data1 + absmin1, xmat1 = xmat1,data2 = data2 + absmin2, xmat2 = xmat2,data3 = data3 + absmin3, xmat3 = xmat3,lik.method = lik.method,tol=0.01)
        
        lik.lambda <- list(par = lik.lambda$minimum, value = lik.lambda$objective, convergence = 0, message = "function optimize used")
      }

    ####DEBUG: 03/24/2010, sometimes return warnings on NA/inf from the optimize function, !checked because one group has variance 0 (e.g, discrete marker zbend)!!
    #if(any(is.na(lik.lambda$par) | is.na(lik.lambda$value) | is.infinite(lik.lambda$par) | is.infinite(lik.lambda$value)))  browser()
    
    lambda.fit <- lik.lambda$par
    
    if (length(lambda.fit) == 1) 
      lambda.fit <- c(lambda.fit, 0)

    data1 <- data1 + lambda.fit[2]
    data2 <- data2 + lambda.fit[2]
    data3 <- data3 + lambda.fit[2]

    ##########02/01/2010, comment the following lines to save time since we don't use the beta and sigmasq estimates and the transformed data
    

    #if (isTRUE(all.equal(unname(lambda.fit[1]), 0)))
    #  {    ###if lambda est=0, log transform (data)
    #    yt1 <- log(data1)
    #    yt2 <- log(data2)
    #    yt3 <- log(data3)
    #  }
    #else
    #  {
    #    yt1 <- ((data1^lambda.fit[1]) - 1)/lambda.fit[1]
    #    yt2 <- ((data2^lambda.fit[1]) - 1)/lambda.fit[1]
    #    yt3 <- ((data3^lambda.fit[1]) - 1)/lambda.fit[1]
    #  }

    #beta1 <- solve(crossprod(xmat1), crossprod(xmat1, yt1))
    #mu1 <- drop(xmat1 %*% beta1)
    #sigmasq1 <- sum((yt1 - mu1)^2)/n1

    #beta2 <- solve(crossprod(xmat2), crossprod(xmat2, yt2))
    #mu2 <- drop(xmat2 %*% beta2)
    #sigmasq2 <- sum((yt2 - mu2)^2)/n2

    #beta3 <- solve(crossprod(xmat3), crossprod(xmat3, yt3))
    #mu3 <- drop(xmat3 %*% beta3)
    #sigmasq3 <- sum((yt3 - mu3)^2)/n3
    
    
    loglik <- negloglik.boxcox.3Grp(lambda.val=lambda.fit,data1=data1, xmat1=xmat1,data2=data2,xmat2=xmat2,data3=data3,xmat3=xmat3, lik.method =lik.method)#negative loglike
    loglik <- -loglik
    
    #beta1 <- drop(beta1)
    #if (length(beta1) > 1)     names(beta1) <- paste("beta1_", 0:(beta1.size - 1), sep = "")
    #beta2 <- drop(beta2)
    #if (length(beta2) > 1)    names(beta2) <- paste("beta2_", 0:(beta2.size - 1), sep = "")
    #beta3 <- drop(beta3)
    #if (length(beta3) > 1) names(beta3) <- paste("beta3_", 0:(beta3.size - 1), sep = "")
        
    #if (length(lik.lambda$par) == 1) lambda.fit <- lambda.fit[1]
    #if (length(lik.lambda$par) == 2) names(lambda.fit) <- c("lambda", "lambda2")
    #res <- list(lambda = lambda.fit, beta1.normal = drop(beta1), sigmasq1.normal = sigmasq1, beta2.normal = drop(beta2), sigmasq2.normal = sigmasq2, beta3.normal = drop(beta3), sigmasq3.normal = sigmasq3, loglik = loglik, optim.results = lik.lambda)
    
    res <- list(lambda=lambda.fit,loglik=loglik)
    #res$call <- call.fc
    #oldClass(res) <- "boxcox.fit"
    return(res)
  }

