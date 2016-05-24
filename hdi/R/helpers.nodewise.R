## This file contains functions that calculate the Z vectors as nodewise lasso
## residuals
## see http://arxiv.org/abs/1110.2563
## and can also calculate the Thetahat matrix as in
## http://arxiv.org/abs/1303.0518

score.nodewiselasso <- function(x,
                                wantTheta = FALSE,
                                verbose = FALSE,
                                lambdaseq = "quantile",
                                parallel = FALSE,
                                ncores = 8,
                                oldschool = FALSE,
                                lambdatuningfactor = 1,
                                cv.verbose = FALSE)
{
  ## Purpose:
  ## This function calculates the score vectors Z OR the matrix of nodewise
  ## regression coefficients Thetahat, depending on the argument wantTheta.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  ## First, a sequence of lambda values over all nodewise regressions is created

  lambdas <-
    switch(lambdaseq,
           "quantile" = nodewise.getlambdasequence(x), ## this is preferred
           "linear"   = nodewise.getlambdasequence.old(x,verbose),
           ## otherwise
           stop("invalid 'lambdaseq': ", lambdaseq))
  
  if(verbose){
    cat("Using the following lambda values:", lambdas, "\n")
  }

  ## 10-fold cv is done over all nodewise regressions to calculate the error
  ## for the different lambda
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x,
                                      parallel = parallel, ncores = ncores,
                                      oldschool = oldschool,
                                      verbose = cv.verbose)
  if(verbose){
    cat(paste("lambda.min is", cvlambdas$lambda.min), "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), "\n")
  }
  
  if(lambdatuningfactor == "lambda.1se"){
    if(verbose)
      cat("lambda.1se used for nodewise tuning\n")
    ## We use lambda.1se for bestlambda now!!!
    bestlambda <- cvlambdas$lambda.1se
    }else{
      if(verbose)
        cat("lambdatuningfactor used is", lambdatuningfactor, "\n")
      
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  
  if(verbose){
    cat("Picked the best lambda:", bestlambda, "\n")
    ##print("with the error ")
    ##print(min(err))
  }

  ## Having picked the 'best lambda', we now generate the final Z or Thetahat
  if(wantTheta){
    out <- score.getThetaforlambda(x = x,
                                   lambda = bestlambda,
                                   parallel = parallel,
                                   ncores = ncores,
                                   oldschool = TRUE,
                                   verbose = verbose)
  }else{
    Z <- score.getZforlambda(x = x, lambda = bestlambda, parallel = parallel,
                             ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out,
                     bestlambda = bestlambda)
  return(return.out)
}

score.getThetaforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                    oldschool = FALSE, verbose = FALSE,
                                    oldtausq = TRUE)
{
  ## Purpose:
  ## This function is for calculating Thetahat once the desired tuning
  ## parameter is known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  message("Calculating Thetahat by doing nodewise regressions and dropping the unpenalized intercept")
  n <- nrow(x)
  p <- ncol(x)
  C <- diag(rep(1,p))
  T2 <- numeric(p)
  
  if(oldschool){
    message("doing getThetaforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i], x[,i])
      coeffs <- as.vector(predict(glmnetfit,x[,-i], type = "coefficients",
                                  s = lambda))[-1]
      ## we just leave out the intercept
      
      C[-i,i] <- -as.vector(coeffs)
      if(oldtausq){
        ## possibly quite incorrect,it ignores the intercept, but intercept is
        ## small anyways. Initial simulations show little difference.
        T2[i] <- as.numeric(crossprod(x[,i]) / n - x[,i] %*% (x[,-i] %*%
                                                              coeffs) / n)
        }else{   
          ##print("now doing the new way of calculating tau^2")
          T2[i] <- as.numeric((x[,i] %*%
                               (x[,i] - predict(glmnetfit,x[,-i],s =lambda)))/n)
        }
    }
  }else{
    stop("not implemented yet!")
  }
  thetahat <- C %*% solve(diag(T2))
  if(verbose){
    cat("1/tau_j^2:", solve(diag(T2)), "\n")
  }
  ##this is thetahat ^ T!!
  thetahat <- t(thetahat)
  
  if(all(thetahat[lower.tri(thetahat)] == 0,
         thetahat[upper.tri(thetahat)] == 0) && verbose)
    cat("Thetahat is a diagonal matrix!\n")
  
  return(thetahat)
}

score.getZforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                oldschool = FALSE)
{
  ## Purpose:
  ## This function is for calculating Z once the desired tuning parameter is
  ## known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version)
    
  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n*p),n)
  
  if(oldschool){
    message("doing getZforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i],x[,i])
      prediction <- predict(glmnetfit,x[,-i],s=lambda)
      Z[,i] <- x[,i] - prediction
    }
  }else{
    ## REPLACING THE FOR LOOP
    if(parallel){
      Z <- mcmapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                    lambda = lambda, mc.cores = ncores)
      
    }else{
      Z <- mapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                  lambda = lambda)
    }
  }
  ## rescale Z such that t(Zj) Xj/n = 1 \-/ j
  Z <- score.rescale(Z,x)
  return(Z)
}

score.getZforlambda.unitfunction <- function(i, x, lambda)
{
  ## Purpose:
  ## Calculate the residuals of a nodewise regression of one column vs the
  ## other columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  glmnetfit  <- glmnet(x[,-i],x[,i])
  prediction <- predict(glmnetfit,x[,-i],s=lambda)
  return(x[,i] - prediction)
}

score.rescale <- function(Z, x)
{
  ## Purpose:
  ## Rescale the Z appropriately such that such that t(Zj) Xj/n = 1 for all j
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  scaleZ <- diag(crossprod(Z,x)) / nrow(x)
  Z      <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z=Z,scaleZ=scaleZ))
}

nodewise.getlambdasequence <- function(x)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## Equidistant quantiles of the complete set of lambda values are returned.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  
  nlambda <- 100 ## use the same value as the glmnet default
  p <- ncol(x)
  
  lambdas <- c()
  for(c in 1:p){
    lambdas <- c(lambdas,glmnet(x[,-c],x[,c])$lambda)
  }
  
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

cv.nodewise.err.unitfunction <- function(c, K, dataselects, x, lambdas,
                                         verbose, p) {
  ## Purpose:
  ## this method returns the K-fold cv error made by the nodewise regression
  ## of the single column c of x vs the other columns for all values of the
  ## tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  
  if(verbose){ ##print some information out about the progress
      ##report every 25%
    interesting.points <- round(c(1/4,2/4,3/4,4/4)*p)
    names(interesting.points) <- c("25%","50%","75%","100%")
    if(c %in% interesting.points){
      cat("The expensive computation is now",
          names(interesting.points)[c == interesting.points],
          "done\n")
    }
  }
  
  ## return  'totalerr' 
  cv.nodewise.totalerr(c = c,
                       K = K,
                       dataselects = dataselects,
                       x = x,
                       lambdas = lambdas)
}

## gets the standard error for one particular lambda
cv.nodewise.stderr <- function(K, x, dataselects, lambda, parallel, ncores)
{
  ## Purpose:
  ## this method returns the standard error 
  ## of the average K-fold cv error made by the nodewise regression
  ## of each column vs the other columns for a single lambda value.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  p <- ncol(x)
  if(parallel){
    totalerr <- mcmapply(cv.nodewise.totalerr,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = lambda,
                         mc.cores = ncores)
  }else{
    totalerr <- mapply(cv.nodewise.totalerr,
                       c = 1:p,
                       K = K,
                       dataselects = list(dataselects = dataselects),
                       x = list(x = x),
                       lambdas = lambda)
  }
  ## get the mean over the variables
  totalerr.varmean <- rowMeans(totalerr)
  
  ## get the stderror over the K;  return stderr.forlambda
  sd(totalerr.varmean) / sqrt(K)
}

cv.nodewise.totalerr <- function(c, K, dataselects, x, lambdas)
{
  ## Purpose:
  ## this method returns the error made for each fold of a K-fold cv
  ## of the nodewise regression of the single column c of x vs the other
  ## columns for all values of the tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  
  for(i in 1:K){ ## loop over the test sets
    whichj <- dataselects == i ##the test part of the data
    
    glmnetfit <- glmnet(x = x[!whichj,-c, drop = FALSE],
                        y = x[!whichj, c, drop = FALSE],
                        lambda = lambdas)
    predictions  <- predict(glmnetfit, newx = x[whichj, -c, drop = FALSE],
                            s = lambdas)
    totalerr[, i] <- apply((x[whichj, c] - predictions)^2, 2, mean)
  }

  totalerr
}


cv.nodewise.bestlambda <- function(lambdas, x, K = 10, parallel = FALSE,
                                   ncores = 8, oldschool = FALSE,
                                   verbose = FALSE)
{
  ## Purpose:
  ## this function finds the optimal tuning parameter value for minimizing
  ## the K-fold cv error of the nodewise regressions.
  ## A second value of the tuning parameter, always bigger or equal to the
  ## former, is returned which is calculated by allowing the cv error to
  ## increase by the amount of
  ## 1 standard error (a similar concept as to what is done in cv.glmnet).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)
  
  ## Based on code from cv.glmnet for sampling the data
  dataselects <- sample(rep(1:K, length = n))
  
  if(oldschool){
    message("doing cv.nodewise.error oldschool")
    totalerr <- numeric(l)
    for(c in 1:p){ ## loop over the nodewise regressions
      for(i in 1:K){ ## loop over the test sets
        whichj <- dataselects == i ## the test part of the data
        
        glmnetfit <- glmnet(x[!whichj,-c,drop=FALSE], x[!whichj,c,drop=FALSE],
                            lambda=lambdas)
        predictions <- predict(glmnetfit,x[whichj, -c, drop = FALSE],
                               s = lambdas)
        totalerr <- totalerr + apply((x[whichj,c]-predictions)^2, 2, mean)
      }
    }
    totalerr <- totalerr / (K * p)
    stop("lambda.1se not implemented for oldschool cv.nodewise.bestlamba")
  }else{
    ## REPLACING THE FOR LOOP

    ##totalerr <- matrix(nrow = l, ncol = p)
    
    if(parallel){
      totalerr <- mcmapply(cv.nodewise.err.unitfunction,
                           c = 1:p,
                           K = K,
                           dataselects = list(dataselects = dataselects),
                           x = list(x = x),
                           lambdas = list(lambdas = lambdas),
                           mc.cores = ncores,
                           SIMPLIFY = FALSE,
                           verbose = verbose,
                           p=p)
    }else{
      totalerr <- mapply(cv.nodewise.err.unitfunction,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = list(lambdas = lambdas),
                         SIMPLIFY = FALSE,
                         verbose = verbose,
                         p = p)
    }
    ## Convert into suitable array (lambda, cv-fold, predictor)
    err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
    err.mean   <- apply(err.array, 1, mean) ## 1 mean for each lambda

    ## calculate mean for every lambda x fold combination (= average over p)
    ## for every lambda then get the standard errors (over folds)
    err.se     <- apply(apply(err.array, c(1, 2), mean), 1, sd) / sqrt(K)
    ##totalerr <- apply(totalerr, 1, mean)
  }

  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]

  stderr.lambda.min <- err.se[pos.min]
##-   stderr.lambda.min <- cv.nodewise.stderr(K = K,
##-                                           x = x,
##-                                           dataselects = dataselects,
##-                                           lambda = lambda.min,
##-                                           parallel = parallel,
##-                                           ncores = ncores)
  list(lambda.min = lambda.min,
       lambda.1se = max(lambdas[err.mean < (min(err.mean) + stderr.lambda.min)]))
}

##DEPRECATED, not really the best choice
nodewise.getlambdasequence.old <- function(x,verbose=FALSE)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## It returns a __linear__ interpolation of lambda values between the max and
  ## min lambda value found.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  nlambda <- 100#take the same value as glmnet does automatically
  p <- ncol(x)
  maxlambda <- 0
  minlambda <- 100
  
  for(c in 1:p){
    lambdas <- glmnet(x[,-c],x[,c])$lambda
    
      ##DEBUG
    if(verbose || is.nan(max(lambdas))){
      cat(paste("c:", c, "\n"))
      cat("lambdas:", lambdas, "\n")
      cat("max(lambdas) max(lambdas,na.rm=TRUE) maxlambda: ",
          max(lambdas), max(lambdas,na.rm=TRUE), maxlambda, "\n")
    }
    if(max(lambdas,na.rm=TRUE) > maxlambda){
      maxlambda <- max(lambdas,na.rm=TRUE)
    }
    if(min(lambdas,na.rm=TRUE) < minlambda){
      minlambda <- min(lambdas, na.rm = TRUE)
    }
  }
  
  lambdas <- seq(minlambda, maxlambda, by = (maxlambda-minlambda)/nlambda)
  ## return
  sort(lambdas, decreasing=TRUE)
}
