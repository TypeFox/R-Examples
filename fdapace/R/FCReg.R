#' Functional Principal Component Analysis Concurrent Regression with Functional dependent variable
#' 
#' Functional concurrent regression for dense or sparse functional data. 
#' 
#' @param expVarScal  A data.frame holding scalar explanatory variables; NAs will be omitted internally.
#' @param expVarFunc  A list holding the FPCA objects for each of the functional explanatory variables.
#' @param depVar  An FPCA object
#' @param regressionType A string defining the type of regression to perform ('dense' or 'sparse'); (default : automatically determined based on 'depVar')
#' @param bwScalar The value of bandwidth to be used for all scalar/functional cross-covariances (default: automatically determined using GCV)
#' @param bwFunct The values of bandwiths to be used for all function/function cross-covariances (default: automatically determined using GCV)
#' @param splineSmooth  Use thin-plate splines during the estimation of the cross-covariance (default: FALSE)
#' @param verbose If TRUE print out the bandwidth used during the GCV procedures of selecting them
#' @param getFitted If TRUE append the fitted values to the return object
#' @param ...  Additional arguments 
#' 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#'
#' \cite{Senturk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)} 

FCReg <- function(depVar,  expVarScal = NULL, expVarFunc = NULL, regressionType = NULL, getFitted = TRUE, 
                  bwScalar = NULL, bwFunct = NULL, splineSmooth = FALSE,  verbose = FALSE){
  
  if ( is.null(regressionType)){
    regressionType = depVar$optns$dataType
  }
  
  if( regressionType == 'Sparse'){
    if (is.null(expVarFunc)){
      stop('You have no functional covariates; just use regular regression functions.')
    }
    
    if (is.null(expVarScal)){
      nZ = 0;
    } else {
      if( is.vector(expVarScal)){
        expVarScal = data.frame(expVarScal)
      }
      nZ =  ncol(expVarScal)
    }
    
    nX = length(expVarFunc)
    ngrid =  length(expVarFunc[[1]]$workGrid)
    
    # Check that all functional objects have the same work-grid
    if (nX >1){
      for (i in 2:nX){
        if ( !all.equal(  expVarFunc[[1]]$workGrid,  expVarFunc[[i]]$workGrid    )){
          stop('All functional covariates need to have the same support. (still)')
        }
      } 
    }
    
    # Define basic arrays
    KovXYZ = array( 0, dim = c(nX+nZ , ngrid) )
    KovXXZZ = array( 0, dim = c(nX + nZ, nX + nZ, ngrid) )
    betaFuncs = array( 0, dim= c(ngrid, nZ + nX) )
    
    # Deal with Z-variables first if they exist
    if( nZ != 0 ){
      
      # Cross-covariance with X-variables
      for (j in 1:nX){ # X-variable counter
        for (i in 1:nZ){ # Z-variable counter
          
          y <- expVarFunc[[j]]$inputData$y 
          t <- expVarFunc[[j]]$inputData$t  
          KovXXZZ[nX+i,j,] = approx(x = expVarFunc[[1]]$obsGrid, xout = expVarFunc[[1]]$workGrid, y = 
                                      GetCrCovYZ(Z = expVarScal[,i], Ly = y, Lt=t, Ymu = expVarFunc[[1]]$mu, 
                                                 bw = bwScalar)$smoothedCC)$y
          KovXXZZ[j,nX+i,] = KovXXZZ[nX+i,j,];
        }
      }
      # Auto-covariance of Z's
      KovXXZZ[(nX+1):(nX+nZ), (nX+1):(nX+nZ),] = cov(expVarScal) 
      # Cross-covariance between Y and Z's
      for (j in 1:nZ){ # Z-variable counter
        y <- depVar$inputData$y 
        t <- depVar$inputData$t  
        KovXYZ[(nX+j),] = approx(x = expVarFunc[[1]]$obsGrid, xout = expVarFunc[[1]]$workGrid, y = 
                                   GetCrCovYZ(Z = expVarScal[,j], Ly = y, Lt=t, Ymu = depVar$mu, bw = bwScalar)$smoothedCC)$y
      }
    } 
    
    # Deal with X-variables next
    # Cross-covariance with other X-variables
    for (j in 1:nX){ # Xj-variable counter
      y1  <- expVarFunc[[j]]$inputData$y 
      t1  <- expVarFunc[[j]]$inputData$t  
      mu1 <- expVarFunc[[j]]$mu
      for (i in j:nX){ # Xi-variable counter
        if (i == j){
          myDiag = diag(expVarFunc[[j]]$fittedCov)
        } else {
          
          y2 <- expVarFunc[[i]]$inputData$y 
          t2 <- expVarFunc[[i]]$inputData$t   
          mu2 <- expVarFunc[[i]]$mu   
          print('x1-x2')
          myDiag = diag(GetCrCovYX( Ly1 = y1, Lt1 = t1, Ly2 = y2, Lt2 = t2, useGAM = splineSmooth, bw1 = bwFunct, bw2 = bwFunct,
                                    Ymu1 =mu1, Ymu2 = mu2)$smoothedCC)
        }
        KovXXZZ[i,j,] = myDiag 
        KovXXZZ[j,i,] = myDiag 
      }
      
      yY <- depVar$inputData$y 
      tY <- depVar$inputData$t   
      print('x-y')
      myDiag = diag(GetCrCovYX( Ly1 = y1, Lt1 = t1, Ly2 = yY, Lt2 = tY, useGAM = splineSmooth, bw1 = bwFunct, bw2 = bwFunct,
                                Ymu1 = mu1, Ymu2 = depVar$mu)$smoothedCC)
      KovXYZ[j,] = myDiag
    }
    
    for (j in 1:length(expVarFunc[[1]]$workGrid)){
      betaFuncs[j,] = solve ( a =  KovXXZZ[,,j], b= KovXYZ[,j] )
    } 
    browser()
    # Get beta0
    beta0 = approx(xout = depVar$workGrid, y = depVar$mu, x = depVar$obsGrid)$y - 0;
    for (j in seq(1, nX, 1)){ 
      beta0 = beta0 - betaFuncs[,j] * approx(xout = depVar$workGrid, y =  expVarFunc[[j]]$mu, x = depVar$obsGrid)$y 
    }
    for (j in seq(1, nZ, 1)){ 
      beta0 = beta0  - betaFuncs[,nX+j] *  mean(expVarScal[,j]);
    }
    
    # Get fitted values
    if (getFitted){
      yhat =  lapply(  depVar$inputData$y , function(x) rep(0, length(x)))
      # For each sample-point
      for (i in seq(1, length(depVar$inputData$t), 1)){
        
        # Get the time-readings and estimate the associated beta vectors 
        ti = depVar$inputData$t[[i]]  
        
        for (j in seq(1, nX, 1)){
          betaij = approx( x = expVarFunc[[1]]$workGrid, y = betaFuncs[,j], xout = ti)$y
          yhat[[i]] = yhat[[i]]  + betaij *  expVarFunc[[j]]$inputData$y[[i]];
        }
        
        for (j in seq(1, nZ, 1)){
          betaij = approx( x = expVarFunc[[1]]$workGrid, y = betaFuncs[,(nX+j)], xout = ti)$y
          yhat[[i]] = yhat[[i]]  + betaij * expVarScal[i,j];
        }
        
        yhat[[i]] =  yhat[[i]] + approx( x =  depVar$workGrid, y = beta0, xout = ti)$y 
        
      }
    }  
    
    FRegObj <- list(betaFunctions = betaFuncs, regMats = list( KovXXZZ = KovXXZZ, KovXYZ = KovXYZ),
                    fittedValues = yhat, beta0 = beta0)
    return(FRegObj)
  }
  
}

