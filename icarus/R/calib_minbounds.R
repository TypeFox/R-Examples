# copyright (C) 2016 A.Rebecq
### Functions solving calibration with minimum bounds (for bounded distances).
### Are all private, used by the main "calibration" function

## Search for tight bounds by simplex algorithm
solveMinBoundsCalib <- function(Xs, d, total, q=NULL,
                           maxIter=500, calibTolerance=1e-06, description=TRUE) {

  if (!requireNamespace("Rglpk", quietly = TRUE) || !requireNamespace("slam", quietly = TRUE)) {
      stop("Package Rglpk needed for this function to work. Please install it.",
            call. = FALSE)
  }

  n <- length(d)
  oneN <- rep(1,n)
  zeroN <- rep(0,n)
  IdentityN <- diag(1,n)
  B <- t(diag(d) %*% Xs)

  a <- c(0,1,rep(0,n))
  
  A1 <- rbind( cbind(-oneN,-oneN,IdentityN) , cbind(-oneN,oneN,IdentityN) )
  b1 <- rep(0,2*n)
  
  A3 <- matrix(cbind(rep(0,nrow(B)+1),rep(0,nrow(B)+1),rbind(B,zeroN)), ncol= n+2, nrow= nrow(B)+1)
  b3 <- c(total,0)
  
  Amat <- rbind(A1,A3)
  bvec <- c(b1,b3)
  const <- c(rep("<=",n), rep(">=",n), rep("==", length(b3)))

  ## Use sparse matrix: slam is loaded along with Rglpk
  Amat_sparse <- slam::as.simple_triplet_matrix(Amat)
  simplexSolution <- Rglpk::Rglpk_solve_LP(obj=a, mat=Amat_sparse, dir=const, rhs=bvec)

  xSol <- simplexSolution$solution
  center <- xSol[1]
  minBounds <- xSol[2]
  gSol <- xSol[3:(n+2)]
  wSol <- gSol * d

  if(description) {
    writeLines("Solution found for calibration on minimal bounds:")
    writeLines(paste("L =",min(gSol)))
    writeLines(paste("U =",max(gSol)))
  }

  return(gSol)

}

## General function for calibration on tight/min bounds
minBoundsCalib <- function(Xs, d, total, q=NULL,
                           maxIter=500, calibTolerance=1e-06, description=TRUE,
                           precisionBounds=1e-4, forceSimplex=FALSE, forceBisection=FALSE) {

  usedSimplex <- FALSE
  
  ## For matrices containing more than 1e8 elements, do not use simplex
  ## (it might cause memory issues)
  if( ( forceSimplex || (nrow(Xs)*ncol(Xs)) <= 1e8 ) && !forceBisection) {

    gSol <- solveMinBoundsCalib(Xs, d, total, q,
                                maxIter, calibTolerance, description)
    
    usedSimplex <- TRUE

    Lmax <- min(gSol)
    Umin <- max(gSol)

    ## As we are sure there is a solution, we can look for it much further
    ## than we normally would. However, there are sometimes numerical issues
    ## on the simplex algorithm, so we stop after a long time anyway, and switch
    ## to bisection
    maxIter <- 5000
    
  } else {

    Lmax <- 1.0
    Umin <- 1.0

  }


  digitsPrec <- abs(log(precisionBounds,10))

#   print("Bornes test :")
  Ltest1 <- round(Lmax - 5*10**(-digitsPrec-1),digitsPrec)
  Utest1 <- round(Umin + 5*10**(-digitsPrec-1),digitsPrec)

  Ltest <- Ltest1
  Utest <- Utest1

#   print(Ltest)
#   print(Utest)

  gFinal <- calib(Xs=Xs, d=d, total=total, method="logit", bounds = c(Ltest,Utest),
                  maxIter=maxIter, calibTolerance=calibTolerance)

  ## If no convergence, bisection to find the true min Bounds
  if(is.null(gFinal)) {

    gTestMin <- calib(Xs=Xs, d=d, total=total, method="linear", maxIter=maxIter, calibTolerance=calibTolerance)
    Llinear <- min(gTestMin)
    Ulinear <- max(gTestMin)

    ## If bounds were found by the simplex method, bisection is made around these bound.
    ## If not, bisection looks for a minimal solution symmetric wrt 1
    if(!usedSimplex) {
      distTo1 <- pmax((1-Llinear), (Ulinear-1))
      LtestMin <- 1-distTo1
      LtestMax <- 1+distTo1
    } else {
      LtestMin <- Llinear
      LtestMax <- Ulinear
    }

    method <- "logit" ## Default bounded method 
    return(bisectMinBounds(c(LtestMin,LtestMax),c(Ltest1,Utest1),gTestMin,
                                       Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description))

    Ltest <- Ltest1
    Utest <- Utest1

  }

  return(gFinal)

}

## Search for tight bounds by bisection (sub-optimal algorithm)
bisectMinBounds <- function(convergentBounds,minBounds,gFinalSauv,
                            Xs,d,total, method,maxIter, calibTolerance, precisionBounds,
                            description=TRUE) {

  newBounds <- (minBounds + convergentBounds) / 2
  Ltest <- newBounds[1]
  Utest <- newBounds[2]

  if(description) {
    writeLines(paste("------ Bisection search : ",newBounds[1],";",newBounds[2]))
  }

  gFinal <- calib(Xs=Xs, d=d, total=total, method="logit", bounds = c(Ltest,Utest),
                   maxIter=maxIter, calibTolerance=calibTolerance)

  if(is.null(gFinal)) {
    return( bisectMinBounds(convergentBounds,c(Ltest,Utest),gFinalSauv,
                            Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description) )
  } else {

    # if( all(abs(gFinal - gFinalSauv) <= rep(precisionBounds, length(gFinal))) ) {
    if( all(abs(c( (max(gFinal) - max(gFinalSauv)) , (min(gFinal) - min(gFinalSauv)))) <= c(precisionBounds,precisionBounds)) ) {
      return(gFinal)
    } else {
      return( bisectMinBounds(c(Ltest,Utest),minBounds,gFinal,
                              Xs,d,total, method,maxIter, calibTolerance, precisionBounds, description) )
    }

  }

}
