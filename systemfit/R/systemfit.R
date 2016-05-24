###   $Id: systemfit.R 1126 2014-08-14 08:49:51Z arne $
###
###            Simultaneous Equation Estimation for R
###
### Copyright 2002-2004 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###                     Arne Henningsen <http://www.arne-henningsen.de>
###
### This file is part of the nlsystemfit library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


systemfit <- function(  formula,
                        method = "OLS",
                        inst=NULL,
                        data=list(),
                        restrict.matrix = NULL,
                        restrict.rhs = NULL,
                        restrict.regMat = NULL,
                        pooled = FALSE,
                        control = systemfit.control( ... ),
                        ... )
{

   ## determine whether we have panel date and thus a panel-like model
   panelLike <- class( data )[1] == "plm.dim"

   ## checking argument 'formula'
   if( panelLike ){
      if( class( formula ) != "formula" ){
         stop( "argument 'formula' must be an object of class 'formula'",
            " for panel-like models" )
      }
   } else {
      # single-equation models
      if( class( formula ) == "formula" ){
         formula <- list( formula )
      } else if( class( formula ) == "list" ){
         if( !all( lapply( formula, class ) == "formula" ) ){
            stop( "the list of argument 'formula' must",
               " contain only objects of class 'formula'" )
         }
      } else {
         stop( "argument 'formula' must be an object of class 'formula'",
            " or a list of objects of class 'formula'" )
      }
   }

   ## checking argument 'method'
   if(!( method %in% c( "OLS", "WLS", "SUR", "2SLS", "W2SLS", "3SLS",
         "LIML", "FIML" ) ) ){
      stop( "The method must be 'OLS', 'WLS', 'SUR',",
         " '2SLS', 'W2SLS', or '3SLS'" )
   }

   ## prepare model and data for panel-like models
   if( panelLike ) {
      if( !is.null( restrict.regMat ) && pooled ){
         stop( "argument 'restrict.regMat' cannot be used for pooled estimation",
            " of panel-like data" )
      }
      result <- .systemfitPanel( formula = formula,
         data = data, pooled = pooled )
      data <- result$wideData
      formula <- result$eqnSystem
      if( pooled ){
         restrict.regMat <- result$restrict.regMat
      }
   }

   ## checking argument 'inst'
   if( method %in% c( "2SLS", "W2SLS", "3SLS" ) ){
      if( is.null( inst ) ) {
         stop( "The methods '2SLS', 'W2SLS', and '3SLS' need instruments" )
      } else if( class( inst ) == "formula" ){
         if( length( inst ) != 2 ){
            stop( "argument 'inst' must be a one-sided formula" )
         }
         inst <- lapply( c( 1:length( formula ) ), function(x) inst )
      } else if( class( inst ) == "list" ){
         if( length( inst ) != length( formula ) ){
            stop( "if different instruments are specified for each equation,",
               " the length of argument 'inst' must be equal to the number",
               " of equations" )
         }
         if( !is.null( names( inst ) ) && !is.null( names( formula ) ) ){
            if( names( inst ) != names( formula ) ){
               warning( "names of formulas for instruments (argument 'inst')",
                  " are not equal to names of equations (argument 'formula')" )
            }
         }
         if( !all( lapply( inst, class ) == "formula" ) ){
            stop( "the list of argument 'inst' must",
               " contain only objects of class 'formula'" )
         }
         if( !all( lapply( inst, length ) == 2 ) ){
            stop( "the list of argument 'inst' must",
               " contain only one-sided formulas" )
         }
      } else {
         stop( "argument 'inst' must be an object of class 'formula'",
            " or a list of objects of class 'formula'" )
      }
   } else {
      if( !is.null( inst ) ) {
         warning( "The methods 'OLS', 'WLS', and 'SUR' do not need instruments;",
            " ignoring argument 'inst'" )
         inst <- NULL
      }
   }

   ## default value of argument 'singleEqSigma'
   if( is.null( control$singleEqSigma ) ) {
      control$singleEqSigma <- ( is.null( restrict.matrix ) & is.null( restrict.regMat ) )
   }

   ## checking argument 'pooled'
   if( !is.logical( pooled ) || length( pooled ) != 1 ){
      stop( "argument 'pooled' must be logical" )
   }

  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  iter    <- NULL                 # number of iterations
  nEq     <- length( formula )       # number of equations
  ssr     <- numeric( nEq ) # sum of squared residuals of each equation
  sigma   <- numeric( nEq ) # estimated sigma (std. dev. of residuals) of each equation

   if( is.null( names( formula ) ) ) {
      eqnLabels <- paste( "eq", c( 1:nEq ), sep = "" )
   } else {
      eqnLabels <- names( formula )
      if( sum( regexpr( " |_", eqnLabels ) != -1 ) > 0 ) {
         stop( "equation labels may not contain blanks (' ') or underscores ('_')" )
      }
   }

   results$call <- match.call() # get the original call

   ## prepare y vectors and X matrices for each equation
   modelFrame <- .prepareData.systemfit( data )
   # list of terms objects of each equation
   termsEq <- list()
   # terms and model frames for the individual equations
   modelFrameEq <- list()
   # list of evaluated model frames of each equation
   evalModelFrameEq <- list()
   # list for vectors of endogenous variables in each equation
   yVecEq  <- list()
   # list for matrices of regressors in each equation
   xMatEq  <- list()
   # number of exogenous variables /(unrestricted) coefficients in each equation
   nCoefEq <- numeric( nEq )
   # names of observations of each equation (with observations with NAs)
   obsNamesNaEq <- list()
   # names of coefficients
   coefNames  <- NULL
   # names of coefficients of each equation
   coefNamesEq <- list()
   # prepare data for individual equations
   for(i in 1:nEq ) {
      modelFrameEq[[ i ]] <- modelFrame
      modelFrameEq[[ i ]]$formula <- formula[[ i ]]
      evalModelFrameEq[[ i ]] <- eval( modelFrameEq[[ i ]] )
      termsEq[[ i ]] <- attr( evalModelFrameEq[[ i ]], "terms" )
      weights <- model.extract( evalModelFrameEq[[ i ]], "weights" )
      yVecEq[[i]] <- model.extract( evalModelFrameEq[[ i ]], "response" )
      xMatEq[[i]] <- model.matrix( termsEq[[ i ]], evalModelFrameEq[[ i ]] )
      obsNamesNaEq[[ i ]] <- rownames( xMatEq[[ i ]] )
      nCoefEq[i] <- ncol(xMatEq[[i]])
      cNamesEq <- NULL
      for(j in 1:nCoefEq[i]) {
         xjName <- colnames( xMatEq[[ i ]] )[ j ]
         if( panelLike && xjName != "(Intercept)" ){
            coefNames <- c( coefNames, xjName )
            cNamesEq <- c( cNamesEq, sub(
               paste( "^", eqnLabels[ i ], "_", sep = "" ), "", xjName ) )
         } else {
            coefNames <- c( coefNames,
               paste( eqnLabels[ i ], xjName, sep = "_" ) )
            cNamesEq <- c( cNamesEq, xjName )
         }
      }
      coefNamesEq[[ i ]] <- cNamesEq
   }
   rm( modelFrameEq, xjName, cNamesEq )

   ## prepare Z matrices of instruments for each equation
   if( !is.null( inst ) ) {
      # list of terms objects of instruments of each equation
      termsInst <- list()
      # model frame of instruments
      modelFrameInst <- list()
      # evaluated model frame of instruments
      evalModelFrameInst <- list()
      # list for matrices of instruments in each equation
      zMatEq  <- list()
      # prepare data for individual equations
      for(i in 1:nEq) {
         modelFrameInst[[ i ]] <- modelFrame
         modelFrameInst[[ i ]]$formula <- inst[[ i ]]
         evalModelFrameInst[[ i ]] <- eval( modelFrameInst[[ i ]] )
         termsInst[[ i ]] <- attr( evalModelFrameInst[[ i ]], "terms" )
         zMatEq[[i]] <- model.matrix( termsInst[[ i ]], evalModelFrameInst[[ i ]] )
      }
   }

   ## check if all endogenous variables, regressors, and instruments
   ## have the same number of observations (including observations with NAs)
   # total number of observations per equation (including NAs)
   nObsWithNa <- length( yVecEq[[ 1 ]] )
   for( i in 1:nEq ){
      if( nObsWithNa != length( yVecEq[[ i ]] ) ) {
         stop( "all equations must have the same number of observations",
            " (including observations with NAs)",
            " but the endogenous variable of equation 1 has ",
            nObsWithNa, " observations",
            " while the endogenous variable of equation ", i, " has ",
            length( yVecEq[[ i ]] ), " observations" )
      }
      if( nObsWithNa != nrow( xMatEq[[ i ]] ) ) {
         stop( "the regressors of each equation must have the same number",
            " of observations as the corresponding endogenous variable",
            " (including observations with NAs)",
            " but the regressors of equation ", i, " have ",
            nrow( xMatEq[[ i ]] ), " observations",
            " while the endogenous variable of this equation has ",
            length( yVecEq[[ i ]] ), " observations" )
      }
      if( !is.null( inst ) ) {
         if( nObsWithNa != nrow( zMatEq[[ i ]] ) ) {
            stop( "the instrumental variables of each equation must have",
               " the same number of observations as the corresponding regressors",
               " (including observations with NAs)",
               " but the instrumental variables of equation ", i, " have ",
               nrow( zMatEq[[ i ]] ), " observations",
               " while the regressors of this equation have ",
               nrow( xMatEq[[ i ]] ), " observations" )
         }
      }
   }

   ## determine valid observations of each equation
   # which observations in each equation have no NAs
   validObsEq <- matrix( NA, nrow = nObsWithNa, ncol = nEq )
   for( i in 1:nEq ){
      validObsEq[ , i ] <- !is.na( yVecEq[[ i ]] ) &
         rowSums( is.na( xMatEq[[ i ]] ) ) == 0
      if( !is.null( inst ) ) {
         validObsEq[ , i ] <- validObsEq[ , i ] &
            rowSums( is.na( zMatEq[[ i ]] ) ) == 0
      }
   }

   ## check if the system of equations is unbalanced
   # which observations have no NAs in all equations
   validObsAll <- rowSums( !validObsEq ) == 0
   unbalanced <- FALSE
   for( i in 1:nEq ) {
      if( any( validObsEq[ !validObsAll, i ] ) ) {
         unbalanced <- TRUE
      }
   }
   if( unbalanced ) {
      warning( "the estimation of systems of equations with unequal numbers",
         " of observations has not been thoroughly tested yet" )
   }
   rm( validObsAll, unbalanced )

   ## remove all observations with NAs
   for( i in 1:nEq ) {
      # vectors of endogenous variables
      yVecEq[[ i ]] <- yVecEq[[ i ]][ validObsEq[ , i ] ]
      # matrices of regressors
      attrAssign <- attributes( xMatEq[[ i ]] )$assign
      xMatEq[[ i ]] <- xMatEq[[ i ]][ validObsEq[ , i ], , drop = FALSE ]
      attributes( xMatEq[[ i ]] )$assign <- attrAssign
      # matrices of instrumental variables
      if( !is.null( inst ) ) {
         zMatEq[[ i ]] <- zMatEq[[ i ]][ validObsEq[ , i ], , drop = FALSE ]
      }
   }
   rm( attrAssign )

   ## prepare matrices for using the Matrix package
   if( control$useMatrix ){
      # attributes of the model matrices
      xMatEqAttr <- list()
      for( i in 1:nEq ) {
         xMatEqAttr[[ i ]] <- attributes( xMatEq[[i]] )
         xMatEq[[ i ]] <- as( xMatEq[[ i ]], "dgeMatrix" )
         if( !is.null( inst ) ) {
            zMatEq[[ i ]] <- as( zMatEq[[ i ]], "dgeMatrix" )
         }
      }
   }

   # number of observations in each equation
   nObsEq  <- numeric( nEq )
   # names of observations of each equation (without observations with NAs)
   obsNamesEq <- list()
   for( i in 1:nEq ){
      obsNamesEq[[ i ]] <- rownames( xMatEq[[ i ]] )
      nObsEq[i] <- length( yVecEq[[i]] )
   }

   # stacked vector of all endogenous variables
   yVecAll <- matrix( 0, 0, 1 )
   for( i in 1:nEq ) {
      yVecAll <- c(yVecAll,yVecEq[[i]])
   }

   # stacked matrices of all regressors
   xMatAll <- .stackMatList( xMatEq, way = "diag",
      useMatrix = control$useMatrix )

   # impose restrictions via argument 'restrict.regMat'
   if( !is.null( restrict.regMat ) ) {
      # checking matrix to modify (post-multiply) the regressor matrix (restrict.regMat)
      if( !is.matrix( restrict.regMat ) ) {
         stop( "argument 'restrict.regMat' must be a matrix" )
      }
      if( nrow( restrict.regMat ) != sum( nCoefEq ) ){
         stop( "argument 'restrict.regMat' must be a matrix with number of rows",
            " equal to the number of all regressors [in this model: ",
            sum( nCoefEq ), "]" )
      }
      # default names for modified regressors and their coefficients
      if( is.null( colnames( restrict.regMat ) ) ){
         colnames( restrict.regMat ) <- paste( "C", c( 1:ncol( restrict.regMat ) ), sep = "" )
      }
      # default rownames for matrix to modify regressors
      if( is.null( rownames( restrict.regMat ) ) ){
         rownames( restrict.regMat ) <- coefNames
      }
      # modify regressor matrix (by restrict.regMat)
      XU <- xMatAll
      xMatAll  <- XU %*% restrict.regMat
      if( control$useMatrix ){
         xMatAll <- as( xMatAll, "dgCMatrix" )
      }
   }

   # fitted values of regressors for IV estimations
   if( !is.null( inst ) ) {
      # stacked matrices of all instruments
      zMatAll <- .stackMatList( zMatEq, way = "diag",
         useMatrix = control$useMatrix )
      # fitted values of all regressors
      xMatHatAll <- calcFittedRegMat( xMatAll = xMatAll, zMatEq = zMatEq, 
         nEq = nEq, nObsEq = nObsEq, 
         useMatrix = control$useMatrix, solvetol = control$solvetol )
   }

   # checking and modifying parameter restrictions
   coefNamesModReg <- if( is.null( restrict.regMat ) ) coefNames else colnames( restrict.regMat )
   if( is.character( restrict.matrix ) ) {
      R.restr <- car::makeHypothesis( coefNamesModReg, restrict.matrix, restrict.rhs )
      if( is.null( dim( R.restr ) ) ){
         R.restr <- t( R.restr )
      }
      q.restr <- R.restr[ , ncol( R.restr ), drop = FALSE ]
      R.restr <- R.restr[ , -ncol( R.restr ), drop = FALSE ]
   } else if( !is.null( restrict.matrix ) ) {
      if( is.null( dim( restrict.matrix ) ) ) {
         R.restr <- t( restrict.matrix )
      } else {
         R.restr <- restrict.matrix
      }
      if( is.null( restrict.rhs ) ) {
         q.restr <- matrix( 0, nrow( restrict.matrix ) ,1 )
      } else {
         if( is.null( dim( restrict.rhs ) ) ) {
            q.restr <- matrix( restrict.rhs, ncol = 1  )
         } else {
            q.restr <- restrict.rhs
         }
      }
   } else {
      R.restr <- NULL
      if( !is.null( restrict.rhs ) ) {
         warning( "ignoring argument 'restrict.rhs',",
            " because argument 'restrict.matrix' is not specified" )
      }
      q.restr <- restrict.rhs
   }

   # row names and column names of restriction matrix and vector
   if( !is.null( R.restr ) ){
      if( is.null( rownames( R.restr ) ) ) {
         rownames( R.restr ) <-
            car::printHypothesis( R.restr, q.restr, coefNamesModReg )
      }
      if( is.null( colnames( R.restr ) ) ) {
         colnames( R.restr ) <- coefNamesModReg
      }
      if( is.null( rownames( q.restr ) ) ) {
         rownames( q.restr ) <-
            car::printHypothesis( R.restr, q.restr, coefNamesModReg )
      }
      if( is.null( colnames( q.restr ) ) ) {
         colnames( q.restr ) <- "*rhs*"
      }
   }

   nObsAll  <- sum( nObsEq )  # total number of observations of all equations
   nCoefAll <- sum( nCoefEq ) # total number of exogenous variables/(unrestricted) coefficients in all equations
   nCoefLiAll <- nCoefAll     # total number of linear independent coefficients in all equations
   nCoefLiEq  <- nCoefEq      # total number of linear independent coefficients in each equation
   if(!is.null(restrict.regMat)) {
      nCoefLiAll <- nCoefLiAll - ( nrow( restrict.regMat ) - ncol( restrict.regMat ) )
      for(i in 1:nEq) {
         nCoefLiEq[i] <- ncol(xMatAll)
         for(j in 1: ncol(xMatAll) ) {
            if(sum(xMatAll[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])),j]^2)==0) {
               nCoefLiEq[i] <- nCoefLiEq[i]-1
            }
         }
      }
   }
   if(!is.null(R.restr)) {
      nCoefLiAll  <- nCoefLiAll - nrow(R.restr)
      if(is.null(restrict.regMat)) {
         for(j in 1:nrow(R.restr)) {
            for(i in 1:nEq) {  # search for restrictions that are NOT cross-equation
               if( sum( R.restr[ j, (1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))]^2) ==
                   sum(R.restr[j,]^2)) {
                  nCoefLiEq[i] <- nCoefLiEq[i]-1
               }
            }
         }
      }
    }
    df <- nObsEq - nCoefLiEq    # degress of freedom of each equation

  ## only for OLS, WLS and SUR estimation
  if( method %in% c( "OLS", "WLS", "SUR" ) ) {
    if(is.null(R.restr)) {
      coef <- solve( crossprod( xMatAll ), crossprod( xMatAll, yVecAll ), tol=control$solvetol )
               # estimated coefficients
    } else {
      W <- .prepareWmatrix( crossprod( xMatAll ), R.restr,
         useMatrix = control$useMatrix )
      V <- c( as.numeric( crossprod( xMatAll, yVecAll ) ), q.restr )
      if( method == "OLS" || control$residCovRestricted ){
         coef <- solve( W, V, tol=control$solvetol )[ 1:ncol(xMatAll) ]
      } else {
         coef <- solve( crossprod( xMatAll ), crossprod( xMatAll, yVecAll ),
            tol = control$solvetol )
      }
    }
  }

  ## only for OLS estimation
  if(method=="OLS") {
    resids <- yVecAll - xMatAll %*% coef                                        # residuals
    if(control$singleEqSigma) {
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )    # residual covariance matrix
      coefCov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
                    # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nCoefLiAll, methodResidCov = control$methodResidCov )
                           # sigma squared
      if(is.null(R.restr)) {
        coefCov <- s2 * solve( crossprod( xMatAll ), tol=control$solvetol )
                          # coefficient covariance matrix
      } else {
        coefCov <- s2 * solve( W, tol=control$solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coefficient covariance matrix
      }
    }
  }

  ## only for WLS estimation
  if( method %in% c( "WLS" ) ||
      ( method %in% c( "SUR" ) && control$residCovWeighted ) ) {
    coefOld  <- coef # coefficients of previous step
    coefDiff <- coef # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(coefDiff^2)/sum(coefOld^2))^0.5>control$tol & iter < control$maxiter^( method == "WLS" ) ) {
      iter  <- iter+1
      coefOld <- coef                # coefficients of previous step
      resids <- yVecAll - xMatAll %*% coef     # residuals
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      coef  <- .calcGLS( xMat = xMatAll, yVec = yVecAll, R.restr = R.restr,
         q.restr = q.restr, sigma = rcov, validObsEq = validObsEq,
         useMatrix = control$useMatrix, solvetol = control$solvetol )
      coefDiff <- coef - coefOld # difference of coefficients between this and previous step
    }
    coefCov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
       solvetol = control$solvetol )
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## only for SUR estimation
  if( method %in% c( "SUR" ) ) {
    coefOld  <- coef # coefficients of previous step
    coefDiff <- coef # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(coefDiff^2)/sum(coefOld^2))^0.5>control$tol & iter < control$maxiter) {
      iter  <- iter+1
      coefOld <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      coef <- .calcGLS( xMat = xMatAll, yVec = yVecAll,
         R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
         solvetol = control$solvetol )     # coefficients
      coefDiff <- coef - coefOld # difference of coefficients between this and previous step
    }
    coefCov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
       solvetol = control$solvetol )
            # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if( method %in% c( "2SLS", "W2SLS", "3SLS" ) ) {
    if(is.null(R.restr)) {
      coef <- solve( crossprod( xMatHatAll ), crossprod( xMatHatAll, yVecAll ), tol=control$solvetol )
         # 2nd stage coefficients
    } else {
      W <- .prepareWmatrix( crossprod(xMatHatAll), R.restr,
         useMatrix = control$useMatrix )
      V <- c( as.numeric( crossprod( xMatHatAll, yVecAll ) ), q.restr )
      if( method == "2SLS" || control$residCovRestricted ){
         coef <- solve( W, V, tol=control$solvetol )[ 1:ncol(xMatAll) ]
      } else {
         coef <- solve( crossprod( xMatHatAll ), crossprod( xMatHatAll, yVecAll ),
            tol = control$solvetol )
      }
    }
    b2 <- coef
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    if(control$singleEqSigma) {
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      coefCov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
         solvetol = control$solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nCoefLiAll, methodResidCov = control$methodResidCov )
                           # sigma squared
      if(is.null(R.restr)) {
        coefCov <- s2 * solve( crossprod( xMatHatAll ), tol=control$solvetol )
                  # coefficient covariance matrix
      } else {
        coefCov <- s2 * solve( W, tol=control$solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if( method %in% c( "W2SLS" ) ||
         ( method %in% c( "3SLS" ) && control$residCovWeighted ) ) {
    coefOld  <- coef # coefficients of previous step
    coefDiff <- coef # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(coefDiff^2)/sum(coefOld^2))^0.5>control$tol & iter < control$maxiter^( method == "W2LS" ) ) {
      iter  <- iter+1
      coefOld <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
         solvetol = control$solvetol )          # (unrestr.) coeffic.
      coefDiff <- coef - coefOld # difference of coefficients between this and previous step
    }
    coefCov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
       solvetol = control$solvetol )  # coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## only for 3SLS estimation
  if( method %in% c( "3SLS" ) ) {
    coefOld  <- coef # coefficients of previous step
    coefDiff <- coef # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(coefDiff^2)/sum(coefOld^2))^0.5>control$tol & iter < control$maxiter) {
      iter  <- iter+1
      coefOld <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
         validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq,
         centered = control$centerResiduals, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      if(control$method3sls=="GLS") {
         coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll,
            R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
            solvetol = control$solvetol )  # (unrestr.) coeffic.
      }
      if(control$method3sls=="IV") {
         coef <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll, yVec = yVecAll,
            R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
            solvetol = control$solvetol )   # (unrestr.) coeffic.
      }
      if(control$method3sls=="GMM") {
        ZtOmega <- .calcXtOmegaInv( xMat = zMatAll, sigma = rcov,
           validObsEq = validObsEq, invertSigma = FALSE, useMatrix = control$useMatrix,
           solvetol = control$solvetol )
        if(is.null(R.restr)) {
          coef <- as.numeric( solve( crossprod( xMatAll, zMatAll ) %*%
            solve( ZtOmega %*% zMatAll, crossprod( zMatAll, xMatAll ),
            tol=control$solvetol ), crossprod( xMatAll, zMatAll ) %*%
            solve( ZtOmega %*% zMatAll, crossprod( zMatAll, yVecAll ),
            tol=control$solvetol ), tol=control$solvetol ) )
        } else {
          W <- .prepareWmatrix( crossprod( xMatAll, zMatAll ) %*%
               solve( ZtOmega %*% zMatAll, crossprod( zMatAll, xMatAll ),
               tol=control$solvetol ), R.restr,
               useMatrix = control$useMatrix )
          V <- c( as.numeric( crossprod( xMatAll, zMatAll ) %*%
            solve( ZtOmega  %*% zMatAll, crossprod( zMatAll, yVecAll ),
            tol = control$solvetol ) ), q.restr )
          Winv <- solve( W, tol=control$solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll),]     # restricted coefficients
        }
      }
      if(control$method3sls=="Schmidt") {
        xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov,
           validObsEq = validObsEq, useMatrix = control$useMatrix,
           solvetol = control$solvetol )
        if(is.null(R.restr)) {
          coef <- as.numeric( solve( crossprod( xMatHatAll, t( xMatHatOmegaInv ) ), 
            xMatHatOmegaInv %*% zMatAll %*%
            solve( crossprod( zMatAll ), crossprod( zMatAll, yVecAll),
            tol=control$solvetol ), tol=control$solvetol ) )
                           # (unrestr.) coeffic.
        } else {
          W <- .prepareWmatrix( crossprod( xMatHatAll, t( xMatHatOmegaInv ) ),
             R.restr, useMatrix = control$useMatrix )
          V <- c( as.numeric( xMatHatOmegaInv %*% zMatAll %*%
            solve( crossprod( zMatAll ), crossprod( zMatAll, yVecAll ),
            tol = control$solvetol ) ), q.restr )
          Winv <- solve( W, tol=control$solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll),]     # restricted coefficients
        }
      }
      if(control$method3sls=="EViews") {
         coef <- b2 + .calcGLS( xMat = xMatHatAll, yVec = ( yVecAll -  xMatAll %*% b2 ),
            R.restr = R.restr, q.restr = q.restr, sigma = rcov,
            validObsEq = validObsEq, useMatrix = control$useMatrix,
            solvetol = control$solvetol )  # (unrestr.) coeffic.
      }
      coefDiff <- coef - coefOld # difference of coefficients between this and previous step
    }
    if(control$method3sls=="GLS") {
       coefCov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
          solvetol = control$solvetol )  # coefficient covariance matrix
    }
    if(control$method3sls=="IV") {
       coefCov <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll,
          R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
          solvetol = control$solvetol )
    }
    if(control$method3sls=="GMM") {
      if(is.null(R.restr)) {
        coefCov <- solve( crossprod( xMatAll, zMatAll ) %*%
           solve( ZtOmega %*% zMatAll, crossprod( zMatAll, xMatAll ),
           tol=control$solvetol ), tol=control$solvetol )
                # final step coefficient covariance matrix
      } else {
        coefCov <- Winv[1:ncol(xMatAll),1:ncol(xMatAll)] # coefficient covariance matrix
      }
    }
    if(control$method3sls=="Schmidt") {
      xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov,
         validObsEq = validObsEq, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      PH <- zMatAll %*%  solve( crossprod( zMatAll ), t( zMatAll ),
         tol=control$solvetol )
      PHOmega <- .calcXtOmegaInv( xMat = t( PH ), sigma = rcov,
         validObsEq = validObsEq, invertSigma = FALSE, useMatrix = control$useMatrix,
         solvetol = control$solvetol )
      if(is.null(R.restr)) {
         coefCov <- solve( xMatHatOmegaInv %*% xMatHatAll,
            xMatHatOmegaInv %*% PHOmega %*% PH %*%
            crossprod( xMatHatOmegaInv, solve( xMatHatOmegaInv %*% xMatHatAll,
            tol=control$solvetol ) ), tol=control$solvetol )
                  # final step coefficient covariance matrix
      } else {
         VV <- .stackMatList(
            list( xMatHatOmegaInv %*% PHOmega %*% PH %*% t( xMatHatOmegaInv ),
            matrix( 0, nrow( R.restr ), nrow( R.restr ) ) ), "diag",
            useMatrix = control$useMatrix )
         coefCov <- ( Winv %*% VV %*% Winv )[ 1:ncol(xMatAll), 1:ncol(xMatAll) ]
                  # coefficient covariance matrix
      }
    }
    if(control$method3sls=="EViews") {
       coefCov <- .calcGLS( xMat = xMatHatAll,
          R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, validObsEq = validObsEq, useMatrix = control$useMatrix,
          solvetol = control$solvetol )  # final step coefficient covariance matrix
    }
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## FIML estimation
  if( method == "FIML" ) {
    fimlResult <- .systemfitFiml( systemfitCall = results$call, nObsEq = nObsEq,
      validObsEq = validObsEq,
      nCoefEq = nCoefLiEq, yVec = yVecAll, xMat = xMatAll, xEq = xMatEq, methodResidCov = control$methodResidCov,
      centerResiduals = control$centerResiduals, solvetol = control$solvetol )
    #print( fimlResult )
    coef <- fimlResult$coefficients
    coefCov <- fimlResult$coefCov
    resids <- fimlResult$resids
  }

  ## for all estimation methods
  fitted.values <- xMatAll %*% coef   # fitted endogenous values
  if(!is.null(restrict.regMat)) {
    coef  <- restrict.regMat %*% coef
    coefCov <- restrict.regMat %*% coefCov %*% t(restrict.regMat)
  }


  ## equation wise results
  for(i in 1:nEq) {
    results$eq[[ i ]] <- list()
    results$eq[[ i ]]$eqnNo    <- i               # equation number
    results$eq[[ i ]]$eqnLabel <- eqnLabels[[i]]
    results$eq[[ i ]]$method   <- method

    results$eq[[ i ]]$residuals <- rep( NA, nObsWithNa )
    results$eq[[ i ]]$residuals[ validObsEq[ , i ] ] <-
      resids[ ( 1 + sum(nObsEq[1:i]) -nObsEq[i] ):( sum(nObsEq[1:i]) ) ]
    names( results$eq[[ i ]]$residuals ) <- obsNamesNaEq[[ i ]]

    results$eq[[ i ]]$coefficients <-
      drop( coef[(1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))] )
              # estimated coefficients of equation i
    results$eq[[ i ]]$coefCov  <- as.matrix(
      coefCov[(1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i])),
        (1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))] )
              # covariance matrix of estimated coefficients of equation i

    # set names
    names( results$eq[[ i ]]$coefficients )  <- coefNamesEq[[ i ]]
    colnames( results$eq[[ i ]]$coefCov ) <- coefNamesEq[[ i ]]
    rownames( results$eq[[ i ]]$coefCov ) <- coefNamesEq[[ i ]]

    results$eq[[ i ]]$fitted.values <- rep( NA, nObsWithNa )
    results$eq[[ i ]]$fitted.values[ validObsEq[ , i ] ] <-
      fitted.values[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
    names( results$eq[[ i ]]$fitted.values ) <- obsNamesNaEq[[ i ]]

    results$eq[[ i ]]$terms    <- termsEq[[ i ]]
    results$eq[[ i ]]$rank     <- nCoefLiEq[i]
      # rank = number of linear independent coefficients
    results$eq[[ i ]]$nCoef.sys    <- nCoefAll
      # total number of coefficients of the entire system
    results$eq[[ i ]]$rank.sys     <- nCoefLiAll
      # rank = number of linear independent coefficients of the entire system
    results$eq[[ i ]]$df.residual  <- df[i]           # degrees of freedom of residuals
    results$eq[[ i ]]$df.residual.sys  <- nObsAll- nCoefLiAll
       # degrees of freedom of residuals of the whole system
    if( control$y ){
      results$eq[[ i ]]$y   <- yVecEq[[i]]     # vector of endogenous variables
      names( results$eq[[ i ]]$y ) <- obsNamesEq[[ i ]]
    }
    if( control$x ){
      results$eq[[ i ]]$x  <- as.matrix( xMatEq[[i]] )
      if( control$useMatrix ){
         attributes( results$eq[[ i ]]$x ) <- xMatEqAttr[[ i ]]
      }
      rownames( results$eq[[ i ]]$x ) <- obsNamesEq[[ i ]]
    }
    if( control$model ){
      results$eq[[ i ]]$model <- evalModelFrameEq[[ i ]]
         # model frame of this equation
      rownames( results$eq[[ i ]]$model ) <- obsNamesNaEq[[ i ]]
      if( !is.null( inst ) ) {
         results$eq[[ i ]]$modelInst <- evalModelFrameInst[[ i ]]
         rownames( results$eq[[ i ]]$modelInst ) <- obsNamesNaEq[[ i ]]
      }
    }
    if( method %in% c( "2SLS", "W2SLS", "3SLS" ) ) {
      results$eq[[ i ]]$inst         <- inst[[i]]
      results$eq[[ i ]]$termsInst <- termsInst[[i]]
      if(  control$z ) {
         results$eq[[ i ]]$z <- zMatEq[[i]]  # matrix of instrumental variables
         rownames( results$eq[[ i ]]$z ) <- obsNamesEq[[ i ]]
      }
    }
    class( results$eq[[ i ]] ) <- "systemfit.equation"
  }

  ## results of the total system
  # all estimated coefficients
  results$coefficients <- as.numeric( drop( coef ) )
  names( results$coefficients ) <- coefNames

  # coefficients covariance matrix
  results$coefCov <- as.matrix( coefCov )
  colnames( results$coefCov ) <- coefNames
  rownames( results$coefCov ) <- coefNames

  # residual covarance matrix used for estimation
  if( method %in% c( "WLS", "W2SLS", "SUR", "3SLS" ) ){
    results$residCovEst <- as.matrix( rcov )
    colnames( results$residCovEst ) <- eqnLabels
    rownames( results$residCovEst ) <- eqnLabels
  }

  # residual covarance matrix
  results$residCov <- .calcResidCov( resids, methodResidCov = control$methodResidCov,
      validObsEq = validObsEq, nCoefEq = nCoefLiEq, xEq = xMatEq,
      centered = control$centerResiduals, solvetol = control$solvetol )
  colnames( results$residCov ) <- eqnLabels
  rownames( results$residCov ) <- eqnLabels

  results$method  <- method
  results$rank    <- nCoefLiAll
     # rank = total number of linear independent coefficients of all equations
  results$df.residual <- nObsAll - nCoefLiAll
     # degrees of freedom of the whole system
  results$iter    <- iter
  results$restrict.matrix <- R.restr
  results$restrict.rhs <- q.restr
  results$restrict.regMat <- restrict.regMat
  results$control <- control
  results$panelLike <- panelLike
  class(results)  <- "systemfit"

  results
}
