## Likelihood function for FIML estimations
.systemfitFimlLik <- function( mlCoef, mlVars ) {

   nObs <- sum( mlVars$nObsEq )
   nEq <- length( mlVars$nObsEq )

   if( length( mlCoef ) != ncol( mlVars$xMat ) ){
      stop( "internal error: argument 'mlCoef' has length ",
         length( mlCoef ), " but must have length ", ncol( mlVars$xMat ) )
   }
   mlResids <- mlVars$yVec - mlVars$xMat %*% mlCoef
   mlSigma <- .calcResidCov( resids = mlResids, methodResidCov = mlVars$methodResidCov,
      validObsEq = mlVars$validObsEq, nCoefEq = mlVars$nCoefEq, xEq = mlVars$xEq,
      centered = mlVars$centerResiduals, diag = FALSE,
      solvetol = mlVars$solvetol )
   likValue <- - nObs * log( 2 * pi ) -
      ( nObs / 2 ) * log( det( mlSigma ) ) -
      ( 1 / 2 ) * .calcXtOmegaInv( mlResids, mlSigma, validObsEq = mlVars$validObsEq,
         solvetol = mlVars$solvetol ) %*% mlResids
   return( likValue )
}

## FIML estimation
.systemfitFiml <- function( systemfitCall, nObsEq, validObsEq, nCoefEq, yVec, xMat, xEq,
      methodResidCov, centerResiduals, solvetol, startCoef = "ITSUR" ) {

   nObs <- sum( nObsEq )
   nEq <- length( nObsEq )

   # starting values
   if( startCoef %in% c( "OLS", "SUR", "ITSUR" ) ) {
      if( startCoef == "OLS" ) {
         systemfitCall$method <- "OLS"
      } else if( startCoef == "SUR" ) {
         systemfitCall$method <- "SUR"
         systemfitCall$maxiter <- 1
      } else if( startCoef == "ITSUR" ) {
         systemfitCall$method <- "SUR"
         systemfitCall$maxiter <- 500
      }
      startResult <- eval( systemfitCall )
      startCoef <- startResult$b
   }

   # variables passed to .systemfitFimlLik
   mlVars <- list()
   mlVars$nObsEq  <- nObsEq
   mlVars$validObsEq <- validObsEq
   mlVars$nCoefEq <- nCoefEq
   mlVars$yVec    <- yVec
   mlVars$xMat    <- xMat
   mlVars$xEq     <- xEq
   mlVars$methodResidCov  <- methodResidCov
   mlVars$centerResiduals <- centerResiduals
   mlVars$sollvetol       <- solvetol

   # list for results of ML estimation
   mlResult <- list()

   mlResult$optim <- optim( startCoef, .systemfitFimlLik, method = "BFGS",
      control = list( fnscale = -1 ), mlVars = mlVars )
   # methods: "BFGS" "Nelder-Mead"
   mlResult$coefficients <- mlResult$optim$par
   mlResult$coefCov <- diag( 1, length( mlResult$coefficients ) )
   mlResult$resids  <- yVec - xMat %*% mlResult$coefficients

   return( mlResult )
}

