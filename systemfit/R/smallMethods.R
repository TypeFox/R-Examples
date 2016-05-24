## this function returns a vector of the
## cross-equation corrlations between eq i and eq j
## from the results set for equation ij
correlation.systemfit <- function( results, eqni, eqnj ) {
  nCoefEq <- NULL
  for( i in 1:length( results$eq ) ) {
     nCoefEq <- c( nCoefEq, length( coef( results$eq[[ i ]] ) ) )
  }
  cij <- vcov( results )[(1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni])),
                      (1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj]))]
  cii <- vcov( results )[(1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni])),
                      (1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni]))]
  cjj <- vcov( results )[(1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj])),
                      (1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj]))]
  rij <- NULL

  for( i in 1:nrow( residuals( results ) ) ) {
    xik    <- model.matrix( results$eq[[eqni]] )[i,]
    xjk    <- model.matrix( results$eq[[eqnj]] )[i,]
    top    <- xik %*% cij %*% xjk
    bottom <- sqrt( ( xik %*% cii %*% xik ) * ( xjk %*% cjj %*% xjk ) )
    rijk   <- top / bottom
    rij    <- rbind( rij, rijk )
  }
  rij
}

## determines the improvement of resultsj (3sls) over
## resultsi (2sls) for equation i and returns a matrix
## of the values, so you can examine the range, mean, etc
se.ratio.systemfit <- function( resultsi, resultsj, eqni ) {
  ratio <- NULL
  for( i in 1:nrow( residuals( resultsi ) ) ) {
    xik    <- model.matrix( resultsi$eq[[eqni]] )[i,]
    top    <- sqrt( xik %*% vcov( resultsi$eq[[eqni]] ) %*% xik )
    bottom <- sqrt( xik %*% vcov( resultsj$eq[[eqni]] ) %*% xik )
    rk     <- top / bottom
    ratio  <- rbind( ratio, rk )
  }
  ratio
}


## return all coefficients
coef.systemfit <- function( object, modified.regMat = FALSE, ... ) {
   if( modified.regMat ){
      if( is.null( object$restrict.regMat ) ){
         stop( "coefficients of the modified regressor matrix are not available,",
            " because argument 'restrict.regMat' has not been used in this estimation." )
      } else {
         return( drop( solve( crossprod( object$restrict.regMat ),
            t( object$restrict.regMat ) %*% coef( object ) ) ) )
      }
   } else {
      return( object$coefficients )
   }
}

## return all coefficients, std.errors, t-values and p-values
coef.summary.systemfit <- function( object, modified.regMat = FALSE, ... ) {
   if( modified.regMat ){
      if( is.null( object$coefModReg ) ){
         stop( "coefficients of the modified regressor matrix are not available,",
            " because argument 'restrict.regMat' has not been used in this estimation." )
      } else {
         return( object$coefModReg )
      }
   } else {
      return( object$coefficients )
   }
}

## return the coefficients of a single equation
coef.systemfit.equation <- function( object, ... ) {
   object$coefficients
}

## return coefficients, std.errors, t-values and p-values of a single equation
coef.summary.systemfit.equation <- function( object, ... ) {
   object$coefficients
}

## return all residuals
residuals.systemfit <- function( object, ... ) {
   result <- data.frame( obsNo = c( 1:length( residuals( object$eq[[1]] ) ) ) )
   for( i in 1:length( object$eq ) ) {
      result[[ object$eq[[i]]$eqnLabel ]] <- residuals( object$eq[[i]] )
   }
   result$obsNo <- NULL
   rownames( result ) <- names( residuals( object$eq[[ 1 ]] ) )
   return( result )
}

## return residuals of a single equation
residuals.systemfit.equation <- function( object, na.rm = FALSE, ... ) {
   if( na.rm ) {
      return( object$residuals[ !is.na( object$residuals ) ] )
   } else {
      return( object$residuals )
   }
}

## return the variance covariance matrix of the coefficients
vcov.systemfit <- function( object, modified.regMat = FALSE, ... ) {
   if( modified.regMat ){
      if( is.null( object$restrict.regMat ) ){
         stop( "coefficients of the modified regressor matrix",
            " and their covariance matrix are not available,",
            " because argument 'restrict.regMat' has not been used in this estimation." )
      } else {
         txtxInv <- solve( crossprod( object$restrict.regMat ) )
         result <- txtxInv %*% t( object$restrict.regMat ) %*% vcov( object ) %*%
            object$restrict.regMat %*% txtxInv
         return( result )
      }
   } else {
      return( object$coefCov )
   }
}

## return the variance covariance matrix of the coefficients of a single equation
vcov.systemfit.equation <- function( object, ... ) {
   object$coefCov
}


## return the fitted values
fitted.systemfit <- function( object, ... ) {
   nEq <- length( object$eq )
   fitted.values <- matrix( NA, length( object$eq[[1]]$fitted.values ), nEq )
   colnames( fitted.values ) <- as.character( 1:ncol( fitted.values ) )
   for(i in 1:nEq )  {
      fitted.values[ , i ]           <- object$eq[[ i ]]$fitted.values
      colnames( fitted.values )[ i ] <- object$eq[[ i ]]$eqnLabel
   }
   rownames( fitted.values ) <- names( fitted( object$eq[[ 1 ]] ) )
   return( as.data.frame( fitted.values ) )
}

## return the fitted values of e single euation
fitted.systemfit.equation <- function( object, na.rm = FALSE, ... ) {
   if( na.rm ) {
      return( object$fitted.values[ !is.na( object$fitted.values ) ] )
   } else {
      return( object$fitted.values )
   }
}


## return model frame of the entire system
model.frame.systemfit <- function( formula, ... ){
   mfColNames <- NULL
   for( i in 1:length( formula$eq ) ) {
      mfi <- model.frame( formula$eq[[ i ]] )
      if( i == 1 ) {
         result <- mfi
      } else {
         for( j in 1:ncol( mfi ) ) {
            if( ! names( mfi )[ j ] %in% names( result ) ) {
               result[[ names( mfi )[ j ] ]] <- mfi[ , j ]
            }
         }
      }
   }
   return( result )
}

## return model frame of a single equation
model.frame.systemfit.equation <- function( formula, ... ){
   if( !is.null( formula$model ) ) {
      result <- formula$model
   } else {
      stop( "returning model frame not possible. Please re-estimate",
         " the system with control variable 'model'",
         " set to TRUE" )
   }
   return( result )
}
