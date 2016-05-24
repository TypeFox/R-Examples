## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, newdata = NULL,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- length( coef( object ) ) != object$rank
         # TRUE if there are restrictions imposed
   }

   for(i in 1:length( object$eq ) )  {
      predicted.i <- predict( object$eq[[ i ]], newdata = newdata,
         se.fit = se.fit, se.pred = se.pred, interval = interval,
         level = level, useDfSys = useDfSys )
      names( predicted.i ) <- paste( object$eq[[ i ]]$eqnLabel, ".",
         names( predicted.i ), sep = "" )
      if( i == 1 ) {
         predicted <- predicted.i
      } else {
         predicted <- cbind( predicted, predicted.i )
      }
   }
   names( predicted ) <- sub( "(?<!\\.se)\\.fit$", ".pred",
         names( predicted ), perl = TRUE )

   return( predicted )
}

## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit.equation <- function( object, newdata = NULL,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- object$nCoef.sys != object$rank.sys
         # TRUE if there are restrictions imposed
   }

   if( is.null( newdata ) ) {
      xMatNoNa <-  model.matrix( object )
      xMat <- matrix( NA, nrow = length( residuals( object ) ),
         ncol = ncol( xMatNoNa ) )
      xMat[ !is.na( residuals( object ) ), ] <- xMatNoNa
      rm( xMatNoNa )
   } else {
      xMat <-  model.matrix( formula( delete.response( object$terms ) ),
         data = newdata )
   }

   # fitted values
   predicted <- data.frame( fit = drop( xMat %*% object$coefficients ) )

   # calculate variance covariance matrices
   if( se.fit | interval == "confidence" ) {
      yCovConf <- drop( xMat %*% object$coefCov %*% t( xMat ) )
   }
   if( se.pred | interval == "prediction" ) {
      sigmaSqr <- sum( residuals( object, na.rm = TRUE )^2 ) / df.residual( object )
      yCovPred <- drop( xMat %*% object$coefCov %*% t( xMat ) + sigmaSqr )
   }
   # standard errors of fitted values
   if( se.fit ) {
      if( length( yCovConf ) == 1 ) {
         predicted[[ "se.fit" ]] <- sqrt( yCovConf )
      } else {
         predicted[[ "se.fit" ]] <- sqrt( diag( yCovConf ) )
      }
   }
   # standard errors of prediction
   if( se.pred ) {
      if( length( yCovPred ) == 1 ) {
         predicted[[ "se.pred" ]] <- sqrt( yCovPred )
      } else {
         predicted[[ "se.pred" ]] <- sqrt( diag( yCovPred ) )
      }
   }

   # confidence intervals
   if( interval == "confidence" ) {
      if( useDfSys ) {
         tval   <- qt( 1 - ( 1- level )/2, object$df.residual.sys )
      } else {
         tval   <- qt( 1 - ( 1- level )/2, object$df.residual )
      }
      if(  length( yCovConf ) == 1 ) {
         stdErConf <- sqrt( yCovConf )
      } else {
         stdErConf <- sqrt( diag( yCovConf ) )
      }
      predicted[[ "lwr" ]] <- predicted$fit - ( tval * stdErConf )
      predicted[[ "upr" ]] <- predicted$fit + ( tval * stdErConf )
   }
   # prediction intervals
   if( interval == "prediction" ) {
      if( useDfSys ) {
         tval   <- qt( 1 - ( 1- level )/2, object$df.residual.sys )
      } else {
         tval   <- qt( 1 - ( 1- level )/2, object$df.residual )
      }
      if( length( yCovPred ) == 1 ) {
         stdErPred <- sqrt( yCovPred )
      } else {
         stdErPred <- sqrt( diag( yCovPred ) )
      }
      predicted[[ "lwr" ]] <- predicted$fit - ( tval * stdErPred )
      predicted[[ "upr" ]] <- predicted$fit + ( tval * stdErPred )
   }

   return( predicted )
}

