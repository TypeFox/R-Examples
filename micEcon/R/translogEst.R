translogEst <- function( yName, xNames, data, shifterNames = NULL,
   dataLogged = FALSE, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( yName, xNames ), varNamesNum = shifterNames )
   }

   result <- quadFuncEst( yName = yName, xNames = xNames, data = logData,
      shifterNames = shifterNames, ... )

   result$r2nonLog <- rSquared( exp( logData[[ yName ]] ),
      exp( logData[[ yName ]] ) - exp( result$fitted ) )

   if( !dataLogged ){
      result$fitted <- exp( result$fitted )
   }

   result$call <- match.call()
   result$yName         <- yName
   result$xNames        <- xNames
   result$shifterNames  <- shifterNames
   result$dataLogged    <- dataLogged
   result$homWeights    <- NULL
   result$regScale      <- NULL

   class( result ) <- "translogEst"
   return( result )
}
