translogCostEst <- function( cName, yName, pNames, data, 
   fNames = NULL, shifterNames = NULL,
   dataLogged = FALSE, homPrice = TRUE, ... ) {

   checkNames( c( cName, yName, pNames, fNames, shifterNames ), names( data ) )

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( cName, yName, pNames, fNames ), 
         varNamesNum = shifterNames )
   }

   qxVars <- yName
   if( homPrice ){
      stop( "imposing linear homogeneity in input prices has not been",
         " implemented yet" )
   } else {
      qxVars <- c( qxVars, pNames )
   }
   qxVars <- c( qxVars, fNames )
      
   result <- quadFuncEst( yName = cName, xNames = qxVars, 
      data = logData, shifterNames = shifterNames, ... )

   result$r2nonLog <- rSquared( exp( logData[[ yName ]] ),
      exp( logData[[ yName ]] ) - exp( result$fitted ) )

   if( !dataLogged ){
      result$fitted <- exp( result$fitted )
   }

   result$call <- match.call()
   result$cName         <- cName
   result$yName         <- yName
   result$pNames        <- pNames
   result$fNames        <- fNames
   result$shifterNames  <- shifterNames
   result$dataLogged    <- dataLogged
   result$homPrice      <- homPrice
   result$xNames        <- NULL
   result$homWeights    <- NULL
   result$regScale      <- NULL

   class( result ) <- "translogCostEst"
   return( result )
}
