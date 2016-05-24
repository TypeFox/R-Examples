## prepare summary results that belong to the whole system
summary.systemfit <- function( object, useDfSys = NULL,
      residCov = TRUE, equations = TRUE, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- length( coef( object ) ) != object$rank
         # TRUE if there are restrictions imposed
   }

   # number of equations
   nEq <- length( object$eq )
   # number of observations per equation
   nObsEq <- rep( NA, nEq )
   validObsEq <- matrix( NA, nrow = nrow( residuals( object ) ), ncol = nEq )
   for( i in 1:nEq ) {
      nObsEq[ i ] <- length( residuals( object$eq[[ i ]], na.rm = TRUE ) )
      validObsEq[ , i ] <- !is.na( residuals( object$eq[[ i ]] ) )
   }
   # total number of observations
   nObs <- sum( nObsEq )

   # preparing objects that will be returned
   result <- list()
   result$call <- object$call
   result$method <- object$method
   result$iter <- object$iter
   result$control <- object$control
   result$residuals <- residuals( object )
   result$residCovEst <- object$residCovEst
   result$residCov <- object$residCov
   if( !is.null( result$residCovEst ) ) {
      dimnames( result$residCovEst ) <- dimnames( result$residCov )
   }
   result$residCor <- cor( residuals( object ), use = "na.or.complete" )
   dimnames( result$residCor ) <- dimnames( result$residCov )
   result$detResidCov <- det( object$residCov, tol = object$control$solvetol )

   # now prepare summury results for the individual equations
   result$eq <- list()
   for( i in 1:length( object$eq ) ) {
       result$eq[[ i ]] <- summary( object$eq[[i]], useDfSys = useDfSys )
   }

   # coefficients, standard errors, ... 
   result$coefCov <- object$coefCov
   coef <- object$coefficients
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( useDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual ) )
   } else {
      pVal <- NULL
      for( i in 1:length( object$eq ) ){
         pVal <- c( pVal, coef( result$eq[[ i ]] )[ , 4 ] )
      }
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( object$rank, object$df.residual )

   # coefficients of the modified regressor matrix
   if( !is.null( object$restrict.regMat ) ) {
      coefModReg <- coef( object, modified.regMat = TRUE )
      stdErModReg <- diag( vcov( object, modified.regMat = TRUE ) )^0.5  # standard errors
      tStatModReg <- coefModReg / stdErModReg    # t-statistic
      if( useDfSys ) {             # p-values
         pValModReg <- 2 * ( 1 - pt( abs( tStatModReg ), object$df.residual ) )
      } else {
         pValModReg <- rep( NA, length( coefModReg ) )
      }
      result$coefModReg <- cbind( coefModReg, stdErModReg, tStatModReg, pValModReg )
      colnames( result$coefModReg ) <- c( "Estimate", "Std. Error",
         "t value", "Pr(>|t|)" )
   }

   # R^2 values
   resid <- NULL
   response <- NULL
   responseMinusMean <- NULL
   for( i in 1:length( object$eq ) ) {
      resid <- c( resid, residuals( object$eq[[ i ]], na.rm = TRUE ) )
      responseEqI <- fitted( object$eq[[ i ]], na.rm = TRUE ) +
         residuals( object$eq[[ i ]], na.rm = TRUE )
      response <- c( response, responseEqI )
      responseMinusMean <- c( responseMinusMean,
         responseEqI - mean( responseEqI ) )
   }

   # OLS R^2 value of the entire system
   rss <- sum( resid^2 )
   tss <- sum( responseMinusMean^2 )
   result$ols.r.squared <- 1 - rss / tss

   # System R^2 value of McElroy (1977)
   # formula from Greene (2003, p. 345 )
   # (first formula, numerator modified to save memory)
   xMat <- matrix( resid, ncol = 1 )
   if( object$control$useMatrix ){
      object$residCov <- as( object$residCov, "dspMatrix" )
      xMat <- as( xMat, "dgCMatrix" )
   }
   rtOmega <- .calcXtOmegaInv( xMat = xMat,
      sigma = object$residCov, validObsEq = validObsEq,
      solvetol = object$control$solvetol,
      useMatrix = object$control$useMatrix )
   yCov <- .calcResidCov( response, methodResidCov = "noDfCor",
      validObsEq = validObsEq, centered = TRUE,
      solvetol = object$control$solvetol )
   residCovInv <- solve( object$residCov, tol = object$control$solvetol )
   denominator <- 0
   for( i in 1:length( object$eq ) ) {
      for( j in 1:length( object$eq ) ) {
         denominator <- denominator + residCovInv[ i, j ] * yCov[ i, j ] *
            ( nObsEq[ i ] * nObsEq[ j ] )^0.5
      }
   }
   result$mcelroy.r.squared <- drop( 1 - ( rtOmega %*% resid ) / denominator )

   result$printEquations <- equations
   result$printResidCov  <- residCov

   class( result ) <- "summary.systemfit"
   return( result )
}

## print summary results of the whole system
print.summary.systemfit <- function( x,
      digits = max( 3, getOption("digits") - 1 ),
      residCov = x$printResidCov, equations = x$printEquations, ... ) {

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("systemfit results \n")
  cat("method: ")
  if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
  cat( paste( x$method, "\n\n"))
  if(!is.null(x$iter)) {
    if(x$iter>1) {
      if(x$iter<x$control$maxiter) {
        cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
      } else {
        cat( paste( "warning: convergence not achieved after", x$iter,
                    "iterations\n\n" ) )
      }
    }
  }

   table.sys <- cbind( round( sum( x$df ),         digits ),
                       round( x$df[2],             digits ),
                       round( sum( sapply( x$eq, function( x ) x$ssr ) ), digits ),
                       round( x$detResidCov,       digits ),
                       round( x$ols.r.squared,     digits ),
                       round( x$mcelroy.r.squared, digits ) )
   rownames( table.sys ) <- c( "system" )
   colnames( table.sys ) <- c( "N", "DF", "SSR", "detRCov", "OLS-R2", "McElroy-R2" )
   print( table.sys, quote = FALSE, right = TRUE, digits = digits )

   cat("\n")

  for(i in 1:length( x$eq ) ) {
    row <- NULL
    row <- cbind( round( sum( x$eq[[i]]$df ),  digits ),
                  round( x$eq[[i]]$df[2], digits ),
                  round( x$eq[[i]]$ssr,   digits ),
                  round( x$eq[[i]]$sigma^2, digits ),
                  round( x$eq[[i]]$sigma,   digits ),
                  round( x$eq[[i]]$r.squared,     digits ),
                  round( x$eq[[i]]$adj.r.squared, digits ))
    table  <- rbind( table, row )
    labels <- rbind( labels, x$eq[[i]]$eqnLabel )
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )
  print(table, quote = FALSE, right = TRUE, digits = digits )

  cat("\n")

   if( residCov ){
      if(!is.null(x$residCovEst)) {
         cat("The covariance matrix of the residuals used for estimation\n")
         print( x$residCovEst, digits = digits )
         cat("\n")
         if( min(eigen( x$residCov, only.values=TRUE)$values) < 0 ) {
            cat("warning: this covariance matrix is NOT positive semidefinit!\n")
            cat("\n")
         }
      }

      cat("The covariance matrix of the residuals\n")
      print( x$residCov, digits = digits )
      cat("\n")

      cat("The correlations of the residuals\n")
      print( x$residCor, digits = digits )
      cat("\n")
   }

   if( equations ){
      ## now print the individual equations
      for(i in 1:length( x$eq ) ) {
         print( x$eq[[i]], digits = digits )
      }
   } else {
      cat( "\nCoefficients:\n" )
      printCoefmat( coef( x ), digits = digits )
   }

  invisible( x )
}


## prepare summary results for a single equation
summary.systemfit.equation <- function( object, useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- object$nCoef.sys != object$rank.sys
         # TRUE if there are restrictions imposed
   }

   # number of observations in this equation
   nObs <- length( residuals( object, na.rm = TRUE ) )

   # preparing objects that will be returned
   result <- list()
   result$eqnLabel <- object$eqnLabel
   result$eqnNo <- object$eqnNo
   result$terms <- object$terms
   result$instruments <- object$inst
   result$method <- object$method
   result$residuals <- object$residuals

   # coefficients, standard errors, ...
   result$coefCov <- object$coefCov
   coef <- object$coefficients
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( useDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual.sys ) )
   } else {
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual ) )
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( length( coef( object ) ), nObs - length( coef( object ) ) )
   result$ssr <- sum( residuals( object, na.rm = TRUE )^2 )
   result$sigma <- sqrt( result$ssr / df.residual( object ) )

   # R^2 values
   response <- fitted( object, na.rm = TRUE ) + residuals( object, na.rm = TRUE )
   rss <- sum( residuals( object, na.rm = TRUE )^2 )
   tss <- sum( ( response - mean( response ) )^2 )
   result$r.squared <- 1 - rss / tss
   result$adj.r.squared <- 1 - ( ( nObs - 1 ) / object$df.residual ) *
      ( 1 - result$r.squared )
   class( result ) <- "summary.systemfit.equation"
   return( result )
}


## print summary results for a single equation
print.summary.systemfit.equation <- function( x,
   digits = max( 3, getOption("digits") - 1 ), ... ) {

  cat("\n")
  cat( x$method, " estimates for '", x$eqnLabel,
         "' (equation ", x$eqnNo, ")\n", sep = "" )

  cat("Model Formula: ")
  print( formula( x$terms ) )
  if(!is.null(x$inst)) {
    cat("Instruments: ")
    print(x$inst)
  }
  cat("\n")

  printCoefmat( x$coefficients, digits = digits )

  cat(paste("\nResidual standard error:", round( x$sigma, digits ),
            "on", x$df[ 2 ], "degrees of freedom\n" ))

  cat( paste( "Number of observations:", round( sum( x$df ), digits ),
              "Degrees of Freedom:", round( x$df[ 2 ], digits ),"\n" ) )

  cat( paste( "SSR:", round( x$ssr, digits ),
              "MSE:", round( x$sigma^2, digits ),
              "Root MSE:", round(x$sigma, digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round( x$r.squared, digits ),
              "Adjusted R-Squared:", round( x$adj.r.squared, digits ),
              "\n" ) )
  cat("\n")
  invisible( x )
}
