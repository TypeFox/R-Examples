cesEst <- function( yName, xNames, data, tName = NULL, vrs = FALSE,
      method = "LM", start = NULL, lower = NULL, upper = NULL,
      multErr = FALSE,
      rho1 = NULL, rho2 = NULL, rho = NULL, returnGridAll = FALSE, 
      random.seed = 123, rhoApprox = c( y = 5e-6, gamma = 5e-6, delta = 5e-6, 
         rho = 1e-3, nu = 5e-6 ), ... ) {

   # y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-nu/rho)
   # s = 1 / ( 1 + rho )

   # check time variable
   if( is.null( tName ) ) {
      withTime <- FALSE
   } else {
      withTime <- TRUE
      if( length( tName ) != 1 || !is.character( tName[1] ) ) {
         stop( "argument 'tName' must be a single character string" )
      }
   }

   # check multErr
   if( length( multErr ) != 1 || !is.logical( multErr )[ 1 ] ) {
      stop( "argument 'multErr' must be a single logical value" )
   }

   checkNames( c( yName, xNames, tName ), names( data ) )

   # abbreviated method
   if( method == "NM" ) {
      method <- "Nelder-Mead"
   }

   # number of exogenous variables
   nExog <- length( xNames )
   if( nExog == 2 ) {
      nested <- FALSE
   } else if( nExog %in% c( 3, 4 ) ) {
      nested <- TRUE
   } else {
      stop( "currently, the CES can be estimated",
         " with two inputs (normal non-nested CES)",
         " or with three or four inputs (nested CES) only" )
   }

   # obtain names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, 
      returnRho1 = is.null( rho1 ), returnRho2 = is.null( rho2 ), 
      returnRho = is.null( rho ), nested = nested, withTime = withTime )

   # check rhoApprox
   if( !nested ) {
      rhoApprox <- cesCheckRhoApprox( rhoApprox = rhoApprox, withY = TRUE,
         withDeriv = TRUE )
   }

   # checking rho1
   if( !is.null( rho1 ) ) {
      if( !nested ) {
         stop( "argument 'rho1' can be used only for estimating",
            " nested CES functions" )
      } else if( !is.numeric( rho1 ) ) {
         stop( "argument 'rho1' must be either 'NULL' or numeric" )
      } else if( min( rho1 ) < -1 ) {
         stop( "the rho1s specified in argument 'rho1'",
            " must not be smaller than '-1'" )
      }
   }

   # checking rho2
   if( !is.null( rho2 ) ) {
      if( !nested ) {
         stop( "argument 'rho2' can be used only for estimating",
            " nested CES functions" )
      } else if( !is.numeric( rho2 ) ) {
         stop( "argument 'rho2' must be either 'NULL' or numeric" )
      } else if( min( rho2 ) < -1 ) {
         stop( "the rho2s specified in argument 'rho2'",
            " must not be smaller than '-1'" )
      }
   }

   # checking "rho"
   if( !is.null( rho ) ) {
      if( !is.numeric( rho ) ) {
         stop( "argument 'rho' must be either 'NULL' or numeric" )
      } else if( min( rho ) < -1 ) {
         stop( "the rhos specified in argument 'rho'",
            " must not be smaller than '-1'" )
      }
   }

   # grid search for rho_1 and rho
   if( length( rho1 ) > 1 || length( rho2 ) > 1  || length( rho ) > 1 ) {
      result <- cesEstGridRho( yName = yName, xNames = xNames, tName = tName,
         data = data, vrs = vrs, method = method, start = start,
         lower = lower, upper = upper, multErr = multErr,
         rho1Values = rho1, rho2Values = rho2, rhoValues = rho, 
         returnAll = returnGridAll,
         random.seed = random.seed, rhoApprox = rhoApprox, ... )
      result$call <- match.call()
      return( result )
   }
   
   # number of parameters
   nParam <- 1 + nExog + vrs + nested * 2 - ( !is.null( rho1 ) ) - 
      ( !is.null( rho2 ) ) - ( !is.null( rho ) ) - ( nested && nExog == 3 ) +
      withTime

   # start values
   start <- cesEstStart( yName = yName, xNames = xNames, data = data,
      tName = tName, vrs = vrs, method = method, start = start, 
      rho1 = rho1, rho2 = rho2, rho = rho, nParam = nParam, nested = nested,
      multErr = multErr )

   # dertermining lower and upper bounds automatically
   if( is.null( lower ) ) {
      lower <- cesCoefBounds( vrs = vrs, returnRho1 = is.null( rho1 ), 
         returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
         method = method, lower = TRUE, nExog = nExog, nested = nested,
         withTime = withTime )
   }
   if( is.null( upper ) ) {
      upper <- cesCoefBounds( vrs = vrs, returnRho1 = is.null( rho1 ), 
         returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
         method = method, lower = FALSE, nExog = nExog, nested = nested,
         withTime = withTime )
   }

   # checking lower and upper bounds
   if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
      if( length( lower ) > 1 && length( lower ) != nParam ) {
         stop( "the lower bound has ", length( lower ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( start < lower ) ) {
         stop( "at least one starting value is smaller than its lower bound" )
      }
      if( length( upper ) > 1 && length( upper ) != nParam ) {
         stop( "the upper bound has ", length( upper ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( start > upper ) ) {
         stop( "at least one starting value is greater than its upper bound" )
      }
      if( length( lower ) == length( upper ) ) {
         if( any( lower > upper ) ) {
            stop( "at least one lower bound is greater than its upper bound" )
         }
      }
   } else if( max( lower ) != -Inf || min( upper ) != Inf ) {
      warning( "lower and upper bounds are ignored in method '", method, "'" )
      lower <- -Inf
      upper <- Inf
   }

   # store the (matched) call
   matchedCall <- match.call()

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by SANN and DE)
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   # prepare list that will be returned
   result <- list()

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      if( nested ) {
         stop( "method 'Kmenta' is not supported for nested models" )
      }
      if( withTime ) {
         stop( "method 'Kmenta' cannot estimate models with time /",
            " technical change yet." )
      }
      if( !is.null( rho ) ) {
         stop( "fixing 'rho' is currently not supported for the",
            " Kmenta approximation" )
      }
      if( multErr ) {
         warning( "method 'Kmenta' ignores argument 'multErr'" )
      }
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data,
         vrs = vrs )
   } else if( method %in% c( "Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B" ) ) {
      if( method %in% c( "Nelder-Mead", "SANN" ) ) {
         result$optim <- optim( par = start, fn = cesRss, data = data,
            method = method, yName = yName, xNames = xNames, vrs = vrs,
            rho1 = rho1, rho2 = rho2, rho = rho, rhoApprox = rhoApprox, 
            nested = nested, tName = tName, multErr = multErr, ... )
      } else {
         result$optim <- optim( par = start, fn = cesRss, gr = cesRssDeriv,
            data = data, method = method, lower = lower, upper = upper, 
            yName = yName, xNames = xNames, vrs = vrs, rho1 = rho1, 
            rho2 = rho2, rho = rho, rhoApprox = rhoApprox, 
            nested = nested, tName = tName, multErr = multErr, ... )
      }
      result$coefficients <- result$optim$par
      result$iter <- result$optim$counts[ !is.na( result$optim$counts ) ]
      if( length( result$iter ) == 1 ) {
         result$iter <- unname( result$iter )
      }
      if( method != "SANN" ) {
         result$convergence <- result$optim$convergence == 0
      }
      result$message <- result$optim$message
      rss <- result$optim$value
   } else if( method == "LM" ) {
      # residual function
      residFun <- function( par, yName, xNames, data, vrs, rho1, rho2, rho,
            rhoApprox, nested, tName, multErr ) {
         # add coefficients rho_1, rho_2, and rho, if they are fixed
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1,
            rho2 = rho2, rho = rho, nExog = length( xNames ), nested = nested )
         if( multErr ) {
            result <- log( data[[ yName ]] ) - 
               log( cesCalc( xNames = xNames, tName = tName,
                  data = data, coef = par, rhoApprox = rhoApprox[ "y" ], 
                  nested = nested ) )
         } else {
            result <- data[[ yName ]] - cesCalc( xNames = xNames, tName = tName,
               data = data, coef = par, rhoApprox = rhoApprox[ "y" ], 
               nested = nested )
         }
         return( result )
      }

      # jacobian function
      jac <- function( par, yName, xNames, data, vrs, rho1, rho2, rho,
            rhoApprox, nested, tName, multErr ) {
         # add coefficients rho_1, rho_2, and rho, if they are fixed
         nPar <- length( par )
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1, 
            rho2 = rho2, rho = rho, nExog = length( xNames ), nested = nested )
         if( multErr ) {
            yHat <- cesCalc( xNames = xNames, tName = tName,
               data = data, coef = par, rhoApprox = rhoApprox[ "y" ], 
               nested = nested )
            return( -c( cesDerivCoef( par = par, xNames = xNames, data = data,
               tName = tName, vrs = vrs, returnRho1 = is.null( rho1 ), 
               returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
               rhoApprox = rhoApprox, nested = nested ) ) /
               rep( yHat, nPar ) )
         } else {
            return( -c( cesDerivCoef( par = par, xNames = xNames, data = data,
               tName = tName, vrs = vrs, returnRho1 = is.null( rho1 ), 
               returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
               rhoApprox = rhoApprox, nested = nested ) ) )
         }
      }

      # perform fit
      result$nls.lm <- nls.lm( par = start, fn = residFun, data = data,
         jac = jac, yName = yName, xNames = xNames, vrs = vrs, rho1 = rho1, 
         rho2 = rho2, rho = rho, rhoApprox = rhoApprox, nested = nested,
         tName = tName, multErr = multErr, ... )
      result$coefficients <- result$nls.lm$par
      result$iter <- result$nls.lm$niter
      result$convergence <- result$nls.lm$info > 0 && result$nls.lm$info < 5
      result$message <- result$nls.lm$message
      rss <- result$nls.lm$deviance
   } else if( method == "Newton" ) {
      cesRss2 <- function( par, yName, xNames, data, vrs, rho1, rho2, rho, 
            rhoApprox, nested, tName, multErr ) {
         result <- cesRss( par = par, yName = yName, xNames = xNames,
            data = data, tName = tName, vrs = vrs, multErr = multErr,
            rho1 = rho1, rho2 = rho2, rho = rho, 
            rhoApprox = rhoApprox[ "y" ], nested = nested )
         attributes( result )$gradient <- cesRssDeriv( par = par, 
            yName = yName, xNames = xNames, tName = tName, data = data, 
            vrs = vrs, multErr = multErr, rho1 = rho1, rho2 = rho2, rho = rho, 
            rhoApprox = rhoApprox, nested = nested )
         return( result )
      }
      # save current setting for warning messages and suppress warning messages
      warnSaved <- options()$warn
      options( warn = -1 )
      # perform fit
      result$nlm <- nlm( f = cesRss2, p = start, data = data,
         yName = yName, xNames = xNames, tName = tName, vrs = vrs, 
         multErr = multErr, rho1 = rho1, rho2 = rho2, rho = rho, 
         rhoApprox = rhoApprox, nested = nested, ... )
      # restore previous setting for warning messages
      options( warn = warnSaved )
      # extract results
      result$coefficients <- result$nlm$estimate
      result$iter <- result$nlm$iterations
      result$convergence <- result$nlm$code <= 2
      rss <- result$nlm$minimum
   } else if( method == "PORT" ) {
      result$nlminb <- nlminb( start = start, objective = cesRss,
         gradient = cesRssDeriv, data = data, yName = yName, xNames = xNames,
         tName = tName, vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
         multErr = multErr, lower = lower, upper = upper, 
         rhoApprox = rhoApprox, nested = nested, ... )
      result$coefficients <- result$nlminb$par
      result$iter <- result$nlminb$iterations
      result$convergence <- result$nlminb$convergence == 0
      result$message <- result$nlminb$message
      rss <- result$nlminb$objective
   } else if( method == "DE" ) {
      result$DEoptim <- DEoptim( fn = cesRss, lower = lower,
         upper = upper, data = data, yName = yName, xNames = xNames,
         tName = tName, vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
         rhoApprox = rhoApprox, nested = nested, multErr = multErr, ... )
      result$coefficients <- result$DEoptim$optim$bestmem
      result$iter <- result$DEoptim$optim$iter
      rss <- result$DEoptim$optim$bestval
   } else if( method == "nls" ) {
      if( !is.null( rho1 ) ) {
         warning( "ignoring argument 'rho1'" )
      }
      if( !is.null( rho2 ) ) {
         warning( "ignoring argument 'rho2'" )
      }
      if( !is.null( rho ) ) {
         warning( "ignoring argument 'rho'" )
      }
      if( nested && nExog == 3 ) {
         nlsFormula <- as.formula( paste( 
            ifelse( multErr, "log( ", "" ), yName, ifelse( multErr, " )", "" ),
            " ~ ", ifelse( multErr, "log( ", "" ), "gamma *", 
            ifelse( withTime, paste( " exp( lambda *", tName, ") *" ), "" ),
            " ( delta * ",
            "( delta_1 * ", xNames[ 1 ], "^(-rho_1)",
            " + ( 1 - delta_1 ) * ", xNames[ 2 ], "^(-rho_1) )",
            "^( rho / rho_1 ) +",
            " ( 1 - delta ) * ", xNames[ 3 ], "^(-rho) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )", 
            ifelse( multErr, " )", "" ),
            sep = "" ) )
      } else if( nested && nExog == 4 ) {
         nlsFormula <- as.formula( paste(
            ifelse( multErr, "log( ", "" ), yName, ifelse( multErr, " )", "" ),
            " ~ ", ifelse( multErr, "log( ", "" ), "gamma *",
            ifelse( withTime, paste( " exp( lambda *", tName, ") *" ), "" ),
            " ( delta * ( delta_1 * ", xNames[ 1 ], "^(-rho_1)",
            " + ( 1 - delta_1 ) * ", xNames[ 2 ], "^(-rho_1) )",
            "^( rho / rho_1 ) +",
            " ( 1 - delta ) * ( delta_2 * ", xNames[ 3 ], "^(-rho_2)",
            " + ( 1 - delta_2 ) * ", xNames[ 4 ], "^(-rho_2) )",
            "^( rho / rho_2 ) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )",
            ifelse( multErr, " )", "" ),
            sep = "" ) )
      } else {
         nlsFormula <- as.formula( paste(
            ifelse( multErr, "log( ", "" ), yName, ifelse( multErr, " )", "" ),
            " ~ ", ifelse( multErr, "log( ", "" ), "gamma *",
            ifelse( withTime, paste( " exp( lambda *", tName, ") *" ), "" ),
            " ( delta * ", xNames[ 1 ], "^(-rho)",
            " + ( 1 - delta ) * ", xNames[ 2 ], "^(-rho) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )",
            ifelse( multErr, " )", "" ),
            sep = "" ) )
      }
      result$nls <- nls( formula = nlsFormula, data = data, start = start,
         ... )
      result$coefficients <- coef( result$nls )
      result$iter <- result$nls$convInfo$finIter
      result$convergence <- result$nls$convInfo$isConv
      if( result$nls$convInfo$stopMessage != "converged" ) {
         result$message <- result$nls$convInfo$stopMessage
      }
      rss <- deviance( result$nls )
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', 'LM', 'Newton', 'PORT',",
         " 'DE', 'nls', or 'Kmenta'" )
   }

   # add names to estimated coefficients
   names( result$coefficients ) <- coefNames

   # add rho_1, rho_2, and rho, if they are fixed
   result$coefficients <- cesCoefAddRho( coef = result$coefficients,
      vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
      nExog = nExog, nested = nested )

   # calculate and return (constant!) elasticities of substitution
   if( !nested ) {
      result$ela <- NA
      if( result$coefficients[ "rho" ] >= -1 ) {
         result$ela <- 1 / ( 1 + result$coefficients[ "rho" ] )
      }
      names( result$ela ) <- "E_1_2 (all)"
   } else {
      if( nExog == 3 ) {
         result$ela <- rep( NA, 2 )
         if( result$coefficients[ "rho_1" ] >= -1 ) {
            result$ela[1] <- 1 / ( 1 + result$coefficients[ "rho_1" ] )
         }
         if( result$coefficients[ "rho" ] >= -1 ) {
            result$ela[2] <- 1 / ( 1 + result$coefficients[ "rho" ] )
         }
         names( result$ela ) <- c( "E_1_2 (HM)", "E_(1,2)_3 (AU)" )
      } else if( nExog == 4 ) {
         result$ela <- rep( NA, 3 )
         if( result$coefficients[ "rho_1" ] >= -1 ) {
            result$ela[1] <- 1 / ( 1 + result$coefficients[ "rho_1" ] )
         }
         if( result$coefficients[ "rho_2" ] >= -1 ) {
            result$ela[2] <- 1 / ( 1 + result$coefficients[ "rho_2" ] )
         }
         if( result$coefficients[ "rho" ] >= -1 ) {
            result$ela[3] <- 1 / ( 1 + result$coefficients[ "rho" ] )
         }
         names( result$ela ) <- 
            c( "E_1_2 (HM)", "E_3_4 (HM)", "E_(1,2)_(3,4) (AU)" )
      }
   }

   # return also the call
   result$call <- matchedCall

   # return the method used for the estimation
   result$method <- method

   # return whether the error term is multiplicative
   result$multErr <- multErr

   # return the starting values
   result$start <- start

   # return lower and upper bounds
   result$lower <- lower
   result$upper <- upper

   # return fixed rho_1
   result$rho1 <- rho1

   # return fixed rho_2
   result$rho2 <- rho2

   # return fixed rho
   result$rho <- rho

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, tName = tName, data = data,
      coef = result$coefficients, rhoApprox = rhoApprox[ "y" ], nested = nested )

   # residuals
   if( multErr ) {
      result$residuals <- log( data[[ yName ]] ) - log( result$fitted.values )
   } else {
      result$residuals <- data[[ yName ]] - result$fitted.values
   }

   # sum of squared residuals
   result$rss <- sum( result$residuals^2 )
   if( method != "Kmenta" ) {
      if( !isTRUE( all.equal( rss, result$rss, tolerance = 1e-3 ) ) ) {
         warning( "internal problem: the minimum of the objective function",
            " returned by the solver (", rss, ") is not equal to the",
            " RSS calculated from the residuals (", result$rss, ")" )
      }
   }

   # unscaled covariance matrix
   gradients <- cesDerivCoef( par = result$coefficients, xNames = xNames,
      tName = tName, data = data, vrs = vrs, rhoApprox = rhoApprox, 
      nested = nested )
   if( multErr ) {
      gradients <- gradients / 
         matrix( rep( result$fitted.values, length( result$coefficients ) ),
            nrow = nrow( gradients ), ncol = ncol( gradients ) )
   }
   result$cov.unscaled <- try( chol2inv( chol( crossprod( gradients ) ) ),
      silent = TRUE )
   if( !is.matrix( result$cov.unscaled ) ) {
      result$cov.unscaled <- matrix( NA, nrow = length( result$coefficients ),
         ncol = length( result$coefficients ) )
   }
   rownames( result$cov.unscaled ) <- names( result$coefficients )
   colnames( result$cov.unscaled ) <- names( result$coefficients )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

