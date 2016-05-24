aidsEst <- function( priceNames, shareNames, totExpName,
      data, method = "LA", priceIndex = "Ls", pxBase = 1,
      hom = TRUE, sym = TRUE,
      shifterNames = NULL, instNames = NULL,
      estMethod = ifelse( is.null( instNames ), "SUR", "3SLS" ),
      ILmaxiter = 50, ILtol = 1e-5, alpha0 = 0, restrict.regMat = FALSE, ... ) {

   if( length( priceNames ) != length( shareNames ) ) {
      stop( "arguments 'priceNames' and 'shareNames' must have the same length" )
   }
   nGoods <- length( priceNames )
   nShifter <- length( shifterNames )
   extractPx <- function( method ) {
      px <- substr( method, 4, nchar( method ) )
      return( px )
   }

   if( ! method %in% c( "LA", "IL", "MK" ) ) {
      if( nchar( method ) >= 4 && substr( method, 3, 3 ) == ":" &&
         substr( method, 1, 2 ) %in% c( "LA", "IL", "MK" ) ) {
            priceIndex <- extractPx( method )
            warning( "using price index specified in argument 'method',",
               " ignoring argument 'priceIndex'" )
         method <- substr( method, 1, 2 )
      } else {
         stop( "argument 'method' must be either",
            " 'LA' (for 'Linear Approximation') or",
            " 'IL' (for 'Iterated Linear Least Squares')" )
      }
   } 

   if( !( priceIndex %in% c( "S", "SL", "P", "L", "Ls", "T" ) ) ) {
      stop( "argument 'priceIndex' that specifies the price index must be either",
         " 'S' (Stone index), 'SL' (Stone index with lagges shares),",
         " 'P' (Paasche index), 'L' (Laspeyres index),",
         " 'Ls' (Laspeyres index, simplified), or",
         " 'T' (Tornqvist index)" )
   }

   if( sym && !hom ) {
      hom <- TRUE  # symmetry implies homogeneity
      warning( "symmetry implies homogeneity: imposing additionally homogeniety" )
   }
   allVarNames <- c( priceNames, shareNames, totExpName, instNames, shifterNames )
   if( sum( is.na( data[ , allVarNames ] ) ) > 0 ) {
      warning( "there are some NAs in the data,",
         " all observations (rows) with NAs are excluded from the analysis" )
      data <- data[ !is.na( rowSums( data[ , allVarNames ] ) ), ]
   }
   nObs   <- nrow( data )      # number of observations
   sample <- if( priceIndex == "SL") c( 2:nObs ) else c( 1:nObs )
   result <- list()
   result$call <- match.call()
   wMeans <- numeric( nGoods )  # mean expenditure shares
   pMeans <- numeric( nGoods )  # mean prices
   xtMean <- mean( data[[ totExpName ]][ sample ] )
   for( i in seq( nGoods ) ) {
      wMeans[ i ] <- mean( data[[ shareNames[ i ] ]][ sample ] )
      pMeans[ i ] <- mean( data[[ priceNames[ i ] ]][ sample ] )
   }
   # log of price index
   lnp  <- aidsPx( priceIndex, priceNames, shareNames = shareNames, data = data,
      base = pxBase )
   # prepare data.frame
   sysData <- data.frame( xt = data[[ totExpName ]],
      lxtr = ( log( data[[ totExpName ]] ) - lnp ) )
   for( i in 1:nGoods ) {
      sysData[[ paste( "w", i, sep = "" ) ]] <- data[[ shareNames[ i ] ]]
      sysData[[ paste( "lp", i, sep = "" ) ]] <- log( data[[ priceNames[ i ] ]] )
   }
   if( is.null( instNames )) {
      ivFormula <- NULL
   } else {
      estMethod <- "3SLS"
      ivFormula <- "~"
      for( i in 1:length( instNames ) ) {
         sysData[[ paste( "i", i, sep = "" ) ]] <- data[[ instNames[ i ] ]]
         ivFormula <- paste( ivFormula, " + i", as.character( i ), sep = "" )
      }
      ivFormula <- as.formula( ivFormula )
   }
   if( is.null( shifterNames )) {
      shifterFormula <- NULL
   } else {
      for( i in 1:length( shifterNames ) ) {
         sysData[[ paste( "s", i, sep = "" ) ]] <- data[[ shifterNames[ i ] ]]
      }
   }
   restr <- .aidsRestr( nGoods = nGoods, hom = hom, sym = sym, restrict.regMat = restrict.regMat, nShifter = nShifter )
      # restrictions for homogeneity and symmetry
   system <- .aidsSystem( nGoods = nGoods, nShifter = nShifter ) # LA-AIDS equation system
   # estimate system
   if( restrict.regMat ) {
      est <- systemfit( system, estMethod, data = sysData, restrict.regMat = restr,
         inst = ivFormula, ... )
   } else {
      est <- systemfit( system, estMethod, data = sysData, restrict.matrix = restr,
         inst = ivFormula, ... )
   }
   if( method == "LA" ) {
      result$coef <- .aidsCoef( coef( est ), nGoods = nGoods, nShifter = nShifter,
         cov = vcov( est ), priceNames = priceNames, shareNames = shareNames,
         shifterNames = shifterNames, df = df.residual( est ) )   # coefficients
      result$wFitted <- aidsCalc( priceNames, totExpName, data = data,
         coef = result$coef, priceIndex = lnp )$shares   # estimated budget shares
      iter <- est$iter
   } else if( method %in% c( "MK", "IL" ) ) {
      b       <- coef( est )# coefficients
      bd      <- b          # difference of coefficients between
                            # this and previous step
      iter    <- est$iter   # iterations of each SUR estimation
      ILiter <- 1          # iterations of IL Loop
      while( ( ( t( bd ) %*% bd ) / ( t( b ) %*% b ) )^0.5 > ILtol &&
            ILiter < ILmaxiter ) {
         ILiter <- ILiter + 1      # iterations of IL Loop
         bl     <- b              # coefficients of previous step
         sysData$lxtr <- log( data[[ totExpName ]] ) -
            aidsPx( "TL", priceNames, shareNames = shareNames, data = data,
            coef = .aidsCoef( coef( est ), nGoods = nGoods, nShifter = nShifter,
               alpha0 = alpha0 ) )
            # real total expenditure using Translog price index
         if( restrict.regMat ) {
            est <- systemfit( system, estMethod, data = sysData, restrict.regMat = restr,
               inst = ivFormula, ... )    # estimate system
         } else {
            est <- systemfit( system, estMethod, data = sysData, restrict.matrix = restr,
               inst = ivFormula, ... )    # estimate system
         }
         iter <- c( iter, est$iter ) # iterations of each estimation
         weightNewCoef <- 1
         b    <- weightNewCoef * coef( est ) + ( 1 - weightNewCoef ) * b # coefficients
         bd   <- b - bl  # difference between coefficients from this
                         # and previous step
      }
      # calculating log of "real" (deflated) total expenditure
      sysData$lxtr <- log( data[[ totExpName ]] ) -
         aidsPx( "TL", priceNames, data = data,
         coef = .aidsCoef( coef( est ), nGoods = nGoods, nShifter = nShifter,
            alpha0 = alpha0 ) )
      # calculating matrix G
      Gmat <- cbind( rep( 1, nObs ), sysData$lxtr )
      for( i in 1:( nGoods ) ) {
         Gmat <- cbind( Gmat, sysData[[ paste( "lp", i, sep = "" ) ]] )
      }
      if( nShifter > 0 ) {
         for( i in 1:nShifter ) {
            Gmat <- cbind( Gmat, sysData[[ paste( "s", i, sep = "" ) ]] )
         }
      }
      # testing matrix G
      if( FALSE ) {
         for( i in 1:( nGoods - 1 ) ) {
            print( fitted( est$eq[[ i ]] ) - Gmat %*%
               coef( est )[ ( ( i - 1 ) * ( nGoods + 2 ) + 1 ):( i * (nGoods + 2 ) ) ] )
         }
      }
      # calculating matrix J
      jacobian <- .aidsJacobian( coef( est ), priceNames, totExpName, data = data,
         shifterNames = shifterNames, alpha0 = alpha0 )
      if( hom ) {
         modRegMat <- .aidsRestr( nGoods = nGoods, nShifter = nShifter,
            hom = hom, sym = sym, restrict.regMat = TRUE )
      } else {
         modRegMat <- diag( ( nGoods - 1 ) * ( nGoods + 2 + nShifter ) )
      }
      # Jmat <- t( modRegMat ) %*% ( diag( nGoods - 1 ) %x% t( Gmat ) ) %*% jacobian
      # JmatInv <- modRegMat %*% solve( Jmat ) %*% t( modRegMat )
      # bcov <- JmatInv  %*% ( est$residCov %x% ( t( Gmat ) %*% Gmat ) ) %*%
      #    t( JmatInv )
      Jmat <- crossprod( modRegMat, ( diag( nGoods - 1 ) %x% t( Gmat ) ) ) %*% jacobian
      JmatInv <- modRegMat %*% solve( Jmat, t( modRegMat ) )
      bcov <- JmatInv  %*% ( est$residCov %x% crossprod( Gmat ) ) %*%
         t( JmatInv )
      result$coef <- .aidsCoef( coef( est ), nGoods = nGoods, nShifter = nShifter,
         cov = bcov, priceNames = priceNames, shareNames = shareNames,
         shifterNames = shifterNames, df = df.residual( est ) )  # coefficients
      result$coef$alpha0 <- alpha0
      result$wFitted <- aidsCalc( priceNames, totExpName, data = data,
         coef = result$coef, priceIndex = "TL" )$shares
         # estimated budget shares
      result$ILiter <- ILiter
   }
   names( result$wFitted ) <- paste( "wFitted", as.character( 1:nGoods ),
      sep = "" )
   result$wResid <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # residuals of shares
   names( result$wResid ) <- paste( "wResid", as.character( 1:nGoods ), sep = "" )
   result$qObs <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qObs ) <- paste( "qObs", as.character( 1:nGoods ), sep = "" )
   result$qFitted <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qFitted ) <- paste( "qFitted", as.character( 1:nGoods ),
      sep = "" )
   result$qResid <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qResid ) <- paste( "qResid", as.character( 1:nGoods ), sep = "" )
   for( i in 1:nGoods ) {
      result$wResid[ , i ] <- data[[ shareNames[ i ] ]] - result$wFitted[ , i ]
      result$qObs[ , i ]   <- data[[ shareNames[ i ] ]] * data[[ totExpName ]] /
         data[[ priceNames[ i ] ]]
      result$qFitted[ , i ] <- result$wFitted[ i ] * data[[ totExpName ]] /
         data[[ priceNames[ i ] ]]
      result$qResid[ , i ] <- result$qObs[ , i ] - result$qFitted[ , i ]
   }
   result$r2 <- numeric( nGoods )
   for( i in 1:( nGoods - 1 ) ) {
      result$r2[ i ] <- summary( est$eq[[ i ]] )$r.squared
   }
   result$r2[ nGoods ] <- rSquared( data[ sample, shareNames[ nGoods ] ],
      result$wResid[ sample, nGoods ] )
   names( result$r2 ) <- shareNames
   result$r2q <- numeric( nGoods ) # R2 values for consumed quantities
   for( i in 1:nGoods ) {
      result$r2q[ i ] <- rSquared( result$qObs[ sample , i ],
         result$qResid[ sample, i ] )
   }
   names( result$r2q ) <- paste( "q_", shareNames, sep = "" )
   result$iter <- iter
   result$est <- est
   result$method <- method
   result$priceIndex <- priceIndex
   result$lnp <- lnp
   result$wMeans <- wMeans
   result$pMeans <- pMeans
   result$xtMean <- xtMean
   result$shareNames <- shareNames
   result$priceNames <- priceNames
   result$totExpName <- totExpName
   result$basePrices <- attributes( result$lnp )$basePrices
   result$baseShares <- attributes( result$lnp )$baseShares

   class( result ) <- "aidsEst"
   return( result )
}
