front41ReadOutput <- function( file = "front41.out" ) {

   result <- list()

   output <- readLines( file )
   line <- 0

   lineSearch <- function( start, regexp, stop = length( output ) ) {
      line2 <- start - 1
      repeat {
         line2 <- line2 + 1
         if( length( grep( regexp, output[ line2 ] ) ) > 0 ) {
            break
         }
         if( line2 >= length( output ) ) {
            stop( "did not find '", regexp, "'" )
         }
         if( line2 >= stop ) {
            line2 <- NULL
            break
         }
      }
      return( line2 )
   }

   rmParts <- function( regexp, string ) {
      for( i in 1:length( regexp ) ) {
         string <- gsub( regexp[[ i ]], "", string )
      }
      return( string )
   }

   getValues <- function( start, name, startNo, stop = length( output ) ) {
      if( is.null( name ) && !is.null( startNo ) ) {
         stop( "if 'name' is NULL, 'startNo' has also to be NULL" )
      }
      if( !is.null( name ) ) {
         line2 <- lineSearch( start,
            paste( "^ *", name, " ", startNo, sep = "" ), stop )
      } else {
         line2 <- start
      }
      values <- NULL
      if( !is.null( line2 ) ) {
         i <- ifelse( is.null( startNo ), 0, startNo )
         repeat {
            strings <- strsplit( output[ line2 + i - ifelse(
               is.null( startNo ), 0, startNo ) ], " +" )[[ 1 ]]
            strings <- strings[ strings != "" ]
            if( !is.null( name ) && regexpr(
               paste( "^", name, "[0-9]+$", sep = "" ), strings[ 1 ] ) > 0 ) {
               strings <- c( name, substr( strings[ 1 ], nchar( name ) + 1,
                  nchar( strings[ 1 ] ) ), strings[ 2:length( strings ) ] )
            }
            if( is.null( name ) && is.null( startNo ) ) {
               savedOptions <- options()
               options( warn = -1 )
               newValues <- as.numeric( strings )
               options( savedOptions )
               if( length( newValues ) > 0 ) {
                  if( is.null( values ) ) {
                     values <- matrix( newValues, nrow = 1 )
                  } else if( length( newValues ) >= ncol( values ) ) {
                     values <- rbind( values, newValues[ 1:ncol( values ) ] )
                  } else {
                     values <- rbind( values, c( newValues,
                        rep( NA, ncol( values ) - length( newValues ) ) ) )
                  }
               } else {
                  break
               }
            } else {
               if( is.null( startNo ) ) {
                  if( strings[ 1 ] == name ) {
                     values <- rbind( values, as.numeric( strings[
                        2:length( strings ) ] ) )
                     rownames( values )[ nrow( values ) ] <- name
                  } else {
                     stop( "did not find '", name )
                  }
                  break
               } else {
                  if( strings[ 1 ] == name && strings[ 2 ] == i ) {
                     values <- rbind( values, as.numeric( strings[
                        3:length( strings ) ] ) )
                     rownames( values )[ nrow( values ) ] <-
                        paste( name, "_", i, sep = "" )
                  } else {
                     break
                  }
               }
            }
            i <- i + 1
         }
      }
      return( values )
   }

   line <- lineSearch( line, "Output from the program FRONTIER" )
   result$version <- rmParts( c( "Output from the program FRONTIER",
      "Version", "[ \\(\\)]" ), output[ line ] )

   line <- lineSearch( line, "instruction file" )
   result$insFile <- rmParts( c( "instruction file", "[ =]" ),
      output[ line ] )

   line <- lineSearch( line, "data file" )
   result$dtaFile <- rmParts( c( "data file", "[ =]" ),
      output[ line ] )

   line <- lineSearch( line, "Frontier .see B.C" )
   result$modelType <- 0
   result$modelTypeName <- rmParts( c( "^ " ), output[ line ] )
   if( length( grep( "Error Components Frontier", result$modelTypeName ) ) > 0 ) {
      result$modelType <- 1
   }
   if( length( grep( "Tech. Eff. Effects Frontier", result$modelTypeName ) ) > 0 ) {
      result$modelType <- 2
   }

   line <- lineSearch( line, "The model is a" )
   result$functionType <- 0
   result$functionTypeName <- rmParts( c( "The model is a ", "^ " ),
      output[ line ] )
   if( length( grep( "production function", result$functionTypeName ) ) > 0 ) {
      result$functionType <- 1
   }
   if( length( grep( "cost function", result$functionTypeName ) ) > 0 ) {
      result$functionType <- 2
   }

   line <- lineSearch( line, "The dependent variable is" )
   if( length( grep( "not logged", output[ line ] ) ) > 0 ) {
      result$logDepVar <- FALSE
   } else {
      result$logDepVar <- TRUE
   }

   line <- lineSearch( line, "the ols estimates are" )
   line <- lineSearch( line, "coefficient[ ]*standard-error[ ]*t-ratio" )
   result$olsResults <- getValues( line, "beta", 0 )
   result$olsResults <- rbind( result$olsResults, rep( NA, ncol(
      result$olsResults ) ) )
   result$olsResults[ nrow( result$olsResults ), 1 ] <- getValues(
      line, "sigma-squared", NULL )
   rownames( result$olsResults )[ nrow( result$olsResults ) ] <- "sigma-squared"
   colnames( result$olsResults ) <- c( "coef", "std.err", "t-ratio" )
   result$nXvars <- nrow( result$olsResults ) - 1

   line <- lineSearch( line, "log likelihood function" )
   result$olsLogl <- as.numeric( rmParts( c( "log likelihood function",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "the estimates after the grid search were" )
   result$gridResults <- getValues( line, "beta", 0 )
   result$gridResults <- rbind( result$gridResults,
      getValues( line, "delta", 0, line + result$nXvars + 10 ) )
   result$gridResults <- rbind( result$gridResults,
      getValues( line, "sigma-squared", NULL ) )
   result$gridResults <- rbind( result$gridResults,
      getValues( line, "gamma", NULL ) )
   if( !is.null( lineSearch( line, "^ *mu *[0-9]", line + result$nXvars + 10 ) ) ) {
      result$gridResults <- rbind( result$gridResults,
         getValues( line, "mu", NULL ) )
   }
   if( !is.null( lineSearch( line, "^ *eta *[0-9]", line + result$nXvars + 10 ) ) ) {
      result$gridResults <- rbind( result$gridResults,
         getValues( line, "eta", NULL ) )
   }
   result$gridResults <- cbind( result$gridResults,
      matrix( NA, nrow = nrow( result$gridResults ), ncol = 2 ) )
   colnames( result$gridResults ) <- c( "coef", "std.err", "t-ratio" )

   line <- lineSearch( line, "the final mle estimates are" )
   result$mleResults <- getValues( line, "beta", 0 )
   if( is.null( lineSearch( line, "delta 0", line + result$nXvars + 10 ) ) ) {
      result$mleResults <- rbind( result$mleResults,
         getValues( line, "delta", 1, line + result$nXvars + 10 ) )
   } else {
      result$mleResults <- rbind( result$mleResults,
         getValues( line, "delta", 0, line + result$nXvars + 10 ) )
   }
   result$mleResults <- rbind( result$mleResults,
      getValues( line, "sigma-squared", NULL ) )
   result$mleResults <- rbind( result$mleResults,
      getValues( line, "gamma", NULL ) )
   if( !is.null( lineSearch( line, "^ *mu *[0-9]", line + result$nXvars + 10 ) ) ) {
      result$mleResults <- rbind( result$mleResults,
         getValues( line, "mu", NULL ) )
   }
   if( !is.null( lineSearch( line, "^ *eta *[0-9]", line + result$nXvars + 10 ) ) ) {
      result$mleResults <- rbind( result$mleResults,
         getValues( line, "eta", NULL ) )
   }
   colnames( result$mleResults ) <- c( "coef", "std.err", "t-ratio" )

   line <- lineSearch( line, "log likelihood function" )
   result$mleLogl <- as.numeric( rmParts( c( "log likelihood function",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "LR test of the one-sided error" )
   result$lrTest <-  as.numeric( rmParts( c( "LR test of the one-sided error",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "with number of restrictions" )
   result$lrTestRestrict <-  as.numeric( rmParts( c( "with number of restrictions",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "number of iterations" )
   result$nIter <-  as.numeric( rmParts( c( "number of iterations",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "\\(maximum number of iterations set at" )
   result$maxIter <-  as.numeric( rmParts( c(
      "\\(maximum number of iterations set at", "[ :\\)]" ), output[ line ] ) )

   line <- lineSearch( line, "number of cross-sections" )
   result$nCross <-  as.numeric( rmParts( c( "number of cross-sections",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "number of time periods" )
   result$nPeriods <-  as.numeric( rmParts( c( "number of time periods",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "total number of observations" )
   result$nObs <-  as.numeric( rmParts( c( "total number of observations",
      "[ =]" ), output[ line ] ) )

   line <- lineSearch( line, "thus there are.*obsns not in the panel" )
   result$nObsMissing <-  as.numeric( rmParts( c( "thus there are",
      "obsns not in the panel", "[ :]" ), output[ line ] ) )

   line <- lineSearch( line, "covariance matrix :" ) + 2
   nCoef <- nrow( result$mleResults )
   result$mleCov <- matrix( NA, nCoef, nCoef )
   myRow <- 1
   myCol <- 1
   while( myRow <= nCoef ) {
      strings <- strsplit( output[ line ], "[ ]+" )[[ 1 ]]
      strings <- strings[ strings != "" ]
      values <- as.numeric( strings )
      nValues <- length( values )
      result$mleCov[ myRow, myCol:( myCol + nValues - 1 ) ] <- values
      if( myCol + nValues > nCoef ) {
         myRow <- myRow + 1
         myCol <- 1
      } else {
         myCol <- myCol + nValues
      }
      line <- line + 1
   }
   rownames( result$mleCov ) <- rownames( result$mleResults )
   colnames( result$mleCov ) <- rownames( result$mleResults )
   line <- lineSearch( line, "efficiency estimates" )
   line <- lineSearch( line, "firm *[year]* *eff\\.-est\\." )
   result$efficiency <- as.data.frame( getValues( line + 2, NULL, NULL ) )
   if( ncol( result$efficiency ) == 3 ) {
      names( result$efficiency ) <- c( "firm", "year", "eff.-est." )
   } else if( ncol( result$efficiency ) == 2 ) {
      names( result$efficiency ) <- c( "firm", "eff.-est." )
   }

   line <- lineSearch( line, "mean efficiency" )
   result$meanEfficiency <-  as.numeric( rmParts( c( "mean efficiency",
      "[ =]" ), output[ line ] ) )

   class( result ) <- "front41Output"
   return( result )
}
