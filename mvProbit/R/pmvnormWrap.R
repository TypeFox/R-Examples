pmvnormWrap <- function( lower = -Inf, upper = Inf, sigma, algorithm, 
   random.seed, nGHK = NULL, ... ) {

   # checking argument 'sigma'
   if( !is.matrix( sigma ) ) {
      stop( "argument 'sigma' must be a matrix" )
   } else if( nrow( sigma ) != ncol( sigma ) ) {
      stop( "argument 'sigma' must be a square matrix" )
   } else if( !all( is.numeric( sigma ) ) ) {
      stop( "argument 'sigma must be a numeric matrix" )
   }

   # checking argument 'lower'
   if( !all( is.numeric( lower ) ) ) {
      stop( "argument 'lower' must be numeric" )
   } else if( length( lower ) == 1 ) {
      lower <- rep( lower, nrow( sigma ) )
   } else if( length( lower ) != nrow( sigma ) ) {
      stop( "argument 'lower' must either be a single numeric value",
         " or a vector with length equal to the number of rows/columns",
         " of argument 'sigma'" )
   }

   # checking argument 'lower'
   if( !all( is.numeric( upper ) ) ) {
      stop( "argument 'upper' must be numeric" )
   } else if( length( upper ) == 1 ) {
      upper <- rep( upper, nrow( sigma ) )
   } else if( length( upper ) != nrow( sigma ) ) {
      stop( "argument 'upper' must either be a single numeric value",
         " or a vector with length equal to the number of rows/columns",
         " of argument 'sigma'" )
   }

   # check argument 'algorithm'
   algOkay <- TRUE
   ghk <- FALSE
   if( !is.list( algorithm ) && length( algorithm ) != 1 ) {
      stop( "argument 'algorithm' must be a single function",
         " or a single character string" )
   } else if( is.function( algorithm ) ) {
      algResult <- do.call( algorithm, list() )
      if( ! class( algResult )[ 1 ] %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
         algOkay <- FALSE
      }
   } else if( is.character( algorithm ) ) {
      if( tolower( algorithm ) == "ghk" ) {
         ghk <- TRUE
      } else if( ! algorithm %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
         algOkay <- FALSE
      }
   } else if( ! class( algorithm ) %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) { 
      algOkay <- FALSE
   }
   if( !algOkay ) {
      stop( "argument 'algorithm' must be either one of the functions",
         " 'GenzBretz()', 'Miwa()', or 'TVPACK()'",
         " or one of the character strings",
         " \"GenzBretz\", \"Miwa\", or \"TVPACK\"" )
   }


   # checking argument 'random.seed'
   if( length( random.seed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical value" )
   } else if( !is.numeric( random.seed ) ) {
      stop( "argument 'random.seed' must be numerical" )
   }

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by pmvnorm)
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   if( ghk ) {
      if( is.null( nGHK ) ) {
         stop( "if the GHK algorithm is used,",
            " argument 'nGHK' must be specified" )
      } else if( length( nGHK ) != 1 ) {
         stop( "argument 'nGHK' must be a single integer value" )
      } else if( !is.numeric( nGHK ) ) {
         stop( "argument 'nGHK' must be numeric" )
      } else if( nGHK <= 0 ) {
         stop( "argument 'nGHK' must be positive" )
      }
      L <- try( t( chol( sigma ) ), silent = TRUE )
      if( class( L ) == "try-error" ) {
         warning( "the correlation matrix is not positive definite" )
         return( NA )
      }
      trunpt <- rep( NA, length( lower ) )
      above <- rep( NA, length( lower ) )
      for( i in 1:length( lower ) ) {
         if( lower[ i ] == -Inf ) {
            trunpt[ i ] <- upper[ i ]
            above[ i ] <- 1
         } else if( upper[ i ] == Inf ) {
            trunpt[ i ] <- lower[ i ]
            above[ i ] <- 0
         } else {
            stop( "if algorithm 'GHK' is used,",
               " either the lower truncation point must be '-Inf'",
               " or the upper truncation point must be 'Inf'" )
         }
      }
      sink(tempfile())
      on.exit( sink(), add = TRUE )
      result <- ghkvec( L = L, trunpt = trunpt, above = above, r = nGHK )
      result <- drop( result )
   } else {
      result <- pmvnorm( lower = lower, upper = upper, sigma = sigma,
         algorithm = algorithm, ... )
   }

   return( result )
}
