mvProbitExp <- function( formula, coef, sigma = NULL, data,
   cond = FALSE, algorithm = "GHK", nGHK = 1000, random.seed = 123, ... ) {

   # checking argument 'formula'
   if( is.list( formula ) ) {
      stop( "using different regressors for the dependent variables",
         " has not been implemented yet. Sorry!" )
   } else if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a data frame" )
   }

   # preparing model matrix
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )

   # preparing model response
   yMat <- model.response( mf )
   if( !is.null( yMat ) ) {
      if( !is.matrix( yMat ) ) {
         stop( "either zero or at least two dependent variables",
            " must be specified in argument 'formula'",
            " (e.g. by 'cbind( y1, y2 ) ~ ...')" )
      } else if( !all( yMat %in% c( 0, 1,TRUE, FALSE ) ) ) {
         stop( "all dependent variables must be either 0, 1, TRUE, or FALSE" )
      }
   }

   # calculate expectations
   result <- mvProbitExpInternal( yMat = yMat, xMat = xMat, 
      coef = coef, sigma = sigma, cond = cond, 
      algorithm = algorithm, nGHK = nGHK, random.seed = random.seed, ... )

   return( result )
}
