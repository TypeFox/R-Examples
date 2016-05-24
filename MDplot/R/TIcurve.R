# load the lambda point information here
load_TIcurve <- function( files,
                          mdEngine = "GROMOS" )
{
  LIST_files <- list()
  for( i in 1:length( files ) )
    LIST_files[[ length( LIST_files ) + 1 ]] <- as.matrix( read.table( files[ i ] ) )
  return( LIST_files )
}

# plot the curve and calculate and plot the integral
TIcurve <- function( lambdas,
                     invertedBackwards = FALSE,
                     energyUnit = "kJ/mol",
                     printValues = TRUE,
                     printErrors = TRUE,
                     errorBarThreshold = 0,
                     barePlot = FALSE,
                     ... )
{
  
  # calculate maximum / minimum and colours / do inversion if necessary
  VEC_errorLimits <- c( min( unlist( lapply( lambdas,
                                             FUN = function( x ) min( x[ , 3 ] ) ) ) ),
                        max( unlist( lapply( lambdas,
                                             FUN = function( x ) max( x[ , 3 ] ) ) ) ) )
  VEC_valueLimits <- c( min( unlist( lapply( lambdas,
                                             FUN = function( x ) min( x[ , 2 ] ) ) ) ) -
                        VEC_errorLimits[ 2 ] * 1.25,
                        max( unlist( lapply( lambdas,
                                             FUN = function( x ) max( x[ , 2 ] ) ) ) ) +
                        VEC_errorLimits[ 2 ] * 1.25 )
  VEC_colours <- c( "black", "red" )
  if( length( lambdas ) > 1 &&
      invertedBackwards )
    lambdas[[ 2 ]][ , 2 ] <- rev( sapply( lambdas[[ 2 ]][ , 2 ],
                                             FUN = function( x ) x * ( -1 ) ) )
  #########
  
  # set proper outer margins and plot it
  if( !barePlot )
    if( length( lambdas ) > 1 && printValues )
      par( oma = c( 3.25, 3.00, 0.45, 0.0 ) )
    else
      par( oma = c( 1.35, 3.00, 0.45, 0.0 ) )
  for( i in 1:length( lambdas ) )
  {
    if( i == 1 )
    {
      TIplot <- plot( lambdas[[ i ]],
                      ylim = VEC_valueLimits,
                      ylab = "",
                      xaxt = ifelse( barePlot, "n", "s" ),
                      yaxt = ifelse( barePlot, "n", "s" ),
                      xaxs = "i", xlab = "",
                      pch = 19, cex = 0.6,
                      col = VEC_colours[ i ],
                      ... )
      if( !barePlot )
      {
        mtext( side = 1, text = expression( lambda ), line = 3, cex = 1.45 )
        mtext( side = 2, 
               text = expression( atop( "<"*frac( paste( partialdiff, "H" ),
                                                  paste( partialdiff, lambda ) )*">", ) ),
               line = 2.4, cex = 1.75, las = 1 )
        mtext( side = 2,
               text = paste( "[",
                             energyUnit,
                             "]",
                             sep = "" ),
               line = 2.4, cex = 1.0, las = 1, padj = 2.25 )
      }
      abline( h = 0, lwd = 1, lty = 3 )
    }
    else
      TIplot <- plot( lambdas[[ i ]],
                      ylim = VEC_valueLimits,
                      ylab = "", yaxt = "n",
                      xaxs = "i", xaxt = "n", xlab = "",
                      pch = 19, cex = 0.6,
                      col = VEC_colours[ i ],
                      ... )
  if( printErrors )
  {
    MAT_positions <- lambdas[[ i ]][ , 1:2 ]
    VEC_errors <- lambdas[[ i ]][ , 3 ]
    for( j in nrow( MAT_positions ):1 )
    {
      if( VEC_errors[ j ] <= errorBarThreshold )
      {
        VEC_errors <- VEC_errors[ -j ]
        MAT_positions <- MAT_positions[ -j, ]
      }
    }
    VEC_curColours <- rep( VEC_colours[ i ], times = nrow( MAT_positions ) )
    plot_segments( MAT_positions,
                   VEC_errors,
                   0.01,
                   col = VEC_curColours )
  }
  par( new = TRUE )
  plot( lambdas[[ i ]],
        ylim = VEC_valueLimits, yaxt = "n", ylab = "",
        xaxs = "i", xaxt = "n", xlab = "", yaxt = "n",
        type = "l", lwd = 1, col = VEC_colours[ i ] )
  par( new = TRUE )
  }
  #########
  
  # integrate over curves
  REAL_forward_integral <- unlist( integrate_curve( lambdas[[ 1 ]] )[ "integral" ] )
  REAL_forward_error <- unlist( integrate_curve( lambdas[[ 1 ]] )[ "error" ] )
  INT_significantForward <- get_sign_digits( REAL_forward_error )
  REAL_forward_integral_rounded <- round( REAL_forward_integral, digits = INT_significantForward )
  REAL_forward_error_rounded <- round( REAL_forward_error, digits = INT_significantForward )
  MAT_integrationResults <- matrix( c( REAL_forward_integral_rounded, REAL_forward_error_rounded ),
                                    ncol = 2 )
  colnames( MAT_integrationResults ) <- c( "deltaG", "error" )
  rownames( MAT_integrationResults ) <- c( "forward" )
  REAL_hysteresis <- NA
  if( !barePlot && printValues )
    mtext( side = 1, line = 4.75, cex = 1.0,
           adj = 1,
           text = substitute( paste( Delta, "G"["forw"], " = ",
                                     REAL_forward_integral_rounded,
                                     #" \u00b1 ",
                                     " +/- ",
                                     REAL_forward_error_rounded,
                                     paste( " [",
                                            energyUnit,
                                            "]",
                                            sep = "" ) ) ) )
  if( length( lambdas ) > 1 && !barePlot && printValues )
  {
    REAL_backward_integral <- unlist( integrate_curve( lambdas[[ 2 ]] )[ "integral" ] )
    REAL_backward_error <- unlist( integrate_curve( lambdas[[ 2 ]] )[ "error" ] )
    INT_significantBackward <- get_sign_digits( REAL_backward_error )
    REAL_backward_integral_rounded <- round( REAL_backward_integral, digits = INT_significantBackward )
    REAL_backward_error_rounded <- round( REAL_backward_error, digits = INT_significantBackward )
    MAT_integrationResults <- rbind( MAT_integrationResults,
                                     c( REAL_backward_integral_rounded, REAL_backward_error_rounded ) )
    rownames( MAT_integrationResults ) <- c( "forward", "backward" )
    mtext( side = 1, line = 6.0, cex = 1.0,
           adj = 1,
           text = substitute( paste( Delta, "G"["back"], " = ",
                                     REAL_backward_integral_rounded,
                                     #" \u00b1 ",
                                     " +/- ",
                                     REAL_backward_error_rounded,
                                     paste( " [",
                                            energyUnit,
                                            "]",
                                            sep = "" ) ) ) )
    REAL_hysteresis <- round( REAL_forward_integral - REAL_backward_integral,
                              digits = min( c( INT_significantForward,
                                               INT_significantBackward ) ) )
    mtext( side = 1, line = 7.25, cex = 1.0,
           adj = 1,
           text = substitute( paste( "hysteresis = ",
                                     REAL_hysteresis,
                                     paste( " [",
                                            energyUnit,
                                            "]",
                                            sep = "" ) ) ) )
  }
  #########
  
  LIST_return <- list( lambdapoints = lambdas,
                       integrationresults = MAT_integrationResults,
                       hysteresis = REAL_hysteresis )
  return( LIST_return )
}