# load timeseries
load_timeseries <- function( files,
                             mdEngine = "GROMOS" )
{
  mdEngine <- toupper( mdEngine )
  if( mdEngine != "GROMOS" &&
      mdEngine != "PLAINTEXT" )
    stop( paste( "The specified 'mdEngine', set to ", mdEngine, " is unknown.", sep = "" ) )
  LIST_return <- list()
  for( i in 1:length( files ) )
  {
    TABLE_input <- NA
    if( mdEngine == "GROMOS" )
    {
      TABLE_input <- read.table( files[ i ] )
    }
    if( mdEngine == "PLAINTEXT" )
    {
      TABLE_input <- read.table( files[ i ] )
    }
    if( length( LIST_return ) == 0 )
      LIST_return <- list( TABLE_input[ , 1 ], TABLE_input[ , 2 ] )
    else
    {
      LIST_return[[ length( LIST_return ) + 1 ]] <- TABLE_input[ , 1 ]
      LIST_return[[ length( LIST_return ) + 1 ]] <- TABLE_input[ , 2 ]
    }
  }
  return( LIST_return )
}

# plot timeseries
timeseries <- function( tsData,
                        printLegend = TRUE,
                        snapshotsPerTimeInt = 1000,
                        timeUnit = "ns",
                        valueName = NA,
                        valueUnit = NA,
                        colours = NA,
                        numberXLabels = 5,
                        names = NA,
                        legendPosition = "bottomright",
                        barePlot = FALSE,
                        ... )
{
  # get boundaries
  REAL_max_value = max( unlist( lapply( tsData[ c( F, T ) ],
                                        FUN = function( x ) max( x ) ) ) )
  REAL_min_value = min( unlist( lapply( tsData[ c( F, T ) ],
                                        FUN = function( x ) min( x ) ) ) )
  VEC_rangeValues <- c()
  if( REAL_max_value < REAL_min_value )
    VEC_rangeValues <- c( REAL_min_value, REAL_max_value )
  else
    VEC_rangeValues <- c( REAL_max_value, REAL_min_value )
  INT_max_snapshot = max( unlist( lapply( tsData[ c( T, F ) ],
                                          FUN = function( x ) max( x ) ) ) )
  INT_min_snapshot = min( unlist( lapply( tsData[ c( T, F ) ],
                                          FUN = function( x ) min( x ) ) ) )
  #########
  
  # set colours and names
  PALETTE_ts_colours <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
  if( is.na( colours ) )
    colours <- PALETTE_ts_colours( length( tsData ) / 2 )
  if( all( is.na( names ) ) )
    names = 1:( length( tsData ) / 2 )
  #########
  
  LIST_return <- list()
  # plot
  for( i in 1:length( tsData ) )
  {
    if( i %% 2 == 1 )
    {
      if( i == 1 )
        plot( tsData[[ i ]], tsData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n",
              yaxt = ifelse( barePlot, "n", "s" ),
              xlab = "", ylab = "",
              ylim = c( VEC_rangeValues[ 1 ] * 0.95, VEC_rangeValues[ 2 ] * 1.05 ),
              xlim = c( INT_min_snapshot, INT_max_snapshot ), ... )
      else
        plot( tsData[[ i ]], tsData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n", yaxt = "n", xlab = "", ylab = "",
              ylim = c( VEC_rangeValues[ 1 ] * 0.95, VEC_rangeValues[ 2 ] * 1.05 ),
              xlim = c( INT_min_snapshot, INT_max_snapshot ) )
      LIST_return[[ length( LIST_return ) + 1 ]] <- list( minValue = min( tsData[[ ( i + 1 ) ]] ),
                                                          maxValue = max( tsData[[ ( i + 1 ) ]] ),
                                                          meanValue = mean( tsData[[ ( i + 1 ) ]] ),
                                                          sd = sd( tsData[[ ( i + 1 ) ]] ) )
      par( new = TRUE )
    }
  }
  #########
  
  # plot the rest
  if( !barePlot )
  {
    axis( 1,
          at = split_equidistant( c( 1, length( tsData[[ 1 ]] ) ), numberXLabels ),
          labels = split_equidistant( c( 1, ( length( tsData[[ 1 ]] ) / snapshotsPerTimeInt ) ), numberXLabels ),
          cex.axis = 1 )
    mtext( side = 1, text = paste( "time [", timeUnit, "]", sep = "" ), line = 3,
           cex = 1 )
    if( !is.na( valueName ) )
    {
      STRING_ylab <- valueName
      if( !is.na( valueUnit ) )
        STRING_ylab <- paste( STRING_ylab,
                              " [",
                              valueUnit,
                              "]",
                              sep = "" )
      mtext( side = 2, text = STRING_ylab, line = 2.75, cex = 1 )
    }
  }
  if( printLegend && !barePlot )
    legend( legendPosition,
            title = "Legend",
            legend = names,
            col = colours,
            lty = 1.0, lwd = 2.0,
            cex = 1 )
  #########
  
  return( LIST_return )
}