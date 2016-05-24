# plot the average RMSD of many runs with spread
rmsd_average <- function( rmsdAvInput,
                          skip = 0,
                          printMeans = FALSE,
                          ... )
{
  
  # initialization
  names <- rep( NA, length( rmsdAvInput ) )
  MAT_result <- NULL
  #########
  
  # calculate the average RMSD and the standard deviation
  for( i in 1:length( rmsdAvInput ) )
  {
    names[ i ] <- rmsdAvInput[[ i ]][[ "name" ]]
    VEC_files <- rmsdAvInput[[ i ]][[ "files" ]]
    VEC_values <- c()
    for( j in 1:length( VEC_files ) )
    {
      VEC_buffer <- read.table( VEC_files[ j ] )[ , 2 ]
      if( skip > 0 )
        VEC_buffer <- VEC_buffer[ skip:length( VEC_buffer ) ]
      VEC_values <- c( VEC_values, VEC_buffer )
    }
    MAT_result <- rbind( MAT_result, c( mean( VEC_values ), sd( VEC_values ) ) )
  }
  #########
  
  # set the column and row names
  colnames( MAT_result ) <- c( "values", "sds" )
  rownames( MAT_result ) <- names
  #########
  
  # plot the bars and the errors
  par( mar = c( 2.5, 4.25, 1.0, 1.5 ) )
  PLOT_positions = barplot( MAT_result[ , 1 ],
                            ylim = c( 0.0, 1.5 * max( MAT_result[ , 1 ] ) ),
                            ylab = "RMSD [nm]",
                            xaxt = "n",
                            ... )
  plot_segments( cbind( PLOT_positions, MAT_result[ , 1 ] ),
                 VEC_spread = MAT_result[ , 2 ] )
  axis( 1,
        at = PLOT_positions,
        labels = rownames( MAT_result ),
        cex.axis = 0.9,
        tick = FALSE )
  #########
  
  # print the means in case specified
  print( paste( "Mean: ",
                round( mean( MAT_result[ , 1 ] ), digits = 2 ),
                " +/- ",
                round( mean( MAT_result[ , 2 ] ), digits = 2 ),
                sep = "" ) )
  #########
}

# load RMSD
load_rmsd <- function( files,
                       mdEngine = "GROMOS" )
{
  mdEngine <- toupper( mdEngine )
  if( mdEngine != "GROMOS" &&
      mdEngine != "GROMACS" )
    stop( paste( "The specified 'mdEngine', set to ", mdEngine, " is unknown.", sep = "" ) )
  LIST_return <- list()
  for( i in 1:length( files ) )
  {
    TABLE_input <- NA
    if( mdEngine == "GROMOS" )
    {
      TABLE_input <- read.table( files[ i ] )
    }
    if( mdEngine == "GROMACS" )
    {
      inputData <- readLines( files[ i ],
                              warn = FALSE )
      inputData <- gsub( "^[#].*", "", inputData )
      inputData <- gsub( "^[@].*", "", inputData )
      inputData <- gsub( "^\\s+|\\s+$", "", inputData )
      inputData <- inputData[ inputData != "" ]
      VEC_input <- as.numeric( unlist( strsplit( inputData, "\\s+" ) ) )
      TABLE_input <- as.table( matrix( VEC_input, byrow = TRUE, ncol = 2 ) )
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

# plot RMSD
rmsd <- function( rmsdData,
                  printLegend = TRUE,
                  snapshotsPerTimeInt = 1000,
                  timeUnit = "ns",
                  rmsdUnit = "nm",
                  colours = NA,
                  names = NA,
                  legendPosition = "bottomright",
                  barePlot = FALSE,
                  ... )
{
  # get boundaries
  REAL_max_RMSD = max( unlist( lapply( rmsdData[ c( F, T ) ],
                                       FUN = function( x ) max( x ) ) ) )
  INT_max_snapshot = max( unlist( lapply( rmsdData[ c( T, F ) ],
                                           FUN = function( x ) max( x ) ) ) )
  #########
  
  # set colours and names
  PALETTE_RMSD_colours <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
  if( is.na( colours ) )
    colours <- PALETTE_RMSD_colours( length( rmsdData ) / 2 )
  if( all( is.na( names ) ) )
    names = 1:( length( rmsdData ) / 2 )
  #########
  
  LIST_return <- list()
  # plot
  for( i in 1:length( rmsdData ) )
  {
    if( i %% 2 == 1 )
    {
      if( i == 1 )
        plot( rmsdData[[ i ]], rmsdData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n",
              yaxt = ifelse( barePlot, "n", "s" ),
              xlab = "", ylab = "",
              ylim = c( 0, REAL_max_RMSD * 1.05 ), xlim = c( 0, INT_max_snapshot ), ... )
      else
        plot( rmsdData[[ i ]], rmsdData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n", yaxt = "n", xlab = "", ylab = "",
              ylim = c( 0, REAL_max_RMSD * 1.05 ), xlim = c( 0, INT_max_snapshot ) )
      LIST_return[[ length( LIST_return ) + 1 ]] <- list( minValue = min( rmsdData[[ ( i + 1 ) ]] ),
                                                          maxValue = max( rmsdData[[ ( i + 1 ) ]] ),
                                                          meanValue = mean( rmsdData[[ ( i + 1 ) ]] ),
                                                          sd = sd( rmsdData[[ ( i + 1 ) ]] ) )
      par( new = TRUE )
    }
  }
  #########
  
  # plot the rest
  if( !barePlot )
  {
    axis( 1,
          at = split_equidistant( c( 1, length( rmsdData[[ 1 ]] ) ), 7 ),
          labels = split_equidistant( c( 1, ( length( rmsdData[[ 1 ]] ) / snapshotsPerTimeInt ) ), 7 ),
          cex.axis = 1 )
    mtext( side = 1, text = paste( "time [", timeUnit, "]", sep = "" ), line = 3,
           cex = 1 )
    mtext( side = 2, text = paste( "RMSD [", rmsdUnit, "]", sep = "" ), line = 2.75,
           cex = 1 )
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