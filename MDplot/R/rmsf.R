# load the RMSF file(s)
load_rmsf <- function( files,
                       mdEngine = "GROMOS" )
{
  LIST_return <- list()
  for( i in 1:length( files ) )
  {
    TABLE_input <- read.table( files[ i ] )
    if( length( LIST_return ) == 0 )
    {
      LIST_return <- list( TABLE_input[ , 1 ], TABLE_input[ , 2 ] )
    }
    else
    {
      LIST_return[[ length( LIST_return ) + 1 ]] <- TABLE_input[ , 1 ]
      LIST_return[[ length( LIST_return ) + 1 ]] <- TABLE_input[ , 2 ]
    }
  }
  return( LIST_return )
}

# plot RMSF of one or multiple files, combined in one list here
rmsf <- function( rmsfData,
                  printLegend = TRUE,
                  rmsfUnit = "nm",
                  colours = NA,
                  residuewise = FALSE,
                  numberXLabels = 7,
                  names = NA,
                  range = NA,
                  legendPosition = "topright",
                  barePlot = FALSE,
                  ... )
{
  # set colours and names
  PALETTE_RMSF_colours <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
  if( all( is.na( colours ) ) )
    colours <- PALETTE_RMSF_colours( length( rmsfData ) / 2 )
  if( all( is.na( names ) ) )
    names <- 1:( length( rmsfData ) / 2 )
  #########
  
  # get boundaries
  REAL_maxRMSF = max( unlist( lapply( rmsfData[ c( F, T ) ],
                                      FUN = function( x ) max( x ) ) ) )
  INT_minAtomnumber = min( unlist( lapply( rmsfData[ c( T, F ) ],
                                           FUN = function( x ) min( x ) ) ) )
  INT_maxAtomnumber = max( unlist( lapply( rmsfData[ c( T, F ) ],
                                           FUN = function( x ) max( x ) ) ) )
  if( all( is.na( range ) ) )
    range <- c( INT_minAtomnumber - 1,
                INT_maxAtomnumber )
  #########
  
  # plot the RMSF for all elements of the list containing the data
  for( i in 1:length( rmsfData ) )
  {
    if( i %% 2 == 1 )
    {
      if( i == 1 )
        plot( rmsfData[[ i ]], rmsfData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n",
              yaxt = ifelse( barePlot, "n", "s" ),
              xlab = "", ylab = "",
              ylim = c( 0, REAL_maxRMSF * 1.05 ), xlim = range,
              ... )
      else
        plot( rmsfData[[ i ]], rmsfData[[ ( i + 1 ) ]], type = "l",
              col = colours[ ceiling( i / 2 ) ], xaxs = "i", yaxs = "i",
              xaxt = "n", yaxt = "n", xlab = "", ylab = "",
              ylim = c( 0, REAL_maxRMSF * 1.05 ), xlim = range )
      par( new = TRUE )
    }
  }
  #########
  
  # plot axis labels and ticks, which are calculated either atom- or residuewise
  if( !barePlot )
  {
    mtext( side = 2, text = paste( "RMSF [", rmsfUnit, "]", sep = "" ), line = 2.4, cex = 1.25 )
    VEC_atomNumbers <- range
    if( !residuewise )
    {
      mtext( side = 1, text = "atom number", line = 3, cex = 1.25 )
      axis( 1, at = split_equidistant( range, numberXLabels ),
            labels = split_equidistant( range, numberXLabels ) )
    }
    else
    {
      mtext( side = 1, text = "residue number", line = 3, cex = 1.25 )
      axis( 1, at = split_equidistant( range, numberXLabels ),
            labels = as.integer( split_equidistant( range, numberXLabels )
                                 / 3 ) )
    }
  }
  #########
  
  # plot the rest
  if( printLegend && !barePlot )
    legend( legendPosition,
            title = "Legend",
            legend = names,
            col = colours,
            lty = 1.0, lwd = 2.0,
            cex = 1.0 )
  #########
  
  return( rmsfData )
}