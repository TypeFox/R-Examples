# load and parse GROMOS hbond timeseries
load_hbond_ts <- function( path,
                           mdEngine = "GROMOS" )
{
  timeseries <- read.table( path )
  return( timeseries )
}

# plot hbond timeseries
# WARNING: "timeRange" does not work
hbond_ts <- function( timeseries,
                      summary,
                      acceptorRange = NA,
                      donorRange = NA,
                      plotOccurences = FALSE,
                      scalingFactorPlot = NA,
                      printNames = FALSE,
                      namesToSingle = FALSE,
                      printAtoms = FALSE,
                      timeUnit = NA,
                      snapshotsPerTimeInt = 1000,
                      timeRange = NA,
                      hbondIndices = NA,
                      barePlot = FALSE,
                      ... )
{
  
  # set all ranges to full, in case they are not specified
  if( all( is.na( acceptorRange ) ) )
    acceptorRange <- c( min( summary[ , 4 ] ), max( summary[ , 4 ] ) )
  if( all( is.na( donorRange ) ) )
    donorRange <- c( min( summary[ , 2 ] ), max( summary[ , 2 ] ) )
  if( all( is.na( hbondIndices ) ) )
    hbondIndices <- list( c( min( summary[ , 1 ] ), max( summary[ , 1 ] ) ) )
  #########

  # select the hbond-IDs of those hbonds, that match all criteria
  VEC_hbondIDs <- c()
  for( i in 1:nrow( summary ) )
  {
    if( summary[ i, 2 ] >= donorRange[ 1 ] &&
        summary[ i, 2 ] <= donorRange[ 2 ] &&
        summary[ i, 4 ] >= acceptorRange[ 1 ] &&
        summary[ i, 4 ] <= acceptorRange[ 2 ] )
    {
      BOOL_add <- FALSE
      for( j in 1:length( hbondIndices ) )
      {
        VEC_hbondIndCur <- hbondIndices[[ j ]]
        if( summary[ i, 1 ] >= VEC_hbondIndCur[ 1 ] &&
            summary[ i, 1 ] <= VEC_hbondIndCur[ 2 ] )
          BOOL_add <- TRUE
      }
      if( BOOL_add )
        VEC_hbondIDs <- c( VEC_hbondIDs, i )
    }
  }
  #########

  # check, that one or more hbonds match the criteria
  # remove all hbonds from the timeseries table, that do not match the criteria
  if( length( VEC_hbondIDs ) == 0 )
    stop( paste( "The selection of acceptor residues ", acceptorRange[ 1 ], ":",
                 acceptorRange[ 2 ], " with donor residues ", donorRange[ 1 ], ":",
                 donorRange[ 2 ], " does not contain any hbonds.", sep = "" ) )
  timeseries <- timeseries[ ( timeseries[ , 2 ] %in% VEC_hbondIDs ), ,
                            drop = FALSE ]
  #########

  # set time axis
  VEC_timeLimits <- c( min( timeseries[ , 1 ] ),
                       max( timeseries[ , 1 ] ) )
  if( !all( is.na( timeRange ) ) )
    VEC_timeLimits <- timeRange
  VEC_timeTicks <- VEC_timeLimits
  if( !is.na( timeUnit ) )
    VEC_timeTicks <- VEC_timeTicks / snapshotsPerTimeInt
  #########

  # get the vector for the names (left plot)
  # in case, names are selected transform the hbond-ID into three or one letter constructs
  VEC_hbondNames <- split_equidistant( VEC_values = c( min( VEC_hbondIDs ),
                                                       max( VEC_hbondIDs ) ),
                                       n = 5 )
  VEC_hbondNamesPositions <- VEC_hbondNames
  if( printNames )
  {
    VEC_hbondNames <- c()
    VEC_hbondNamesPositions <- VEC_hbondIDs
    for( i in 1:length( VEC_hbondIDs ) )
    {
      LIST_tableLine <- summary[ summary[ , 1 ] == VEC_hbondIDs[ i ], ]
      if( namesToSingle )
      {
        LIST_tableLine[[ "resDonorName" ]] <- translate_aminoacids( LIST_tableLine[[ "resDonorName" ]],
                                                                    switchMode = 1 )
        LIST_tableLine[[ "resAcceptorName" ]] <- translate_aminoacids( LIST_tableLine[[ "resAcceptorName" ]],
                                                                       switchMode = 1 ) 
      }
      VEC_hbondNames <- c( VEC_hbondNames,
                           paste( paste( LIST_tableLine[[ "resDonorName" ]],
                                         LIST_tableLine[[ "resDonor" ]],
                                         ifelse( printAtoms,
                                                 paste( ":",
                                                        LIST_tableLine[[ "atomDonorName" ]],
                                                        sep = "" ),
                                                 "" ),
                                         sep = "" ),
                                  "->",
                                  paste( LIST_tableLine[[ "resAcceptorName" ]],
                                         LIST_tableLine[[ "resAcceptor" ]],
                                         ifelse( printAtoms,
                                                 paste( ":",
                                                        LIST_tableLine[[ "atomAcceptorName" ]],
                                                        sep = "" ),
                                                 "" ),
                                         sep = "" ) ) )
    }
  }
  #########

  # calculate a scaling factor in dependence of the number of hbonds to be plotted
  if( is.na( scalingFactorPlot ) )
    scalingFactorPlot <- 0.75 / log( length( VEC_hbondIDs ) )
  #########

  # in case the occurences are to be plotted, make some space
  # do make sure, that the y axis has enough space for long labels
  INT_additionForAtomNumbers <- ifelse( printAtoms,
                                        3.0,
                                        0.0 )
  if( !barePlot )
    par( mar = c( 4.0,
                  ifelse( !namesToSingle && printNames,
                          7.95 + INT_additionForAtomNumbers,
                          5.95 + INT_additionForAtomNumbers ),
                  3.25,
                  ifelse( plotOccurences,
                          0.0,
                          2.0 ) ) )
  else
    par( mar = c( 4.0,
                  3.25,
                  3.25,
                  ifelse( plotOccurences,
                          0.0,
                          2.0 ) ) )
  if( plotOccurences )
    layout( matrix( 1:2, ncol = 2 ), widths = c( 2, 1 ), heights = c( 1, 1 ) )
  #########
  
  # plot (left graph)
  plot( timeseries,
        xlim = VEC_timeLimits, xaxs = "i", xaxt = "n", xlab = "",
        ylim = c( min( VEC_hbondIDs ), max( VEC_hbondIDs ) ), yaxt = "n", ylab = "",
        type = "n" )
  LIST_ellipsis <- list( ... )
  if( !barePlot )
  {
    mtext( side = 3, line = 1.25, cex = 1.45,
           text = ifelse( is.null( LIST_ellipsis[[ "main" ]] ),
                          "Hbond timeseries",
                          LIST_ellipsis[[ "main" ]] ),
           adj = ifelse( plotOccurences,
                         0.795,
                         0.5 ) )
    mtext( side = 1, line = 2.45, cex = 1.0,
           text = paste( "time",
                         ifelse( !is.na( timeUnit ),
                         paste( "[",
                                timeUnit,
                                "]",
                                sep = "" ),
                         "[snapshots]" ) ) )
  }
  if( !printNames && !barePlot )
    mtext( side = 2, line = 2.45, cex = 1, text = "hbonds [#]" )
  if( !barePlot )
  {
    axis( 1,
          at = split_equidistant( VEC_values = VEC_timeLimits,
                                  n = 6,
                                  BOOL_removeLast = TRUE,
                                  BOOL_roundDown = TRUE ),
          labels = split_equidistant( VEC_values = VEC_timeTicks,
                                      n = 6,
                                      BOOL_removeLast = TRUE,
                                      BOOL_roundDown = TRUE ) )
    axis( 2,
          at = VEC_hbondNamesPositions,
          labels = VEC_hbondNames,
          las = ifelse( printNames,
                        1,
                        0 ),
          cex.axis = ifelse( printNames,
                             scalingFactorPlot * 2.15,
                             1 ) )
  }
  segments( timeseries[ , 1 ],
            timeseries[ , 2 ] - scalingFactorPlot,
            timeseries[ , 1 ],
            timeseries[ , 2 ] + scalingFactorPlot,
            lwd = 0.15 )
  #########

  # in case the occurences are to be plotted, add them on the right hand side
  VEC_occurences <- c()
  for( i in 1:length( VEC_hbondIDs ) )
  {
    INT_numberAppearances <- sum( timeseries[ , 2 ] == VEC_hbondIDs[ i ] )
    
    # calculate percent from the time limits
    # CAUTION: might differ from the overall score out of the summary file
    VEC_occurences <- c( VEC_occurences, ( INT_numberAppearances / 
                                             ( VEC_timeLimits[ 2 ] - VEC_timeLimits[ 1 ] ) * 
                                             100 ) )
  }
  if( plotOccurences )
  {
    par( mar = c( 4.0, 0.0, 3.25, 3.0 ) )
    
    # plot
    plot( VEC_hbondIDs~VEC_occurences,
          ylim = c( min( VEC_hbondIDs ),
                    max( VEC_hbondIDs ) ),
          yaxt = "n", ylab = "",
          xaxs = "i", xlim = c( 0, 100 ), xlab = "",
          xaxt = ifelse( barePlot, "n", "s" ),
          type = "n" )
    if( !barePlot )
      mtext( side = 1, line = 2.45, cex = 1, text = "occurence [%]" )
    segments( rep( 0, length( VEC_hbondIDs ) ),
              VEC_hbondIDs,
              VEC_occurences,
              VEC_hbondIDs,
              lwd = 2.0 )
    segments( VEC_occurences,
              VEC_hbondIDs - scalingFactorPlot * 0.65,
              VEC_occurences,
              VEC_hbondIDs + scalingFactorPlot * 0.65,
              lwd = 2.0 )
  }
  #########
  MAT_return <- matrix( c( VEC_hbondIDs,
                           VEC_occurences ),
                        ncol = 2, byrow = FALSE )
  colnames( MAT_return ) <- c( "hbondID", "occurence [%]" )
  return( MAT_return )
}

# load and parse GROMOS hbond output
load_hbond <- function( path,
                        mdEngine = "GROMOS" )
{
  CON_input <- file( path, open = "r" )
  
  #skip first comments
  readLines( CON_input, n = 23, warn = FALSE )
  LIST_buffer <- list()
  INT_line <- 1
  #########
  
  # read line by line and split by one or multiple spaces
  # second block (three-centered hbonds) is skipped
  while( length( STRING_theLine <- readLines( CON_input, n = 1, warn = FALSE ) ) > 0 )
  {
    VEC_buffer <- unlist( strsplit( STRING_theLine, split = " +" ) )
    if( length( VEC_buffer ) == 0 )
      break
    LIST_buffer[[ INT_line ]] <- VEC_buffer
    INT_line <- INT_line + 1
  }
  close( CON_input )
  #########
  
  # transform list of lists into table using temporary file reads / writes and select
  # the important columns only
  DF_buffer <- data.frame( matrix( unlist( LIST_buffer ), ncol = 21, byrow = TRUE ),
                           stringsAsFactors = FALSE )
  input <- DF_buffer[ , c( -20, -19, -18, -15, -14, -12, -7, -6, -3, -1 ) ]
  STRING_tempFile <- tempfile( "MDplot_Hbonds" )
  on.exit( unlink( STRING_tempFile, recursive = FALSE, force = FALSE ) )
  write.table( input, file = STRING_tempFile )
  input <- read.table( STRING_tempFile )
  #########
  
  # update column names
  colnames( input ) <- c( "hbondID", "resDonor", "resDonorName",
                          "resAcceptor", "resAcceptorName", "atomDonor",
                          "atomDonorName", "atomH", "atomAcceptor",
                          "atomAcceptorName", "percentage" )
  #########
  return( input )
}

# plot the hbond information
hbond <- function( hbonds,
                   plotMethod = "residue-wise",
                   acceptorRange = NA,
                   donorRange = NA,
                   printLegend = TRUE,
                   barePlot = FALSE,
                   ... )
{
  
  # residue-wise: sum all residue-to-residue hbonds up (per-atom contributions)
  if( plotMethod == "residue-wise" )
  {
    MAT_result <- matrix()[ -1, -1, drop = FALSE ]
    VEC_boundariesDonor <- c( min( hbonds[ , 2 ] ), max( hbonds[ , 2 ] ) )
    VEC_boundariesAcceptor <- c( min( hbonds[ , 4 ] ), max( hbonds[ , 4 ] ) )
    
    # identify, whether hbond belongs to already identified hbond and add
    # the percentage there
    for( i in 1:nrow( hbonds ) )
    {
      BOOL_found = FALSE
      if( nrow( MAT_result ) != 0 )
      {
        for( j in 1:nrow( MAT_result ) )
        {
          if( ( MAT_result[ j, 2 ] == hbonds[ i, 2 ] ) &&
              ( MAT_result[ j, 4 ] == hbonds[ i, 4 ] ) )
          {
            MAT_result[ j, 9 ] <- MAT_result[ j, 9 ] + hbonds[ i, 9 ]
            BOOL_found = TRUE
            break
          }
        }
      }
      if( !BOOL_found )
      {
        MAT_result <- rbind( MAT_result, hbonds[ i, ] )
      }
    }
    #########
    
    # remove unused columns and define zoom area if specified
    MAT_result <- MAT_result[ , c( -8, -7, -6, -5, -3, -1 ) ]
    MAT_result <- MAT_result[ order( MAT_result[ , 1 ], MAT_result[ , 2 ] ), ]
    PALETTE_residuewise <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
    if( all( !is.na( acceptorRange ) ) )
    {
      MAT_result <- MAT_result[ MAT_result[ , 2 ] >= acceptorRange[ 1 ], , drop = FALSE ]
      MAT_result <- MAT_result[ MAT_result[ , 2 ] <= acceptorRange[ 2 ], , drop = FALSE ]
    }
    if( all( !is.na( donorRange ) ) )
    {
      MAT_result <- MAT_result[ MAT_result[ , 1 ] >= donorRange[ 1 ], , drop = FALSE ]
      MAT_result <- MAT_result[ MAT_result[ , 1 ] <= donorRange[ 2 ], , drop = FALSE ]
    }
    #########
    
    # normalize values to 0 to 100% and plot graph and legend (in case)
    VEC_normalized <- lapply( as.list( MAT_result[ , 3 ] ), function( x ) x / 100 )
    if( printLegend && !barePlot )
    {
      par( mar = c( 4.25, 4.25, 3.25, 0.25 ) )
      layout( matrix( 1:2, ncol = 2 ), widths = c( 2.5, 0.5 ), heights = c( 1.0, 1.0 ) )
    }
    PALETTE_colors <- colorRampPalette( brewer.pal( 11, 'Spectral' ) )
    PALETTE_colors_rev <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
    plot( MAT_result[ , 1:2 ],
          col = PALETTE_colors_rev( 10 )[ as.numeric( cut( as.numeric( VEC_normalized ), breaks = 10 ) )  ],
          pch = 19,
          cex = 0.75,
          xlab = ifelse( barePlot, "", "Donor [residue number]" ),
          ylab = ifelse( barePlot, "", "Acceptor [residue number]" ),
          xaxt = ifelse( barePlot, "n", "s" ),
          yaxt = ifelse( barePlot, "n", "s" ),
          ... )
    if( printLegend && !barePlot )
    {
      par( mar = c( 4.25, 0.25, 3.25, 1.0 ) )
      legend_image <- as.raster( matrix( PALETTE_colors( 10 ), ncol = 1 ) )
      plot( c( 0, 2 ), c( 0, 1 ), type = 'n', axes = F, xlab = '', ylab = '', main = '[%]' )
      text( x = 1.5, y = seq( 0, 1, l = 5 ),
            labels = seq( 0, 100, l = 5 ) )
      rasterImage( legend_image, 0, 0, 1, 1 )
    }
    #########
  }
  else
  {
    stop( paste( "Error: plot method ", plotMethod, " is unknown. Process aborted!" ) )
  }
  return( hbonds )
}