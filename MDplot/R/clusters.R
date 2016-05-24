# load function for "MDplot_clusters_timeseries"
load_clusters_ts <- function( path,
                              lengths,
                              names = NA,
                              mdEngine = "GROMOS" )
{
  
  # read input and get rid of column 2, which is unnecessary and pack everything in a list
  # of lists, where the first element is the name of the trajectory and the second is the
  # list of clusterIDs between the boundaries specified by the lengths of every trajectory
  TABLE_buf <- read.table( path )[ , -2 ]
  LIST_return <- list()
  INT_startLine <- 1  
  if( all( is.na( names ) ) ||
      length( lengths ) != length( names ) )
    names <- sapply( seq( 1, length( lengths ) ),
                         function( x ) paste( "trajectory",
                                              x,
                                              sep = "" ) )
  for( i in 1:length( lengths ) )
  {
    LIST_return[[ length( LIST_return ) + 1 ]] <- list( names[ i ], TABLE_buf[ INT_startLine:( INT_startLine + ( lengths[ i ] - 1 ) ),
                                                                               2 ] )
    INT_startLine <- INT_startLine + lengths[ i ]
  }
  #########
  return( LIST_return )
}

# plot timeseries of the clusters
clusters_ts <- function( clustersDataTS,
                         clustersNumber = NA,
                         selectTraj = NA,
                         selectTime = NA,
                         timeUnit = NA,
                         snapshotsPerTimeInt = 1000,
                         ... )
{
  
  # check all the input specified and select from data if necessary
  # in addition, generate the fake plotting matrix for the timeseries plot to generate the properly spanned area
  if( is.na( clustersNumber ) )
    clustersNumber <- max( unlist( lapply( clustersDataTS, FUN = function( x ) unlist( x[[ 2 ]] ) ) ) )
  if( all( !is.na( selectTraj ) ) )
    clustersDataTS <- clustersDataTS[ selectTraj ] 
  if( all( !is.na( selectTime ) ) )
    for( i in 1:length( clustersDataTS ) )
      clustersDataTS[[ i ]][[ 2 ]] <- clustersDataTS[[ i ]][[ 2 ]][ selectTime[ 1 ]:selectTime[ 2 ] ]
  INT_maxSnapshots <- max( unlist( lapply( clustersDataTS, FUN = function( x ) length( unlist( x[[ 2 ]] ) ) ) ) )
  MAT_plotSpan <- matrix( nrow = length( clustersDataTS ),
                          ncol = INT_maxSnapshots )
  #########
  
  # do the colors
  PALETTE_colours <- colorRampPalette( brewer.pal( 11, "Spectral" ) )
  COLOURS_clusters <- PALETTE_colours( clustersNumber )
  #########

  # occurences and plot device division indeed
  par( mar = c( 3.0, 7.0, 4.5, 2.0 ) )
  REAL_widthOfOccurences <- ifelse( ( clustersNumber / 10 ) > 1.0,
                                    1.0,
                                    clustersNumber / 10 )
  if( !is.null( list( ... )[[ "main"]] ) )
  {
    layout( matrix( c( 1, 1, 2, 0, 3, 3 ), nrow = 3, byrow = TRUE ),
            widths = c( 1.0,
                        1 - REAL_widthOfOccurences,
                        1.0 ),
            heights = c( 0.4, 1.0, 1.5 ) )
    plot.new()
    mtext( side = 3, padj = 0, cex = 1.45, text = list( ... )[[ "main" ]] )
  }
  else
  {
    layout( matrix( c( 1, 0, 2, 2 ), nrow = 2, byrow = TRUE ),
            widths = c( REAL_widthOfOccurences,
                        1.0 - REAL_widthOfOccurences,
                        1.0 ),
            heights = c( 1.5, 1.5, 1.5 ) )
  }
  #########
  
  # calculate the percentages for the clusters
  VEC_occurences <- c()
  VEC_allClusterIDs <- unlist( lapply( clustersDataTS, function( x ) x[[ 2 ]] ) )
  for( i in 1:clustersNumber )
    VEC_occurences <- c( VEC_occurences,
                         sum( VEC_allClusterIDs == i ) / length( VEC_allClusterIDs ) * 100 )
  MAT_printResults <- matrix( setNumberDigits( VEC_occurences,
                                               2 ),
                              ncol = length( VEC_occurences ) )
  VEC_colNames <- c()
  for( i in 1:length( VEC_occurences ) )
    VEC_colNames <- c( VEC_colNames,
                       paste( "cluster",
                              i,
                              sep = "" ) )
  colnames( MAT_printResults ) <- VEC_colNames
  #########
  
  # plot the top plot: all the occurences in percents
  PLOT_bp <- barplot( VEC_occurences,
                      xaxs = "i", yaxt = "n", xlab = "",
                      yaxs = "i", yaxt = "n", ylab = "populations",
                      bty = "n",
                      col = COLOURS_clusters, cex.lab = 1.45 )
  mtext( sapply( VEC_occurences,
                 FUN = function( x ) paste( round( x, digits = 1 ),
                                            "%" ), 
                 simplify = TRUE ),
         side = 3,
         las = 3,
         at = PLOT_bp,
         line = 0.45 )
  axis( 1,
        labels = 1:length( VEC_occurences ),
        at = PLOT_bp,
        tick = FALSE,
        line = -0.45 )
  #########
  
  # reset margins
  # plot the bottom plot: set the framework for the plotting later on
  par( mar = c( 4.0, 7.0, 0.0, 2.0 ) )
  plot( MAT_plotSpan,
        xlim = c( 1, INT_maxSnapshots ),
        xaxs = "i", xaxt = "n", xlab = "",
        ylim = c( 0.575, length( clustersDataTS ) + 0.425 ),
        yaxt = "n", ylab = "", yaxs = "i",
        bty = "n", type = "n" )
  axis( 1,
        at = split_equidistant( VEC_values = c( 0, INT_maxSnapshots ),
                                n = 5 ),
        labels = split_equidistant( VEC_values = c( 0, INT_maxSnapshots ),
                                    n = 5 ) /
                 ifelse( !is.na( timeUnit ),
                         snapshotsPerTimeInt,
                         1 ),
        tick = FALSE,
        line = -0.45,
        cex.axis = 1.25 )
  axis( 2,
        at = 1:length( clustersDataTS ),
        labels = unlist( lapply( clustersDataTS,
                                function( x ) x[[ 1 ]] ) ),
        tick = FALSE,
        las = 1,
        cex.axis = 1.25 )
  mtext( side = 1, line = 2.25, cex = 1,
         text = paste( "time [",
                       ifelse( is.na( timeUnit ),
                               "snapshots",
                               timeUnit ),
                       "]",
                       sep = "" ) )
  #########
  
  # plot the coloured lines representing the cluster occurences over time here
  if( length( clustersDataTS ) > 1 )
  {
    for( i in 1:length( clustersDataTS ) )
    {
      VEC_trajOccurences <- c()
      VEC_clusterIDs <- unlist( clustersDataTS[[ i ]][[ 2 ]] )
      VEC_xTicks <- seq( 1, length( clustersDataTS[[ i ]][[ 2 ]] ) )
      for( j in 1:clustersNumber )
      {
        segments( VEC_xTicks[ VEC_clusterIDs == j ],
                  rep( i - 0.425, length( VEC_xTicks[ VEC_clusterIDs == j ] ) ),
                  VEC_xTicks[ VEC_clusterIDs == j ],
                  rep( i + 0.425, length( VEC_xTicks[ VEC_clusterIDs == j ] ) ),
                  lwd = 0.65,
                  col = COLOURS_clusters[ j ] )
        VEC_trajOccurences <- c( VEC_trajOccurences,
                                 length( VEC_clusterIDs[ VEC_clusterIDs == j ] ) / length( VEC_clusterIDs ) * 100 )
      }
      MAT_printResults <- rbind( MAT_printResults,
                                 setNumberDigits( VEC_trajOccurences,
                                                  2 ) )
    }
    rownames( MAT_printResults ) <- c( "overall",
                                       unlist( lapply( clustersDataTS, function( x ) x[[ 1 ]] ) ) )
  }
  #########
  
  return( MAT_printResults )
}

# load function for "MDplot_clusters"
load_clusters <- function( path,
                           names = NA,
                           mdEngine = "GROMOS" )
{
  
  # load and transpose matrix
  MAT_pre <- as.matrix( read.table( path ) )[ , -1  ]
  MAT_pre <- MAT_pre[ , ( ( ncol( MAT_pre ) / 2 ) + 1 ):ncol( MAT_pre ) ]
  MAT_pre <- t( MAT_pre )
  if( all( !is.na( names ) ) &&
      length( names ) == nrow( MAT_pre ) )
  {
    rownames( MAT_pre ) <- names
  }
  else
  {
    rownames( MAT_pre ) <- 1:nrow( MAT_pre )
  }
  #########
  return( MAT_pre )
}

# plot the clusters
clusters <- function( clusters,
                      clustersNumber = NA,
                      legendTitle = "trajectories",
                      barePlot = FALSE,
                      ... )
{
  # reduce number of clusters, in case specified and take care of the trajectory names
  if( !is.na( clustersNumber ) )
    clusters <- clusters[ , 1:clustersNumber ]
  colnames( clusters ) <- 1:ncol( clusters )
  #########
  
  # plot clusters
  PALETTE_clusters <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
  COLOURS_CLUSTERS <- PALETTE_clusters( nrow( clusters ) )
  defaultArguments <- list( xlab = ifelse( barePlot, "", "clusters" ),
                            main = "",
                            col = COLOURS_CLUSTERS )
  ellipsis <- list( ... )
  defaultArguments[ names( ellipsis ) ] <- ellipsis
  ellipsis[ names( defaultArguments ) ] <- defaultArguments
  do.call( what = barplot,
           c( list( height = clusters,
                    xaxt = ifelse( barePlot, "n", "s" ),
                    yaxt = ifelse( barePlot, "n", "s" ) ),
              ellipsis ) )
  if( !barePlot )
    legend( "topright", inset = 0.045, legend = rownames( clusters ),
            title = legendTitle, box.lty = 0, box.lwd = 0, 
            col = COLOURS_CLUSTERS, pch = 19, cex = 1.25 )
  #########
  
  VEC_clusterNames <- c()
  for( i in 1:ncol( clusters ) )
    VEC_clusterNames <- c( VEC_clusterNames,
                           paste( "cluster",
                                  i,
                                  sep = "" ) )
  colnames( clusters ) <- VEC_clusterNames
  return( clusters )
}