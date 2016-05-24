# load ramachandran data into a matrix, select columns
load_ramachandran <- function( path,
                               angleColumns = c( 1, 2 ),
                               shiftAngles = NA,
                               mdEngine = "GROMOS" )
{
  
  # load and parse matrix, return result
  MAT_buffer <- as.matrix( read.table( path ) )
  if( ncol( MAT_buffer ) < max( angleColumns ) )
    stop( paste( "Error while loading and parsing file '",
                 path, "' since the number of columns is less than ",
                 "the maximum column number specified.",
                 sep = "" ) )
  MAT_input <- MAT_buffer[ , angleColumns ]
  if( !is.na( shiftAngles ) )
  {
    MAT_input[ , 1 ] <- MAT_input[ , 1 ] + shiftAngles
    MAT_input[ , 2 ] <- MAT_input[ , 2 ] + shiftAngles
  }
  return( MAT_input )
}


# plot the angles on a x = ( -180, 180 ) to y = ( -180, 180 ) area
ramachandran <- function( dihedrals,
                          xBins = 150,
                          yBins = 150,
                          heatFun = "norm", 
                          structureAreas = c(),
                          plotType = "sparse",
                          printLegend = FALSE,
                          heatUnits = NA,
                          plotContour = FALSE,
                          barePlot = FALSE,
                          ... )
{
  
  # settings (small offset for label printing required)
  VEC_xTicks  <- c( -135,  -90,  -45,    0,   45,   90,  135 )
  VEC_xLabels <- c( -135,  -90,  -45,    0,   45,   90,  135 )
  PALETTE_sparse <- colorRampPalette( rev( brewer.pal( 11, 'Spectral' ) ) )
  PALETTE_comic <- colorRampPalette( c( "white", rev( brewer.pal( 11, 'Spectral' ) ) ) )
  PALETTE_fancy <- colorRampPalette( c( "lightgrey", rev( brewer.pal( 11, 'Spectral' ) ) ) )
  #########
  
  # determine function for plotting (logarithm is useful for lots of points)
  FUN_heatFun = function( x ) log( x )
  if( heatFun == "log" )
  {}
  else if( heatFun == "norm" )
    FUN_heatFun = function( x ) x
  else
    stop( paste( "Error: the function '", heatFun, "' for parameter 'heatFun' is not defined.", sep = "" ) )
  #########
  
  # plotting
  LIST_filled <- fill_bins( dihedrals,
                            INT_xbins = xBins,
                            INT_ybins = yBins,
                            VEC_xLim = c( -180, 180 ),
                            VEC_yLim = c( -180, 180 ) )
  VEC_heatValues <- c()
  VEC_palette <- NA
  if( printLegend && !barePlot )
  {
    layout( matrix( 1:2, ncol = 2 ), widths = c( 0.75, 0.25 ), heights = c( 1, 1 ) )
    par( mar = c( 4.0, 4.0, 4.0, 0.0 ) )
    VEC_heatValues <- LIST_filled[[ "freq2D" ]]
    VEC_heatValues[ is.na( VEC_heatValues ) ] <- 0
  }
  if( plotType == "sparse" )
  {
    VEC_palette <- PALETTE_sparse( 21 )
    hist2d( dihedrals, nbins = c( xBins, yBins ), same.scale = FALSE, na.rm = TRUE, 
            show = TRUE, col = VEC_palette, xlab = "", ylab = "",
            xaxs = "i", xaxt = "n", yaxs = "i", yaxt = "n",
            ann = TRUE, xaxt = "n", yaxt = "n", FUN = function( x ) FUN_heatFun( length( x ) ),
            xlim = c( -180, 180 ), ylim = c( -180, 180 ), ... )
    if( !barePlot )
    {
      axis( 1, at = VEC_xTicks, labels = VEC_xLabels, cex.axis = 1.0 )
      axis( 2, at = VEC_xTicks, labels = VEC_xLabels, cex.axis = 1.0 )
      mtext( side = 1, text = expression( paste( phi, " [", degree, "]" ) ), line = 2.8, cex = 1.45 )
      mtext( side = 2, text = expression( paste( psi, " [", degree, "]" ) ), line = 2.3, cex = 1.45 )
    }
    
    # contour
    # TODO: separate this functionality somehow
    if( plotContour && ncol( dihedrals ) < 3 )
    {
      DF_frequencies <- as.data.frame( table( findInterval( dihedrals[ , 1 ],
                                                            LIST_filled[[ "xBins" ]] ),
                                              findInterval( dihedrals[ , 2 ],
                                                            LIST_filled[[ "yBins" ]] ) ) )
      DF_frequencies[ , 1 ] <- as.numeric( DF_frequencies[ , 1 ] )
      DF_frequencies[ , 2 ] <- as.numeric( DF_frequencies[ , 2 ] )
      freq2D <- diag( LIST_filled[[ "xBins" ]] ) * 0
      freq2D[ cbind( DF_frequencies[ , 1 ], DF_frequencies[ , 2 ] ) ] <- DF_frequencies[ , 3 ]
      contour( LIST_filled[[ "xBins" ]],
               LIST_filled[[ "yBins" ]],
               freq2D,
               add = TRUE,
               drawlabels = FALSE,
               xaxs = "i", yaxs = "i" )
    }
    #########
  }
  if( plotType == "comic" )
  {
    # heatmap
    VEC_palette <- PALETTE_comic( 100 )
    image( LIST_filled[[ "xBins" ]],
           LIST_filled[[ "yBins" ]],
           FUN_heatFun( LIST_filled[[ "freq2D" ]] ),
           col = VEC_palette,
            xaxt = "n", xlab = "",
            yaxt = "n", ylab = "",
           ... )
    ##########
    
    # contour
    # TODO: separate this functionality somehow
    if( plotContour && ncol( dihedrals ) < 3 )
    {
      DF_frequencies <- as.data.frame( table( findInterval( dihedrals[ , 1 ],
                                                            LIST_filled[[ "xBins" ]] ),
                                              findInterval( dihedrals[ , 2 ],
                                                            LIST_filled[[ "yBins" ]] ) ) )
      DF_frequencies[ , 1 ] <- as.numeric( DF_frequencies[ , 1 ] )
      DF_frequencies[ , 2 ] <- as.numeric( DF_frequencies[ , 2 ] )
      freq2D <- diag( LIST_filled[[ "xBins" ]] ) * 0
      freq2D[ cbind( DF_frequencies[ , 1 ], DF_frequencies[ , 2 ] ) ] <- DF_frequencies[ , 3 ]
      contour( LIST_filled[[ "xBins" ]],
               LIST_filled[[ "yBins" ]],
               freq2D,
               add = TRUE,
               drawlabels = FALSE,
               xaxs = "i", yaxs = "i" )
    }
    #########
    
    if( !barePlot )
    {
      axis( 1, at = VEC_xTicks, labels = VEC_xLabels, cex.axis = 1.0 )
      axis( 2, at = VEC_xTicks, labels = VEC_xLabels, cex.axis = 1.0 )
      mtext( side = 1, text = expression( paste( phi, " [", degree, "]" ) ), line = 2.8, cex = 1.45 )
      mtext( side = 2, text = expression( paste( psi, " [", degree, "]" ) ), line = 2.3, cex = 1.45 )
    }
  }
  if( plotType == "fancy" )
  {
    #DF_frequencies <- as.data.frame( table( findInterval( dihedrals[ , 1 ],
    #                                                      LIST_filled[[ "xBins" ]] ),
    #                                        findInterval( dihedrals[ , 2 ],
    #                                                      LIST_filled[[ "yBins" ]] ) ) )
    #DF_frequencies[ , 1 ] <- as.numeric( DF_frequencies[ , 1 ] )
    #DF_frequencies[ , 2 ] <- as.numeric( DF_frequencies[ , 2 ] )
    #freq2D <- diag( LIST_filled[[ "xBins" ]] ) * 0
    #freq2D[ cbind( DF_frequencies[ , 1 ], DF_frequencies[ , 2 ] ) ] <- DF_frequencies[ , 3 ]
    freq2D <- LIST_filled[[ "freq2D" ]]
    freq2D[ is.na( freq2D ) ] <- 0
    INT_numberColours <- 100
    VEC_palette <- PALETTE_fancy( INT_numberColours )
    zFacetValue <- freq2D[ -1, -1 ] + 
                   freq2D[ -1, -ncol( freq2D ) ] +
                   freq2D[ -nrow( freq2D ), -1 ] + 
                   freq2D[ -nrow( freq2D ), -ncol( freq2D ) ]
    zFacetCol <- cut( zFacetValue, INT_numberColours )
    par( oma = c( 1.0, 0.0, 0.0, 0.0 ), mar = c( 0.0, 0.0, 4.5, 0.0 ) )
    perspMatrix <- persp( freq2D,
                          col = VEC_palette[ zFacetCol ],
                          box = FALSE,
                          axes = FALSE,
                          theta = 15,
                          r = 3,
                          d = 0.75,
                          ... )
    lines( trans3d( seq( 0, 1, by = 0.25 ), 0, 0, perspMatrix ), col = "black" )
    #lines( trans3d( 0, seq( 0, 1, by = 0.25 ), 0, perspMatrix ), col = "black" )
    #lines( trans3d( 0, 0, seq( 0, max( DF_frequencies[ , 3 ] ), length.out = 4 ), perspMatrix ),
    #       col = "black" )
    
    # x-axis
    tick.start <- trans3d( seq( 1 / 8, 7 / 8, by = 1 / 8 ), 0, 0, perspMatrix )
    tick.end <- trans3d( seq( 1 / 8, 7 / 8, by = 1 / 8 ), -0.05, 0, perspMatrix )
    segments( tick.start$x, tick.start$y, tick.end$x, tick.end$y )
    labels <- c( -135, -90, -45, 0, 45, 90, 135 )
    label.pos <- trans3d( seq( 1 / 8, 7 / 8, by = 1 / 8 ), -0.09, 0, perspMatrix )
    if( !barePlot )
    {
      text( label.pos$x, label.pos$y, labels = labels, adj = c( 0, NA ), cex = 0.9 )
      text( label.pos$x[ 4 ], label.pos$y[ 4 ] - 0.0225,
            labels = c( expression( paste( phi, " [", degree, "]" ) ) ),
            cex = 1.25 )
    }
    
    # y-axis
    tick.start <- trans3d( 0, seq( 1 / 8, 7 / 8, by = 1 / 8 ), 0, perspMatrix )
    tick.end <- trans3d( -0.0175, seq( 1 / 8, 7 / 8, by = 1 / 8 ), 0, perspMatrix )
    segments( tick.start$x, tick.start$y, tick.end$x, tick.end$y )
    labels <- c( -135, -90, -45, 0, 45, 90, 135 )
    label.pos <- trans3d( -0.11, seq( 1 / 8, 7 / 8, by = 1 / 8 ), 0, perspMatrix )
    if( !barePlot )
    {
      text( label.pos$x, label.pos$y, labels = labels, adj = c( 0, NA ), cex = 0.9 )
      text( label.pos$x[ 4 ] - 0.025, label.pos$y[ 4 ] + 0.015,
            labels = c( expression( paste( psi, " [", degree, "]" ) ) ),
            cex = 1.25 )
      #label.pos <- trans3d( seq( 0, 1, by = 0.25 ), ( 0 - 0.1 ), 0, perspMatrix )
      #text( label.pos$x, label.pos$y, labels = labels, adj = c( 0, NA ), cex = 1.0 )
    }
  }
  #########
  
  # if specified, print all Ramachandran regions and their labels
  if( length( structureAreas ) > 0 )
    for( i in 1:length( structureAreas ) )
    {
      if( length( structureAreas[[ i ]] ) < 3 )
      {
        next
      }
      for( j in 2:length( structureAreas[[ i ]] ) )
      {
        if( length( structureAreas[[ i ]] ) > j )
        {
          segments( structureAreas[[ i ]][[ j ]][[ 1 ]], structureAreas[[ i ]][[ j ]][[ 2 ]], 
                    structureAreas[[ i ]][[ j + 1 ]][[ 1 ]], structureAreas[[ i ]][[ j + 1 ]][[ 2 ]],
                    col = "black", lwd = 2 )
        }
        else
        {
          segments( structureAreas[[ i ]][[ j ]][[ 1 ]], structureAreas[[ i ]][[ j ]][[ 2 ]], 
                    structureAreas[[ i ]][[ 2 ]][[ 1 ]], structureAreas[[ i ]][[ 2 ]][[ 2 ]],
                    col = "black", lwd = 2 )
        }
      }
      cen <- calculate_mid( structureAreas[[ i ]][ -1 ] )
      text( cen[[ 1 ]], cen[[ 2 ]], labels = structureAreas[[ i ]][[ 1 ]], cex = 2.15 )
    }
  #########
  
  # print legend if specified to do so
  if( printLegend && !barePlot )
  {
    if( ncol( dihedrals ) == 3 )
      VEC_heatValues <- dihedrals[ , 3 ]
    if( plotType != "fancy" )
      VEC_palette <- c( "white",
                        VEC_palette )
    legend_image <- as.raster( matrix( rev( VEC_palette ) ), ncol = 1 )
    par( mar = c( 4.0, 1.0, 4.0, 2.5 ) )
    STRING_legendCaption <- ifelse( is.na( heatUnits ),
                                    "Legend",
                                    paste( "Legend ", heatUnits, sep = "" ) )
    plot( c( 0, 2 ), c( 0, 1 ), type = 'n',
          axes = F, xlab = '', ylab = '',
          main = STRING_legendCaption )
    text( x = 1.5, y = seq( 0, 1, l = 5 ),
          labels = round( seq( min( VEC_heatValues ),
                               max( VEC_heatValues ),
                               l = 5 ), digits = 0 ) )
    rasterImage( legend_image, 0, 0, 1, 1 )
  }
  #########
  
  return( LIST_filled )
}