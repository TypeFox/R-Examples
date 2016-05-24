#' Get colors for plots
#' 
#' Get gradient default colors for plots
#' 
#' @param col_vec chr string of plot colors to use. Any color palette from RColorBrewer can be used as a named input. Palettes from grDevices must be supplied as the returned string of colors for each palette.
#' 
#' @details This is a convenience function for retrieving a color palette.  Palettes from RColorBrewer will use the maximum number of colors.  The default palette is 'Spectral'. 
#' 
#' @return A character vector of colors in hexidecimal notation.
gradcols <- function(col_vec = NULL){

  cols <- RColorBrewer::brewer.pal(11, 'Spectral')
    
  # color ramp for pts
  if(!is.null(col_vec)){
 
    # get color palette if provided, otherwise user-supplied
    chk_cols <- row.names(RColorBrewer::brewer.pal.info)
    
    if(any(chk_cols %in% col_vec)){
      
      col_vec <- chk_cols[which(chk_cols %in% col_vec)][1]
      max_cols <- RColorBrewer::brewer.pal.info[col_vec, 'maxcolors']
      cols <- RColorBrewer::brewer.pal(max_cols, col_vec)
      
    } else {
      
      cols <- col_vec
      
    }
   
  }
  
  return(cols)

}