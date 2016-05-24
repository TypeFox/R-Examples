

#' Conversion of vector/array/matrix to row-matrix form
#'
#' This function manipulates a numeric data structure of vector/array/matrix
#' type in order to obtain a matrix representation.
#' For 1D data structures and column/row arrays and matrices the output is
#' turned in a matrix format with just one row.
#' If the input structure is rectangualar, instead, it is only converted in
#' matrix format.
#'
#' @section Warning:
#' The function is \bold{not} supposed to work with arbitrary N-dimensional
#' arrays.
#'
#' @param D a generic array, matrix or vector to be converted in
#' row-matrix format.
#'
#' @examples
#'
#' toRowMatrixForm( 1 : 10 )
#'
#' toRowMatrixForm( array( 1 : 10, dim = c(1,10 ) ) )
#'
#' toRowMatrixForm( array( 1 : 10, dim = c( 10, 1 ) ) )
#'
#' toRowMatrixForm( matrix( 1 : 10, ncol = 10, nrow = 1 ) )
#'
#' toRowMatrixForm( matrix( 1 : 10, ncol = 1, nrow = 10 ) )
#'
#' toRowMatrixForm( matrix( 1 : 12, ncol = 3, nrow = 4 ) )
#'
#' toRowMatrixForm( matrix( 1 : 12, ncol = 4, nrow = 3 ) )
#'
#' @export
toRowMatrixForm = function( D )
{
  if( is.null( dim( D ) ) |
      is.array( D ) &
      length( dim( D ) ) == 1 )
  {
    D = t( as.matrix( D ) )

  } else if( is.matrix( D ) ) {

    if( ncol( D ) == 1 ){
      D = t( D )
    }
  } else {
    stop( 'Error: unsupported value provided to toRowMatrixForm' )
  }

  return( D )
}

#' Function to setup alpha value for a set of colours
#'
#' \code{set_alpha} manipulates a vector of colour representations in order
#' to setup the alpha value, and get the desired transparency level.
#'
#' @param col a vector of colours
#' @param alpha the value(s) of alpha for (each of) the colors.
#'
#' @seealso \code{\link{fDColorPalette}}
#'
#' @examples
#'
#' original_col = c( 'blue', 'red', 'green', 'yellow' )
#'
#' alpha_col = set_alpha( original_col, 0.5 )
#'
#' alpha_col = set_alpha( original_col, c(0.5, 0.5, 0.2, 0.1 ) )
#'
#' dev.new()
#' par( mfrow = c( 1, 2 ) )
#'
#' plot( seq_along( original_col ),
#'       seq_along( original_col ),
#'       col = original_col,
#'       pch = 16,
#'       cex = 2,
#'       main = 'Original colours' )
#'
#' plot( seq_along( alpha_col ),
#'       seq_along( alpha_col ),
#'       col = alpha_col,
#'       pch = 16,
#'       cex = 2,
#'       main = 'Alpha colours' )
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @export
set_alpha = function( col, alpha )
{
  alpha = alpha * 255

  rgb_colors = rbind( col2rgb( col ), alpha = alpha )

  return( apply( rgb_colors, 2, function( x )( rgb( x[ 1 ],
                                                    x[ 2 ],
                                                    x[ 3 ],
                                                    x[ 4 ],
                                                    maxColorValue = 255 ) ) ) )
}

#' A set of fancy color to plot functional datasets
#'
#' This function can be used to generate a palette of colors useful to plot
#' functional datasets with the \code{plot} methods.
#'
#' The function, built around \code{scales::hue_pal}, allows to set up the
#' HCL parameters of the set of colors desired, and besides to set up the
#' alpha channel value.
#'
#' @param N number of different colors (ideally, functional observations).
#' @param hue_range the range of hues in the HCL scheme.
#' @param alpha the alpha channel parameter(s) of the colors (transparency).
#' @param ... additional parameters to be passed to \code{scales::hue_pal}
#'
#' @seealso \code{\link{plot.fData}}, \code{\link{plot.mfData}}
#'
#' @examples
#'
#' N = 1e2
#' angular_grid = seq( 0, 359, length.out = N )
#'
#' dev.new()
#' plot( angular_grid, angular_grid,
#'       col = fDColorPalette( N, hue_range = c( 0, 359 ), alpha = 1 ),
#'       pch = 16, cex = 3 )
#'
#' @export
fDColorPalette = function( N, hue_range = c( 0, 360 ), alpha = 0.8, ... )
{
  return( set_alpha( scales::hue_pal( h = hue_range, ... )( N ),
                     alpha ) )
}
