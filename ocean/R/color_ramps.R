#' Generate a sequence of colors for plotting bathymetric data.
#' 
#' \code{bathy.colors(n)} generates a sequence of \eqn{n} colors along a
#' linear scale from light grey to pure blue. 
#'
#' @param n The number of colors to return.
#' @param alpha Alpha values to be passed to \code{rgb()}.
#' @return A vector of blue scale colors.
#'
#' @examples {
#' # Plot a colorbar using bathy.colors
#' image(matrix(seq(100), 100), col=bathy.colors(100))
#' }
bathy.colors <- function(n, alpha=1)
    return(rgb(seq(0.9,0,len=n), seq(0.9,0,len=n), 1, alpha))

## TODO: Use hsv colors or rgb to support alpha values
#' Generate a sequence of colors alog the jet colormap.
#'
#' \code{jet.colors(n)} generates a sequence of \eqn{n} colors from dark blue
#' to cyan to yellow to dark red. It is similar to the default color schemes
#' in Python's matplotlib or MATLAB.
#'
#' @param n The number of colors to return.
#' @param alpha The transparency value of the colors. See \code{?rgb} for
#'              details.
#' @return A vecotr of colors along the jet colorramp.
#'
#' @examples {
#' # Plot a colorbar with jet.colors
#' image(matrix(seq(100), 100), col=jet.colors(100))
#' }
jet.colors <- function(n, alpha=1) {
    if(n > 0) {
        if(length(alpha) != 1 & length(alpha) != n) {
            print('Warning: using only first alpha value')
            alpha <- alpha[1]
        }
        if(length(alpha) == 1) {
            alpha <- rep(alpha, n)
        }
        ## TODO Include alpha values
        return(colorRampPalette(c('#000066', 'blue', 'cyan', 'yellow',
                                  'red', '#660000'))(n))
    } else {
        ## Return an empty character string if they requested nothing.
        character()
    }
}

