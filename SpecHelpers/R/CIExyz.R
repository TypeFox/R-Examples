

##' Spectral Locus for the 1931 CIE chromaticity diagram
##'
##' This data set gives wavelengths every 1.0 nm, along
##' with the associated CIE xyz values for the spectral locus of
##' the 1931 CIE chromaticity diagram.  They are called xyz values
##' here as they are called that in the original source, but they
##' are also known as xyY or XYZ values.
##'
##' @seealso \code{\link{plotCIEchrom}} for examples of this data in use.
##' 
##' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
##'
##'
##' @keywords data
##'
##' @source Color Vision Research Lab.  \url{cvrl.ioo.ucl.ac.uk/index.htm}
##' Go to this URL, then choose 'NEW CIE XYZ...' In the new page that opens,
##' go to 'New physiologically-relevant CIE x,y chromaticity coordinates
##' (proposed)' and get the 2-deg coordinates at 1.0 nm resolution
##'
##' @format A data frame with 4400 observations each with the
##' following 4 variables:
##' \describe{
##'   \item{wavelength}{wavelength in nm}
##'   \item{x}{x values}
##'   \item{y}{y values}
##'   \item{z}{z values}
##' }
##'
##' @examples
##' data(CIExyz)
"CIExyz"