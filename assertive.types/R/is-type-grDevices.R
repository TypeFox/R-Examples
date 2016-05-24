#' Is the input a raster?
#'
#' Checks to see if the input is a raster.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_raster} wraps \code{is.raster}, providing more 
#' information on failure. \code{is_a_raster} returns \code{TRUE} if the 
#' input is raster and scalar.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[grDevices]{is.raster}}.
#' @examples
#' m <- matrix(hcl(0, 80, seq(50, 80, 10)), nrow=4, ncol=5)
#' assert_is_raster(as.raster(m))
#' \dontrun{
#' #These examples should fail.
#' assert_is_raster(m)
#' }
#' @importFrom assertive.base is2
#' @export
is_raster <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "raster", .xname)
}

