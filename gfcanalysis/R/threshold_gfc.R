#' Threshold the GFC product
#'
#' Uses the GFC data output from \code{\link{extract_gfc}} to make an 
#' thresholded layer stack with five layers: treecover2000, loss, gain,
#' lossyear, and datamask layers. See Details for the coding used in each 
#' layer. Note that the output file format is fixed as GeoTIFF with LZW 
#' compression.
#' 
#' The output uses the following codes to describe forest change at each pixel:
#'
#' \bold{Band 1 (forest2000)}
#' \tabular{lc}{
#'     Non-forest           \tab 0 \cr
#'     Forest               \tab 1 \cr
#' }
#'
#' \bold{Band 2 (lossyear)}
#' \tabular{lc}{
#'     No loss      \tab 0  \cr
#'     Loss in 2001 \tab 1  \cr
#'     Loss in 2002 \tab 2  \cr
#'     Loss in 2003 \tab 3  \cr
#'     Loss in 2004 \tab 4  \cr
#'     Loss in 2005 \tab 5  \cr
#'     Loss in 2006 \tab 6  \cr
#'     Loss in 2007 \tab 7  \cr
#'     Loss in 2008 \tab 8  \cr
#'     Loss in 2009 \tab 9  \cr
#'     Loss in 2010 \tab 10 \cr
#'     Loss in 2011 \tab 11 \cr
#'     Loss in 2012 \tab 12 \cr
#'     Loss in 2013 \tab 13 \cr
#'     Loss in 2014 \tab 14 \cr
#' }
#' Note that lossyear is zero for pixels that were not forested in 2000, and 
#' that the 2013 and 2014 loss layers are not available in the original 2013 
#' Hansen dataset (the 2013 loss layer is available in the 2014 and 2015 
#' updates, while the 2014 loss layer is available in the 2015 update only).
#'
#' \bold{Band 3 (gain)}
#' \tabular{lc}{
#'     No gain \tab 0 \cr
#'     Gain    \tab 1 \cr
#' }
#' Note that gain is zero for pixels that were forested in 2000
#'
#' \bold{Band 4 (lossgain)}
#' \tabular{lc}{
#'     No loss and gain \tab 0 \cr
#'     Loss and gain    \tab 1 \cr
#' }
#' Note that loss and gain is difficult to interpret from the thresholded 
#' product, as the original GFC product does not provide information on the 
#' sequence (loss then gain, or gain then loss), or the levels of canopy cover 
#' reached prior to loss (for gain then loss) or after loss (for loss then gain 
#' pixels). The layer is calculated here as: \code{lossgain <- gain & (lossyear 
#' != 0)}, where loss year and gain are the original GFC gain and lossyear 
#' layers, respectively.
#'
#' \bold{Band 5 (datamask)}
#' \tabular{lc}{
#'     No data \tab 0 \cr
#'     Land    \tab 1 \cr
#'     Water   \tab 2 \cr
#' }
#'
#' @seealso \code{\link{extract_gfc}}
#'
#' @export
#' @import raster
#' @param gfc extract of GFC product for a given AOI (see 
#' \code{\link{extract_gfc}})
#' @param forest_threshold percent woody vegetation to use as a threshold for 
#' mapping forest/non-forest
#' @param ... additional arguments as for writeRaster, such as \code{filename} 
#' or \code{overwrite}
#' @return \code{RasterBrick} with thresholded GFC product (see details above)
threshold_gfc <- function(gfc, forest_threshold=25, ...) {
    names(gfc) <- c('treecover2000', 'loss', 'gain', 'lossyear', 'datamask')
    recode_gfc <- function(treecover2000, loss, gain, lossyear, datamask) {
        forest2000 <- treecover2000 > forest_threshold
        # Note forest2000 is a binary variable
        lossyear_recode <- lossyear * forest2000
        gain_recode <- gain & (!forest2000)
        lossgain <- gain & (lossyear != 0)
        thresholded <- matrix(c(forest2000, lossyear_recode, gain_recode, 
                                lossgain, datamask),
                             nrow=length(forest2000),
                             ncol=5)
        return(thresholded)
    }
    thresholded <- overlay(gfc, fun=recode_gfc, datatype='INT1U', 
                           format='GTiff', options="COMPRESS=LZW", ...)
    names(thresholded) <- c('forest2000', 'lossyear', 'gain', 'lossgain', 
                            'datamask')
    return(thresholded)
}
