#' Generate an annual stack of forest change from GFC product
#'
#' Uses thresholded GFC data as output by \code{\link{threshold_gfc}} to make 
#' an annualized layer stack of forest change. See Details for the class codes 
#' used in the annual raster stack. The \code{\link{animate_annual}} function 
#' can be used to produce an animation of forest change from the generated 
#' layer stack.
#' 
#' The output raster stack uses the following codes to describe forest change 
#' at each pixel:
#' \tabular{lc}{
#'     Nodata               \tab 0 \cr
#'     Forest               \tab 1 \cr
#'     Non-forest           \tab 2 \cr
#'     Forest loss          \tab 3 \cr
#'     Forest gain          \tab 4 \cr
#'     Forest loss and gain \tab 5 \cr
#'     Water                \tab 6 \cr
#' }
#'
#' @seealso  \code{\link{threshold_gfc}}, \code{\link{animate_annual}}
#'
#' @export
#' @import raster
#' @param gfc thresholded extract of GFC product for a given AOI (see 
#' \code{\link{threshold_gfc}})
#' @param data_year which version of the Hansen data was used when
annual_stack <- function(gfc, data_year=2015) {
    names(gfc) <- c('forest2000', 'lossyear', 'gain', 'lossgain', 'datamask')
    out <- raster(gfc)
    layer_names <- paste0('y', gen_year_list(data_year))
    for (n in 1:length(layer_names)) {
        if (n == 1) {
            # Code forest as 1, non-forest as 2
            this_year <- gfc$forest2000
            this_year[this_year == 0] <- 2 # non-forest
        } else {
            this_year <- raster(out, layer=(n-1))
            # Code forest loss (can't have loss in the first year, 2000)
            this_year[(gfc$lossyear == n) & !(gfc$gain)] <- 3 # loss
        }
        # Code gain (no years are attributed to gain). Gain can only occur 
        # where loss does not also occur (loss == 0), as gain and loss is coded 
        # separately below.
        this_year[gfc$gain & (gfc$lossyear == 0)] <- 4 # gain
        this_year[gfc$lossgain] <- 5 # loss and gain
        this_year[gfc$datamask == 2] <- 6 # water
        names(this_year) <- layer_names[n]
        out <- addLayer(out, this_year)
    }
    out[gfc$datamask == 0] <- 0 # missing
    out <- setZ(out, as.Date(paste0(gen_year_list(data_year), '-1-1')))
    names(out) <- layer_names
    return(out)
}
