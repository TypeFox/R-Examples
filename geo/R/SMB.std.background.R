#' Standard background for the SMB.
#' 
#' Plots standard background for the SMB survey, aka Icelandic Groundfish
#' Survey.
#' 
#' 
#' @param depth Vector of isobaths to draw, default none
#' @param depthcol Colors for isobaths
#' @param depthlty Line types for isobaths
#' @param depthlwd Line widths for isobaths
#' @param eyjar Additional data to plot, e.g. islands
#' @param depthlab Label depths? Defaults to FALSE
#' @param depthlabcsi Depthlabel size
#' @param \dots Additional arguemnts to \code{geoplot} and \code{geolines}
#' @return No value returned, plots standard layout SMB/IGFS plot.
#' @note Depth labelling etc might be cleaned up, use of \code{eyjar} argument
#' could change, calls data object \code{depthloc} with a few positions for
#' isobath labels.
#' @seealso \code{\link{geoplot}}, \code{\link{gbplot}},
#' \code{\link{geolines}}, \code{\link{eyjar}}
#' @keywords hplot
#' @export SMB.std.background
SMB.std.background <-
function(depth, depthcol = 1, depthlty = 1, depthlwd = 1, eyjar, depthlab,
        depthlabcsi = 0.12, ...)
{
        SMB.limits <- list(lat = c(62.85, 67.5), lon = c(-27.8, -9.8))
        geoplot(xlim = SMB.limits, ...)
        if(!missing(depth))
                gbplot(depth, depthcol, depthlty, depthlwd, depthlab,
                        depthlabcsi)
        if(!missing(eyjar))
                geolines(eyjar, ...)
}

