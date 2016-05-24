#' Plots Legend 
#' 
#' Plots the legend for resistivity values.
#' 
#' @param .Object either a single Profile or a ProfileSet.
#' @param legend.lab label of legend (default: expression(paste("Resistivity [", Omega, "]"))).
#' @param minData minimum value.
#' @param maxData maximum value.
#' @param breaks Break points in sorted order to indicate the intervals for assigning the colors. 
#' Note that if there are nlevel colors there should be (nlevel+1) breakpoints. 
#' If breaks is not specified (nlevel+1) equally spaced breaks are created where the first and last bin have their midpoints at the minimum and maximum values in z or at zlim.
#' @param lab.breaks number of breaks.
#' @param legend.line distance in units of character height (as in mtext) of the legend label from the color bar. 
#' Make this larger if the label collides with the color axis labels.
#' @param nlevel number of color levels.
#' @param col vector of colors.
#' @param trafo transformation to be done on data (default: log).
#' For linear scale: function(x) x.
#' @param backtrafo back transformation to plot correct labels (default: exp).
#' For linear scale: function(x) x.
#' @param horizontal If false legend will be a vertical strip on the right side. If true (default) the legend strip will be along the bottom.
#' @param ... image.plot arguments.
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{ProfileSet-class}},
#' \code{\link{plot3dXyz}},
#' @examples 
#' # data(sinkhole)
#' 
#' # plotLegend(sinkhole)
#' 
#' # for linear scale:
#' # plotLegend(sinkhole@profiles[[1]], 
#' #            trafo=function(x) x, 
#' #            backtrafo=function(x) x,
#' #            minData=100, maxData=50000)
setGeneric("plotLegend", function(.Object, 
                                  legend.lab=expression(paste("Resistivity [", Omega, " m]")),
                                  minData=0, maxData=999999,
                                  breaks = NULL, legend.line=2.2,
                                  nlevel=18,
                                  lab.breaks=c(),                                   
                                  horizontal=T,
                                  col=colors, trafo=log, backtrafo=exp, ...) { 
  
  standardGeneric("plotLegend")
  if(length(lab.breaks) > 0)
    nlevel <- length(lab.breaks)
  if(length(lab.breaks) == 0)
    lab.breaks <- round(backtrafo(seq(trafo(minData),
                                      trafo(maxData), 
                                      length.out=nlevel)))
  image.plot(legend.only=TRUE, add=F, breaks=breaks,
             zlim= c(minData, maxData),
             legend.line = legend.line,
             legend.lab = legend.lab,
             nlevel=nlevel, 
             col=colorRampPalette(col)(nlevel-1),
             lab.breaks=lab.breaks, 
             horizontal=horizontal, ...)
})

#' @rdname plotLegend
#' @export
setMethod("plotLegend", signature(.Object="ProfileSet"),
          function(.Object, legend.lab,
                   minData=.Object@minData, maxData=.Object@maxData) {
          })

#' @rdname plotLegend
#' @export
setMethod("plotLegend", signature(.Object="Profile"),
          function(.Object, legend.lab,
                   minData=.Object@xyzData@minData, 
                   maxData=.Object@xyzData@maxData) {
          })