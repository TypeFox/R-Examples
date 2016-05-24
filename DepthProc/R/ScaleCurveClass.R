#' @rdname getPlot-methods
#' @aliases getPlot,ScaleCurveList
#' @export
setMethod("getPlot", "ScaleCurveList", function(object)
{
  p = .getPlot(object)
  p = p + ggtitle(object[[1]]@title)
  p = p + ylab("Volume")
  p = p + xlab("Alpha")
  return(p)
})


