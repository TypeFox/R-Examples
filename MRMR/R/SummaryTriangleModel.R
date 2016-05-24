#' summaryTriangleModel
#' 
#' @export summaryTriangleModel
#' @include TriangleModel.R
#' 
#' @param objTriangleModel TriangleModel object
#' 
#' @return A vector of intervals
#' 
#' @seealso \code{\link{CreateCumulative}}, \code{\link{CreatePriors}}
#' 
summaryTriangleModel = function(objTriangleModel)
{
  print(objTriangleModel@BP)
  print(objTriangleModel@SW)
  print(objTriangleModel@DW)
  
}

#setMethod("summary", signature = "objTriangleModel", definition = summaryTriangleModel)
#setMethod("summary", definition = summaryTriangleModel)