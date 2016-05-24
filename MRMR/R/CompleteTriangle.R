#' @title CompleteTriangle
#' 
#' @description
#' This function will bind the projected values to the base triangle data for a "complete" triangle. This facilitates comparison of ultimates between multiple TriangleModels.
#' 
#' @export CompleteTriangle
#' @include TriangleProjection.R
#' @param objProjection A TriangleProjection object
#' @return A data frame with the sample data (the "upper triangle") bound with the projected data (the "lower triangle").
#' 
CompleteTriangle = function(objProjection)
{
  objTriangle = objProjection@TriangleModel@Triangle
  dfBase = objTriangle@TriangleData
  
  dfProjection = objProjection@ProjectionData
  
  commonCols = intersect(colnames(dfBase), colnames(dfProjection))
  
  dfAll = rbind(dfBase[,commonCols], dfProjection[,commonCols]) 
  
  dfAll
}
