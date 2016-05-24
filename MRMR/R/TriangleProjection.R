#' TriangleProjection class
#' 
#' @description
#' TriangleProjection is an S4 class used to project values.
#' 
#' @seealso \code{\link{newTriangle}}
#' 
#' @name TriangleProjection-class
#' @rdname TriangleProjection-class
#' @exportClass TriangleProjection
#' 
setClass("TriangleProjection"
         , representation(TriangleModel = "TriangleModel"
                          , ProjectToDev = "logical"
                          , MaxDev = "numeric"
                          , AsOfDate = "POSIXct"
                          , ProjectionData = "data.frame"))

#' @title TriangleProjection
#' 
#' @description
#' This will construct a TriangleProjection object
#' 
#' @export 
#' 
#' @param objTriangleModel A TriangleModel object
#' @param ProjectToDev Boolean indicating whether one is projecting to a maximum 
#'        development interval. If this paramter is FALSE, there must be an argument for AsOfDate
#' @param MaxDev The maximum development interval to which to project.
#' @param AsOfDate The date to which one wants to project.
#' 
#' @include Triangle.R
#' @include TriangleModel.R
TriangleProjection = function(objTriangleModel
                              , ProjectToDev = TRUE
                              , MaxDev = 10
                              , AsOfDate = NULL)
{
  dfLatest = LatestDiagonal(objTriangleModel@Triangle)
  
  lOriginYear = dlply(dfLatest, "OriginPeriodStart")
  
  if (ProjectToDev){
    dfProjection = ProjectToDev(objTriangleModel, lOriginYear, MaxDev)
  } else {
    dfProjection = ProjectToDate(objTriangleModel, lOriginYear, AsOfDate)
  }
  
  if (objTriangleModel@FitCategory == "DevInteger"){
    dfProjection$FitCategory = TailFunction(dfProjection$DevInteger, objTriangleModel@Tail)
    dfProjection$FitCategory = factor(dfProjection$FitCategory)
  }
  
  objTriangle = objTriangleModel@Triangle
  strResponse = objTriangleModel@Response
  
  if (is.StochasticMeasure(objTriangle, strResponse)){
    dfProjection = ProjectStochasticValues(objTriangleModel, dfProjection)
  } else {
    dfProjection[objTriangleModel@Response] = ProjectStaticValues(objTriangleModel, dfProjection)
  }
  
  row.names(dfProjection) = NULL
  
  if(is.null(AsOfDate)) AsOfDate = mdy("01-01-1900")
  
  proj = new("TriangleProjection"
             , TriangleModel = objTriangleModel
             , ProjectToDev = ProjectToDev
             , MaxDev = MaxDev
             , AsOfDate = AsOfDate
             , ProjectionData = dfProjection)
  
  proj
}