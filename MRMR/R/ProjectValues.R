#' @include TriangleModel.R

ProjectStaticValues = function(objTriangleModel, dfProjection)
{
  fit = objTriangleModel@Fit
  dfProjection = predict(fit, newdata = dfProjection)
  
  dfProjection
}

ProjectStochasticValues = function(objTriangleModel, dfProjection)
{
  fit = objTriangleModel@Fit
  
  theCols = GetStochasticColumnNames(objTriangleModel@Response)
  priorCol = theCols[1]
  incrCol = theCols[2]
  cumulCol = theCols[3]
  theResponse = objTriangleModel@Response
  
  evalDates = unique(dfProjection$EvaluationDate)
  evalDates = sort(evalDates)
  
  maxDev = max(dfProjection$DevInteger)
  
  for (i in seq_along(evalDates)) {
    whichRows = (dfProjection$EvaluationDate == evalDates[i])
    dfProjection[whichRows, theResponse] = predict(fit, newdata = dfProjection[whichRows,])
          
    if (theResponse == incrCol) {
      dfProjection[whichRows, cumulCol] = dfProjection[whichRows, incrCol] + dfProjection[whichRows, priorCol]
    }else {
      dfProjection[whichRows, incrCol] = dfProjection[whichRows, cumulCol] - dfProjection[whichRows, priorCol]
    }
    
    if(i != length(evalDates)){
      nextRows = (dfProjection$EvaluationDate == evalDates[i+1])
      currentRows = (dfProjection$DevInteger %in% (dfProjection$DevInteger[nextRows] - 1)
                     & dfProjection$EvaluationDate == evalDates[i])
      dfProjection[nextRows, priorCol] = dfProjection[currentRows, cumulCol]
    }
  }
  
#  dfProjection = predict(fit, newdata = dfProjection)
  
  dfProjection
}