#' PlotModelGoF
#' @export PlotModelGoF
#' @include TriangleModel.R
#' 
#' @description
#' This function will plot the F distribution associated with the TriangleModel, along with a vertical line indicating the 
#' F statistic for this model.
#' 
#' @param 
#' objTriangleModel A TriangleModel object
#' 
#' @return A vector of intervals
#' 
#' @seealso \code{\link{PlotModelFactors}}
#' 
PlotModelGoF = function(objTriangleModel)
{
  fit = objTriangleModel@Fit
  
  fstat = summary(fit)$fstatistic
  df1 = fstat[2]
  df2 = fstat[3]
  
  xlow = qf(0.025, df1, df2)
  xhigh = qf(0.975, df1, df2)
  x = seq(xlow, xhigh, length.out = 200)
   
  df = data.frame(x = x, y = df(x, df1, df2))
  y = NULL
  plt = ggplot(df, aes(x, y)) + geom_line() + geom_vline(xintercept = fstat[1])
  plt
  
}