#' plot.Triangle
#' 
#' @export 
#' @include Triangle.R
#' 
#' @param objTriangle A triangle object
#' @param Response The measure being plotted
#' @param Predictor The variable used to predict the response
#' @param Group The name of the group column used to group the data. The default is OriginPeriodStart
#' @param Lines Draw lines to connect the observations?
#' @param FitLines Draw a line of best fit? Note that fit lines will have an intercept
#' 
plotTriangle = function(objTriangle, Response
                        , Predictor, Group = "OriginPeriodStart"
                        , Lines = TRUE, FitLines = FALSE)
{
  # TODO: allow for null parameters
  df = objTriangle@TriangleData

  df$WhichGroup = factor(df[, Group])
  df$WhichResponse = df[, Response]
  df$WhichPredictor = df[,Predictor]
  
  df = df[!is.na(df$WhichPredictor), ]
  df = df[!is.na(df$WhichResponse), ]
  df = df[!is.na(df$WhichGroup), ]
  
  plotTitle = paste0(Response, " by ", Predictor, " grouped by ", Group)
  
  WhichPredictor = NULL
  WhichResponse = NULL
  WhichGroup = NULL
  plt = ggplot(df, aes(x = WhichPredictor, y = WhichResponse, group = WhichGroup, colour = WhichGroup)) 
  plt = plt + geom_point() + labs(title = plotTitle) + xlab(Predictor) + ylab(Response)
  
  if (Lines) plt = plt + geom_line() 
  if (FitLines) plt = plt + stat_smooth(method = lm, se = FALSE)
  
  plt
  return (plt)
}