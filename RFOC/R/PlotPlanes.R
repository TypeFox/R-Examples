`PlotPlanes` <-
function(MEC, col1=1, col2=3)
{
  #  planes are already in strike and dip format
  if(is.null(MEC$UP)) { MEC$UP = FALSE }
  if(missing(col1)) col1 = 1
  if(missing(col2)) col2 = 3

  
LP1 = lowplane( MEC$az1, MEC$dip1, col=col1, UP=MEC$UP)
LP2 = lowplane( MEC$az2, MEC$dip2, col=col2, UP=MEC$UP)
invisible(list(LP1=LP1, col1=col1, LP2=LP2, col2=col2))
}

