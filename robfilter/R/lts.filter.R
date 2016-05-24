lts.filter <- function(y, width, h= floor(width/2)+1, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, method="LTS", online=online, h=h, extrapolate=extrapolate))}
