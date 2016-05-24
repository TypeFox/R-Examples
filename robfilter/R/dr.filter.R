dr.filter <- function(y, width, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, method="DR", online=online, extrapolate=extrapolate))}
