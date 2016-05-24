lqd.filter <- function(y, width, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, method="LQD", online=online, extrapolate=extrapolate))}
