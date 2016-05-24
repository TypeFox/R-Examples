rm.filter <- function(y, width, minNonNAs=5, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, method="RM", online=online, extrapolate=extrapolate, minNonNAs=minNonNAs))}
