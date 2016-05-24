med.filter <- function(y, width, minNonNAs = 5, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, minNonNAs=minNonNAs, method="MED", online=online, extrapolate=extrapolate))}
