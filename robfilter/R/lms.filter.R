lms.filter <- function(y, width, online=FALSE, extrapolate=TRUE)
{return(robreg.filter(y=y, width=width, method="LMS", online=online, extrapolate=extrapolate))}
