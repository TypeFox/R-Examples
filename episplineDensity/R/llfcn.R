llfcn <- function (x, epiparameters, softinfo) 
{
#
# Compute the reference density. The name of the function to compute
# that (e.g. "dnorm") had better be found in softinfo$KLDensity, with
# parameters softinfo$KLDensityParams
#
if (!any (names (softinfo) == "KLDensity"))
    stop ("No KL density supplied. How did we get here?")
if (!exists (softinfo$KLDensity))
    stop (paste ("KL density", softinfo$KLDensity, "not found"))

param.list <- softinfo$KLDensityParams
param.list$x <- x # data value is always "x" in dnorm, dexp, etc.
h <- do.call (softinfo$KLDensity, param.list)

#
# h holds the density; now compute and return the integrand 
#
if (h == 0) return (h)
return (h * log(h))
}
