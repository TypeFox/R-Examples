# Rearranges the arguments of the *.stat functions
#    Used in conjunction with optimize
MetricI <- function(i,fun=ksStat,...)
{
   return(fun(i=round(i),...))
}

