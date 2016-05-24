"boa.plot.open" <- function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   do.call(options()$device, args = list())
   created <- dev.cur()
   boa.par(dev.list = intersect(c(boa.par("dev.list"), created), dev.list()))

   return(created)
}
