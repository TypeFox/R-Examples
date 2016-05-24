"boa.plot.close" <- function(which = dev.cur())
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   shutdown <- NULL
   current <- boa.par("dev.list")
   idx <- is.element(current, which)
   for(i in intersect(current[idx], dev.list())) {
      shutdown <- dev.off(i)
   }
   boa.par(dev.list = current[!idx])

   return(shutdown)
}
