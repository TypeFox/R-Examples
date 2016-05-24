"boa.gewekePwr" <- function(link)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   spec <- NULL
   pnames <- boa.pnames(link)
   if(nrow(link) > 10) {
      for(i in pnames) {
         spec <- c(spec, spectrum0(link[, i])$spec)
      }
   } else {
      spec <- rep(NA, length(pnames))
   }

   return(structure(spec, names = pnames))
}
