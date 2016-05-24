"boa.getparms" <-
function(link, pnames)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   idx <- is.element(dimnames(link)[[2]], pnames)
   if(any(idx))  result <- link[, idx, drop = FALSE]

   return(result)
}
