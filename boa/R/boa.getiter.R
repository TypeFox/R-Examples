"boa.getiter" <-
function(link, iter)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   idx <- is.element(dimnames(link)[[1]], iter)
   if(any(idx))  result <- link[idx, , drop = FALSE]

   return(result)
}
