"boa.batchMeans" <-
function(link, size)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   result <- NULL
   riter <- range(boa.iter(link))
   iter.first <- riter[1]
   iter.last <- iter.first + size - 1
   while(iter.last <= riter[2]) {
      result <- cbind(result, colMeans(boa.getiter(link, iter.first:iter.last)))
      iter.first <- iter.last + 1
      iter.last <- iter.last + size
   }
   if(is.null(result))  result <- matrix(NA, nrow = ncol(link), ncol = 1)
   dimnames(result) <- list(boa.pnames(link), 1:ncol(result))

   return(result)
}
