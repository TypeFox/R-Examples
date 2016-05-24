"boa.transform" <-
function(x, support, inv = FALSE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   transform <- x
   support <- range(support)
   bounded <- is.finite(support)
   if(!inv) {
      if(all(bounded)) {
         p <- (x - support[1]) / abs(diff(support))
         transform <- log(p / (1 - p))
      } else if(bounded[1]) {
         transform <- log(x - support[1])
      } else if(bounded[2]) {
         transform <- log(support[2] - x)
      }
   } else {
      if(all(bounded)) {
         transform <- abs(diff(support)) / (1 + exp(-x)) + support[1]
      } else if(bounded[1]) {
         transform <- exp(x) + support[1]
      } else if(bounded[2]) {
         transform <- support[2] - exp(x)
      }
   }

   return(transform)
}
