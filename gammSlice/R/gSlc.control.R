gSlc.control <-
function(nBurnin = 5000, nIter = 5000, nThin = 5,fixedEffPriorVar = 1e10,sdPriorScale = 1e5)
{  
   parameters <- c(nIter,nBurnin,nThin,fixedEffPriorVar,sdPriorScale)
   if (!is.numeric(parameters)) stop("Parameters must be numbers.")

   if (nIter != as.integer(nIter)) stop("Iteration number must be integer.")
   if (nBurnin != as.integer(nBurnin)) stop("Burnin number must be integer.")
   if (nThin != as.integer(nThin)) stop("Thin number must be integer.")
   
   if ( !(fixedEffPriorVar >0)) stop("Variance must be positive.")   
   if ( !(sdPriorScale > 0)) stop("Standard Scale of Half Cauchy distribution must be positive." )

   if (nIter < 1) stop ("Iteration number must be positive.")
   if (nThin < 1) stop ("Thin number  must be positive.")
   if (nBurnin < 0) stop ("Burnin number cannot be negative.")
   
   if (nThin > nIter) stop("Iteration number must be greater than Thin number.")
  
   ret <- list(nIter=nIter,nBurnin = nBurnin, nThin = nThin, fixedEffPriorVar = fixedEffPriorVar, sdPriorScale = sdPriorScale)
   return(ret)
}
