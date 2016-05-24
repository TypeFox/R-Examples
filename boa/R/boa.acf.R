"boa.acf" <-
function(link, lags)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   pnames <- boa.pnames(link)
   result <- matrix(NA, nrow = ncol(link), ncol = length(lags),
                    dimnames = list(pnames, paste("Lag", lags)))
   lags <- lags[lags <= (nrow(link) - 1)]
   n.lags <- length(lags)
   if(n.lags > 0) {
      idx <- 1:n.lags
      lag.max <- max(lags)
      for(i in pnames) {
         result[i, idx] <- acf(link[, i], lag.max = lag.max,
                               plot = FALSE)$acf[lags + 1]
      }
   }

   return(result)
}
