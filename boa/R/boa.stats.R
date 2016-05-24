"boa.stats" <-
function(link, probs, batch.size)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   niter <- nrow(link)
   nparms <- ncol(link)
   sd <- sqrt(colVars(link))
   batch <- t(boa.batchMeans(link, batch.size))
   result <- cbind(colMeans(link),
                   sd,
                   sd / sqrt(niter),
                   sqrt(boa.gewekePwr(link) / niter),
                   sqrt(colVars(batch) / nrow(batch)),
                   boa.acf(batch, 1),
                   t(apply(link, 2, quantile, probs = probs)),
                   matrix(rep(range(boa.iter(link)), nparms), nrow = nparms,
                          ncol = 2, byrow = TRUE),
                   rep(niter, nparms))
   dimnames(result)[[2]] <- c("Mean", "SD", "Naive SE", "MC Error", "Batch SE",
                              "Batch ACF", probs, "MinIter", "MaxIter",
                              "Sample")

   return(result)
}
