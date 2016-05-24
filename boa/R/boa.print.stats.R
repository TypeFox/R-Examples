"boa.print.stats" <-
function(probs = boa.par("quantiles"),
                            batch.size = boa.par("batch.size"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "SUMMARY STATISTICS:\n",
       "===================\n\n",
       "Bin size for calculating Batch SE and (Lag 1) ACF = ", batch.size,
       "\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      print(boa.stats(work[[i]], probs, batch.size))
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
