"boa.print.acf" <-
function(lags = boa.par("acf.lags"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "LAGS AND AUTOCORRELATIONS:\n",
       "==========================\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      print(boa.acf(work[[i]], lags))
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
