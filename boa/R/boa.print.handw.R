"boa.print.handw" <-
function(error = boa.par("handw.error"),
                            alpha = boa.par("alpha"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "HEIDLEBERGER AND WELCH STATIONARITY AND INTERVAL HALFWIDTH TESTS:\n",
       "=================================================================\n\n",
       "Halfwidth test accuracy = ", error, "\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      print(boa.handw(work[[i]], error, alpha))
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
