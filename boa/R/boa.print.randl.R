"boa.print.randl" <-
function(q = boa.par("randl.q"),
                            error = boa.par("randl.error"),
                            prob = 1 - boa.par("alpha"),
                            delta = boa.par("randl.delta"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "RAFTERY AND LEWIS CONVERGENCE DIAGNOSTIC:\n",
       "=========================================\n\n",
       "Quantile = ", q, "\n",
       "Accuracy = +/- ", error, "\n",
       "Probability = ", prob, "\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      randl <- boa.randl(work[[i]], q, error, prob, delta)
      if(is.matrix(randl)) {
         print(randl)
      } else {
         cat("******* Warning *******\n",
             "Available chain length is ", nrow(work[[i]]), ".\n",
             "Re-run simulation for at least ", randl, " iterations\n",
             "OR reduce the quantile, accuracy, or probability to be ",
             "estimated.\n", sep = "")
      }
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
