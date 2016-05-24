"boa.print.geweke" <-
function(p.first = boa.par("geweke.first"),
                             p.last = boa.par("geweke.last"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "GEWEKE CONVERGENCE DIAGNOSTIC:\n",
       "==============================\n\n",
       "Fraction in first window = ", p.first, "\n",
       "Fraction in last window = ", p.last, "\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      print(t(boa.geweke(work[[i]], p.first, p.last)))
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
