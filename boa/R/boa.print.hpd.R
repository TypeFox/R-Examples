"boa.print.hpd" <-
function(alpha = boa.par("alpha"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "HIGHEST PROBABILITY DENSITY INTERVALS:\n",
       "======================================\n\n",
       "Alpha level = ", alpha, "\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      link <- work[[i]]
      pnames <- boa.pnames(link)
      R <- NULL
      for(j in pnames) {
         R <- rbind(R, boa.hpd(link[, j], alpha))
      }
      dimnames(R)[[1]] <- pnames
      print(R)
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
