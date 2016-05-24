"boa.print.cor" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   work <- boa.chain("work")
   cat("\n",
       "CROSS-CORRELATION MATRIX:\n",
       "=========================\n", sep = "")
   for(i in names(work)) {
      header <- paste("\nChain: ", i, "\n", sep = "")
      cat(header, rep("-", nchar(header) - 2), "\n\n", sep = "")
      corr <- round(cor(work[[i]]), digits = options()$digits)
      lt <- lower.tri(corr, diag = TRUE)
      corr[!lt] <- ""
      print(corr, quote = FALSE)
      cat("\nPress <ENTER> to continue")
      readline()
   }
   invisible()
}
