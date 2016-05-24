"boa.print.gandr" <-
function(alpha = boa.par("alpha"),
                            win = boa.par("gandr.win"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{    work <- boa.chain("work")
   work.support <- boa.chain("work.support")
   cat("\n",
       "BROOKS, GELMAN, AND RUBIN CONVERGENCE DIAGNOSTICS:\n",
       "==================================================\n\n", sep = "")
   gandr <- boa.chain.gandr(work, work.support, alpha, window = win)
   cat("Iterations used = ", paste(gandr$window, collapse=":"), "\n\n",
       sep = "")
   cat("Potential Scale Reduction Factors\n",
       "---------------------------------\n\n", sep = "")
   print(gandr$psrf)
   cat("\n",
       "Multivariate Potential Scale Reduction Factor = ",
       round(gandr$mpsrf, digits = options()$digits), "\n", sep = "")
   cat("\n",
       "Corrected Scale Reduction Factors\n",
       "---------------------------------\n\n", sep = "")
   print(gandr$csrf)
   cat("\nPress <ENTER> to continue")
   readline()
   invisible()
}
