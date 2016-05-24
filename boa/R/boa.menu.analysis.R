"boa.menu.analysis" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nANALYSIS MENU\n============="
   choices <- c("Back",
                "---------------------------+",
                "Descriptive Statistics  >> |",
                "Convergence Diagnostics >> |",
                "---------------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = idx <- boa.menu.stats(),
         "4" = idx <- boa.menu.coda(),
         "5" = NULL
      )
   }

   return(abs(idx))
}
