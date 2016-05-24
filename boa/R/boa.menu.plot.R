"boa.menu.plot" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nPLOT MENU\n========="
   choices <- c("Back",
                "---------------------------+",
                "Descriptive             >> |",
                "Convergence Diagnostics >> |",
                "---------------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = idx <- boa.menu.plotdesc(),
         "4" = idx <- boa.menu.plotcoda(),
         "5" = NULL
      )
   }

   return(abs(idx))
}
