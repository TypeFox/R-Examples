"boa.menu.par" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nGLOBAL OPTIONS MENU\n==================="
   choices <- c("Back",
                "----------------+",
                "Analysis...     |",
                "Import Files... |",
                "Plot...         |",
                "All...          |",
                "----------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = boa.menu.setpar(c("Descriptive", "Convergence")),
         "4" = boa.menu.setpar("Import"),
         "5" = boa.menu.setpar(c("Descriptive Plot", "Convergence Plot", "Plot")),
         "6" = boa.menu.setpar(),
         "7" = NULL
      )
   }

   return(abs(idx))
}
