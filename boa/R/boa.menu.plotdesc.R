"boa.menu.plotdesc" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nDESCRIPTIVE PLOT MENU\n---------------------"
   choices <- c("Back",
                "-----------------+",
                "Autocorrelations |",
                "Density          |",
                "Running Mean     |",
                "Trace            |",
                "Options...       |",
                "-----------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = boa.plot("acf"),
         "4" = boa.plot("density"),
         "5" = boa.plot("history"),
         "6" = boa.plot("trace"),
         "7" = boa.menu.setpar(c("Descriptive Plot", "Plot")),
         "8" = NULL
      )
   }

   return(abs(idx))
}
