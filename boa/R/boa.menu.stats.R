"boa.menu.stats" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nDESCRIPTIVE STATISTICS MENU\n---------------------------"
   choices <- c("Back",
                "---------------------------------------+",
                "Autocorrelations                       |",
                "Correlation Matrix                     |",
                "Highest Probability Density Intervals  |",
                "Summary Statistics                     |",
                "Options...                             |",
                "---------------------------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = boa.print.acf(),
         "4" = boa.print.cor(),
         "5" = boa.print.hpd(),
         "6" = boa.print.stats(),
         "7" = boa.menu.setpar("Descriptive"),
         "8" = NULL
      )
   }

   return(abs(idx))
}
