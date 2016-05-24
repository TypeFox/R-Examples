"boa.menu.plotcoda" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- paste("\nCONVERGENCE DIAGNOSTICS PLOT MENU",
                   "\n---------------------------------", sep = "")
   choices <- c("Back",
                "----------------+",
                "Brooks & Gelman |",
                "Gelman & Rubin  |",
                "Geweke          |",
                "Options...      |",
                "----------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = boa.plot("bandg"),
         "4" = boa.plot("gandr"),
         "5" = boa.plot("geweke"),
         "6" = boa.menu.setpar(c("Convergence Plot", "Plot")),
         "7" = NULL
      )
   }

   return(abs(idx))
}
