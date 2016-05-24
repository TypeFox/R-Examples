"boa.menu.coda" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nCONVERGENCE DIAGNOSTICS MENU\n----------------------------"
   choices <- c("Back",
                "------------------------+",
                "Brooks, Gelman, & Rubin |",
                "Geweke                  |",
                "Heidelberger & Welch    |",
                "Raftery & Lewis         |",
                "Options...              |",
                "------------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = boa.print.gandr(),
         "4" = boa.print.geweke(),
         "5" = boa.print.handw(),
         "6" = boa.print.randl(),
         "7" = boa.menu.setpar("Convergence"),
         "6" = NULL
      )
   }

   return(abs(idx))
}
