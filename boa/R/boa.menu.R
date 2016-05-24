"boa.menu" <-
function(recover = FALSE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   boa.init(recover)
   if(!recover)
      cat("NOTE: if the event of a menu system crash, type\n",
          "\"boa.menu(recover = TRUE)\" to restart and recover your work.\n",
          sep = "")

   mtitle <- "\nBOA MAIN MENU\n*************"
   choices <- c("File            >>", "Data Management >>",
                "Analysis        >>", "Plot            >>",
                "Options         >>", "Window          >>")
   idx <- 1
   while(idx != 99) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- boa.menu.file(),
         "2" = boa.menu.data(),
         "3" = boa.menu.analysis(),
         "4" = boa.menu.plot(),
         "5" = boa.menu.par(),
         "6" = boa.menu.window()
      )
   }
   boa.quit()
   invisible()
}
