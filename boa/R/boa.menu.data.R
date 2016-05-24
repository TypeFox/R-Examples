"boa.menu.data" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nDATA MANAGEMENT MENU\n===================="
   choices <- c("Back",
                "---------------------------+",
                "Chains                  >> |",
                "Parameters              >> |",
                "Display Working Dataset    |",
                "Display Master Dataset     |",
                "Reset                      |",
                "---------------------------+")
   idx <- 1
   while(idx > 0) {
      sync <- boa.chain("work.sync")
      if(sync) {
         choices[7] <- "*****                      |"
      } else {
         choices[7] <- "Reset                      |"
      }
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = idx <- boa.menu.chains(),
         "4" = idx <- boa.menu.parms(),
         "5" = { boa.print.info()
                 cat("\nPress <ENTER> to continue")
                 readline()
               },
         "6" = { boa.print.info("master")
                 cat("\nPress <ENTER> to continue")
                 readline()
               },
         "7" = { if(!sync) {
                    boa.chain.reset()
                    cat("+++ Master Dataset successfully copied to Working",
                        "Dataset +++\n")
                 }
               },
         "8" = NULL
      )
   }

   return(abs(idx))
}
