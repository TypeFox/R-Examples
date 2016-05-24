"boa.menu.file" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   mtitle <- "\nFILE MENU\n========="
   choices <- c("Back",
                "-----------------------+",
                "Import Data         >> |",
                "Save Session           |",
                "Load Session           |",
                "Exit BOA               |",
                "-----------------------+")
   idx <- 1
   while(idx > 0) {
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = idx <- boa.menu.import(),
         "4" = { saveas <- "y"
                 cat("\nEnter name of object to which to save the session",
                     "data [none]\n")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(length(value)) {
                    if(!boa.save(value)) {
                       cat("Object already exists.  Overwrite (y/n) [n]?\n")
                       saveas <- scan(what = "", n = 1, strip.white = TRUE)
                    }
                    if(length(saveas) > 0 && saveas == "y") {
                       boa.save(value, replace = TRUE)
                       cat("+++ Data successfully saved +++\n")
						  }
                 }
               },
         "5" = { cat("\nEnter name of object to load [none]\n")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(length(value) && boa.load(value))
                    cat("+++ Data successfully loaded +++\n")
               },
         "6" = { cat("\nDo you really want to EXIT (y/n) [n]?\n")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(length(value) > 0 && value == "y")  idx <- -99
               },
         "7" = NULL
      )
   }

   return(abs(idx))
}
