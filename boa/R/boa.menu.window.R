"boa.menu.window" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   choices <- c("Back",
                "------------------------+",
                "Previous                |",
                "Next                    |",
                "Save to Postscript File |",
                "Close                   |",
                "Close All               |",
                "------------------------+")
   idx <- 1
   while(idx > 0) {
      mtitle <- dev.cur()
      mtitle <- paste("\nWINDOW ", mtitle, " MENU\n============",
                      paste(rep("=", nchar(mtitle)), collapse = ""), sep = "")
      idx <- menu(choices, title = mtitle)
      switch(idx,
         "1" = idx <- -1,
         "2" = NULL,
         "3" = dev.set(dev.prev()),
         "4" = dev.set(dev.next()),
         "5" = { cat("\nEnter name of file to which to save the plot [none]\n")
                 value <- scan(what = "", n = 1, strip.white = TRUE)
                 if(length(value) > 0) {
                    dev.print(device = postscript,
                              file = paste(boa.par("path"), value, sep = "/"))
                 }
               },
         "6" = boa.plot.close(),
         "7" = boa.plot.close(boa.par("dev.list")),
         "8" = NULL
      )
   }

   return(abs(idx))
}
