"boa.menu.setpar" <-
function(group)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   if(missing(group)) {
      value <- boa.print.par()
   } else {
      value <- boa.print.par(group)
   }
   par.names <- value[, "par"]
   par.notes <- value[, "note"]
   cat("\n")
   idx <- ""
   while((length(idx) > 0) && !is.element(idx, seq(par.names))) {
      cat("Specify parameter to change or press <ENTER> to continue\n")
      idx <- scan(what = "", n = 1, strip.white = TRUE)
   }
   if(length(idx) > 0) {
      cat("\n")
      idx <- as.numeric(idx)
      if(nchar(par.notes[idx]))  cat("DESCRIPTION:", par.notes[idx], "\n")
      switch(data.class(boa.par(par.names[idx])),
         "numeric"   = value <- boa.getinput("\nEnter new numeric value\n"),
         "character" = { cat("\nEnter new character string\n")
                         value <- scan(what = "", n = 1, sep = "\n")
                       },
         "function"  = value <- boa.getinput("\nEnter new function followed by a blank line\n",
                                             n = -1),
         "list" = value <- boa.getinput("\nEnter new list of values followed by a blank line\n",
                                        n = -1),
         "logical"   = value <- boa.getinput("\nEnter new logical value\n"),
         value <- NULL
      )
      boa.par(structure(list(value), names = par.names[idx]))
   }
   invisible()
}
