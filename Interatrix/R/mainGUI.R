InteratrixGUI <- function() {
  
  chooseMethod <- tktoplevel()
  tkwm.title(chooseMethod, "Choose method")
    
  lab <- tklabel(chooseMethod, text = "Select the method to use:")
  
  methVar <- tclVar(1)
  rb1 <- tkradiobutton(chooseMethod, text = "corrected chi-square", value = 1, variable = methVar)
  rb2 <- tkradiobutton(chooseMethod, text = "chi-square corrected also for the cumulative effect of age", value = 2, variable = methVar)

  meth <- tclvalue(methVar)
  b1 <- tkbutton(chooseMethod, text = "Cancel", command = function() tkdestroy(chooseMethod))
	b2 <- tkbutton(chooseMethod, text = "Choose", default = "active", command = function() {
    meth <- tclvalue(methVar)
    if(meth == 1) {
        tkdestroy(chooseMethod)
        chi2CorrGUI()
    }
    if(meth == 2) {
        tkdestroy(chooseMethod)
        chi2CorrAgeGUI()
    }    
  })

	tkgrid(lab, row = 0, column = 0, columnspan = 2, pady = c(10, 10))
  tkgrid(rb1, row = 1, column = 0, columnspan = 2, sticky = "w")
  tkgrid(rb2, row = 2, column = 0, columnspan = 2, sticky = "w", padx = c(0, 10))
  tkgrid(b1, row = 3, column = 0, sticky = "e", pady = c(10, 10))
  tkgrid(b2, row = 3, column = 1, sticky = "w", pady = c(10, 10))
}
