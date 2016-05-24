## Helper function for 'setData()'
## Define 'data' window components
## See 'drawTables()' for drawing window components

drawWindow <-
function(window, outcomeName, this){
  tkwm.title(window,
             paste("Data:",
                   ifelse(DALYtclvalue("diseaseName") != "",
                          DALYtclvalue("diseaseName"),
                          "disease"),
                   ">",
                   ifelse(tclvalue(outcomeName) != "",
                          tclvalue(outcomeName),
                          paste("outcome", this))))
  tkwm.geometry(window, "+0+0")

  DALYassign("Top", tkframe(window))
  DALYassign("Left", tkframe(DALYget("Top"), padx = 15, pady = 10))
  DALYassign("Right", tkframe(DALYget("Top"), padx = 15, pady = 10))
  DALYassign("Bottom", tkframe(window, padx = 0, pady = 10))

  for (i in seq(1, 4))
    DALYassign(paste("Frame", i, sep = ""),
               tkframe(DALYget("Left"),
                       relief = "groove", borderwidth = 2, bg = "white"))
  for (i in seq(5, 8))
    DALYassign(paste("Frame", i, sep = ""),
               tkframe(DALYget("Right"),
                       relief = "groove", borderwidth = 2, bg = "white"))

  lbl <- lapply(DALYget("txtlbl"), paste, this, sep = "")
  Lbl <- lapply(DALYget("txtLbl"), paste, this, sep = "")
  save <- c(paste(".", lbl, sep = ""),
          paste(".dist", Lbl, sep = ""),
		  paste(".strat", Lbl, sep = ""))
  chck <- paste(".", lbl, sep = "")

  DALYassign("save.but",
             tkbutton(DALYget("Bottom"), width = 10,
                      text = "OK",
                      command = function(){ saveWindow(window, save, chck) }))
  DALYassign("cancel.but",
             tkbutton(DALYget("Bottom"), width = 10,
                      text = "Cancel",
                      command = function(){ cancelWindow(window) }))
}