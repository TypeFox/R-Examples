## Create 'Life Expectancy' window

setLifeExp.startup <-
function(){
  ## Define active window state
  DALYassign("active.windows", TRUE, item = "LE.win")
  
  ## Read 'LE' from 'DALY' database
  DALYupdate(".LE")
  
  ## Create 'Life Expectancy' window
  DALYassign("LE.win", tktoplevel(padx = 15, pady = 10))
  tkwm.title(DALYget("LE.win"), "Life Expectancy")

  ## Catch closing of main window
  tkbind(DALYget("LE.win"), "<Destroy>", function(){ DALYdestroy("LE.win") })
  
  LE.top <- tkframe(DALYget("LE.win"), padx = 5, pady = 10)
  LE.bottom <- tkframe(DALYget("LE.win"), padx = 5, pady = 10)
  reset.but <- tkbutton(LE.bottom, width = 25,
                        text = "reset standard life expectancy",
                        command = setStdLE)

  LE.table <- tkwidget(LE.top, "table", resizeborders = "none",
                       selectmode = "extended", multiline = "0",
                       rowseparator = "\n", colseparator = "\t",
                       rows = 22, cols = 3, colwidth = 12,
                       titlerows = 1, titlecols = 1, background = "white",
                       variable = DALYget(".LE"))

  save.but  <- tkbutton(LE.bottom, width = 10,
                        text = "OK",
                        command = function(){
                                    saveWindow(DALYget("LE.win"),
                                    list(".LE"), list(".LE"),
									allow.null = FALSE) })
  cancel.but  <- tkbutton(LE.bottom, width = 10,
                          text = "Cancel",
                          command = function(){
						              cancelWindow(DALYget("LE.win")) })

  tkgrid(LE.table)
  tkgrid(LE.top)

  tkgrid(save.but, cancel.but)
  tkgrid(tklabel(LE.bottom, text = ""))
  tkgrid(reset.but, columnspan = 2)
  tkgrid(LE.bottom)
}