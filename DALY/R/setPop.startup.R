## Create 'population' window

setPop.startup <-
function(){
  ## Define active window state
  DALYassign("active.windows", TRUE, item = "pop.win")

  ## Define main window
  DALYassign("pop.win", tktoplevel(padx = 15, pady = 15))
  tkwm.title(DALYget("pop.win"), "Population")
  pop.top <- tkframe(DALYget("pop.win"), padx = 5, pady = 5)
  pop.bottom <- tkframe(DALYget("pop.win"), padx = 5, pady = 5)
  
  ## Catch closing of main window
  tkbind(DALYget("pop.win"), "<Destroy>",
         function(){ DALYdestroy("pop.win") })
  
  ## Read 'pop' from 'DALY' database
  DALYupdate(".pop")

  pop.table <- tkwidget(pop.top, "table", resizeborders = "none",
                        selectmode = "extended", multiline = "0",
                        rowseparator = "\n", colseparator = "\t",
                        rows = 6, cols = 3, titlerows = 1, titlecols = 1,
                        colwidth = 15, background = "white",
                        variable = DALYget(".pop"))

  tkgrid(pop.table)

  save.but  <- tkbutton(pop.bottom, width = 10,
                        text = "OK",
                        command = function(){
                                    saveWindow(DALYget("pop.win"),
                                    list(".pop"), list(".pop")) })
  cancel.but  <- tkbutton(pop.bottom, width = 10,
                          text = "Cancel",
                          command = function(){ 
						              cancelWindow(DALYget("pop.win")) })
  tkgrid(save.but, cancel.but)

  tkgrid(pop.top)
  tkgrid(pop.bottom)
}