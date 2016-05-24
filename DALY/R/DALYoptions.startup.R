## Create 'DALY options' window

DALYoptions.startup <-
function(){
  ## Define active window state
  DALYassign("active.windows", TRUE, item = "opt.win")
  
  ## Create 'DALY options' window
  DALYassign("opt.win", tktoplevel(padx = 15, pady = 15))
  tkwm.title(DALYget("opt.win"), "Options")

  opt.top <- tkframe(DALYget("opt.win"), padx = 5, pady = 5,
                  relief = "groove", borderwidth = 2)
  opt.mid <- tkframe(DALYget("opt.win"), padx = 5, pady = 5,
                  relief = "groove", borderwidth = 2)
  opt.bottom <- tkframe(DALYget("opt.win"), padx = 5, pady = 5)

  ## Catch closing of main window
  tkbind(DALYget("opt.win"), "<Destroy>",
         function(){ DALYdestroy("opt.win") })
  
  ## Read 'it' from 'DALY' database
  DALYupdate(".it")
  itEntry <- tkentry(opt.top, width = 8, textvariable = DALYget(".it"))

  ## Read options from 'DALY' database
  DALYupdate(".optOP")
  DALYupdate(".optOC")
  DALYupdate(".optRA")
  DALYupdate(".optHist")
  
  ## Define comboboxes & check button
  cbOPval <- c("Summed over age/sex classes", "Per age/sex class")
  cbOCval <- c("Summed over outcomes", "Per outcome")
  cbRAval <- c("Absolute", "Relative (per 1000 pop)")
  cbOutput <- ttkcombobox(opt.mid, state = "readonly", width = 30,
                          values = cbOPval, textvariable = DALYget(".optOP"))
  cbOutcomes <- ttkcombobox(opt.mid, state = "readonly", width = 30,
                          values = cbOCval, textvariable = DALYget(".optOC"))
  cbRelative <- ttkcombobox(opt.mid, state = "readonly", width = 30,
                          values = cbRAval, textvariable = DALYget(".optRA"))


  checkHist <- tkcheckbutton(opt.mid, variable = DALYget(".optHist"),
                             text = "DALY histogram")

  ## Define buttons
  save.but  <- tkbutton(opt.bottom, width = 10,
                        text = "OK",
                        command = function(){
						            saveOptWin(
									  list(".it", ".optHist",
									       ".optOP", ".optOC", ".optRA"),
                                      list(".it")) })
  cancel.but  <- tkbutton(opt.bottom, width = 10,
                          text = "Cancel",
                          command = function(){
						              cancelWindow(DALYget("opt.win")) })
									  
  ## Lay out GUI grids
  tkgrid(tklabel(opt.top, text = "Iterations   "), itEntry)
  tkgrid(opt.top, sticky = "ew")

  tkgrid(tklabel(DALYget("opt.win"), text = ""))
  
  tkgrid(tklabel(opt.mid, text = "Output",
                 font = tkfont.create(weight = "bold", size = "9")),
         columnspan = 2, sticky = "w")
  tkgrid(cbOutput, sticky = "w")
  tkgrid(cbOutcomes, sticky = "w")
  tkgrid(cbRelative, sticky = "w")
  tkgrid(checkHist, sticky = "w")
  tkgrid(opt.mid, sticky = "ew")

  tkgrid(tklabel(DALYget("opt.win"), text = ""))

  tkgrid(save.but, cancel.but)
  tkgrid(opt.bottom)
}