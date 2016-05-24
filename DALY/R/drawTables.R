## Helper function for 'setData()'
## Draw 'data' window
## See 'drawWindow()' for definition of window components

drawTables <-
function(tblInc, tblOns, tblDWt, tblMrt, tblTrt, tblDur, tblDWn, tblLxp, 
         cbdInc, cbdOns, cbdDWt, cbdMrt, cbdTrt, cbdDur, cbdDWn, cbdLxp, 
         cbsInc, cbsOns, cbsDWt, cbsMrt, cbsTrt, cbsDur, cbsDWn, cbsLxp){
  small <- tkfont.create(size = "2")
  label <- tklabel(DALYget("Frame1"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  incidence: cases/1000 persons/year")
  tkgrid(label, cbdInc, cbsInc)
  tkgrid.configure(cbdInc, sticky = "e")
  tkgrid.configure(cbsInc, sticky = "e")
  tkgrid(tblInc, columnspan = 3)
  tkgrid(DALYget("Frame1"))
  tkgrid(tklabel(DALYget("Left"), text = "", font = small))

  label <- tklabel(DALYget("Frame2"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  onset: age in years")
  tkgrid(label, cbdOns, cbsOns)
  tkgrid.configure(cbdOns, sticky = "e")
  tkgrid.configure(cbsOns, sticky = "e")
  tkgrid(tblOns, columnspan = 3)
  tkgrid(DALYget("Frame2"))
  tkgrid(tklabel(DALYget("Left"), text = "", font = small))

  label <- tklabel(DALYget("Frame3"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  DW treated: range [0-1]")
  tkgrid(label, cbdDWt, cbsDWt)
  tkgrid.configure(cbdDWt, sticky = "e")
  tkgrid.configure(cbsDWt, sticky = "e")
  tkgrid(tblDWt, columnspan = 3)
  tkgrid(DALYget("Frame3"))
  tkgrid(tklabel(DALYget("Left"), text = "", font = small))

  label <- tklabel(DALYget("Frame4"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  mortality: deaths/1000 persons/year")
  tkgrid(label, cbdMrt, cbsMrt)
  tkgrid.configure(cbdMrt, sticky = "e")
  tkgrid.configure(cbsMrt, sticky = "e")
  tkgrid(tblMrt, columnspan = 3)
  tkgrid(DALYget("Frame4"))

  label <- tklabel(DALYget("Frame5"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  treatment: proportion [0-1]")
  tkgrid(label, cbdTrt, cbsTrt)
  tkgrid.configure(cbdTrt, sticky = "e")
  tkgrid.configure(cbsTrt, sticky = "e")
  tkgrid(tblTrt, columnspan = 3)
  tkgrid(DALYget("Frame5"))
  tkgrid(tklabel(DALYget("Right"), text = "", font = small))

  label <- tklabel(DALYget("Frame6"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  duration: years")
  tkgrid(label, cbdDur, cbsDur)
  tkgrid.configure(cbdDur, sticky = "e")
  tkgrid.configure(cbsDur, sticky = "e")
  tkgrid(tblDur, columnspan = 3)
  tkgrid(DALYget("Frame6"))
  tkgrid(tklabel(DALYget("Right"), text = "", font = small))

  label <- tklabel(DALYget("Frame7"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  DW untreated: range [0-1]")
  tkgrid(label, cbdDWn, cbsDWn)
  tkgrid.configure(cbdDWn, sticky = "e")
  tkgrid.configure(cbsDWn, sticky = "e")
  tkgrid(tblDWn, columnspan = 3)
  tkgrid(DALYget("Frame7"))
  tkgrid(tklabel(DALYget("Right"), text = "", font = small))

  label <- tklabel(DALYget("Frame8"),
                   bg = "white", width = "36", anchor = "w",
                   text = "  average age at death: age in years")
  tkgrid(label, cbdLxp, cbsLxp)
  tkgrid.configure(cbdLxp, sticky = "e")
  tkgrid.configure(cbsLxp, sticky = "e")
  tkgrid(tblLxp, columnspan = 3)
  tkgrid(DALYget("Frame8"))

  tkgrid(DALYget("Left"), DALYget("Right"))
  tkgrid(DALYget("save.but"), DALYget("cancel.but"))
  tkgrid(DALYget("Top"))
  tkgrid(DALYget("Bottom"))
}