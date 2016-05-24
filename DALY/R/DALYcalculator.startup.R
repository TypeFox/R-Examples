## Create DALY Calculator main window

DALYcalculator.startup <-
function(){
  ## Define active window state
  DALYassign("active.windows", TRUE, item = "main")

  ## Set some variables: 'data' & 'options'
  defineVariables()
  setStdLE()

  ## Info message
  info <- paste("DALY Calculator 1.4.0 (2014-10-27)",
                "\nhttp://users.ugent.be/~bdvleess/DALYcalculator",
                "\nDeveloped and maintained by:",
                "  Brecht Devleesschauwer <Brecht.Devleesschauwer@UGent.be>",
                "\nWith contributions from:",
				"  Scott McDonald",
                "  Juanita Haagsma",
                "  Nicolas Praet",
                "  Arie Havelaar",
                "  Niko Speybroeck",
                sep = "\n")

  ## Create main DALY Calculator window
  DALYassign("main", tktoplevel(padx = 10, pady = 10))
  tkwm.title(DALYget("main"), "DALY Calculator")

  Frame1 <- tkframe(DALYget("main"), padx = 80, pady = 5,
                    relief = "groove", borderwidth = 2)
  Frame2 <- tkframe(DALYget("main"), padx = 15, pady = 10,
                    relief = "groove", borderwidth = 2)
  Frame3 <- tkframe(DALYget("main"), padx = 62, pady = 5,
                    relief = "groove", borderwidth = 2)
  Frame4 <- tkframe(DALYget("main"), pady = 10)
  
  ## Catch closing of main window
  tkbind(DALYget("main"), "<Destroy>", function(){ DALYdestroy.main() })

  ## Define menu structure
  DALYmenu <- tkmenu(DALYget("main"))
  tkconfigure(DALYget("main"), menu = DALYmenu)
  fileMenu <- tkmenu(DALYmenu, tearoff = FALSE)
  setsMenu <- tkmenu(DALYmenu, tearoff = FALSE)
  helpMenu <- tkmenu(DALYmenu, tearoff = FALSE)
  examplesMenu <- tkmenu(helpMenu, tearoff = FALSE)

  #(1) File
  tkadd(fileMenu, "command",
        label = "Load DALY data from file...",
        command = function(){ readDALYdata(NULL) })
  tkadd(fileMenu, "command",
        label = "Save DALY data to file...",
        command = function(){ saveDALYdata() })
  tkadd(fileMenu, "separator")

  tkadd(fileMenu, "command",
        label = "Reset DALY Calculator",
        command = reset)
  tkadd(fileMenu, "separator")

  tkadd(fileMenu, "command",
        label = "Exit",
        command = function(){ DALYdestroy.main() })
  tkadd(DALYmenu, "cascade",
        label = "File", menu = fileMenu)

  #(2) Settings
  tkadd(setsMenu, "command",
        label = "Life Expectancy table...",
        command = setLifeExp)
  tkadd(setsMenu, "command",
        label = "Options...",
        command = DALYoptions)
  tkadd(DALYmenu, "cascade",
        label = "Settings", menu = setsMenu)

  #(3) Help
  tkadd(helpMenu, "cascade",
        label = "Load examples", menu = examplesMenu)
  tkadd(examplesMenu, "command",
        label = "Neurocysticercosis, Cameroon (Praet et al., 2009)",
        command = function(){ setDALYexample(1) })
  tkadd(examplesMenu, "command",
        label = "Toxoplasmosis, the Netherlands (Kortbeek et al., 2009)",
        command = function(){ setDALYexample(2) })
  tkadd(helpMenu, "separator")

  tkadd(helpMenu, "command",
        label = "Html help",
        command = function(){ openHelpFile("DALYcalculator") })
  tkadd(helpMenu, "command",
        label = "DALY Calculator manual (PDF)",
        command = DALYmanual)
  tkadd(helpMenu, "separator")

  tkadd(helpMenu, "command",
        label = "Package description",
        command = function(){ openHelpFile("DALY-package") })
  tkadd(helpMenu, "command",
        label = "DALY Calculator Info",
        command = function(){ tkmessageBox(message = info,
                                           title = "Information") })
  tkadd(DALYmenu, "cascade",
        label = "Help", menu = helpMenu)

  ## Define buttons
  setPop.but <- tkbutton(Frame1, text = "set population",
                         width = 14, padx = 5, pady = 5,
                         command = setPop)

  outcome1.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(1) })
  outcome2.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(2) })
  outcome3.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(3) })
  outcome4.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(4) })
  outcome5.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(5) })
  outcome6.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(6) })
  outcome7.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(7) })
  outcome8.but <- tkbutton(Frame2, text = "set data", padx = 2,
                           command = function(){ setData(8) })

  calculate.but <- tkbutton(Frame4, text = "CALCULATE DALYs",
                            width = 37, padx = 5, pady = 10,
                            command = function() getDALY(button.call = TRUE))

  ## Define Disease & Outcome entry boxes
  diseaseEntry <- tkentry(Frame2, width = 20,
                          textvariable = DALYget("diseaseName"))
  for (i in seq(8))
    DALYassign(paste("outcome", i, "Entry", sep = ""),
           tkentry(Frame2, width = 20,
              textvariable = DALYget(paste("outcome", i, "Name", sep = ""))))

  ## Define AW & DR widgets
  awList <- c("Yes", "No")
  awCombo <- ttkcombobox(Frame3, state = "readonly", width = 5,
                         values = awList, textvariable = DALYget(".aw"))
  drEntry <- tkentry(Frame3, width = 7,
                     textvariable = DALYget(".dr"))

  ## Lay out GUI grid
  tkgrid(setPop.but)
  tkgrid(Frame1)

  bold <- tkfont.create(weight = "bold", size = "11")
  tkgrid(tklabel(Frame2, text = "disease ", font = bold), diseaseEntry)
  for (i in seq(8))
    tkgrid(tklabel(Frame2, text = paste("outcome ", i, " ", sep = "")),
	       DALYget(paste("outcome", i, "Entry", sep = "")),
		   get(paste("outcome", i, ".but", sep = "")))
  tkgrid(Frame2, pady = 10)

  tkgrid(tklabel(Frame3, text = "age weighting "), awCombo)
  tkgrid(tklabel(Frame3, text = "discount rate (%) "), drEntry)
  tkgrid(Frame3)

  tkgrid(calculate.but)
  tkgrid(Frame4)

  ## Clean-up 'DALY' database
  rm(list = paste("outcome", seq(8), "Entry", sep = ""),
     envir = DALYenv())
}