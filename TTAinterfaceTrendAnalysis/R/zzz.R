.onAttach <- function(lib,pkg) {

        themes <- tk2theme.list()
		if ("aqua" %in% themes) { # This must be aquaTk on a Mac
			try(tk2theme("aqua"), silent = TRUE)
		} else if ("vista" %in% themes) { # This must be Vista or Win 7
			try(tk2theme("vista"), silent = TRUE)
		} else if ("xpnative" %in% themes) { # This must be XP
			try(tk2theme("xpnative"), silent = TRUE)
		} else if ("winnative" %in% themes) { # This must be a pre-XP windows
			try(tk2theme("winnative"), silent = TRUE)
		} else if ("radiance" %in% themes) {
			try(tk2theme("radiance"), silent = TRUE)
		} else { # A modern "default" theme that fit not too bad in many situations
			try(tk2theme("clearlooks"), silent = TRUE) }	
			
  startup <- tktoplevel()
  tkwm.geometry(startup, "230x110")
  tkwm.resizable(startup, 0,0)
  tktitle(startup) <- "Start Panel"
  LOAD1 <- function() {  tclRequire("BWidget")
  Envir$sel.theme <- tclvalue(Theme)                                                       # exporte la valeur du theme dans l'interface
  TTAinterface()   }
  imgStart <- tclVar()                                                                                                
  tcl("image","create","photo",imgStart,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgStart.gif",fsep=.Platform$file.sep))
  LOAD1.but <- tk2button(startup, image=imgStart, text=" Ready to Start ! ", compound="right", command=LOAD1)
  
  label.theme <- tklabel(startup, text="Choose a theme:")
  themes.list <- c(tk2theme.list(), "clearlooks", "Auto")
  Theme <- tclVar("Auto")
  cb.theme <- tk2combobox(startup, values=themes.list, textvariable=Theme, state="readonly") 
  
  tkpack(tklabel(startup, text=""), side="top")
  tkpack(LOAD1.but, side="top")
  tkpack(label.theme, cb.theme, side="top")
  tkpack(tk2separator(startup), fill="x", side="top")

  Envir$pversion <- utils::packageVersion("TTAinterfaceTrendAnalysis")
  subtext <- tklabel(startup,text= paste("TTAinterface v", Envir$pversion, " launch panel", sep=""))
  tkconfigure(subtext, font=tkfont.create(size=7))
  tkpack(subtext, side="top", fill="both")

  tk2tip(subtext,
  if(getRversion() >= "3.0.0") {
   text= paste("Need R v3.0.0+ ; Your R version:", getRversion(), "(Good)") } else {
   text= paste("Need R v3.0.0+ ; Your R version:", getRversion(), "(Need update)") } )

}

if(getRversion() >= "3.0.0") globalVariables(names=c("Category","Salinity","Depth","MONTHS", "param", "depth","sal","site","npsu"), package="TTAinterfaceTrendAnalysis")