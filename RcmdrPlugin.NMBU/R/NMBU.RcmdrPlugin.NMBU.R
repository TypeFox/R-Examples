# $Id: NMBU.RcmdrPlugin.NMBU.R 35 2014-01-10 21:17:26Z khliland $

#
# RcmdrPlugin.NMBU - Extensions of the R Commander for statistical teaching
#                       at the Norwegian University of Life Sciences.
#
# Organisation:
# NMBU.Utilities.R  - small functions and look-ups
# NMBU.GUI.Data.R   - graphical user interface functions
# NMBU.GUI.Statistics.R - graphical user interface functions
# NMBU.GUI.Models.R - graphical user interface functions
# NMBU.GUI.File.R   - graphical user interface functions
# NMBU.GUI.Graphs.R - graphical user interface functions
# NMBU.Graphics.R   - plot functions
# NMBU.Statistics.R - statistical functions

.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
            Commander()
        }
    }
	putRcmdr("allFileName",NULL) # Prepare variable on package load
#	options(contrasts=c('contr.sum','contr.poly'))
	tkbind(CommanderWindow(), "<Control-e>", onExportE)
	tkbind(CommanderWindow(), "<Control-E>", onExportN)
	if(!("package:nmbu"%in%search())){
	  packageStartupMessage(paste("\n----------------------------------------------------------"), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("Re-open a closed R Commander with the command: Commander()\nPlease report bugs to kristian.liland@nmbu.no\n"), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("Currently installed:\n", R.Version()$version.string, "", sep=""), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("R Commander version ", packageDescription("Rcmdr")[["Version"]], "", sep=""), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("mixlm version ", packageDescription("mixlm")[["Version"]], "", sep=""), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("RcmdrPlugin.NMBU version ", packageDescription("RcmdrPlugin.NMBU")[["Version"]], "", sep=""), domain="R-RcmdrPlugin.NMBU")
	  packageStartupMessage(paste("----------------------------------------------------------\n"), domain="R-RcmdrPlugin.NMBU")
	}
}

.onLoad <- function(...){
	packagesAvailable <- function(packages){
		sapply(sapply(packages, find.package, quiet=TRUE),
				function(x) length(x) != 0)
	}
	if (!interactive()) return()
	save.options <- options(warn=-1)
	on.exit(options(save.options))
	required.packages <- rev(c("abind", "e1071",  "lattice", "leaps", "xtable",
					 "multcomp", "lme4", "gmodels", "pbkrtest", "vcd", "mvtnorm", "nnet"))
	if (WindowsP()) required.packages <- c(required.packages, c("RODBC", "XLConnect"))
	check <- options("Rcmdr")[[1]]$check.packages
	if (length(check) > 0 && !check) return()
	packages.to.check <- required.packages
	available.packages <- packagesAvailable(packages.to.check)
	missing.packages <- packages.to.check[!available.packages]
	if (any(!available.packages)) {
		response <- tkmessageBox(message=paste(gettext("The following packages used by RcmdrPlugin.NMBU are missing:\n", domain="R-Rcmdr"),
						paste(missing.packages, collapse=", "),
						gettext("\nWithout these packages, some features will not be available.", domain="R-Rcmdr"),
						gettext("\nInstall these packages?", domain="R-Rcmdr")),
				icon="error", type="yesno")
		if (tclvalue(response) == "yes") {
			top <- tktoplevel(borderwidth=10)
#			tkwm.iconbitmap(top, system.file("etc", "R-logo.ico", package="Rcmdr"))
			tkwm.title(top, gettext("Install Missing Packages", domain="R-Rcmdr"))
			locationFrame <- tkframe(top)
			locationVariable <- tclVar("CRAN")
			CRANbutton <- ttkradiobutton(locationFrame, variable=locationVariable, value="CRAN")
#         Note: Bioconductor code not currently necessary
#            BioconductorButton <- ttkradiobutton(locationFrame, variable=locationVariable, value="Bioconductor")
			localButton <- ttkradiobutton(locationFrame, variable=locationVariable, value="local")
			directoryVariable <- tclVar("")
			directoryFrame <- tkframe(locationFrame)
			onBrowse <- function(){
				tclvalue(directoryVariable) <- tclvalue(tkchooseDirectory(parent=top))
			}
			browseButton <- buttonRcmdr(directoryFrame, text=gettext("Browse...", domain="R-Rcmdr"), width="12", command=onBrowse, borderwidth=3)
			locationField <- ttkentry(directoryFrame, width="20", textvariable=directoryVariable)
			locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
					command=function(...) tkxview(locationField, ...))
			tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
			tkgrid(labelRcmdr(top, text=gettext("Install Packages From:", domain="R-Rcmdr"), fg="blue"), sticky="nw")
			tkgrid(labelRcmdr(directoryFrame, text=gettext("Specify package  \ndirectory:", domain="R-Rcmdr"), justify="left"),
					locationField, sticky="w")
			tkgrid(browseButton, locationScroll, sticky="w")
			tkgrid(locationScroll, sticky="ew")
			tkgrid(labelRcmdr(locationFrame, text="CRAN"), CRANbutton, sticky="w")
#            tkgrid(labelRcmdr(locationFrame, text="Bioconductor"), BioconductorButton, sticky="w")
			tkgrid(labelRcmdr(locationFrame, text=gettext("Local package directory\n(must include PACKAGES index file)", domain="R-Rcmdr"),
							justify="left"), localButton, directoryFrame, sticky="nw")
			tkgrid(locationFrame, sticky="w")
			tkgrid(labelRcmdr(top, text=""))
			onOK <- function(){
				errorMessage <- function() tkmessageBox(message=paste(
									gettext("The following packages were not found at the specified location:\n", domain="R-Rcmdr"),
									paste(missing.packages[!present], collapse=", ")),  icon="warning", type="ok")
				tkgrab.release(top)
				tkdestroy(top)
				location <- tclvalue(locationVariable)
				if (location == "CRAN") {
					packages <- available.packages()[,1]
					present <- missing.packages %in% packages
					if (!all(present)) errorMessage()
					if (!any(present)) return()
					install.packages(missing.packages[present], lib=.libPaths()[1])		
				}
#                else if (location == "Bioconductor") {
#                    packages <- CRAN.packages(CRAN=getOption("BIOC"))[,1]
#                    present <- missing.packages %in% packages
#                    if (!all(present)) errorMessage()
#                    install.packages(missing.packages[present], lib=.libPaths()[1],
#                        CRAN=getOption("BIOC"))
#                    }
				else {
					directory <- paste("file:", tclvalue(directoryVariable), sep="")
					packages <- available.packages(contriburl=directory)[,1]
					present <- missing.packages %in% packages
					if (!all(present)) errorMessage()
					if (!any(present)) return()
					install.packages(missing.packages[present], contriburl=directory,
							lib=.libPaths()[1])
				}
			}
			onCancel <- function(){
				tkgrab.release(top)
				tkdestroy(top)
				return()
			}
			onHelp <- function() help("install.packages")
			buttonsFrame <- tkframe(top)
			OKbutton <- buttonRcmdr(buttonsFrame, text="OK", foreground="darkgreen", width="12", command=onOK, default="active",
					borderwidth=3)
			cancelButton <- buttonRcmdr(buttonsFrame, text=gettext("Cancel", domain="R-Rcmdr"), foreground="red", width="12", command=onCancel,
					borderwidth=3)
			helpButton <- buttonRcmdr(buttonsFrame, text=gettext("Help", domain="R-Rcmdr"), width="12", command=onHelp, borderwidth=3)
			tkgrid(OKbutton, labelRcmdr(buttonsFrame, text="  "), cancelButton,
					labelRcmdr(buttonsFrame, text="            "),
					helpButton, sticky="w")
			tkgrid(buttonsFrame, sticky="w")
			for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
			tkgrid.columnconfigure(top, 0, weight=0)
			.Tcl("update idletasks")
			tkwm.resizable(top, 0, 0)
			tkbind(top, "<Return>", onOK)
			tkwm.deiconify(top)
			tkgrab.set(top)
		#	tkfocus(top)
			tkwait.window(top)
		}
	}
}
