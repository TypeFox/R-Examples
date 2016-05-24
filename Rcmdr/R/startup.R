# last modified 2016-02-22 by J. Fox
#  applied patch to improve window behaviour supplied by Milan Bouchet-Valat 2011-09-22

.onAttach <- function(...){
    if (interactive()){
        if (!(exists(".RcmdrEnv") && is.environment(RcmdrEnv()) &&
                exists("commanderWindow", RcmdrEnv()) &&
                !is.null(get("commanderWindow", RcmdrEnv())))){
            Commander()
            packageStartupMessage(gettext("\nRcmdr Version", domain="R-Rcmdr"), " ",
                getRcmdr("RcmdrVersion"), "\n")
            return()
        }
    }
    else {
        packageStartupMessage("The Commander GUI is launched only in interactive sessions",
            domain="R-Rcmdr")
        return()
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
    if (MacOSXP()){
        #       PATH <- system2("/usr/libexec/path_helper", "-s", stdout=TRUE)
        #       PATH <- sub("\"; export PATH;$", "", sub("^PATH=\\\"", "", PATH))
        #       Sys.setenv(PATH=PATH)
        PATH <- Sys.getenv("PATH")
        putRcmdr("PATH", PATH)
        PATH <- unlist(strsplit(PATH, .Platform$path.sep, fixed=TRUE))
        if (MacOSXP("15.0.0")){
            if (length(grep("^/Library/TeX/texbin$", PATH)) == 0) {
                PATH[length(PATH) + 1] <- "/Library/TeX/texbin"
                Sys.setenv(PATH=paste(PATH, collapse=.Platform$path.sep))
            }
        }
        else{
            if (length(grep("^/usr/texbin$", PATH)) == 0) {
                PATH[length(PATH) + 1] <- "/usr/texbin"
                Sys.setenv(PATH=paste(PATH, collapse=.Platform$path.sep))
            }
        }
    }
    required.packages <- rev(c("abind", "aplpack", "car", "colorspace", 
        "effects", "e1071", "foreign", "Hmisc", "knitr", "lattice", "leaps", "lmtest",
        "markdown", "MASS", "mgcv", "multcomp", "nlme", "nnet", "RcmdrMisc", "relimp", "rgl",
        "rglwidget", "rmarkdown", "sandwich", "sem", "XLConnect"))
#    if (WindowsP()) required.packages <- c(required.packages, "RODBC")
    check <- options("Rcmdr")[[1]]$check.packages
    if (length(check) > 0 && !check) return()
    packages.to.check <- required.packages
    available.packages <- packagesAvailable(packages.to.check)
    missing.packages <- packages.to.check[!available.packages]
    if (any(!available.packages)) {
        response <- tkmessageBox(message=paste(gettext("The following packages used by Rcmdr are missing:\n", domain="R-Rcmdr"),
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

.onUnload <- function(...){
  if (MacOSXP()){
    PATH <- getRcmdr("PATH", fail=FALSE)
    if (!is.null(PATH)) Sys.setenv(PATH=PATH)
  }
}
