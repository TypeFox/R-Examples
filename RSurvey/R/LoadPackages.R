# This function loads R packages suggested by RSurvey. If a suggested
# package is unavailable on the local computer an attempt is made to
# acquire the package from CRAN using an existing network connection.

LoadPackages <- function() {

  ## Additional functions

  # Install missing packages from CRAN mirror
  InstallPackages <- function() {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- which(cran.mirrors$Name %in% as.character(tclvalue(repos.var)))
    repo <- cran.mirrors$URL[idx]
    contriburl <- contrib.url(repos=repo, type=getOption("pkgType"))
    cran.pkgs <- available.packages(contriburl)

    idxs <- as.integer(tkcurselection(frame1.lst.2.2)) + 1L
    missing.pkgs <- missing.pkgs[idxs]

    is.on.cran <- missing.pkgs %in% cran.pkgs
    available.pkgs <- missing.pkgs[is.on.cran]
    unavailable.pkgs <- missing.pkgs[!is.on.cran]
    if (length(unavailable.pkgs) > 0) {
      msg <- paste0("The following packages are unavailable on this ",
                    "CRAN mirror:\n\n",
                    paste(paste0("\'", unavailable.pkgs, "\'"), collapse=", "),
                    "\n\nWould you like to try a different CRAN mirror?")
      ans <- tkmessageBox(icon="question", message=msg, title="CRAN",
                          type="yesno", parent=tt)
      if (tolower(substr(ans, 1, 1)) == "y")
        return()
    }
    if (length(available.pkgs) > 0)
      install.packages(available.pkgs, repos=repo, verbose=TRUE)

    # Load packages into current session
    for (pkg in available.pkgs) {
      is.loaded <- suppressWarnings(require(pkg, character.only=TRUE,
                                            warn.conflicts=FALSE, quietly=TRUE))
      if (!is.loaded)
        warning(paste("unable to load suggested package:", pkg), call.=FALSE)
    }

    tclvalue(tt.done.var) <- 1
  }

  ## Main program

  # Suggested packages
  txt <- readLines(system.file("DESCRIPTION", package="RSurvey"))
  pkgs <- sub(",$", "", strsplit(txt[grep("^Suggests:", txt)], " ")[[1]][-1])

  # Account for missing packages

  is.pkg.missing <- !pkgs %in% .packages(all.available=TRUE)
  if (any(is.pkg.missing)) {
    missing.pkgs <- pkgs[is.pkg.missing]
    cran.mirrors <- getCRANmirrors(all=FALSE, local.only=FALSE)
    default.repo <- getOption("repos")
    idx <- which(sub("/$", "", cran.mirrors$URL) %in%
                 sub("/$", "", default.repo["CRAN"]))
    idx <- if (length(idx) > 0) idx[1] else 1
    repos.var <- tclVar(cran.mirrors$Name[idx])
    rlogo.var <- tclVar()
    tt.done.var <- tclVar(0)
    pkgs.var <- tclVar()
    for (i in seq_along(missing.pkgs))
      tcl("lappend", pkgs.var, missing.pkgs[i])

    # Open GUI
    tclServiceMode(FALSE)
    tt <- tktoplevel()
    tktitle(tt) <- "Manage Packages"
    tkwm.resizable(tt, 0, 0)

    # Frame 0, ok and cancel buttons
    frame0 <- tkframe(tt, relief="flat")
    frame0.but.2 <- ttkbutton(frame0, width=12, text="OK", default="active",
                              command=InstallPackages)
    frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                              command=function() tclvalue(tt.done.var) <- 1)
    tkgrid("x", frame0.but.2, frame0.but.3, sticky="se", pady=10)
    tkgrid.columnconfigure(frame0, 0, weight=1)
    tkgrid.configure(frame0.but.2, padx=c(10, 2))
    tkgrid.configure(frame0.but.3, padx=c(2, 10))
    tkpack(frame0, fill="x", side="bottom", anchor="e")

    # Frame 1, message and mirror selection
    frame1 <- tkframe(tt, relief="flat", background="white")
    if ("RSurvey" %in% .packages(all.available=TRUE))
      f <- system.file("images/rlogo.gif", package="RSurvey")
    else
      f <- file.path(getwd(), "inst", "images", "rlogo.gif")
    tkimage.create("photo", rlogo.var, format="GIF", file=f)

    frame1.lab.1.1 <- ttklabel(frame1, image=rlogo.var, background="white")

    txt <- "The following suggested (not required) packages are missing:"
    frame1.lab.1.2 <- ttklabel(frame1, text=txt, justify="left",
                               background="white")

    frame1.lst.2.2 <- tklistbox(frame1, selectmode="extended",
                                activestyle="none", relief="groove", height=4, width=25,
                                exportselection=FALSE, listvariable=pkgs.var,
                                highlightthickness=0)
    frame1.ysc.2.3 <- ttkscrollbar(frame1, orient="vertical")
    tkconfigure(frame1.lst.2.2, background="white",
                yscrollcommand=paste(.Tk.ID(frame1.ysc.2.3), "set"))
    tkconfigure(frame1.ysc.2.3, command=paste(.Tk.ID(frame1.lst.2.2), "yview"))
    tkselection.set(frame1.lst.2.2, 0, length(missing.pkgs))

    txt <- paste("Some features will not be available without these packages.",
                 "By default all missing packages are selected for installation.",
                 sep="\n")
    frame1.lab.3.2 <- ttklabel(frame1, text=txt, justify="left",
                               background="white")

    txt <- paste("The selected packages will be installed from the",
                 "Comprehensive R Archive Network (CRAN).", sep="\n")
    frame1.lab.4.2 <- ttklabel(frame1, text=txt, justify="left",
                               background="white")

    frame1.lab.5.2 <- ttklabel(frame1, text="CRAN mirror",
                               justify="left", background="white")
    frame1.box.5.3 <- ttkcombobox(frame1, textvariable=repos.var, width=25,
                                  values=cran.mirrors$Name, state="readonly")
    tcl(frame1.box.5.3, "current", 0)

    tkgrid(frame1.lab.1.1, frame1.lab.1.2, "x", pady=c(30, 8))
    tkgrid("x", frame1.lst.2.2, frame1.ysc.2.3, pady=c(0, 8))
    tkgrid("x", frame1.lab.3.2, "x", pady=c(0, 8))
    tkgrid("x", frame1.lab.4.2, "x", pady=c(0, 8))
    tkgrid("x", frame1.lab.5.2, frame1.box.5.3, pady=c(0, 30))

    tkgrid.configure(frame1.lab.1.1, padx=c(40, 20), rowspan=4, sticky="n")
    tkgrid.configure(frame1.lab.1.2, padx=c(0, 40), columnspan=2, sticky="w")
    tkgrid.configure(frame1.lst.2.2, padx=c(0, 0), sticky="e")
    tkgrid.configure(frame1.ysc.2.3, padx=c(0, 40), sticky="nsw")
    tkgrid.configure(frame1.lab.3.2, padx=c(0, 40), columnspan=2, sticky="w")
    tkgrid.configure(frame1.lab.4.2, padx=c(0, 40), columnspan=2, sticky="w")
    tkgrid.configure(frame1.lab.5.2, padx=c( 0,  4), sticky="e")
    tkgrid.configure(frame1.box.5.3, padx=c( 0, 40), sticky="w")

    tkpack(frame1)

    # Binds events
    tclServiceMode(TRUE)
    tkbind(tt, "<Return>", InstallPackages)
    tkbind(tt, "<Key-space>", InstallPackages)

    # GUI control
    tkfocus(tt)
    tkgrab(tt)
    tkwait.variable(tt.done.var)
    tclServiceMode(FALSE)
    tkgrab.release(tt)
    tkdestroy(tt)
    tclServiceMode(TRUE)
  }

  # Warn if Tktable is unavailable
  tcl.pkg <- tryCatch(tcl("package", "require", "Tktable"), error=identity)
  if (inherits(tcl.pkg, "error")) {
    msg <- paste("Tcl package Tktable is missing and is strongly recommended",
                 "for full functionality of RSurvey.\n\n ",
                 "http://tktable.sourceforge.net")
    tkmessageBox(icon="warning", message=msg, title="Missing Tktable",
                 type="ok")
  }
}
