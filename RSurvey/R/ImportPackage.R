# A GUI for importing data sets from R packages.

ImportPackage <- function(classes=NULL, parent=NULL) {

  ## Additional functions

  # Load data set
  LoadDataset <- function() {
    tkconfigure(tt, cursor="watch")
    tclServiceMode(FALSE)
    on.exit(tkconfigure(tt, cursor="arrow"))
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- as.integer(tkcurselection(frame1.lst.2.1)) + 1L
    pkg.name <- pkg.names[idx]
    idx <- as.integer(tkcurselection(frame1.lst.2.4))
    pkg.item <- paste(as.character(tkget(frame1.lst.2.4, idx)), collapse=" ")
    e <- environment(LoadDataset)
    txt <- paste0("data(", pkg.item, ", package=\"", pkg.name, "\", envir=e)")
    ds.name <- eval(parse(text=txt))
    rtn <<- list(d=eval(parse(text=paste("(", ds.name, ")")), envir=e),
                 src=c(dataset=ds.name, package=pkg.name,
                       accessed=format(Sys.time())))
    tclvalue(tt.done.var) <- 1
  }

  # Check if package(s) are loaded
  IsPackageLoaded <- function(pkg.names) {
    vapply(pkg.names, function(i) paste("package", i, sep=":") %in% search(),
           TRUE)
  }

  # Load package
  LoadPackage <- function(pkg.name) {
    tkconfigure(tt, cursor="watch")
    tclServiceMode(FALSE)
    on.exit(tkconfigure(tt, cursor="arrow"))
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- as.integer(tkcurselection(frame1.lst.2.1)) + 1L
    pkg.name <- pkg.names[idx]
    lib <- paste("package", pkg.name, sep=":")
    if (!lib %in% search())
      suppressPackageStartupMessages(require(pkg.name, quietly=TRUE,
                                             warn.conflicts=FALSE,
                                             character.only=TRUE))
    if (lib %in% search()) {
      idx <- as.integer(tcl(frame1.box.3.1, "current"))
      if (idx == 2L) {
        tcl(frame1.box.3.1, "current", 0)
        SelectPackageType()
        tkselection.clear(frame1.lst.2.1, 0, "end")
        idx <- which(pkg.names %in% pkg.name) - 1L
        tkselection.set(frame1.lst.2.1, idx)
        tksee(frame1.lst.2.1, idx)
      }
      SelectPackage()
    } else {
      tkmessageBox(icon="error", message="Unable to load package.",
                   title="Error", type="ok", parent=tt)
    }
    tkfocus(frame1.lst.2.1)
  }

  # Summarize data sets in package
  SummarizePackage <- function() {
    idx <- as.integer(tkcurselection(frame1.lst.2.1)) + 1L
    pkg.name <- pkg.names[idx]
    pkg.datasets <- ds.list[[pkg.name]]
    pkg.datasets <- pkg.datasets[order(pkg.datasets[, 1]), , drop=FALSE]
    nmax <- max(nchar(pkg.datasets[, "Item"]))
    if (nmax < 20)
      nmax <- 20
    items <- sprintf(paste0("%-", nmax, "s"), pkg.datasets[, "Item"])
    txt <- paste0("Data sets in package ", sQuote(pkg.name), ":\n")
    txt <- c(txt, paste(items, pkg.datasets[, "Title"], sep="  "))
    EditText(txt, read.only=TRUE, win.title="Summary of Data Sets",
             is.fixed.width.font=TRUE, parent=tt)
    tkfocus(frame1.lst.2.1)
  }

  # Describe data set
  DescribeDataset <- function() {
    idx <- as.integer(tkcurselection(frame1.lst.2.1)) + 1L
    pkg.name <- pkg.names[idx]
    idx <- as.integer(tkcurselection(frame1.lst.2.4))
    pkg.item <- paste(as.character(tkget(frame1.lst.2.4, idx)), collapse=" ")
    ans <- try(help(pkg.item, package=(pkg.name)), silent=TRUE)
    if (inherits(ans, "try-error"))
      tkmessageBox(icon="error", message="Problem with dataset documentation.",
                   title="Error", type="ok", parent=tt)
    else
      print(ans)
    tkfocus(frame1.lst.2.4)
  }

  # Get class of data set
  GetClass <- function(pkg.name, pkg.item) {
    e <- environment(GetClass)
    txt <- paste0("data(", pkg.item, ", package=\"", pkg.name, "\", envir=e)")
    suppressWarnings(eval(parse(text=txt)))
    txt <- paste("(", pkg.item, ")")
    pkg.item.class <- class(try(eval(parse(text=txt), envir=e), silent=TRUE))
    return(pkg.item.class[1])
  }

  # GUI control for select package
  SelectPackage <- function() {
    tkconfigure(tt, cursor="watch")
    tclServiceMode(FALSE)
    on.exit(tkconfigure(tt, cursor="arrow"))
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- as.integer(tkcurselection(frame1.lst.2.1)) + 1L
    pkg.name <- pkg.names[idx]
    if (IsPackageLoaded(pkg.name)) {
      tkconfigure(frame1.but.4.2, state="disabled", default="disabled")
      pkg.datasets <- ds.list[[pkg.name]]
      all.pkg.items <- pkg.datasets[, "Item"]
      if (is.null(ds.class[[pkg.name]]))
        ds.class[[pkg.name]] <<- vapply(all.pkg.items,
                                        function(i) GetClass(pkg.name, i), "")

      ds.classes <- vapply(ds.class[[pkg.name]],
                           function(i) paste(i, collapse=" "), "")
      ds.in.error <- ds.classes %in% "try-error"
      ds.class.vals <- sort(unique(ds.classes[!ds.in.error]))

      if (as.logical(as.integer(tclvalue(fit.for.loading.var))))
        ds.class.vals <- ds.class.vals[ds.class.vals %in% classes]

      if (length(ds.class.vals) > 0)
        ds.class.vals <- c("Show all classes", ds.class.vals)
      else
        ds.class.vals <- "{Show all classes}"

      old.class <- paste(as.character(tkget(frame1.box.3.4)), collapse=" ")
      tkconfigure(frame1.box.3.4, value=ds.class.vals, state="readonly")

      idx <- which(ds.class.vals %in% old.class)
      idx <- if (length(idx) > 0) idx - 1L else 0L
      tcl(frame1.box.3.4, "current", idx)

      if (idx > 0)
        pkg.items <- all.pkg.items[ds.classes %in% ds.class.vals[idx + 1L]]
      else
        pkg.items <- all.pkg.items[ds.classes %in% ds.class.vals]

      if (length(pkg.items) > 0)
        pkg.items <- sort(pkg.items)
      tkselection.clear(frame1.lst.2.4, 0, "end")
      tclvalue(dataset.var) <- ""
      for (i in seq_along(pkg.items))
        tcl("lappend", dataset.var, pkg.items[i])
      if (length(pkg.items) > 0) {
        tkselection.set(frame1.lst.2.4, 0)
        tkconfigure(frame1.but.4.4, state="normal")
      } else {
        tkconfigure(frame1.but.4.4, state="disabled")
      }
    } else {
      tkconfigure(frame1.box.3.4, value="{}")
      tcl(frame1.box.3.4, "current", 0)
      tkconfigure(frame1.box.3.4, state="disabled")

      tkconfigure(frame0.but.1.2, state="disabled")
      tkconfigure(frame1.but.4.2, state="normal", default="active")
      tkconfigure(frame1.but.4.4, state="disabled")
      tkselection.clear(frame1.lst.2.4, 0, "end")
      tclvalue(dataset.var) <- ""
    }
    SelectDataset()
  }

  # GUI control for select data set
  SelectDataset <- function() {
    idx <- as.integer(tkcurselection(frame1.lst.2.4))
    if (length(idx) == 0) {
      tkconfigure(frame0.but.1.2, state="disabled", default="disabled")
      return()
    }
    pkg.item <- paste(as.character(tkget(frame1.lst.2.4, idx)), collapse=" ")
    pkg.name <- pkg.names[as.integer(tkcurselection(frame1.lst.2.1)) + 1L]
    idx <- which(ds.list[[pkg.name]][, "Item"] %in% pkg.item)
    is.valid <- ds.class[[pkg.name]][idx] %in% classes
    if (is.null(classes) | is.valid)
      tkconfigure(frame0.but.1.2, state="normal", default="active")
    else
      tkconfigure(frame0.but.1.2, state="disabled", default="disabled")
  }

  # GUI control for select package type
  SelectPackageType <- function() {
    tkconfigure(tt, cursor="watch")
    tclServiceMode(FALSE)
    on.exit(tkconfigure(tt, cursor="arrow"))
    on.exit(tclServiceMode(TRUE), add=TRUE)
    idx <- as.integer(tcl(frame1.box.3.1, "current"))
    if (idx == 0L) {
      pkg.names <<- all.pkgs
    } else {
      is.pkg.loaded <- IsPackageLoaded(all.pkgs)
      if (idx == 1L)
        pkg.names <<- all.pkgs[is.pkg.loaded]
      else
        pkg.names <<- all.pkgs[!is.pkg.loaded]
    }
    tkselection.clear(frame1.lst.2.1, 0, "end")
    tclvalue(package.var) <- ""
    for (i in seq_along(pkg.names))
      tcl("lappend", package.var, pkg.names[i])
    tkselection.set(frame1.lst.2.1, 0)
    SelectPackage()
    tkfocus(frame1.lst.2.1)
  }

  # GUI control for select class type
  SelectClassType <- function() {
    SelectPackage()
    if (length(as.integer(tkcurselection(frame1.lst.2.4))) == 0)
      tkfocus(frame1.lst.2.1)
    else
      tkfocus(frame1.lst.2.4)
  }

  ## Main program

  # Initialize values

  all.pkgs <- .packages(all.available=TRUE, lib.loc=.libPaths())
  all.pkgs <- all.pkgs[!all.pkgs %in% c("Rcmdr")]

  all.ds <- suppressWarnings(data(package=all.pkgs)$results)
  all.pkgs <- sort(unique(all.ds[, "Package"]))

  Fun <- function(i) all.ds[all.ds[, "Package"] == i, c("Item", "Title"),
                            drop=FALSE]
  ds.list <- sapply(all.pkgs, Fun, simplify=FALSE)

  ds.class <- list()

  pkg.type.vals <- paste(c("Show all", "Loaded", "Unloaded"), "packages")
  ds.class.vals <- "{}"

  pkg.names <- NULL
  rtn <- NULL

  # Assign variables linked to Tk widgets

  package.var <- tclVar()
  dataset.var <- tclVar()

  fit.for.loading.var <- tclVar(1)

  tt.done.var <- tclVar(0)

  # Open GUI

  tclServiceMode(FALSE)

  tt <- tktoplevel()

  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }

  tktitle(tt) <- "Import Data From Package"

  # Frame 0 contains load, cancel, and help buttons, and size grip

  frame0 <- ttkframe(tt, relief="flat")

  frame0.but.1.2 <- ttkbutton(frame0, width=12, text="Load",
                              command=LoadDataset)
  frame0.but.1.3 <- ttkbutton(frame0, width=12, text="Cancel",
                              command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.1.4 <- ttkbutton(frame0, width=12, text="Help",
                              command=function() {
                                print(help("ImportPackage",
                                           package="RSurvey"))
                              })
  frame0.grp.1.5 <- ttksizegrip(frame0)

  tkgrid("x", frame0.but.1.2, frame0.but.1.3, frame0.but.1.4, frame0.grp.1.5)

  tkgrid.configure(frame0.but.1.2, padx=c(10, 0))
  tkgrid.configure(frame0.but.1.3, padx=4)
  tkgrid.configure(frame0.but.1.4, padx=c(0, 10), columnspan=2)

  tkgrid.configure(frame0.but.1.2, frame0.but.1.3, frame0.but.1.4,
                   pady=c(15, 10))

  tkgrid.configure(frame0.grp.1.5, sticky="se")
  tkraise(frame0.but.1.4, frame0.grp.1.5)
  tkgrid.columnconfigure(frame0, 0, weight=1)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  tkconfigure(frame0.but.1.2, state="disabled")

  # Frame 1, package and dataset

  frame1 <- ttkframe(tt, relief="flat", padding=0, borderwidth=0)

  frame1.lab.1.1 <- ttklabel(frame1, text="Package", foreground="#141414")
  frame1.lab.1.4 <- ttklabel(frame1, text="Data set", foreground="#141414")

  frame1.lst.2.1 <- tklistbox(frame1, selectmode="browse", activestyle="none",
                              relief="flat", borderwidth=5, width=30, height=8,
                              exportselection=FALSE, listvariable=package.var,
                              highlightthickness=0)
  frame1.ysc.2.3 <- ttkscrollbar(frame1, orient="vertical")
  tkconfigure(frame1.lst.2.1, background="white",
              yscrollcommand=paste(.Tk.ID(frame1.ysc.2.3), "set"))
  tkconfigure(frame1.ysc.2.3, command=paste(.Tk.ID(frame1.lst.2.1), "yview"))

  frame1.lst.2.4 <- tklistbox(frame1, selectmode="browse", activestyle="none",
                              relief="flat", borderwidth=5, width=30, height=8,
                              exportselection=FALSE, listvariable=dataset.var,
                              highlightthickness=0)
  frame1.ysc.2.6 <- ttkscrollbar(frame1, orient="vertical")
  tkconfigure(frame1.lst.2.4, background="white",
              yscrollcommand=paste(.Tk.ID(frame1.ysc.2.6), "set"))
  tkconfigure(frame1.ysc.2.6, command=paste(.Tk.ID(frame1.lst.2.4), "yview"))

  frame1.box.3.1 <- ttkcombobox(frame1, state="readonly", value=pkg.type.vals)
  frame1.box.3.4 <- ttkcombobox(frame1, state="readonly", value=ds.class.vals)

  frame1.but.4.1 <- ttkbutton(frame1, width=10, text="Summary",
                              command=SummarizePackage)
  frame1.but.4.2 <- ttkbutton(frame1, width=10, text="Load",
                              command=LoadPackage)
  frame1.but.4.4 <- ttkbutton(frame1, width=10, text="Describe",
                              command=DescribeDataset)

  if (is.null(classes))
    frame1.chk.4.5 <- "x"
  else
    frame1.chk.4.5 <- ttkcheckbutton(frame1, text="Fit for loading",
                                     variable=fit.for.loading.var,
                                     command=function() {
                                               SelectPackage()
                                               tkfocus(frame1.lst.2.1)
                                             })

  tkgrid(frame1.lab.1.1, "x", "x", frame1.lab.1.4, "x", "x", pady=c(10, 0))
  tkgrid(frame1.lst.2.1, "x", frame1.ysc.2.3, frame1.lst.2.4, "x",
         frame1.ysc.2.6)
  tkgrid(frame1.box.3.1, "x", "x", frame1.box.3.4, "x", "x", pady=4)
  tkgrid(frame1.but.4.1, frame1.but.4.2, "x", frame1.but.4.4, frame1.chk.4.5,
         "x")

  tkgrid.configure(frame1.lab.1.1, columnspan=3)
  tkgrid.configure(frame1.lab.1.4, columnspan=3)
  tkgrid.configure(frame1.lst.2.1, columnspan=2)
  tkgrid.configure(frame1.lst.2.4, columnspan=2)
  tkgrid.configure(frame1.box.3.1, columnspan=2)
  tkgrid.configure(frame1.box.3.4, columnspan=2)

  tkgrid.configure(frame1.lab.1.1, frame1.lab.1.4, sticky="w")
  tkgrid.configure(frame1.lst.2.1, frame1.lst.2.4, sticky="nswe")
  tkgrid.configure(frame1.ysc.2.3, frame1.ysc.2.6, sticky="ns")
  tkgrid.configure(frame1.box.3.1, frame1.box.3.4, sticky="we")
  tkgrid.configure(frame1.but.4.2, frame1.but.4.4, frame1.chk.4.5, sticky="w")

  tkgrid.configure(frame1.ysc.2.3, padx=c(0, 25))
  tkgrid.configure(frame1.but.4.1, padx=c(0, 4))
  tkgrid.configure(frame1.but.4.4, padx=c(0, 8))

  tkgrid.columnconfigure(frame1, 1, minsize=85, weight=1)
  tkgrid.columnconfigure(frame1, 4, minsize=85, weight=1)
  tkgrid.rowconfigure(frame1, 1, weight=1)

  tkpack(frame1, fill="both", expand=TRUE, anchor="nw", padx=10)

  tkselection.set(frame1.lst.2.1, 0)
  tcl(frame1.box.3.1, "current", 0)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(frame1.lst.2.1, "<<ListboxSelect>>", SelectPackage)
  tkbind(frame1.lst.2.4, "<<ListboxSelect>>", SelectDataset)

  tkbind(frame1.box.3.1, "<<ComboboxSelected>>", SelectPackageType)
  tkbind(frame1.box.3.4, "<<ComboboxSelected>>", SelectClassType)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  # GUI control

  SelectPackageType()

  tkfocus(frame1.lst.2.1)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  invisible(rtn)
}
