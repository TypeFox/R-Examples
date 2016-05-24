# A GUI for establishing find and replace arguments.

Search <- function(is.replace=FALSE, defaults=NULL, parent=NULL) {

  ## Additional functions

  # Return find and replace parameters
  ReturnParameters <- function(is.replace.first=FALSE) {
    find.what <- as.character(tclvalue(find.what.var))
    if (is.replace)
      replace.with <- as.character(tclvalue(replace.with.var))
    else
      replace.with <- NULL
    is.match.word <- as.logical(as.integer(tclvalue(match.word.var)))
    is.match.case <- as.logical(as.integer(tclvalue(match.case.var)))
    is.reg.exps   <- as.logical(as.integer(tclvalue(reg.exps.var)))
    is.search.col <- as.logical(as.integer(tclvalue(search.col.var)))
    is.perl       <- as.logical(as.integer(tclvalue(perl.var)))
    is.search.sel <- as.logical(as.integer(tclvalue(search.sel.var)))
    rtn <<- list(find.what=find.what, replace.with=replace.with,
                 is.match.word=is.match.word, is.match.case=is.match.case,
                 is.reg.exps=is.reg.exps, is.search.col=is.search.col,
                 is.perl=is.perl, is.replace.first=is.replace.first,
                 is.search.sel=is.search.sel)
    tclvalue(tt.done.var) <- 1
  }

  # Toggle match word
  ToggleMatchWord <- function() {
    is.match.word <- as.logical(as.integer(tclvalue(match.word.var)))
    if (is.match.word) {
      tclvalue(reg.exps.var) <- FALSE
      tkconfigure(frame3.chk.4.1, state="disabled")
    } else {
      tkconfigure(frame3.chk.4.1, state="normal")
    }
    ToggleRegExps()
  }

  # Toggle regular expression
  ToggleRegExps <- function() {
    is.reg.exps <- as.logical(as.integer(tclvalue(reg.exps.var)))
    if (is.reg.exps) {
      tkconfigure(frame3.rad.4.2, state="normal")
      tkconfigure(frame3.rad.4.3, state="normal")
    } else {
      tkconfigure(frame3.rad.4.2, state="disabled")
      tkconfigure(frame3.rad.4.3, state="disabled")
    }
  }

  ## Main program

  # Assigin global variables
  rtn <- NULL

  # Assign variables linked to Tk widgets
  find.what.var     <- tclVar()
  replace.with.var  <- tclVar()
  match.word.var    <- tclVar(0)
  match.case.var    <- tclVar(0)
  reg.exps.var      <- tclVar(0)
  search.col.var    <- tclVar(0)
  search.sel.var    <- tclVar(0)
  perl.var          <- tclVar(0)
  tt.done.var       <- tclVar(0)

  # Set default values
  replace.with <- ""
  if (!is.null(defaults) && is.list(defaults)) {
    if (!is.null(defaults$find.what) && is.character(defaults$find.what))
      tclvalue(find.what.var) <- defaults$find.what
    if (!is.null(defaults$replace.with) && is.character(defaults$replace.with))
      tclvalue(replace.with.var) <- defaults$replace.with
    if (!is.null(defaults$is.match.word) && is.logical(defaults$is.match.word))
      tclvalue(match.word.var) <- defaults$is.match.word
    if (!is.null(defaults$is.match.case) && is.logical(defaults$is.match.case))
      tclvalue(match.case.var) <- defaults$is.match.case
    if (!is.null(defaults$is.reg.exps) && is.logical(defaults$is.reg.exps))
      tclvalue(reg.exps.var) <- defaults$is.reg.exps
    if (!is.null(defaults$is.search.col) && is.logical(defaults$is.search.col))
      tclvalue(search.col.var) <- defaults$is.search.col
    if (!is.null(defaults$is.perl) && is.logical(defaults$is.perl))
      tclvalue(perl.var) <- defaults$is.perl
    if (!is.null(defaults$is.search.sel) && is.logical(defaults$is.search.sel))
      tclvalue(search.sel.var) <- defaults$is.search.sel
  }

  # Open GUI
  tclServiceMode(FALSE)
  tt <- tktoplevel()
  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- "Search"
  tkwm.resizable(tt, 1, 0)

  # Frame 0

  frame0 <- ttkframe(tt, relief="flat")

  if (is.replace) {
    frame0.but.1.2 <- ttkbutton(frame0, width=12, text="Replace",
                                command=function() ReturnParameters(TRUE))
    frame0.but.1.3 <- ttkbutton(frame0, width=12, text="Replace All",
                                command=function() ReturnParameters(FALSE))
    frame0.but.1.4 <- ttkbutton(frame0, width=12, text="Cancel",
                                command=function() tclvalue(tt.done.var) <- 1)
  } else {
    frame0.but.1.2 <- "x"
    frame0.but.1.3 <- ttkbutton(frame0, width=12, text="Find",
                                command=function() ReturnParameters())
    frame0.but.1.4 <- ttkbutton(frame0, width=12, text="Cancel",
                                command=function() tclvalue(tt.done.var) <- 1)
  }
  frame0.but.1.5 <- ttkbutton(frame0, width=12, text="Help",
                              command=function() {
                                print(help("Search", package="RSurvey"))
                              })

  tkgrid("x", frame0.but.1.2, frame0.but.1.3, frame0.but.1.4, frame0.but.1.5)

  tkgrid.columnconfigure(frame0, 0, weight=1)

  tkgrid.configure(frame0.but.1.3, frame0.but.1.4, padx=c(4, 0))
  if (is.replace)
    tkgrid.configure(frame0.but.1.2, padx=c(10, 0))
  else
    tkgrid.configure(frame0.but.1.3, padx=c(10, 0))
  tkgrid.configure(frame0.but.1.5, padx=c(4, 10), pady=10)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1

  frame1 <- ttkframe(tt, relief="flat")

  frame1.lab.1.1 <- ttklabel(frame1, text="Find what:", foreground="#141414")
  frame1.ent.2.1 <- ttkentry(frame1, width=10, font="TkFixedFont",
                             textvariable=find.what.var)
  tkgrid(frame1.lab.1.1, sticky="w")
  tkgrid(frame1.ent.2.1, sticky="we")
  tkgrid.columnconfigure(frame1, 0, weight=1)

  if (is.replace) {
    frame1.lab.3.1 <- ttklabel(frame1, text="Replace with:",
                               foreground="#141414")
    frame1.ent.4.1 <- ttkentry(frame1, width=10, font="TkFixedFont",
                               textvariable=replace.with.var)
    tkgrid(frame1.lab.3.1, sticky="w", pady=c(10, 0))
    tkgrid(frame1.ent.4.1, sticky="we")
    tkgrid.rowconfigure(frame1, 3, weight=1)
  }

  tkpack(frame1, fill="x", expand="yes", padx=10, pady=10)

  # Frame 3

  frame3 <- ttkframe(tt, relief="flat")

  frame3.chk.1.1 <- ttkcheckbutton(frame3, text="Match case",
                                   variable=match.case.var)
  frame3.chk.2.1 <- ttkcheckbutton(frame3, text="Match entire cell contents",
                                   variable=match.word.var,
                                   command=function() ToggleMatchWord())
  frame3.chk.3.1 <- ttkcheckbutton(frame3, text="Search in selected cells",
                                   variable=search.sel.var)
  frame3.chk.4.1 <- ttkcheckbutton(frame3, text="Search using regular expressions",
                                   variable=reg.exps.var,
                                   command=function() ToggleRegExps())
  frame3.rad.4.2 <- ttkradiobutton(frame3, variable=perl.var, value=FALSE,
                                   text="Unix")
  frame3.rad.4.3 <- ttkradiobutton(frame3, variable=perl.var, value=TRUE,
                                   text="Perl")

  tkgrid(frame3.chk.1.1, "x", "x", sticky="w")
  tkgrid(frame3.chk.2.1, "x", "x", sticky="w", pady=2)

  tkgrid(frame3.chk.3.1, "x", "x", sticky="w")
  tkgrid(frame3.chk.4.1, frame3.rad.4.2, frame3.rad.4.3, sticky="w", pady=c(2, 10))

  tkgrid.configure(frame3.rad.4.2, padx=c(4, 4))

  tkpack(frame3, fill="x", padx=10)

  # GUI control

  tclServiceMode(TRUE)

  ToggleRegExps()

  tkfocus(frame1.ent.2.1)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(rtn)
}
