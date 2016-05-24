# A GUI for renaming values in a vector of character strings.

Rename <- function(names=NULL, cur.name=NULL, win.title=NULL, parent=NULL) {

  ## Additional functions

  # Update entry
  UpdateEntry <- function() {
    if (tclvalue(cur.var) != "" && !(tclvalue(new.var) %in% new.names))
      new.names[names %in% tclvalue(cur.var)] <<- tclvalue(new.var)
    tclvalue(new.var) <- new.names[names %in% tclvalue(old.var)]
    tclvalue(cur.var) <- tclvalue(old.var)
  }

  # Save renamed values
  SaveNames <- function() {
    UpdateEntry()
    rtn.names <<- new.names
    tclvalue(tt.done.var) <- 1
  }

  ## Main program

  if (is.null(names))
    return(NULL)

  rtn.names <- new.names <- names

  # Assign the variables linked to Tk widgets

  old.var <- tclVar("")
  new.var <- tclVar("")
  cur.var <- tclVar("")

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

  if (!is.null(win.title))
    tktitle(tt) <- win.title

  tkwm.resizable(tt, 1, 0)

  # Frame 0

  frame0 <- ttkframe(tt, relief="flat")

  frame0.but.2 <- ttkbutton(frame0, width=12, text="OK",
                            command=SaveNames)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("Rename", package="RSurvey"))
                            })

  tkgrid("x", frame0.but.2, frame0.but.3, frame0.but.4,
         sticky="se", pady=c(15, 10), padx=c(4, 0))
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkgrid.configure(frame0.but.2, padx=c(40, 0))
  tkgrid.configure(frame0.but.4, padx=c(4, 10))
  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1

  frame1 <- ttkframe(tt, relief="flat")

  frame1.lab.1 <- ttklabel(frame1, text="Old name")
  frame1.lab.2 <- ttklabel(frame1, text="New name")

  if (length(names) == 1)
    prep.names <- paste0("{", names, "}")
  else
    prep.names <- names

  frame1.box.1 <- ttkcombobox(frame1, state="readonly", values=prep.names,
                              textvariable=old.var)

  frame1.ent.1 <- ttkentry(frame1, textvariable=new.var)

  if (!is.null(cur.name) && cur.name %in% names)
    tcl(frame1.box.1, "current", match(cur.name, names) - 1)

  tkgrid(frame1.lab.1, frame1.box.1)
  tkgrid(frame1.lab.2, frame1.ent.1, pady=c(10, 0))

  tkgrid.configure(frame1.lab.1, frame1.lab.2, sticky="w", padx=c(0, 2))
  tkgrid.configure(frame1.box.1, frame1.ent.1, sticky="we")

  tkgrid.columnconfigure(frame1, 1, weight=1, minsize=25)

  tkpack(frame1, fill="x", expand=TRUE, padx=10, pady=c(10, 0))

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame1.box.1, "<<ComboboxSelected>>", UpdateEntry)

  # GUI control

  UpdateEntry()

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(rtn.names)
}
