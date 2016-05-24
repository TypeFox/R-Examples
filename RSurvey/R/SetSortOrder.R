# A GUI for specifying the variable used to sort the data set.

SetSortOrder <- function(col.ids, sort.on=NULL, parent=NULL) {

  ## Additional functions

  # Save sort order
  SaveSortOrder <- function() {
    col.id <- as.character(tclvalue(col.id.var))
    decreasing <- as.logical(as.integer(tclvalue(decreasing.var)))
    na.last <- as.integer(tclvalue(na.last.var))
    if (na.last %in% 0:1)
      na.last <- as.logical(na.last)
    else
      na.last <- NA
    if (col.id == "") {
      sort.on <- NULL
    } else {
      sort.on <- which(col.ids == col.id)
      attr(sort.on, "decreasing") <- decreasing
      attr(sort.on, "na.last") <- na.last
    }
    rtn <<- sort.on
    tclvalue(tt.done.var) <- 1
  }

  ## Main program

  rtn <- sort.on

  # Assign variables linked to Tk widgets
  col.id.var     <- tclVar()
  decreasing.var <- tclVar(0)
  na.last.var    <- tclVar(1)
  tt.done.var    <- tclVar(0)

  # Set variables
  idx <- 0L
  if (is.integer(sort.on) && idx %in% seq_along(col.ids)) {
    idx <- as.integer(sort.on)
    decreasing <- attr(sort.on, "decreasing")
    if (!is.null(decreasing))
      tclvalue(decreasing.var) <- as.logical(decreasing)
    na.last <- attr(sort.on, "na.last")
    if (!is.null(na.last)) {
      if (is.logical(na.last))
        tclvalue(na.last.var) <- as.integer(na.last)
      if (is.na(na.last))
        tclvalue(na.last.var) <- 2
    }
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

  tktitle(tt) <- "Sort Order"

  tkwm.resizable(tt, 1, 0)

  # Frame 0

  frame0 <- tkframe(tt, relief="flat")

  frame0.but.2 <- ttkbutton(frame0, width=12, text="OK",
                            command=SaveSortOrder)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("SetSortOrder", package="RSurvey"))
                            })

  tkgrid("x", frame0.but.2, frame0.but.3, frame0.but.4,
         pady=c(15, 10), padx=c(4, 0))
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkgrid.configure(frame0.but.4, padx=c(4, 10))
  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1

  frame1 <- ttkframe(tt, relief="flat")

  frame1.lab.1.1 <- tklabel(frame1, text="Variable to sort on")

  vals <- c("", col.ids)
  if (length(vals) == 1)
    vals <- paste0("{", vals, "}")
  frame1.box.1.2 <- ttkcombobox(frame1, state="readonly",
                                textvariable=col.id.var, values=vals)
  tcl(frame1.box.1.2, "current", idx)

  frame1.lab.2.2 <- ttklabel(frame1, text="Order")
  frame1.rad.2.3 <- ttkradiobutton(frame1, variable=decreasing.var, value=FALSE,
                                   text="increasing", width=10)
  frame1.rad.3.3 <- ttkradiobutton(frame1, variable=decreasing.var, value=TRUE,
                                   text="decreasing", width=10)

  frame1.lab.2.4 <- ttklabel(frame1, text="NAs")
  frame1.rad.2.5 <- ttkradiobutton(frame1, variable=na.last.var, value=1,
                                   text="place last")
  frame1.rad.3.5 <- ttkradiobutton(frame1, variable=na.last.var, value=0,
                                   text="place first")
  frame1.rad.4.5 <- ttkradiobutton(frame1, variable=na.last.var, value=2,
                                   text="remove")

  tkgrid(frame1.lab.1.1, frame1.box.1.2, pady=c(15, 5))
  tkgrid("x", frame1.lab.2.2, frame1.rad.2.3, frame1.lab.2.4, frame1.rad.2.5,
         "x")
  tkgrid("x", "x", frame1.rad.3.3, "x", frame1.rad.3.5, "x")
  tkgrid("x", "x", "x", "x", frame1.rad.4.5, "x")

  tkgrid.configure(frame1.box.1.2, sticky="ew", columnspan=6)
  tkgrid.configure(frame1.lab.2.2, padx=c(0, 4))
  tkgrid.configure(frame1.lab.2.4, padx=c(20, 4))
  tkgrid.configure(frame1.rad.2.5, frame1.rad.3.5, frame1.rad.4.5, sticky="w")

  tkgrid.columnconfigure(frame1, 6, weight=1, minsize=0)

  tkpack(frame1, fill="x", padx=10)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  # GUI control

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  return(rtn)
}
