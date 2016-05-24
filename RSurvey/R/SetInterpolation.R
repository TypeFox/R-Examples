# A GUI for specifying the interpolation parameters.

SetInterpolation <- function(parent=NULL) {

  ## Additional functions

  # Update parameter values

  UpdatePar <- function() {
    vars <- c("grid.res", "grid.mba")

    old <- sapply(vars, function(i) Data(i))

    grid.res <- list(x=as.numeric(tclvalue(grid.dx.var)),
                     y=as.numeric(tclvalue(grid.dy.var)))
    if (all(vapply(grid.res, function(i) is.na(i), TRUE)))
      Data("grid.res", NULL)
    else
      Data("grid.res", grid.res)

    grid.mba <- list(n=as.integer(tclvalue(mba.n.var)),
                     m=as.integer(tclvalue(mba.m.var)),
                     h=as.integer(tclvalue(mba.h.var)))
    if (all(vapply(grid.mba, function(i) is.na(i), TRUE)))
      Data("grid.mba", NULL)
    else
      Data("grid.mba", grid.mba)

    new <- sapply(vars, function(i) Data(i))

    if (!identical(old, new))
      Data("data.grd", NULL)

    tclvalue(tt.done.var) <- 1
  }

  ## Main program

  # Assign the variables linked to Tk widgets

  grid.dx.var <- tclVar()
  grid.dy.var <- tclVar()
  mba.n.var   <- tclVar()
  mba.m.var   <- tclVar()
  mba.h.var   <- tclVar()

  grid.res <- Data("grid.res")
  if (!is.na(grid.res$x))
    tclvalue(grid.dx.var) <- as.numeric(grid.res$x)
  if (!is.na(grid.res$y))
    tclvalue(grid.dy.var) <- as.numeric(grid.res$y)

  grid.mba <- Data("grid.mba")
  if (!is.na(grid.mba$n))
    tclvalue(mba.n.var) <- as.integer(grid.mba$n)
  if (!is.na(grid.mba$m))
    tclvalue(mba.m.var) <- as.integer(grid.mba$m)
  if (!is.na(grid.mba$h))
    tclvalue(mba.h.var) <- as.integer(grid.mba$h)

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
  tktitle(tt) <- "Set Interpolation Method"

  tkwm.resizable(tt, 1, 0)

  # Frame 0, ok and cancel buttons

  frame0 <- tkframe(tt, relief="flat")

  frame0.but.2 <- ttkbutton(frame0, width=12, text="OK",
                            command=UpdatePar)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)

  frame0.but.4 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("SetInterpolation", package="RSurvey"))
                            })
  tkgrid("x", frame0.but.2, frame0.but.3, frame0.but.4,
         pady=c(15, 10), padx=c(4, 0))
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkgrid.configure(frame0.but.4, padx=c(4, 10))
  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1, interpolation parameteres

  frame1 <- ttkframe(tt, relief="flat")

  txt <- "Interpolated-grid spacing along the x-axis"
  frame1.lab.1.1 <- ttklabel(frame1, text=txt)
  txt <- "Interpolated-grid spacing along the y-axis"
  frame1.lab.2.1 <- ttklabel(frame1, text=txt)

  frame1.ent.1.2 <- ttkentry(frame1, width=15, textvariable=grid.dx.var)
  frame1.ent.2.2 <- ttkentry(frame1, width=15, textvariable=grid.dy.var)

  tkgrid(frame1.lab.1.1, frame1.ent.1.2, pady=c(20, 4))
  tkgrid(frame1.lab.2.1, frame1.ent.2.2, pady=c(0, 10))

  tkgrid.configure(frame1.lab.1.1, frame1.lab.2.1, sticky="w", padx=c(0, 2))
  tkgrid.configure(frame1.ent.1.2, frame1.ent.2.2, sticky="we")

  tkgrid.columnconfigure(frame1, 1, weight=1, minsize=20)

  tkpack(frame1, fill="both", expand="yes", padx=30)

  # Frame 2, MBA input parameters

  frame2 <- ttklabelframe(tt, relief="flat", borderwidth=10, padding=0,
                          text="Multilevel B-spline approximation")

  txt <- "Initial size of the spline space along the x-axis"
  frame2.lab.1.1 <- ttklabel(frame2, text=txt)
  txt <- "Initial size of the spline space along the y-axis"
  frame2.lab.2.1 <- ttklabel(frame2, text=txt)
  txt <- "Number of levels in the hierarchical construction"
  frame2.lab.3.1 <- ttklabel(frame2, text=txt)

  frame2.ent.1.2 <- ttkentry(frame2, width=15, textvariable=mba.n.var)
  frame2.ent.2.2 <- ttkentry(frame2, width=15, textvariable=mba.m.var)
  frame2.ent.3.2 <- ttkentry(frame2, width=15, textvariable=mba.h.var)

  tkgrid(frame2.lab.1.1, frame2.ent.1.2, pady=c(0, 4))
  tkgrid(frame2.lab.2.1, frame2.ent.2.2, pady=c(0, 4))
  tkgrid(frame2.lab.3.1, frame2.ent.3.2)

  tkgrid.configure(frame2.lab.1.1, frame2.lab.2.1, frame2.lab.3.1, sticky="w",
                   padx=c(0, 2))
  tkgrid.configure(frame2.ent.1.2, frame2.ent.2.2, frame2.ent.3.2, sticky="we")

  tkgrid.columnconfigure(frame2, 1, weight=1, minsize=20)

  tkpack(frame2, fill="both", expand="yes", padx=10, pady=c(0, 0))

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame1.ent.1.2, "<KeyRelease>",
         function() {
           tclvalue(grid.dx.var) <- CheckEntry("numeric", tclvalue(grid.dx.var))
         })
  tkbind(frame1.ent.1.2, "<KeyRelease>",
         function() {
           tclvalue(grid.dy.var) <- CheckEntry("numeric", tclvalue(grid.dy.var))
         })

  # GUI control

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)
}
