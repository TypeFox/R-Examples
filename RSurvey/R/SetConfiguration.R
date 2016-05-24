# A GUI for specifying window geometry and universal plotting parameters.

SetConfiguration <- function(parent=NULL) {

  ## Additional functions

  UpdatePar <- function() {
    val <- as.numeric(tclvalue(width.var))
    Data("width", if (is.na(val)) NULL else val)

    val <- as.integer(tclvalue(nlevels.var))
    Data("nlevels", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(cex.pts.var))
    Data("cex.pts", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(asp.yx.var))
    Data("asp.yx", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(asp.zx.var))
    Data("asp.zx", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(vmax.var))
    Data("vmax", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(vxby.var))
    Data("vxby", if (is.na(val)) NULL else val)

    val <- as.numeric(tclvalue(vyby.var))
    Data("vyby", if (is.na(val)) NULL else val)

    Data("rkey", as.integer(tclvalue(rkey.var)))
    Data("img.contour", as.integer(tclvalue(img.contour.var)))
    Data("show.lines", as.integer(tclvalue(show.lines.var)))
    Data("show.points", as.integer(tclvalue(show.points.var)))
    Data("show.poly", as.integer(tclvalue(show.poly.var)))
    Data("vuni", as.integer(tclvalue(vuni.var)))
    Data("show.2.axes", as.integer(tclvalue(show.2.axes.var)))
    Data("minor.ticks", as.integer(tclvalue(minor.ticks.var)))
    Data("ticks.inside", as.integer(tclvalue(ticks.inside.var)))
    Data("rm.pnt.line", as.integer(tclvalue(rm.pnt.line.var)))

    tclvalue(tt.done.var) <- 1
  }

  ## Main program

  # Assign variables linked to Tk widgets

  width.var        <- tclVar()
  nlevels.var      <- tclVar()
  cex.pts.var      <- tclVar()
  asp.yx.var       <- tclVar()
  asp.zx.var       <- tclVar()
  vmax.var         <- tclVar()
  vxby.var         <- tclVar()
  vyby.var         <- tclVar()
  rkey.var         <- tclVar()
  show.poly.var    <- tclVar()
  img.contour.var  <- tclVar()
  show.lines.var   <- tclVar()
  show.points.var  <- tclVar()
  vuni.var         <- tclVar()
  show.2.axes.var  <- tclVar()
  minor.ticks.var  <- tclVar()
  ticks.inside.var <- tclVar()
  rm.pnt.line.var  <- tclVar()

  if (!is.null(Data("width")))
    tclvalue(width.var) <- Data("width")
  if (!is.null(Data("nlevels")))
    tclvalue(nlevels.var) <- Data("nlevels")
  if (!is.null(Data("cex.pts")))
    tclvalue(cex.pts.var) <- Data("cex.pts")
  if (!is.null(Data("asp.yx")))
    tclvalue(asp.yx.var) <- Data("asp.yx")
  if (!is.null(Data("asp.zx")))
    tclvalue(asp.zx.var) <- Data("asp.zx")
  if (!is.null(Data("vmax")))
    tclvalue(vmax.var) <- Data("vmax")
  if (!is.null(Data("vxby")))
    tclvalue(vxby.var) <- Data("vxby")
  if (!is.null(Data("vyby")))
    tclvalue(vyby.var) <- Data("vyby")
  if (!is.null(Data("rkey")))
    tclvalue(rkey.var) <- Data("rkey")
  if (!is.null(Data("show.poly")))
    tclvalue(show.poly.var) <- Data("show.poly")
  if (!is.null(Data("img.contour")))
    tclvalue(img.contour.var) <- Data("img.contour")
  if (!is.null(Data("show.lines")))
    tclvalue(show.lines.var) <- Data("show.lines")
  if (!is.null(Data("show.points")))
    tclvalue(show.points.var) <- Data("show.points")
  if (!is.null(Data("vuni")))
    tclvalue(vuni.var) <- Data("vuni")
  if (!is.null(Data("show.2.axes")))
    tclvalue(show.2.axes.var) <- Data("show.2.axes")
  if (!is.null(Data("minor.ticks")))
    tclvalue(minor.ticks.var) <- Data("minor.ticks")
  if (!is.null(Data("ticks.inside")))
    tclvalue(ticks.inside.var) <- Data("ticks.inside")
  if (!is.null(Data("rm.pnt.line")))
    tclvalue(rm.pnt.line.var) <- Data("rm.pnt.line")

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
  tktitle(tt) <- "Configuration"

  tkwm.resizable(tt, 1, 0)

  # Frame 0 contains ok and cancel buttons

  frame0 <- ttkframe(tt, relief="flat")
  frame0.but.2 <- ttkbutton(frame0, width=12, text="OK", command=UpdatePar)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("SetConfiguration", package="RSurvey"))
                            })
  tkgrid("x", frame0.but.2, frame0.but.3, frame0.but.4,
         sticky="se", pady=10, padx=c(4, 0))
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkgrid.configure(frame0.but.4, padx=c(4, 10))
  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Paned window

  pw <- ttkpanedwindow(tt, orient="horizontal")

  # Frame 1 contains parameters

  frame1 <- ttkframe(pw, relief="flat", borderwidth=0, padding=10)

  txt <- "Width of plotting window canvas, in inches"
  frame1.lab.1.1 <- ttklabel(frame1, text=txt)
  txt <- "Approximate number of contour levels"
  frame1.lab.2.1 <- ttklabel(frame1, text=txt)
  txt <- "Scaling for point symbols"
  frame1.lab.3.1 <- ttklabel(frame1, text=txt)
  txt <- "Horizontal aspect ratio"
  frame1.lab.4.1 <- ttklabel(frame1, text=txt)
  txt <- "Vertical aspect ratio"
  frame1.lab.5.1 <- ttklabel(frame1, text=txt)
  txt <- "Maximum arrow length, in inches"
  frame1.lab.6.1 <- ttklabel(frame1, text=txt)
  txt <- "Increment for sequence of arrows in x direction"
  frame1.lab.7.1 <- ttklabel(frame1, text=txt)
  txt <- "Increment for sequence of arrows in y direction"
  frame1.lab.8.1 <- ttklabel(frame1, text=txt)

  frame1.ent.1.2 <- ttkentry(frame1, width=8, textvariable=width.var)
  frame1.ent.2.2 <- ttkentry(frame1, width=8, textvariable=nlevels.var)
  frame1.ent.3.2 <- ttkentry(frame1, width=8, textvariable=cex.pts.var)
  frame1.ent.4.2 <- ttkentry(frame1, width=8, textvariable=asp.yx.var)
  frame1.ent.5.2 <- ttkentry(frame1, width=8, textvariable=asp.zx.var)
  frame1.ent.6.2 <- ttkentry(frame1, width=8, textvariable=vmax.var)
  frame1.ent.7.2 <- ttkentry(frame1, width=8, textvariable=vxby.var)
  frame1.ent.8.2 <- ttkentry(frame1, width=8, textvariable=vyby.var)

  tkgrid(frame1.lab.1.1, frame1.ent.1.2, pady=c(15, 4))
  tkgrid(frame1.lab.2.1, frame1.ent.2.2, pady=c(0, 4))
  tkgrid(frame1.lab.3.1, frame1.ent.3.2, pady=c(0, 4))
  tkgrid(frame1.lab.4.1, frame1.ent.4.2, pady=c(0, 4))
  tkgrid(frame1.lab.5.1, frame1.ent.5.2, pady=c(0, 4))
  tkgrid(frame1.lab.6.1, frame1.ent.6.2, pady=c(0, 4))
  tkgrid(frame1.lab.7.1, frame1.ent.7.2, pady=c(0, 4))
  tkgrid(frame1.lab.8.1, frame1.ent.8.2)

  tkgrid.configure(frame1.lab.1.1, frame1.lab.2.1, frame1.lab.3.1,
                   frame1.lab.4.1, frame1.lab.5.1, frame1.lab.6.1,
                   frame1.lab.7.1, frame1.lab.8.1,
                   sticky="w")
  tkgrid.configure(frame1.ent.1.2, frame1.ent.2.2, frame1.ent.3.2,
                   frame1.ent.4.2, frame1.ent.5.2, frame1.ent.6.2,
                   frame1.ent.7.2, frame1.ent.8.2, padx=c(2, 15),
                   sticky="we")

  tkgrid.columnconfigure(frame1, 1, weight=1, minsize=6)

  # Frame 2 contains plot features

  frame2 <- ttkframe(pw, relief="flat", borderwidth=0, padding=10)

  txt <- "Reverse legend"
  frame2.chk.01.1 <- ttkcheckbutton(frame2, text=txt, variable=rkey.var)
  txt <- "Show polygons"
  frame2.chk.02.1 <- ttkcheckbutton(frame2, text=txt, variable=show.poly.var)
  txt <- "Use image contour"
  frame2.chk.03.1 <- ttkcheckbutton(frame2, text=txt, variable=img.contour.var)
  txt <- "Show contour lines"
  frame2.chk.04.1 <- ttkcheckbutton(frame2, text=txt, variable=show.lines.var)
  txt <- "Show points on maps"
  frame2.chk.05.1 <- ttkcheckbutton(frame2, text=txt, variable=show.points.var)
  txt <- "Use uniform arrow lengths"
  frame2.chk.06.1 <- ttkcheckbutton(frame2, text=txt, variable=vuni.var)

  txt <- "Show tickmarks on second axis"
  frame2.chk.07.1 <- ttkcheckbutton(frame2, text=txt, variable=show.2.axes.var)
  txt <- "Add minor tickmarks"
  frame2.chk.08.1 <- ttkcheckbutton(frame2, text=txt, variable=minor.ticks.var)
  txt <- "Place tickmarks inside plot region"
  frame2.chk.09.1 <- ttkcheckbutton(frame2, text=txt, variable=ticks.inside.var)
  txt <- "Remove point symbol boundary line"
  frame2.chk.10.1 <- ttkcheckbutton(frame2, text=txt, variable=rm.pnt.line.var)

  tkgrid(frame2.chk.01.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.02.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.03.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.04.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.05.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.06.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.07.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.08.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.09.1, sticky="w", pady=c(0, 2))
  tkgrid(frame2.chk.10.1, sticky="w")

  # Final layout
  tkgrid(frame1, frame2, sticky="nswe")
  tkgrid.columnconfigure(pw, 0, weight=2)
  tkpack(pw, fill="x", expand=TRUE)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame1.ent.1.2, "<KeyRelease>",
         function() {
           tclvalue(width.var) <- CheckEntry("numeric", tclvalue(width.var))
         })
  tkbind(frame1.ent.2.2, "<KeyRelease>",
         function() {
           tclvalue(nlevels.var) <- CheckEntry("integer", tclvalue(nlevels.var))
         })
  tkbind(frame1.ent.3.2, "<KeyRelease>",
         function() {
           tclvalue(cex.pts.var) <- CheckEntry("numeric", tclvalue(cex.pts.var))
         })
  tkbind(frame1.ent.4.2, "<KeyRelease>",
         function() {
           tclvalue(asp.yx.var) <- CheckEntry("numeric", tclvalue(asp.yx.var))
         })
  tkbind(frame1.ent.5.2, "<KeyRelease>",
         function() {
           tclvalue(asp.zx.var) <- CheckEntry("numeric", tclvalue(asp.zx.var))
         })
  tkbind(frame1.ent.6.2, "<KeyRelease>",
         function() {
           tclvalue(vmax.var) <- CheckEntry("numeric", tclvalue(vmax.var))
         })
  tkbind(frame1.ent.7.2, "<KeyRelease>",
         function() {
           tclvalue(vxby.var) <- CheckEntry("integer", tclvalue(vxby.var))
         })
  tkbind(frame1.ent.8.2, "<KeyRelease>",
         function() {
           tclvalue(vyby.var) <- CheckEntry("integer", tclvalue(vyby.var))
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
