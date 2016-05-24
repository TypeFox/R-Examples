# A GUI for selecting a graphic color.

ChooseColor <- function(col, parent=NULL) {

  ## Additional functions

  # Save color and quit
  SaveColor <- function() {
    col.hex <- Txt2hex(tclvalue(col.var))
    if (col.hex == "")
      col.hex <- NA
    rtn.col <<- col.hex
    tclvalue(tt.done.var) <- 1
  }

  # Draw polygon on color sample
  DrawSamplePolygon <- function(fill) {
    tcl(frame0.cvs.1, "delete", "col")
    pts <- .Tcl.args(c(0, 0, dx, 0, dx, dy, 0, dy))
    tkcreate(frame0.cvs.1, "polygon", pts, fill=fill, outline="", tag="col")
  }

  # Draw polygon on color chart
  DrawChartPolygon <- function(i, j, fill, outline, tag) {
    x1 <- j * dx - dx
    y1 <- i * dy - dy
    x2 <- j * dx
    y2 <- i * dy
    pts <- .Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2) - 0.5)
    tkcreate(frame1.cvs, "polygon", pts, fill=fill, outline=outline, tag=tag)
  }

  # Build color chart
  BuildColorChart <- function() {
    for (i in seq_len(m)) {
      for (j in seq_len(n)) {
        DrawChartPolygon(i, j, fill=d[i, j], outline="black", tag="")
      }
    }
  }

  # Update color ramp
  UpdateColorRamp <- function(col.hex) {
    col.ramp <<- colorRampPalette(c("#FFFFFF", col.hex, "#000000"),
                                  space="Lab")(n + 4)[3:(n + 2)]
    tcl(frame2.cvs, "delete", "ramp")
    dx <- (w - 1) / n
    x2 <- 1
    y1 <- 1
    y2 <- dy / 2
    for (i in col.ramp) {
      x1 <- x2
      x2 <- x1 + dx
      pts <- .Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2))
      tkcreate(frame2.cvs, "polygon", pts, fill=i, tag="ramp")
    }
  }

  # Change color

  ChangeColor <- function(col.hex, is.hsva=FALSE, is.color=FALSE) {
    col.hex <- substr(col.hex, 1, 7)
    if (is.na(col.hex) || col.hex == "")
      col.hex <- "#000000"

    tcl(frame1.cvs, "delete", "browse")
    if (col.hex %in% d) {
      ij <- which(d == col.hex, arr.ind=TRUE)[1, ]
      i <- ij[1]
      j <- ij[2]
      DrawChartPolygon(i, j, fill="", outline="white", tag="browse")
    }

    DrawSamplePolygon(col.hex)
    UpdateColorRamp(col.hex)

    if (!is.hsva) {
      col.rgb <- col2rgb(col.hex, alpha=FALSE)
      col.hsv <- rgb2hsv(col.rgb, maxColorValue=255)
      nh <<- col.hsv[1]
      ns <<- col.hsv[2]
      nv <<- col.hsv[3]
      tclvalue(h.scl.var) <- nh
      tclvalue(s.scl.var) <- ns
      tclvalue(v.scl.var) <- nv
      tclvalue(a.scl.var) <- na
      tclvalue(h.ent.var) <- sprintf("%.2f", nh)
      tclvalue(s.ent.var) <- sprintf("%.2f", ns)
      tclvalue(v.ent.var) <- sprintf("%.2f", nv)
      tclvalue(a.ent.var) <- sprintf("%.2f", na)
    }
    if (!is.color)
      tclvalue(col.var) <- Hsv2hex()
  }

  # Select ramp color
  SelectRampColor <- function(x) {
    i <- ceiling((as.numeric(x)) / ((w - 1) / n))
    col.hex <- col.ramp[i]
    ChangeColor(col.hex)
  }

  # Select chart color
  SelectChartColor <- function(x, y) {
    tcl(frame1.cvs, "delete", "browse")
    i <- ceiling((as.numeric(y)) / dy)
    j <- ceiling((as.numeric(x)) / dx)
    if (i == 0)
      i <- 1
    if (j == 0)
      j <- 1
    DrawChartPolygon(i, j, fill="", outline="white", tag="browse")
    ChangeColor(d[i, j])
  }

  # Coerce text string to hexadecimal color
  Txt2hex <- function(txt) {
    txt <- CheckColorStr(as.character(txt))
    if (substr(txt, 1, 1) != "#")
      txt <- paste0("#", txt)
    if (is.transparent)
      fmt <- "%08s"
    else
      fmt <- "%06s"
    fmt.txt <- gsub(" ", "0", sprintf(fmt, substr(txt, 2, nchar(txt))))
    txt <- paste0("#", fmt.txt)
    if (inherits(try(col2rgb(txt), silent=TRUE), "try-error")) {
      if (is.transparent)
        txt <- "#000000FF"
      else
        txt <- "#000000"
    }
    return(txt)
  }

  # Coerce numeric HSV values to hexadecimal color
  Hsv2hex <- function() {
    if (is.transparent)
      col.hex <- hsv(h=nh, s=ns, v=nv, alpha=na)
    else
      col.hex <- hsv(h=nh, s=ns, v=nv)
    return(col.hex)
  }

  # Check range of numeric color attribute
  CheckColorNum <- function(...) {
    num <- suppressWarnings(as.numeric(...))
    if (is.na(num) || num < 0) {
      num <- 0
    } else if (num > 1) {
      num <- 1
    }
    return(num)
  }

  # Check hexadecimal character string
  CheckColorStr <- function(txt) {
    txt <- as.character(txt)
    if (is.transparent)
      txt <- substring(txt, 1, 9)
    else
      txt <- substring(txt, 1, 7)
    sep.txt <- strsplit(txt, "")[[1]]
    idxs <- which(vapply(sep.txt, function(i) i %in% hex.digits, TRUE))
    txt <- paste(sep.txt[idxs], collapse="")
    return(txt)
  }

  # Update based on change in scales

  ScaleH <- function(...) {
    nh <<- as.numeric(...)
    tclvalue(h.scl.var) <- nh
    tclvalue(h.ent.var) <- sprintf("%.2f", nh)
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  ScaleS <- function(...) {
    ns <<- as.numeric(...)
    tclvalue(s.scl.var) <- ns
    tclvalue(s.ent.var) <- sprintf("%.2f", ns)
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  ScaleV <- function(...) {
    nv <<- as.numeric(...)
    tclvalue(v.scl.var) <- nv
    tclvalue(v.ent.var) <- sprintf("%.2f", nv)
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  ScaleA <- function(...) {
    na <<- as.numeric(...)
    tclvalue(a.ent.var) <- sprintf("%.2f", na)
    tclvalue(col.var) <- Hsv2hex()
  }

  # Update based on change in numeric color attributes

  EntryH <- function() {
    txt <- tclvalue(h.ent.var)
    nh <<- CheckColorNum(txt)
    tclvalue(h.scl.var) <- nh
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  EntryS <- function() {
    txt <- tclvalue(s.ent.var)
    ns <<- CheckColorNum(txt)
    tclvalue(s.scl.var) <- ns
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  EntryV <- function() {
    txt <- tclvalue(v.ent.var)
    nv <<- CheckColorNum(txt)
    tclvalue(v.scl.var) <- nv
    ChangeColor(Hsv2hex(), is.hsva=TRUE)
  }

  EntryA <- function() {
    txt <- tclvalue(a.ent.var)
    na <<- CheckColorNum(txt)
    tclvalue(a.scl.var) <- na
    tclvalue(col.var) <- Hsv2hex()
  }

  # Toggle transparency check-box

  ToggleTransparency <- function() {
    is.transparent <<- as.logical(as.integer(tclvalue(trans.var)))
    if (is.transparent) {
      tkconfigure(frame3.lab.4.1, state="noraml")
      tcl(frame3.scl.4.2, "state", "!disabled")
      tkconfigure(frame3.ent.4.3, state="noraml")
    } else {
      tkconfigure(frame3.lab.4.1, state="disabled")
      tcl(frame3.scl.4.2, "state", "disabled")
      tkconfigure(frame3.ent.4.3, state="disabled")
    }
    tclvalue(col.var) <- Hsv2hex()
  }

  # Edit color entry

  EditColorEntry <- function() {
    txt <- CheckColorStr(tclvalue(col.var))
    tclvalue(col.var) <- txt
    col.hex <- Txt2hex(txt)
    ChangeColor(col.hex, is.color=TRUE)
  }

  ## Main program

  if (!requireNamespace("colorspace", quietly=TRUE))
    stop()

  # Color chart information
  m <- 12
  dx <- dy <- 20
  d1 <- cbind(colorspace::rainbow_hcl(m), colorspace::heat_hcl(m),
              colorspace::terrain_hcl(m),
              rev(gray.colors(m, start=0.1, end=0.9, gamma=1.0)))
  d2 <- c("#000000", "#000033", "#000066", "#000099", "#0000CC", "#0000FF",
          "#990000", "#990033", "#990066", "#990099", "#9900CC", "#9900FF",
          "#003300", "#003333", "#003366", "#003399", "#0033CC", "#0033FF",
          "#993300", "#993333", "#993366", "#993399", "#9933CC", "#9933FF",
          "#006600", "#006633", "#006666", "#006699", "#0066CC", "#0066FF",
          "#996600", "#996633", "#996666", "#996699", "#9966CC", "#9966FF",
          "#009900", "#009933", "#009966", "#009999", "#0099CC", "#0099FF",
          "#999900", "#999933", "#999966", "#999999", "#9999CC", "#9999FF",
          "#00CC00", "#00CC33", "#00CC66", "#00CC99", "#00CCCC", "#00CCFF",
          "#99CC00", "#99CC33", "#99CC66", "#99CC99", "#99CCCC", "#99CCFF",
          "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC", "#00FFFF",
          "#99FF00", "#99FF33", "#99FF66", "#99FF99", "#99FFCC", "#99FFFF",
          "#330000", "#330033", "#330066", "#330099", "#3300CC", "#3300FF",
          "#CC0000", "#CC0033", "#CC0066", "#CC0099", "#CC00CC", "#CC00FF",
          "#333300", "#333333", "#333366", "#333399", "#3333CC", "#3333FF",
          "#CC3300", "#CC3333", "#CC3366", "#CC3399", "#CC33CC", "#CC33FF",
          "#336600", "#336633", "#336666", "#336699", "#3366CC", "#3366FF",
          "#CC6600", "#CC6633", "#CC6666", "#CC6699", "#CC66CC", "#CC66FF",
          "#339900", "#339933", "#339966", "#339999", "#3399CC", "#3399FF",
          "#CC9900", "#CC9933", "#CC9966", "#CC9999", "#CC99CC", "#CC99FF",
          "#33CC00", "#33CC33", "#33CC66", "#33CC99", "#33CCCC", "#33CCFF",
          "#CCCC00", "#CCCC33", "#CCCC66", "#CCCC99", "#CCCCCC", "#CCCCFF",
          "#33FF00", "#33FF33", "#33FF66", "#33FF99", "#33FFCC", "#33FFFF",
          "#CCFF00", "#CCFF33", "#CCFF66", "#CCFF99", "#CCFFCC", "#CCFFFF",
          "#660000", "#660033", "#660066", "#660099", "#6600CC", "#6600FF",
          "#FF0000", "#FF0033", "#FF0066", "#FF0099", "#FF00CC", "#FF00FF",
          "#663300", "#663333", "#663366", "#663399", "#6633CC", "#6633FF",
          "#FF3300", "#FF3333", "#FF3366", "#FF3399", "#FF33CC", "#FF33FF",
          "#666600", "#666633", "#666666", "#666699", "#6666CC", "#6666FF",
          "#FF6600", "#FF6633", "#FF6666", "#FF6699", "#FF66CC", "#FF66FF",
          "#669900", "#669933", "#669966", "#669999", "#6699CC", "#6699FF",
          "#FF9900", "#FF9933", "#FF9966", "#FF9999", "#FF99CC", "#FF99FF",
          "#66CC00", "#66CC33", "#66CC66", "#66CC99", "#66CCCC", "#66CCFF",
          "#FFCC00", "#FFCC33", "#FFCC66", "#FFCC99", "#FFCCCC", "#FFCCFF",
          "#66FF00", "#66FF33", "#66FF66", "#66FF99", "#66FFCC", "#66FFFF",
          "#FFFF00", "#FFFF33", "#FFFF66", "#FFFF99", "#FFFFCC", "#FFFFFF")
  d2 <- matrix(d2, nrow=m, ncol=length(d2) / m)
  d <- cbind(d1, d2)
  n <- ncol(d)
  w <- dx * n
  h <- dy * m

  # All possible digits in color character string
  hex.digits <- list(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                     "a", "b", "c", "d", "e", "f",
                     "A", "B", "C", "D", "E", "F", "#")

  # Initialize return color
  rtn.col <- NULL

  # Initialize color ramp palette
  col.ramp <- NULL

  # Account for improper color argument
  if (missing(col) || !inherits(col, "character"))
    col <- "#000000"

 # Transparency status
  if (nchar(col) == 9)
    is.transparent <- TRUE
  else
    is.transparent <- FALSE

  # Color is intially required to be in hexadecimal format
  col.hex <- Txt2hex(col)

  # Initialize hue, saturation, value, and alpha color components
  col.rgb <- col2rgb(col.hex, alpha=TRUE)
  col.hsv <- rgb2hsv(col.rgb[seq_len(3)], maxColorValue=255)
  nh <- col.hsv[1]
  ns <- col.hsv[2]
  nv <- col.hsv[3]
  na <- col.rgb[4] / 255

  # Assign additional variables linked to Tk widgets

  col.var <- tclVar(col.hex)

  trans.var <- tclVar(is.transparent)

  h.scl.var <- tclVar(nh)
  s.scl.var <- tclVar(ns)
  v.scl.var <- tclVar(nv)
  a.scl.var <- tclVar(na)
  h.ent.var <- tclVar(sprintf("%.2f", nh))
  s.ent.var <- tclVar(sprintf("%.2f", ns))
  v.ent.var <- tclVar(sprintf("%.2f", nv))
  a.ent.var <- tclVar(sprintf("%.2f", na))

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
  tkwm.resizable(tt, 0, 0)
  tktitle(tt) <- "Choose Color"

  # Frame 0, ok and cancel buttons, and sample color

  frame0 <- ttkframe(tt, relief="flat")

  frame0.cvs.1 <- tkcanvas(frame0, relief="flat", width=dx - 1, height=dy - 1,
                           background="white", confine=TRUE, closeenough=0,
                           borderwidth=0, highlightthickness=0)
  frame0.ent.2 <- ttkentry(frame0, textvariable=col.var, width=12)
  frame0.chk.3 <- ttkcheckbutton(frame0, variable=trans.var,
                                 text="Transparency",
                                 command=ToggleTransparency)

  frame0.but.5 <- ttkbutton(frame0, width=12, text="OK", command=SaveColor)
  frame0.but.6 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() {
                              rtn.col <<- NULL
                              tclvalue(tt.done.var) <- 1
                            })

  tkgrid(frame0.cvs.1, frame0.ent.2,  frame0.chk.3, "x", frame0.but.5,
         frame0.but.6, pady=c(10, 10))
  tkgrid.columnconfigure(frame0, 3, weight=1)

  tkgrid.configure(frame0.cvs.1, sticky="w", padx=c(11, 1), pady=11)
  tkgrid.configure(frame0.ent.2, padx=c(4, 4))

  tkgrid.configure(frame0.but.5, sticky="e", padx=c(0, 4))
  tkgrid.configure(frame0.but.6, sticky="w", padx=c(0, 10), rowspan=2)

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1, color chart
  frame1 <- ttkframe(tt, relief="flat")
  frame1.cvs <- tkcanvas(frame1, relief="flat", width=w + 1, height=h + 1,
                         background="black", confine=TRUE, closeenough=0,
                         borderwidth=0, highlightthickness=0)
  tkgrid(frame1.cvs, padx=10, pady=10)
  tkpack(frame1)

  # Frame 2, color ramp
  frame2 <- ttkframe(tt, relief="flat")
  frame2.cvs <- tkcanvas(frame2, relief="flat", width=w + 1, height=dy / 2 + 1,
                         background="black", confine=TRUE, closeenough=0,
                         borderwidth=0, highlightthickness=0)

  tkgrid(frame2.cvs, padx=10, pady=c(0, 10))
  tkpack(frame2)

  # Frame 3, color attriubte sliders

  frame3 <- ttkframe(tt, relief="flat")

  frame3.lab.1.1 <- ttklabel(frame3, text="H", justify="center")
  frame3.lab.2.1 <- ttklabel(frame3, text="S", justify="center")
  frame3.lab.3.1 <- ttklabel(frame3, text="V", justify="center")
  frame3.lab.4.1 <- ttklabel(frame3, text="A", justify="center")

  frame3.scl.1.2 <- tkwidget(frame3, "ttk::scale", from=0, to=1,
                             orient="horizontal", value=nh, variable=h.scl.var,
                             command=function(...) ScaleH(...))
  frame3.scl.2.2 <- tkwidget(frame3, "ttk::scale", from=0, to=1,
                             orient="horizontal", value=ns, variable=s.scl.var,
                             command=function(...) ScaleS(...))
  frame3.scl.3.2 <- tkwidget(frame3, "ttk::scale", from=0, to=1,
                             orient="horizontal", value=nv, variable=v.scl.var,
                             command=function(...) ScaleV(...))
  frame3.scl.4.2 <- tkwidget(frame3, "ttk::scale", from=0, to=1,
                             orient="horizontal", value=na, variable=a.scl.var,
                             command=function(...) ScaleA(...))

  frame3.ent.1.3 <- ttkentry(frame3, textvariable=h.ent.var, width=4)
  frame3.ent.2.3 <- ttkentry(frame3, textvariable=s.ent.var, width=4)
  frame3.ent.3.3 <- ttkentry(frame3, textvariable=v.ent.var, width=4)
  frame3.ent.4.3 <- ttkentry(frame3, textvariable=a.ent.var, width=4)

  ToggleTransparency()

  tkgrid(frame3.lab.1.1, frame3.scl.1.2, frame3.ent.1.3, pady=c(0, 5))
  tkgrid(frame3.lab.2.1, frame3.scl.2.2, frame3.ent.2.3, pady=c(0, 5))
  tkgrid(frame3.lab.3.1, frame3.scl.3.2, frame3.ent.3.3, pady=c(0, 5))
  tkgrid(frame3.lab.4.1, frame3.scl.4.2, frame3.ent.4.3)

  tkgrid.configure(frame3.lab.1.1, frame3.lab.2.1, frame3.lab.3.1,
                   frame3.lab.4.1, padx=c(10, 2))
  tkgrid.configure(frame3.scl.1.2, frame3.scl.2.2, frame3.scl.3.2,
                   frame3.scl.4.2, sticky="we", padx=2)
  tkgrid.configure(frame3.ent.1.3, frame3.ent.2.3, frame3.ent.3.3,
                   frame3.ent.4.3, padx=c(10, 10))

  tkgrid.columnconfigure(frame3, 1, weight=1)

  tkpack(frame3, fill="x")

  # Initial commands
  BuildColorChart()
  ChangeColor(col.hex)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame0.ent.2, "<KeyRelease>", EditColorEntry)
  tkbind(frame0.ent.2, "<Return>", SaveColor)

  tkbind(frame1.cvs, "<ButtonPress>", function(x, y) SelectChartColor(x, y))
  tkbind(frame2.cvs, "<ButtonPress>", function(x) SelectRampColor(x))

  tkbind(frame3.ent.1.3, "<KeyRelease>", EntryH)
  tkbind(frame3.ent.2.3, "<KeyRelease>", EntryS)
  tkbind(frame3.ent.3.3, "<KeyRelease>", EntryV)
  tkbind(frame3.ent.4.3, "<KeyRelease>", EntryA)

  # GUI control

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  invisible(rtn.col)
}
