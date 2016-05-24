# A GUI for building histograms.

BuildHistogram <- function(d, var.names=NULL, var.default=1L, parent=NULL) {

  ## Additional functions

  # Calculate and plot histogram

  CalcHist <- function(draw.plot=TRUE) {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))

    idx <- as.integer(tcl(frame1.box.1.2, "current")) + 1L
    xlab <- as.character(var.names[idx])

    type <- as.integer(tclvalue(breaks.var))
    if (type == 1L) {
      breaks <- fun.names[as.integer(tcl(frame2.box.2.1, "current")) + 1L]
    } else if (type == 2L) {
      breaks <- as.integer(tclvalue(single.var))
    } else if (type == 3L) {
      s <- as.character(tclvalue(vector.var))
      str.split <- unlist(strsplit(s, "[[:space:]]"))
      num.split <- suppressWarnings(as.numeric(str.split))
      breaks <- num.split[!is.na(num.split)]
      if (length(breaks) == 0)
        return()
    }

    right <- as.logical(as.integer(tclvalue(right.var)))
    obs <- as.logical(as.integer(tclvalue(obs.var)))
    freq <- as.logical(as.integer(tclvalue(freq.var)))
    obj <- try(hist(d[[idx]], breaks=breaks, right=right, plot=FALSE), silent=TRUE)
    if (inherits(obj, "try-error")) {
      msg <- "Unable to construct historgram."
      tkmessageBox(icon="error", message=msg, detail=obj, title="Error",
                   type="ok", parent=tt)
      return()
    }

    if (draw.plot) {
      if (dev.cur() == dev) {
        dev.new()
        par(mar=c(5, 5, 2, 2) + 0.1)
      }
      plot(obj, col="light grey", freq=freq, main=NULL, xlab=xlab)
      if (obs)
        rug(d[[idx]], quiet=TRUE)
      if (!freq) {
        bandwidth <- as.numeric(tclvalue(bandwidth.var))
        if (!is.na(bandwidth) && bandwidth > 0) {
          dens <- density(d[[idx]], adjust=bandwidth, na.rm=TRUE)
          lines(dens, col="blue")
        }
      }
    } else {
      obj$xname <- xlab
      txt <- paste(c(capture.output(obj), ""), collapse="\n")
      EditText(txt, read.only=TRUE, win.title="Histogram Description",
               is.fixed.width.font=TRUE, parent=tt)
    }
  }

  # Adjust scale for number of cells
  AdjustScaleSingle <- function(x) {
    idx <- as.integer(tcl(frame1.box.1.2, "current")) + 1L
    breaks <- as.integer(x * (maxs[idx] - 1) + 1)
    if (breaks != as.integer(tclvalue(single.var))) {
      tclvalue(single.var) <- breaks
      PlotHist()
    }
  }

  # Adjust scale for bandwidth in density estimate
  AdjustScaleBandwidth <- function(x) {
    tclvalue(bandwidth.var) <- round(as.numeric(x), digits=1)
    PlotHist()
  }

  # Plot Histogram
  PlotHist <- function() {
    if (dev.cur() > dev)
      CalcHist()
  }

  # Toggle state on break options
  ToggleStateBreaks <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    type <- as.integer(tclvalue(breaks.var))
    states <- rep(FALSE, 4)
    states[type] <- TRUE
    s <- if (states[1]) "readonly" else "disabled"
    tkconfigure(frame2.box.2.1, state=s)
    s <- if (states[2]) "!disabled" else "disabled"
    tcl(frame2.scl.4.1, "state", s)
    s <- if (states[2]) "normal" else "disabled"
    tkconfigure(frame2.ent.4.3, state=s)
    s <- if (states[3]) "normal" else "disabled"
    tkconfigure(frame2.ent.6.1, state=s)
    if (states[1]) {
      tkfocus(frame2.box.2.1)
    } else if (states[2]) {
      tkfocus(frame2.ent.4.3)
    } else if (states[3]) {
      tkfocus(frame2.ent.6.1)
    }
    PlotHist()
  }

  # Toggle state on bandwidth
  ToggleStateBandwidth <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    freq <- as.logical(as.integer(tclvalue(freq.var)))
    s <- if (freq) "disabled" else "!disabled"
    tcl(frame4.scl.2.1, "state", s)
    s <- if (freq) "disabled" else "normal"
    tkconfigure(frame4.ent.2.3, state=s)
    if (freq)
      tkfocus(frame4.ent.2.3)
    PlotHist()
  }

  ## Main program

  # Check input arguments

  if (inherits(d, c("data.frame", "matrix"))) {
    d <- as.list(d)
  } else if (!is.list(d) && is.null(dim(d))) {
    d <- list(d)
  }
  if (!is.list(d) || length(d) == 0L)
    stop()

  if (!is.character(var.names) || length(var.names) != length(d)) {
    var.names <- names(d)
    if (is.null(var.names))
      var.names <- paste0("Unknown (", length(d), ")")
  }

  d <- lapply(d, function(i) try(as.numeric(i), silent=TRUE))

  Fun <- function(i) !(inherits(i, "try-error") | all(is.na(i)))
  is.num <- vapply(d, Fun, TRUE)
  d <- d[is.num]
  if (length(d) == 0L) {
    msg <- "None of the variables can be converted to class 'numeric'."
    tkmessageBox(icon="error", message=msg, title="Error", type="ok", parent=tt)
    return()
  }
  var.names <- var.names[is.num]

  if (inherits(var.default, c("character", "integer"))) {
    if (is.character(var.default))
      var.default <- which(var.default[1] == var.names)
    if (length(var.default) != 1L || !var.default %in% seq_len(length(d)))
      var.default <- 1L
  } else {
    var.default <- 1L
  }

  # Set limits and default value
  maxs <- vapply(d, function(i) length(unique(i)), 0L)
  maxs[maxs > 100] <- 100
  maxs[maxs <  10] <-  10
  defs <- vapply(d, function(i) length(hist(i, plot=FALSE)$breaks), 0L)
  xdef <- (defs[var.default] - 1L) / (maxs[var.default] - 1L)

  # Initialize device
  dev <- dev.cur()

  # Assign the variables linked to Tk widgets

  fun.names <- c("sturges", "scott", "freedman-diaconis")

  breaks.var <- tclVar(1L)
  single.var <- tclVar(defs[var.default])
  scale.sgl.var <- tclVar(xdef)
  vector.var <- tclVar()
  right.var <- tclVar(TRUE)
  obs.var <- tclVar(FALSE)
  freq.var <- tclVar(TRUE)
  bandwidth.var <- tclVar(1)
  scale.bw.var <- tclVar(1)

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
  tktitle(tt) <- "Build Histogram"
  tkwm.resizable(tt, 1, 0)

  # Frame 0

  frame0 <- ttkframe(tt, relief="flat")

  frame0.but.1 <- ttkbutton(frame0, width=12, text="View",
                            command=function() CalcHist(draw.plot=FALSE))
  frame0.but.2 <- ttkbutton(frame0, width=12, text="Plot",
                            command=function() CalcHist())
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Close",
                            command=function() tclvalue(tt.done.var) <- 1)
  frame0.but.5 <- ttkbutton(frame0, width=12, text="Help",
                            command=function() {
                              print(help("BuildHistogram", package="RSurvey"))
                            })
  tkgrid(frame0.but.1, frame0.but.2, "x", frame0.but.4, frame0.but.5, pady=10)
  tkgrid.configure(frame0.but.1, padx=c(10, 4))
  tkgrid.configure(frame0.but.5, padx=c(4, 10))

  tkgrid.columnconfigure(frame0, 2, weight=1, minsize=15)
  tkpack(frame0, fill="x", expand=TRUE, side="bottom", anchor="e")

  # Frame 1

  frame1 <- ttkframe(tt, relief="flat")

  frame1.lab.1.1 <- ttklabel(frame1, text="Variable")
  frame1.box.1.2 <- ttkcombobox(frame1, state="readonly")
  tkgrid(frame1.lab.1.1, frame1.box.1.2, pady=c(10, 0))

  tkgrid.configure(frame1.lab.1.1, sticky="e", padx=c(10, 2))
  tkgrid.configure(frame1.box.1.2, sticky="we", padx=c(0, 10))

  tkconfigure(frame1.box.1.2, value=var.names)
  tcl(frame1.box.1.2, "current", var.default - 1L)

  tkgrid.columnconfigure(frame1, 1, weight=1, minsize=25)
  tkpack(frame1, fill="x", expand=TRUE, padx=10, pady=5)

  # Frame 2

  frame2 <- ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Breaks")

  txt <- "A function to compute the number of cells for the histogram"
  frame2.rbt.1.1 <- ttkradiobutton(frame2, variable=breaks.var, value=1L,
                                  text=txt, command=ToggleStateBreaks)
  frame2.box.2.1 <- ttkcombobox(frame2, state="readonly")
  tkconfigure(frame2.box.2.1, value=fun.names)
  tcl(frame2.box.2.1, "current", 0)

  txt <- "A single number giving the suggested number of cells"
  frame2.rbt.3.1 <- ttkradiobutton(frame2, variable=breaks.var, value=2L,
                                  text=txt, command=ToggleStateBreaks)
  frame2.scl.4.1 <- tkwidget(frame2, "ttk::scale", from=0, to=1,
                             orient="horizontal", variable=scale.sgl.var,
                             command=function(...) {
                               AdjustScaleSingle(x=as.numeric(...))
                             })
  frame2.ent.4.3 <- ttkentry(frame2, width=4, textvariable=single.var)

  txt <- "A vector giving the breakpoints between cells"
  frame2.rbt.5.1 <- ttkradiobutton(frame2, variable=breaks.var, value=3L,
                                  text=txt, command=ToggleStateBreaks)
  frame2.ent.6.1 <- ttkentry(frame2, width=15, textvariable=vector.var)

  tkgrid(frame2.rbt.1.1, sticky="w", columnspan=3)
  tkgrid(frame2.box.2.1, padx=c(20, 0), pady=c(0, 10), sticky="we",
         columnspan=3)
  tkgrid(frame2.rbt.3.1, sticky="w", columnspan=3)
  tkgrid(frame2.scl.4.1, "x", frame2.ent.4.3, pady=c(0, 5), sticky="we")
  tkgrid(frame2.rbt.5.1, sticky="w", columnspan=3)
  tkgrid(frame2.ent.6.1, padx=c(20, 0), pady=c(0, 5), sticky="we",
         columnspan=3)

  tkgrid.configure(frame2.scl.4.1, columnspan=2, padx=c(20, 4))
  tkgrid.columnconfigure(frame2, 1, weight=1, minsize=50)
  tkpack(frame2, fill="x", expand=TRUE, padx=10, pady=5)

  # Frame 3
  frame3 <- ttkframe(tt, relief="flat")
  txt <- "The histogram cells are right-closed (left-open) intervals"
  frame3.chk.1.1 <- ttkcheckbutton(frame3, text=txt, variable=right.var,
                                   command=PlotHist)
  txt <- "Show individual observations"
  frame3.chk.2.1 <- ttkcheckbutton(frame3, text=txt, variable=obs.var,
                                   command=PlotHist)
  tkgrid(frame3.chk.1.1, sticky="w")
  tkgrid(frame3.chk.2.1, sticky="w")
  tkpack(frame3, fill="x", expand=TRUE, padx=20, pady=5)

  # Frame 4

  frame4 <- ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Axis scaling and density estimate")
  frame4.rbt.1.1 <- ttkradiobutton(frame4, variable=freq.var, value=TRUE,
                                   text="Frequences",
                                   command=ToggleStateBandwidth)
  frame4.rbt.1.2 <- ttkradiobutton(frame4, variable=freq.var, value=FALSE,
                                   text="Probability densities (adjust bandwidth)",
                                   command=ToggleStateBandwidth)

  frame4.scl.2.1 <- tkwidget(frame4, "ttk::scale", from=0, to=2,
                             orient="horizontal", variable=scale.bw.var,
                             command=function(...) {
                               AdjustScaleBandwidth(x=as.numeric(...))
                             })
  frame4.ent.2.3 <- ttkentry(frame4, width=4, textvariable=bandwidth.var)

  tkgrid(frame4.rbt.1.1, frame4.rbt.1.2, "x")
  tkgrid(frame4.scl.2.1, "x", frame4.ent.2.3, pady=c(0, 5), sticky="we")

  tkgrid.configure(frame4.rbt.1.1, padx=c(0, 10))
  tkgrid.configure(frame4.scl.2.1, columnspan=2, padx=c(20, 4))

  tkpack(frame4, fill="x", expand=TRUE, padx=10, pady=5)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(frame1.box.1.2, "<<ComboboxSelected>>",
         function() {
           idx <- as.integer(tcl(frame1.box.1.2, "current")) + 1L
           x <- as.numeric(tclvalue(scale.sgl.var))
           tclvalue(single.var) <- as.integer(x * (maxs[idx] - 1) + 1)
           PlotHist()
         })
  tkbind(frame2.box.2.1, "<<ComboboxSelected>>", PlotHist)

  tkbind(frame2.ent.4.3, "<Return>",
         function() {
           idx <- as.integer(tcl(frame1.box.1.2, "current")) + 1L
           ent <- as.integer(CheckEntry("integer", tclvalue(single.var)))
           if (is.na(ent)) {
             ent <- 1L
           } else if (ent > maxs[idx]) {
             ent <- maxs[idx]
           }
           tclvalue(single.var) <- ent
           tclvalue(scale.sgl.var) <- (ent - 1) / (maxs[idx] - 1)
           PlotHist()
         })

  tkbind(frame4.ent.2.3, "<Return>",
         function() {
           ent <- as.numeric(CheckEntry("numeric", tclvalue(bandwidth.var)))
           if (is.na(ent)) {
             ent <- 1L
           } else if (ent < 0) {
             ent <- 0
           } else if (ent > 2) {
             ent <- 2
           }
           tclvalue(bandwidth.var) <- ent
           tclvalue(scale.bw.var) <- ent
           PlotHist()
         })

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  # GUI control

  ToggleStateBreaks()
  ToggleStateBandwidth()

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)

  if (!is.null(parent))
    tkfocus(parent)

  tclServiceMode(TRUE)
  invisible()
}
