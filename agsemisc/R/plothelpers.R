plotfit <- function(fit, ...) {
  # plots an lm object on a 2 by 3 panel
  oldpar <- par(mfrow=c(2,3))
  plot(fit, which=1:6, ...)
  invisible(par(oldpar))
}


plotf <- function(theplot, file="plot%03d", type=c("active"),
                 size=c(6,4,10,96), prepare=a.resetplotparams, ...) {
  # theplot is either an expression (a call to an old plot function)
  # or a Trellis object.
  # This function performs the plotting of either to possibly multiple devices.
  if (length(type) == 1) type = strsplit(type, ",")[[1]]
  size = or.else(size[1:4], c(6,4,10,96))  # apply default where missing
  wi = size[1]      # width in inches
  hi = size[2]      # height in inches
  points = size[3]  # base size of text (in points) for vector formats
  dpi = size[4]     # resolution (in dots per inch) for bitmaps
  ppoints = size[3]*dpi/72  # base size of text for bitmap formats
  do.plot = function(theplot, prep=TRUE, devoff=TRUE) {  # the workhorse
    if (devoff) {
      on.exit(dev.off())  # in case of an error during plotting
    }
    if (prep & is.function(prepare)) prepare()
    if (mode(theplot) == "expression")
      eval(theplot)
    else
      print(theplot)
    if (devoff) { on.exit(); dev.off() }
  }
  if (("active" %in% type) | ("current" %in% type)) {  # current device
    # plot to current device
    do.plot(theplot, prep=FALSE, devoff=FALSE)
    if ("eps" %in% type) # encapsulated postscript
      dev.copy2eps(file = paste(file, ".eps", sep=""),
                   width=wi, height=hi, pointsize=points, ...)
  }
  if ("ps" %in% type) {  # postscript
    postscript(file = paste(file, ".ps", sep=""),
               print.it=FALSE, onefile=FALSE, 
               width=wi, height=hi, pointsize=points, ...)
    do.plot(theplot)
  }
  if ("pdf" %in% type) {  # PDF
    pdf(file = paste(file, ".pdf", sep=""),
        onefile=FALSE, 
        width=wi, height=hi, pointsize=points, ...)
    do.plot(theplot)
  }
  if ("png" %in% type) {  # PNG
    png(filename = paste(file, ".png", sep=""),
        width=wi*dpi, height=hi*dpi, pointsize=ppoints, res=dpi, ...)
    do.plot(theplot)
  }
  if (("jpeg" %in% type) | ("jpg" %in% type)) {  # JPEG
    jpeg(filename = paste(file, ".jpg", sep=""),
         width=wi*dpi, height=hi*dpi, pointsize=ppoints, res=dpi, ...)
    do.plot(theplot)
  }
  if ("wmf" %in% type) {  # Windows Metafile
    if (exists("win.metafile")) {
      do.call("win.metafile", list(filename = paste(file, ".wmf", sep=""),
                                   width=wi, height=hi, pointsize=points, ...))
      do.plot(theplot)
    }
    else {
      warning("plotf: ignoring type 'wmf' (not supported on this machine)")
    }
  }
  invisible()  # return nothing
}


prepanel.0 <- function(x,y) {
  # forces that zero is included on the axes
  valid <- !is.na(x) & !is.na(y)  # ignore pairs for which x or y is NA
  result = list()
  if (is.numeric(x)) result$xlim = range(x[valid], 0)
  if (is.numeric(y)) result$ylim = range(y[valid], 0)
  result
}


a.resetplotparams = function() {
  # we use a colored-on-white color scheme. Should be turned into a theme.
  #-- normal graphics: background and margins
  par(bg="white", mar=c(4,4,1,0.2)+0.2)
  #-- Lattice colors and lines
  bg = trellis.par.get("background")
  bg$col = "white"
  trellis.par.set("background", bg)
  b.r = trellis.par.get("reference.line")
  b.r$col = "gray"
  b.r$lwd = 1
  b.r$lty = "dotted"
  trellis.par.set("reference.line", b.r)
  r.l = trellis.par.get("regression.line") # nonstandard extension for panel.xy
  if (is.null(r.l)) r.l = list()
  r.l$col.lm   = "blue";    r.l$lty.lm   = 1; r.l$lwd.lm   = 2
  r.l$col.conf = "blue";    r.l$lty.conf = 3; r.l$lwd.conf = 1
  r.l$col.pred = "blue";    r.l$lty.pred = 3; r.l$lwd.pred = 1
  r.l$col.lqs  = "red";     r.l$lty.lqs  = 2; r.l$lwd.lqs  = 2
  r.l$col.loess= "violet";  r.l$lty.loess= 1; r.l$lwd.loess= 2
  trellis.par.set("regression.line", r.l)
  a.l = trellis.par.get("ab.line") # nonstandard extension for panel.xy
  if (is.null(a.l)) a.l = list()
  a.l$col  = "darkblue"
  a.l$lty  = 1
  a.l$lwd  = 2
  trellis.par.set("ab.line", a.l)
  #-- Lattice boxplot box and whiskers
  b.r = trellis.par.get("box.rectangle")
  b.r$col = "blue"
  b.r$fill = "yellow"
  b.r$lwd = 2
  b.r$lty = 1
  trellis.par.set("box.rectangle", b.r)
  b.u = trellis.par.get("box.umbrella")
  b.u$col = "blue"
  b.u$lwd = 1
  b.u$lty = 1
  trellis.par.set("box.umbrella", b.u)
  p.l = trellis.par.get("plot.line")
  p.l$col = "black"
  trellis.par.set("plot.line", p.l)
  #-- Lattice plot symbols
  s.s = trellis.par.get("superpose.symbol")
  s.s$pch = c(1, 2, 5, 6, 0, 3, 4)
  s.s$col = c("red","blue","darkgreen","orange","violet","black","green")
  s.s$cex = 0.7
  trellis.par.set("superpose.symbol", s.s)
  s.l = trellis.par.get("superpose.line")
  s.l$col = c("red","blue","darkgreen","orange","violet","black","green")
  s.l$lwd = 1
  s.l$lty = c("solid","dashed","dotted","dotdash","longdash","twodash","AAAAAADB")
  trellis.par.set("superpose.line", s.l)
  p.s = trellis.par.get("plot.symbol")
  p.s$pch = s.s$pch[1]
  p.s$col = s.s$col[1]
  p.s$cex = s.s$cex
  trellis.par.set("plot.symbol", p.s)
}
