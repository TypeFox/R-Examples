choose_palette <- function(pal = diverge_hcl, n = 7L, parent = NULL) {
# A GUI for selecting a color palette.

  # Additional functions (subroutines)
  
  # Choose a file interactively
  ChooseFile <- function(cmd, win.title, initialfile=NULL, 
                         defaultextension=NULL) {
    filetypes <- "{{R Source Files} {.R}} {{All files} {*}}"
    if (cmd == "Open") {
      args <- list("tk_getOpenFile")
    } else {
      args <- list("tk_getSaveFile")
      if (defaultextension == ".txt") 
        filetypes <- "{{Text Files} {.txt}} {{All files} {*}}"
    }
    args[["title"]] <- win.title
    args[["parent"]] <- tt
    args[["initialdir"]] <- initialdir
    args[["filetypes"]] <- filetypes
    
    if (!is.null(initialfile))
      args[["initialfile"]] <- initialfile
    if (!is.null(defaultextension))
      args[["defaultextension"]] <- defaultextension
    
    f <- tcltk::tclvalue(do.call(tcltk::tcl, args))
    if (!nzchar(f))
      return()
    initialdir <<- dirname(f)
    f
  }

  # Open palette from file
  OpenPaletteFromFile <- function() {
    f <- ChooseFile(cmd="Open", win.title="Open Palette File")
    if (is.null(f))
      return()
    pal <- dget(file=f)
    ConvertPaletteToAttributes(pal)
    AssignAttributesToWidgets()
    UpdateDataType()
  }

  # Save palette to file
  SavePaletteToFile <- function() {
    f <- ChooseFile(cmd="Save As", win.title="Save Palette As",
                    initialfile="color_palette", defaultextension=".R")
    if (is.null(f))
      return()
    pal <- GetPalette(h1, h2, c1, c2, l1, l2, p1, p2, fixup)
    dput(pal, file=f)
  }

  # Save colors to file
  SaveColorsToFile <- function(type) {
    pal <- GetPalette(h1, h2, c1, c2, l1, l2, p1, p2, fixup)
    
    cols <- try(hex2RGB(pal(n)), silent=TRUE)
    if (inherits(cols, "try-error")) {
      msg <- "Palette results in invaild hexadecimal colors."
      tcltk::tkmessageBox(icon="error", message=msg, title="Color Error",
                   parent=tt)
      return()
    }
    
    f <- ChooseFile(cmd="Save As", win.title="Save Colors As",
                    initialfile=paste("colors_", type, sep=""),
                    defaultextension=".txt")
    if (is.null(f))
      return()
    
    if (type == "HEX") {
      writehex(cols, file=f)
    } else {
      if (type == "sRGB") {
        cols <- as(cols, "sRGB")@coords
      } else if (type == "HSV") {
        cols <- as(cols, "HSV")@coords
      } else if (type == "HCL") {
        cols <- as(cols, "polarLUV")@coords
      } else if (type == "CMYK") {
        cols <- as(cols, "RGB")@coords
        red   <- cols[, "R"]
        green <- cols[, "G"]
        blue  <- cols[, "B"]
        black <- sapply(1:n, function(i) min(c(1 - red[i], 1 - green[i],
                                               1 - blue[i])))
        cyan <- (1 - red - black) / (1 - black)
        magenta <- (1 - green - black) / (1 - black)
        yellow <- (1 - blue - black) / (1 - black)
        cols <- as.matrix(as.data.frame(list(C=cyan, M=black, Y=yellow,
                                             K=black)))
      }
      write.table(cols, file=f, quote=FALSE, row.names=FALSE, sep="\t")
    }
  }

  # Save palette and quit
  SavePalette <- function() {
    pal.rtn <<- GetPalette(h1, h2, c1, c2, l1, l2, p1, p2, fixup)
    tcltk::tclvalue(tt.done.var) <- 1
  }

  # Scale change
  ScaleChange <- function(x, v, x.ent.var) {
    if (x == get(v))
      return
    assign(v, x, inherits=TRUE)
    fmt <- ifelse(v %in% c("p1", "p2"), "%.1f", "%.0f")
    tcltk::tclvalue(x.ent.var) <- sprintf(fmt, x)
    DrawPalette(v == "n")
  }

  # Entry change
  EntryChange <- function(v, x.lim, x.ent.var, x.scl.var) {
    x <- suppressWarnings(as.numeric(tcltk::tclvalue(x.ent.var)))
    if (is.na(x))
      return()
    if (x < x.lim[1]) {
      tcltk::tclvalue(x.ent.var) <- x.lim[1]
      x <- x.lim[1]
    } else if (x > x.lim[2]) {
      tcltk::tclvalue(x.ent.var) <- x.lim[2]
      x <- x.lim[2]
    }
    assign(v, x, inherits=TRUE)
    tcltk::tclvalue(x.scl.var) <- x
    DrawPalette(v == "n")
  }

  # Get color palette as function of n
  GetPalette <- function(h1, h2, c1, c2, l1, l2, p1, p2, fixup) {
    fixup <- as.logical(fixup)
    type <- as.character(tcltk::tclvalue(nature.var))
    if (type == "Qualitative") {
      f <- rainbow_hcl
      formals(f) <- eval(substitute(alist(n=, c=d1, l=d2, start=d3, end=d4,
                                          fixup=d5, gamma=NULL, alpha=1, ...=),
                                    list(d1=c1, d2=l1, d3=h1, d4=h2, d5=fixup)))
    } else if (type == "Sequential (single hue)") {
      f <- sequential_hcl
      formals(f) <- eval(substitute(alist(n=, h=d1, c.=d2, l=d3, power=d4,
                                          fixup=d5, gamma=NULL, alpha=1, ...=),
                                    list(d1=h1, d2=c(c1, c2), d3=c(l1, l2),
                                         d4=p1, d5=fixup)))
    } else if (type == "Sequential (multiple hues)") {
      f <- heat_hcl
      formals(f) <- eval(substitute(alist(n=, h=d1, c.=d2, l=d3, power=d4,
                                          fixup=d5, gamma=NULL, alpha=1, ...=),
                                    list(d1=c(h1, h2), d2=c(c1, c2),
                                         d3=c(l1, l2), d4=c(p1, p2), d5=fixup)))
    } else if (type == "Diverging") {
      f <- diverge_hcl
      formals(f) <- eval(substitute(alist(n=, h=d1, c=d2, l=d3, power=d4,
                                          fixup=d5, gamma=NULL, alpha=1, ...=),
                                    list(d1=c(h1, h2), d2=c1, d3=c(l1, l2),
                                         d4=p1, d5=fixup)))
    }
    f
  }

  # Draw palette
  DrawPalette <- function(is.n=FALSE) {
    pal <- GetPalette(h1, h2, c1, c2, l1, l2, p1, p2, fixup)
    if (!is.n)
      tcltk::tcl(frame2.cvs, "delete", "browse")
    tcltk::tcl(frame7.cvs, "delete", "pal")
    pal.cols <- pal(n)
    pal.cols[is.na(pal.cols)] <- "#FFFFFF"
    if (as.logical(as.integer(tcltk::tclvalue(desaturation.var))))
      pal.cols <- desaturate(pal.cols)
    if (as.logical(as.integer(tcltk::tclvalue(colorblind.var)))) {
      type <- as.character(tcltk::tclvalue(colorblind.type.var))
      pal.cols <- dichromat::dichromat(pal.cols, type=type)
    }
    dx <- (cvs.width - 1) / n
    x2 <- 1
    y1 <- 1
    y2 <- cvs.height
    for (i in pal.cols) {
      x1 <- x2
      x2 <- x1 + dx
      pts <- tcltk::.Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2))
      tcltk::tkcreate(frame7.cvs, "polygon", pts, fill=i, tag="pal")
    }
    RegenExample(pal.cols)
  }

  # Update data type
  UpdateDataType <- function() {
    type <- as.character(tcltk::tclvalue(nature.var))
    if (type == "Qualitative") {
      is.normal <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
      default.pals <<- qual.pals
    } else if (type == "Sequential (single hue)") {
      is.normal <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
      default.pals <<- seqs.pals
    } else if (type == "Sequential (multiple hues)") {
      is.normal <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
      default.pals <<- seqm.pals
    } else if (type == "Diverging") {
      is.normal <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
      default.pals <<- dive.pals
    }

    # Default palettes
    tcltk::tcl(frame2.cvs, "delete", "default")
    x1 <- 10
    for (i in 1:length(default.pals)) {
      pal <- do.call(GetPalette, args=as.list(default.pals[[i]]))
      pal.cols <- pal(5)
      pal.cols[is.na(pal.cols)] <- "#FFFFFF"
      y2 <- 10
      for (j in pal.cols) {
        x2 <- x1 + 20
        y1 <- y2
        y2 <- y1 + 10
        pts <- tcltk::.Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2))
        tcltk::tkcreate(frame2.cvs, "polygon", pts, fill=j, tag="default")
      }
      x1 <- x1 + 30
    }

    s <- ifelse(is.normal, "normal", "disabled")
    tcltk::tkconfigure(frame3.lab.2.1, state=s[1])
    tcltk::tkconfigure(frame3.lab.4.1, state=s[2])
    tcltk::tkconfigure(frame3.lab.6.1, state=s[3])
    tcltk::tkconfigure(frame3.lab.7.1, state=s[4])
    tcltk::tkconfigure(frame3.lab.8.1, state=s[5])
    tcltk::tkconfigure(frame3.ent.2.3, state=s[1])
    tcltk::tkconfigure(frame3.ent.4.3, state=s[2])
    tcltk::tkconfigure(frame3.ent.6.3, state=s[3])
    tcltk::tkconfigure(frame3.ent.7.3, state=s[4])
    tcltk::tkconfigure(frame3.ent.8.3, state=s[5])
    s <- ifelse(is.normal, "!disabled", "disabled")
    tcltk::tcl(frame3.scl.2.2, "state", s[1])
    tcltk::tcl(frame3.scl.4.2, "state", s[2])
    tcltk::tcl(frame3.scl.6.2, "state", s[3])
    tcltk::tcl(frame3.scl.7.2, "state", s[4])
    tcltk::tcl(frame3.scl.8.2, "state", s[5])

    DrawPalette()
  }

  # Select default palette
  SelectDefaultPalette <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (is.na(x) | is.na(y))
      return()
    y1 <- 5
    y2 <- 65
    if (y < y1 | y > y2)
      return()
    max.x <- length(default.pals) * 30 + 10
    if (x < 5 | x > max.x)
      return()
    x.seq <- seq(5, max.x, by=30)
    i <- findInterval(x, x.seq, rightmost.closed=TRUE)
    x1 <- x.seq[i]
    x2 <- x.seq[i + 1]
    for (j in 1:length(vars)) {
      val <- default.pals[[i]][j]
      if (is.na(val))
        val <- 0
      assign(vars[j], val, inherits=TRUE)
    }
    AssignAttributesToWidgets()
    DrawPalette()
    pts <- tcltk::.Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2) - 0.5)
    tcltk::tkcreate(frame2.cvs, "polygon", pts, fill="", outline="black", tag="browse")
  }

  # Convert palette to attributes
  ConvertPaletteToAttributes <- function(pal) {
    pal.attributes <- NULL
    if (inherits(pal, "function")) {
      what <- c("numeric", "integer")
      q.args <- c("c", "l", "start", "end")
      d.args <- c("h", "c",  "l", "power")
      s.args <- c("h", "c.", "l", "power")
      arg <- sapply(formals(pal), function(i) {if (is.call(i)) eval(i) else i})
      if (!is.null(arg$fixup) && is.logical(arg$fixup))
        fix.up <- as.integer(arg$fixup)
      else
        fix.up <- 1
      if (all(sapply(q.args, function(i) inherits(arg[[i]], what)))) {
        tcltk::tclvalue(nature.var) <- "Qualitative"
        pal.attributes <- c(arg$start, arg$end, arg$c, NA, arg$l, NA, NA, NA, fix.up)
      } else if (all(sapply(s.args, function(i) inherits(arg[[i]], what)))) {
        if (length(arg$h) == 1 && length(arg$p) == 1) {
          tcltk::tclvalue(nature.var) <- "Sequential (single hue)"
          pal.attributes <- c(arg$h, NA, arg$c., arg$l, arg$power, NA, fix.up)
        } else {
          tcltk::tclvalue(nature.var) <- "Sequential (multiple hues)"
          pal.attributes <- c(arg$h, arg$c., arg$l, arg$power, fix.up)
        }
      } else if (all(sapply(d.args, function(i) inherits(arg[[i]], what)))) {
        tcltk::tclvalue(nature.var) <- "Diverging"
        pal.attributes <- c(arg$h, arg$c, NA, arg$l, arg$power, NA, fix.up)
      }
    }
    if (is.null(pal.attributes)) {
      tcltk::tclvalue(nature.var) <- "Sequential (multiple hues)"
      pal.attributes <- seqm.pals[[4]]
    }
    for (i in 1:length(vars)) {
      if (is.na(pal.attributes[i]))
        assign(vars[i], 0, inherits=TRUE)
      else
        assign(vars[i], pal.attributes[i], inherits=TRUE)
    }
    AssignAttributesToWidgets()
  }

  # Assign attributes to widgets
  AssignAttributesToWidgets <- function() {
    tcltk::tclvalue(h1.ent.var) <- sprintf("%.0f", h1)
    tcltk::tclvalue(h2.ent.var) <- sprintf("%.0f", h2)
    tcltk::tclvalue(c1.ent.var) <- sprintf("%.0f", c1)
    tcltk::tclvalue(c2.ent.var) <- sprintf("%.0f", c2)
    tcltk::tclvalue(l1.ent.var) <- sprintf("%.0f", l1)
    tcltk::tclvalue(l2.ent.var) <- sprintf("%.0f", l2)
    tcltk::tclvalue(p1.ent.var) <- sprintf("%.1f", p1)
    tcltk::tclvalue(p2.ent.var) <- sprintf("%.1f", p2)
    tcltk::tclvalue(h1.scl.var) <- h1
    tcltk::tclvalue(h2.scl.var) <- h2
    tcltk::tclvalue(c1.scl.var) <- c1
    tcltk::tclvalue(c2.scl.var) <- c2
    tcltk::tclvalue(l1.scl.var) <- l1
    tcltk::tclvalue(l2.scl.var) <- l2
    tcltk::tclvalue(p1.scl.var) <- p1
    tcltk::tclvalue(p2.scl.var) <- p2
    tcltk::tclvalue(fixup.var)  <- fixup
  }

  # Show example plot
  ShowExample <- function() {
    if (!dev.example %in% dev.list()) {
      dev.new(width=7L, height=7L)
      dev.example <<- dev.cur()
    }
    par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
    DrawPalette(is.n=TRUE)
  }

  # Regenerate example plot
  RegenExample <- function(pal.cols) {
    if (dev.example %in% dev.list())
      dev.set(which=dev.example)
    else
      return()
    PlotExample <- get(paste("Plot", tcltk::tclvalue(example.var), sep=""))
    if (as.logical(as.integer(tcltk::tclvalue(reverse.var))))
      pal.cols <- rev(pal.cols)
    PlotExample(pal.cols)
  }
  
  # Plot map example 
  PlotMap <- function(pal.cols) {
    n <- length(pal.cols)
    plot(0, 0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
         xlim=c(-88.5, -78.6), ylim=c(30.2, 35.2), asp=1)
     polygon(colorspace::USSouthPolygon,
             col=pal.cols[cut(na.omit(colorspace::USSouthPolygon$z), 
             breaks=0:n / n)])
  }
  
  # Plot heatmap example
  PlotHeatmap <- function(pal.cols) {
    image(volcano, col=rev(pal.cols), xaxt="n", yaxt="n", useRaster=TRUE)
  }
  
  # Plot scatter example
  PlotScatter <- function(pal.cols) {
    
    # Generate artificial data 
    if (is.null(xyhclust)) {
      set.seed(1071)
      x0 <- sin(pi * 1:60 / 30) / 5
      y0 <- cos(pi * 1:60 / 30) / 5
      xr <- c(0.1, -0.6, -0.7, -0.9,  0.4,  1.3, 1.0)
      yr <- c(0.3,  1.0,  0.1, -0.9, -0.8, -0.4, 0.6)
      dat <- data.frame(
        x=c(x0 + xr[1], x0 + xr[2], x0 + xr[3], x0 + xr[4], x0 + xr[5], 
            x0 + xr[6], x0 + xr[7]),
        y=c(y0 + yr[1], y0 + yr[2], y0 + yr[3], y0 + yr[4], y0 + yr[5], 
            y0 + yr[6], y0 + yr[7])
      )
      attr(dat, "hclust") <- hclust(dist(dat), method="ward")
      dat$xerror <- rnorm(nrow(dat), sd=runif(nrow(dat), 0.05, 0.45))
      dat$yerror <- rnorm(nrow(dat), sd=runif(nrow(dat), 0.05, 0.45))
      xyhclust <<- dat
    }
    plot(xyhclust$x + xyhclust$xerror, xyhclust$y + xyhclust$yerror,
         col="black", bg=pal.cols[cutree(attr(xyhclust, "hclust"), length(pal.cols))],
         xlab="", ylab="", axes=FALSE, pch=21, cex=1.3)
  }
  
  # Plot spine example
  PlotSpine <- function(pal.cols) {
    n <- length(pal.cols)
    
    # Rectangle dimensions
    off <- 0.015
    widths <- c(0.05, 0.1, 0.15, 0.1, 0.2, 0.08, 0.12, 0.16, 0.04)
    k <- length(widths)
    heights <- sapply(
      c(2.5, 1.2, 2.7, 1, 1.3, 0.7, 0.4, 0.2, 1.7),
      function(p) (0:n / n)^(1 / p)
    )
   
    # Rectangle coordinates
    xleft0 <- c(0, head(cumsum(widths + off), -1))
    xleft <- rep(xleft0, each=n)
    xright <- xleft + rep(widths, each=n)
    ybottom <- as.vector(heights[-(n + 1), ])
    ytop <- as.vector(heights[-1, ])
   
    # Draw rectangles, borders, and annotation
    plot(0, 0, xlim=c(0, sum(widths) + off * (k - 1)), ylim=c(0, 1),
         xaxs="i", yaxs="i", main="", xlab="", ylab="",
         type="n", axes=FALSE)
    rect(xleft, ybottom, xright, ytop, col=rep(pal.cols, k),
         border=if(n < 10) "black" else "transparent")
    if(n >= 10) rect(xleft0, 0, xleft0 + widths, 1, border="black")
  }
  
  # Plot bar example
  PlotBar <- function(pal.cols) {
    barplot(cbind(1.1 + abs(sin(0.5 + seq_along(pal.cols))) / 3,
            1.9 + abs(cos(1.1 + seq_along(pal.cols))) / 3,
            0.7 + abs(sin(1.5 + seq_along(pal.cols))) / 3,
            0.3 + abs(cos(0.8 + seq_along(pal.cols))) / 3),
            beside=TRUE, col=pal.cols, axes=FALSE)
  }

  # Plot pie example
  PlotPie <- function(pal.cols) {
    pie(0.01 + abs(sin(0.5 + seq_along(pal.cols))), labels="", col=pal.cols, 
        radius=1)
  }
  
  # Plot perspective example
  PlotPerspective <- function(pal.cols) {
    
    # Mixture of bivariate normals
    n <- 31
    x1 <- x2 <- seq(-3, 3, length=n)
    y <- outer(x1, x2, 
               function(x, y) {
                 0.5 * dnorm(x, mean=-1, sd=0.80) * dnorm(y, mean=-1, sd=0.80) +
                 0.5 * dnorm(x, mean= 1, sd=0.72) * dnorm(y, mean= 1, sd=0.72)
               })
  
    # Compute color based on density
    if (length(pal.cols) > 1) {
      facet <- cut(y[-1, -1] + y[-1, -n] + y[-n, -1] + y[-n, -n], 
                   length(pal.cols))
      cols <- rev(pal.cols)[facet]
    } else {
      cols <- pal.cols
    }
  
    # Perspective plot coding z-axis with color
    persp(x1, x2, y, col=cols, phi=28, theta=20, r=5, xlab="", ylab="", zlab="")
  }
  
  # Plot mosaic example
  PlotMosaic <- function(pal.cols) {
    if (is.null(msc.matrix)) {
      set.seed(1071)
      mat <- list()
      for (i in 1:50) {
        mat[[i]] <- matrix(runif(i * 10, min=-1, max=1), nrow=10, ncol=i)
      }
      msc.matrix <<- mat
    }
    image(msc.matrix[[n]], col=pal.cols, xaxt="n", yaxt="n")
  }
  
  # Plot lines example
  PlotLines <- function(pal.cols) {
    n <- length(pal.cols)
    plot(NULL, xlab="", ylab="", xaxt="n", yaxt="n", type="n", 
         xlim=c(0, 6), ylim=c(1.5, n + 1.5))
    s <- 2:(n + 1)
    rev.s <- rev(s)
    rev.pal.cols <- rev(pal.cols)
    lwd <- 6
    if (n > 5)
     lwd <- lwd -1
    if (n > 15)
     lwd <- lwd -1
    if (n > 25)
     lwd <- lwd -1
    segments(1 / s, s, 2 + 1 / rev.s, rev.s, pal.cols, lwd=lwd)
    segments(2 + 1 / s, s, 4 - 1 / s, s, rev.pal.cols, lwd=lwd)
    segments(4 - 1 / s, s, 6 - 1 / s, rev.s, rev.pal.cols, lwd=lwd)
  }


  # Main program
  
  # Initialize directory
  initialdir <- getwd()
  
  # Initialize return palette
  pal.rtn <- NULL

  # Initialize default palettes
  default.pals <- NULL
  
  # Initialize data for scatter plot example
  xyhclust <- NULL
  
  # Initialize data for mosaic plot example
  msc.matrix <- NULL
  
  # Flag graphics device
  dev.example <- 1

  # Set default and initial palettes
  h1 <- h2 <- c1 <- c2 <- l1 <- l2 <- p1 <- p2 <- 0
  fixup <- 1
  vars <- c("h1", "h2", "c1", "c2", "l1", "l2", "p1", "p2", "fixup")

  qual.pals <- list()
  qual.pals[[1]]  <- c(  0,  288,  35, NA, 85, NA,  NA,  NA, 1) # ColorBrewer.org: Pastel1
  qual.pals[[2]]  <- c( 10,  320,  50, NA, 80, NA,  NA,  NA, 1) # ColorBrewer.org: Set3
  qual.pals[[3]]  <- c(  0,  288,  60, NA, 70, NA,  NA,  NA, 1) # ColorBrewer.org: Set2
  qual.pals[[4]]  <- c(  0,  288,  50, NA, 60, NA,  NA,  NA, 1) # ColorBrewer.org: Dark2
  qual.pals[[5]]  <- c(  0,  300,  80, NA, 60, NA,  NA,  NA, 1) # JCF/Z: Even Darker
  qual.pals[[6]]  <- c( 90,  -30,  50, NA, 70, NA,  NA,  NA, 1) # Z+KH+PM-09, Fig.4: Warm (based on Ihaka-03)
  qual.pals[[7]]  <- c(270,  150,  50, NA, 70, NA,  NA,  NA, 1) # Z+KH+PM-09, Fig.4: Cold (based on Ihaka-03)
  qual.pals[[8]]  <- c( 60,  240,  50, NA, 70, NA,  NA,  NA, 1) # Z+KH+PM-09, Fig.4: Harmonic (based on Ihaka-03)
  qual.pals[[9]]  <- c( 30,  300,  50, NA, 70, NA,  NA,  NA, 1) # Z+KH+PM-09, Fig.4: Dynamic (based on Ihaka-03)

  seqs.pals <- list()
  seqs.pals[[1]]  <- c(  0,   NA,   0,  0, 15, 95, 1.3,  NA, 1) # ColorBrewer.org: Greys
  seqs.pals[[2]]  <- c(280,   NA,  60,  5, 20, 95, 0.7, 1.3, 1) # ColorBrewer.org: Purples
  seqs.pals[[3]]  <- c(260,   NA,  80, 10, 30, 95, 0.7, 1.3, 1) # ColorBrewer.org: Blues
  seqs.pals[[4]]  <- c(135,   NA,  50, 10, 40, 95, 0.4, 1.3, 1) # ColorBrewer.org: Greens
  seqs.pals[[5]]  <- c( 20,   NA,  80,  5, 35, 95, 0.6, 1.3, 1) # ColorBrewer.org: Oranges
  seqs.pals[[6]]  <- c( 10,   NA,  80, 10, 30, 95, 0.7, 1.3, 1) # JCF/Z: Reds
  seqs.pals[[7]]  <- c(260,   NA,  80,  0, 30, 90, 1.5,  NA, 1) # Z+KH+PM-09, Fig.5: Blues
  seqs.pals[[8]]  <- c(260,   NA,   0,  0, 30, 90, 1.5,  NA, 1) # Z+KH+PM-09, Fig.5: Light Grays
  
  seqm.pals <- list()
  seqm.pals[[1]]  <- c(300,  200,  60,  0, 25, 95, 0.7, 1.3, 1) # ColorBrewer.org: BuPu
  seqm.pals[[2]]  <- c(370,  280,  80,  5, 25, 95, 0.7, 1.3, 1) # ColorBrewer.org: PuRd
  seqm.pals[[3]]  <- c(140,   80,  40, 10, 35, 95, 0.7, 1.7, 1) # ColorBrewer.org: YlGn
  seqm.pals[[4]]  <- c(265,   80,  60, 10, 25, 95, 0.7, 2.0, 1) # ColorBrewer.org: YlGnBu
  seqm.pals[[5]]  <- c( 10,   85,  80, 10, 25, 95, 0.4, 1.3, 1) # ColorBrewer.org: YlOrRd
  seqm.pals[[6]]  <- c(  0,   90,  80, 30, 30, 90, 0.2, 2.0, 1) # JCF/Z: alternative to heat_hcl
  seqm.pals[[7]]  <- c(  0,   90, 100, 30, 50, 90, 0.2, 1.0, 1) # Z+KH+PM-09, Fig.5: heat_hcl
  seqm.pals[[8]]  <- c(130,   30,  65,  0, 45, 90, 0.5, 1.5, 1) # JCF/Z: alternative to terrain_hcl
  seqm.pals[[9]]  <- c(130,    0,  80,  0, 60, 95, 0.1, 1.0, 1) # Z+KH+PM-09, Fig.5: terrain_hcl
  seqm.pals[[10]] <- c(  0, -100,  80, 40, 40, 75, 1.0, 1.0, 1) # Z+KH+PM-09, Fig.5: Red-Blue

  dive.pals <- list()
  dive.pals[[1]]  <- c(340,  128,  45, NA, 35, 95, 0.7, 1.3, 1) # ColorBrewer.org: PiYG
  dive.pals[[2]]  <- c(300,  128,  45, NA, 30, 95, 0.7, 1.3, 1) # ColorBrewer.org: PRGn
  dive.pals[[3]]  <- c( 55,  160,  30, NA, 35, 95, 0.7, 1.3, 1) # ColorBrewer.org: BrBG
  dive.pals[[4]]  <- c( 40,  270,  45, NA, 30, 95, 0.7, 1.3, 1) # ColorBrewer.org: PuOr
  dive.pals[[5]]  <- c( 12,  265,  80, NA, 25, 95, 0.7, 1.3, 1) # ColorBrewer.org: RdBu
  dive.pals[[6]]  <- c(260,    0,  80, NA, 30, 90, 1.5,  NA, 1) # Z+KH+PM-09, Fig.6: Blue-Red (high luminance contrast)
  dive.pals[[7]]  <- c(260,    0, 100, NA, 50, 90, 1.0,  NA, 1) # Z+KH+PM-09, Fig.6: Blue-Red (medium luminance contrast)
  dive.pals[[8]]  <- c(130,   43, 100, NA, 70, 90, 1.0,  NA, 1) # Z+KH+PM-09, Fig.6: Green-Orange (low luminance contrast)
  dive.pals[[9]]  <- c(180,  330,  59, NA, 75, 95, 1.5,  NA, 1) # Z+KH+PM-09, Fig.6: Blue-Red (low luminance contrast)

  # Set limits for palette attributes
  n.lim <- c(   1,  50)
  h.lim <- c(-360, 360)
  c.lim <- c(   0, 100)
  l.lim <- c(   0, 100)
  p.lim <- c(   0,   3)

  # Set dimensions on palette canvas
  cvs.width <- 328 # 30 * 10 + 10 + 18
  cvs.height <- 25

  # Assign additional variables linked to Tk widgets
  example.var <- tcltk::tclVar()
  nature.var <- tcltk::tclVar()

  n.scl.var <- tcltk::tclVar(n)
  n.ent.var <- tcltk::tclVar(n)

  h1.scl.var <- tcltk::tclVar()
  h1.ent.var <- tcltk::tclVar()
  h2.scl.var <- tcltk::tclVar()
  h2.ent.var <- tcltk::tclVar()
  c1.scl.var <- tcltk::tclVar()
  c1.ent.var <- tcltk::tclVar()
  c2.scl.var <- tcltk::tclVar()
  c2.ent.var <- tcltk::tclVar()
  l1.scl.var <- tcltk::tclVar()
  l1.ent.var <- tcltk::tclVar()
  l2.scl.var <- tcltk::tclVar()
  l2.ent.var <- tcltk::tclVar()
  p1.scl.var <- tcltk::tclVar()
  p1.ent.var <- tcltk::tclVar()
  p2.scl.var <- tcltk::tclVar()
  p2.ent.var <- tcltk::tclVar()
  
  fixup.var <- tcltk::tclVar(fixup)
  reverse.var <- tcltk::tclVar(FALSE)
  desaturation.var <- tcltk::tclVar(FALSE)
  colorblind.var <- tcltk::tclVar(FALSE)
  colorblind.type.var <- tcltk::tclVar("deutan")

  tt.done.var <- tcltk::tclVar(0)

  # Open GUI
  tcltk::tclServiceMode(FALSE)

  tt <- tcltk::tktoplevel()
  if (!is.null(parent)) {
    tcltk::tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tcltk::tkwm.geometry(parent)), "\\+"))
    tcltk::tkwm.geometry(tt, paste("+", as.integer(geo[2]) + 25,
                            "+", as.integer(geo[3]) + 25, sep=""))
  }
  tcltk::tkwm.resizable(tt, 0, 0)
  tcltk::tktitle(tt) <- "Choose Color Palette"

  # Top file menu
  top.menu <- tcltk::tkmenu(tt, tearoff=0)
  menu.file <- tcltk::tkmenu(tt, tearoff=0)
  tcltk::tkadd(top.menu, "cascade", label="File", menu=menu.file, underline=0)

  tcltk::tkadd(menu.file, "command", label="Open palette", accelerator="Ctrl+O",
        command=OpenPaletteFromFile)
  tcltk::tkadd(menu.file, "command", label="Save palette as",
        accelerator="Shift+Ctrl+S", command=SavePaletteToFile)

  menu.file.colors <- tcltk::tkmenu(tt, tearoff=0)
  tcltk::tkadd(menu.file.colors, "command", label="HEX",
        command=function() SaveColorsToFile("HEX"))
  tcltk::tkadd(menu.file.colors, "command", label="sRGB",
        command=function()  SaveColorsToFile("sRGB"))
  tcltk::tkadd(menu.file.colors, "command", label="HSV",
        command=function()  SaveColorsToFile("HSV"))
  tcltk::tkadd(menu.file.colors, "command", label="HCL",
        command=function()  SaveColorsToFile("HCL"))
  tcltk::tkadd(menu.file.colors, "command", label="CMYK",
        command=function()  SaveColorsToFile("CMYK"))
  tcltk::tkadd(menu.file, "cascade", label="Save colors as", menu=menu.file.colors)

  tcltk::tkconfigure(tt, menu=top.menu)

  # Frame 0, ok and cancel buttons
  frame0 <- tcltk::ttkframe(tt, relief="flat")
  frame0.but.3 <- tcltk::ttkbutton(frame0, width=12, text="OK", command=SavePalette)
  frame0.but.4 <- tcltk::ttkbutton(frame0, width=12, text="Cancel",
                            command=function() {
                              pal.rtn <<- NULL
                              tcltk::tclvalue(tt.done.var) <- 1
                            })
  tcltk::tkgrid("x", frame0.but.3, frame0.but.4, pady=c(10, 10))
  tcltk::tkgrid.configure(frame0.but.3, sticky="e")
  tcltk::tkgrid.configure(frame0.but.4, sticky="w", padx=c(4, 10))
  tcltk::tkgrid.columnconfigure(frame0, 0, weight=1)

  tcltk::tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Frame 1, choose nature of data
  frame1 <- tcltk::ttkframe(tt, relief="flat")
  frame1.lab.1 <- tcltk::ttklabel(frame1, text="The nature of your data")
  frame1.box.2 <- tcltk::ttkcombobox(frame1, state="readonly", textvariable=nature.var,
                              values=c("Qualitative", "Sequential (single hue)",
                                       "Sequential (multiple hues)",
                                       "Diverging"))

  tcltk::tkgrid(frame1.lab.1, frame1.box.2, pady=c(10, 0))
  tcltk::tkgrid.configure(frame1.lab.1, padx=c(10, 2))
  tcltk::tkgrid.configure(frame1.box.2, padx=c(0, 10), sticky="we")

  tcltk::tkgrid.columnconfigure(frame1, 1, weight=1)

  tcltk::tkpack(frame1, fill="x")

  # Frame 2, default color schemes
  frame2 <- tcltk::ttklabelframe(tt, relief="flat", borderwidth=5, padding=5,
                          text="Default color schemes")
  frame2.cvs <- tcltk::tkcanvas(frame2, relief="flat", width=30 * 10 + 10, height=70,
                         background="white", confine=TRUE, closeenough=0,
                         borderwidth=0, highlightthickness=0)
  tcltk::tkgrid(frame2.cvs, sticky="we")
  tcltk::tkgrid.columnconfigure(frame2, 0, weight=1)
  tcltk::tkpack(frame2, fill="x", padx=10, pady=10)

  # Frame 3, color description
  txt <- "Palette description: Hue, Chroma, Luminance, Power"
  frame3 <- tcltk::ttklabelframe(tt, relief="flat", borderwidth=5, padding=5, text=txt)

  frame3.lab.1.1 <- tcltk::ttklabel(frame3, text="H1", width=2)
  frame3.lab.2.1 <- tcltk::ttklabel(frame3, text="H2", width=2)
  frame3.lab.3.1 <- tcltk::ttklabel(frame3, text="C1", width=2)
  frame3.lab.4.1 <- tcltk::ttklabel(frame3, text="C2", width=2)
  frame3.lab.5.1 <- tcltk::ttklabel(frame3, text="L1", width=2)
  frame3.lab.6.1 <- tcltk::ttklabel(frame3, text="L2", width=2)
  frame3.lab.7.1 <- tcltk::ttklabel(frame3, text="P1", width=2)
  frame3.lab.8.1 <- tcltk::ttklabel(frame3, text="P2", width=2)

  frame3.ent.1.3 <- tcltk::ttkentry(frame3, textvariable=h1.ent.var, width=4)
  frame3.ent.2.3 <- tcltk::ttkentry(frame3, textvariable=h2.ent.var, width=4)
  frame3.ent.3.3 <- tcltk::ttkentry(frame3, textvariable=c1.ent.var, width=4)
  frame3.ent.4.3 <- tcltk::ttkentry(frame3, textvariable=c2.ent.var, width=4)
  frame3.ent.5.3 <- tcltk::ttkentry(frame3, textvariable=l1.ent.var, width=4)
  frame3.ent.6.3 <- tcltk::ttkentry(frame3, textvariable=l2.ent.var, width=4)
  frame3.ent.7.3 <- tcltk::ttkentry(frame3, textvariable=p1.ent.var, width=4)
  frame3.ent.8.3 <- tcltk::ttkentry(frame3, textvariable=p2.ent.var, width=4)

  frame3.scl.1.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=h.lim[1], to=h.lim[2],
                             orient="horizontal", value=h1, variable=h1.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="h1",
                                           x.ent.var=h1.ent.var)
                             })
  frame3.scl.2.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=h.lim[1], to=h.lim[2],
                             orient="horizontal", value=h2, variable=h2.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="h2",
                                           x.ent.var=h2.ent.var)
                             })
  frame3.scl.3.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=c.lim[1], to=c.lim[2],
                             orient="horizontal", value=c1, variable=c1.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="c1",
                                           x.ent.var=c1.ent.var)
                             })
  frame3.scl.4.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=c.lim[1], to=c.lim[2],
                             orient="horizontal", value=c2, variable=c2.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="c2",
                                           x.ent.var=c2.ent.var)
                             })
  frame3.scl.5.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=l.lim[1], to=l.lim[2],
                             orient="horizontal", value=l1, variable=l1.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="l1",
                                           x.ent.var=l1.ent.var)
                             })
  frame3.scl.6.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=l.lim[1], to=l.lim[2],
                             orient="horizontal", value=l2, variable=l2.scl.var,
                             command=function(...) {
                               ScaleChange(x=round(as.numeric(...)), v="l2",
                                           x.ent.var=l2.ent.var)
                             })
  frame3.scl.7.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=p.lim[1], to=p.lim[2],
                             orient="horizontal", value=p1, variable=p1.scl.var,
                             command=function(...) {
                               ScaleChange(x=as.numeric(...), v="p1",
                                           x.ent.var=p1.ent.var)
                             })
  frame3.scl.8.2 <- tcltk::tkwidget(frame3, "ttk::scale", from=p.lim[1], to=p.lim[2],
                             orient="horizontal", value=p2, variable=p2.scl.var,
                             command=function(...) {
                               ScaleChange(x=as.numeric(...), v="p2",
                                           x.ent.var=p2.ent.var)
                             })

  tcltk::tkgrid(frame3.lab.1.1, frame3.scl.1.2, frame3.ent.1.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.2.1, frame3.scl.2.2, frame3.ent.2.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.3.1, frame3.scl.3.2, frame3.ent.3.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.4.1, frame3.scl.4.2, frame3.ent.4.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.5.1, frame3.scl.5.2, frame3.ent.5.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.6.1, frame3.scl.6.2, frame3.ent.6.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.7.1, frame3.scl.7.2, frame3.ent.7.3, pady=c(0, 5))
  tcltk::tkgrid(frame3.lab.8.1, frame3.scl.8.2, frame3.ent.8.3)

  tcltk::tkgrid.configure(frame3.scl.1.2, frame3.scl.2.2, frame3.scl.3.2,
                   frame3.scl.4.2, frame3.scl.5.2, frame3.scl.6.2,
                   frame3.scl.7.2, frame3.scl.8.2,
                   sticky="we", padx=c(4, 10))

  tcltk::tkgrid.columnconfigure(frame3, 1, weight=1)

  tcltk::tkpack(frame3, fill="x", padx=10, pady=0)
  
  # Frame 4, color palette fixup
  frame4 <- tcltk::ttkframe(tt, relief="flat")
  txt <- "Correct all colors to valid RGB color model values"
  frame4.chk.1 <- tcltk::ttkcheckbutton(frame4, text=txt, variable=fixup.var,
                                 command=function() {
                                   fixup <<- as.integer(tcltk::tclvalue(fixup.var))
                                   DrawPalette(is.n=TRUE)
                                 })
  tcltk::tkgrid.configure(frame4.chk.1, padx=c(12, 0), pady=c(2, 0))
  tcltk::tkpack(frame4, fill="x")
  
  # Frame 5, number of colors in palette
  txt <- "Number of colors in palette"
  frame5 <- tcltk::ttklabelframe(tt, relief="flat", borderwidth=5, padding=5, text=txt)

  frame5.lab.1 <- tcltk::ttklabel(frame5, text="n", width=2)
  frame5.ent.3 <- tcltk::ttkentry(frame5, textvariable=n.ent.var, width=4)
  frame5.scl.2 <- tcltk::tkwidget(frame5, "ttk::scale", from=n.lim[1], to=n.lim[2],
                           orient="horizontal", value=n, variable=n.scl.var,
                           command=function(...) {
                             ScaleChange(x=round(as.numeric(...)), v="n",
                                         x.ent.var=n.ent.var)
                           })

  tcltk::tkgrid(frame5.lab.1, frame5.scl.2, frame5.ent.3)
  tcltk::tkgrid.configure(frame5.scl.2, sticky="we", padx=c(4, 10))
  tcltk::tkgrid.columnconfigure(frame5, 1, weight=1)

  tcltk::tkpack(frame5, fill="x", padx=10, pady=10)

  # Frame 6, example plots and reverse colors
  frame6 <- tcltk::ttklabelframe(tt, relief="flat", borderwidth=5, padding=5, 
                          text="Show example")
  frame6.lab.1 <- tcltk::ttklabel(frame6, text="Plot type")
  frame6.box.2 <- tcltk::ttkcombobox(frame6, state="readonly", 
                              textvariable=example.var,
                              values=c("Map", "Heatmap", "Scatter", "Spine", 
                                       "Bar", "Pie", "Perspective", "Mosaic", 
                                       "Lines"))
  frame6.chk.3 <- tcltk::ttkcheckbutton(frame6, text="Reverse colors", 
                                 variable=reverse.var, command=ShowExample)
  tcltk::tkgrid(frame6.lab.1, frame6.box.2, frame6.chk.3)
  tcltk::tkgrid.configure(frame6.box.2, padx=c(2, 10), sticky="we")
  tcltk::tkgrid.columnconfigure(frame6, 1, weight=1)
  tcltk::tkpack(frame6, fill="x", padx=10, pady=0)
  
 # Frame 7, color palette and robustness checks
  frame7 <- tcltk::ttkframe(tt, relief="flat")
  frame7.cvs <- tcltk::tkcanvas(frame7, relief="flat",
                         width=cvs.width + 1, height=cvs.height + 1,
                         background="black", confine=TRUE, closeenough=0,
                         borderwidth=0, highlightthickness=0)
  tcltk::tkgrid(frame7.cvs, padx=10, pady=c(12,0))

  frame7.chk.1 <- tcltk::ttkcheckbutton(frame7, text="Desaturation",
                                 variable=desaturation.var,
                                 command=function() DrawPalette(is.n=TRUE))

  is.pkg <- requireNamespace("dichromat", quietly=FALSE)
  if (is.pkg) {
    frame7.chk.2 <- tcltk::ttkcheckbutton(frame7, text="Color blindness:",
                                   variable=colorblind.var,
                                   command=function() DrawPalette(is.n=TRUE))
    frame7.rb.3 <- tcltk::ttkradiobutton(frame7, variable=colorblind.type.var,
                                  value="deutan", text="deutan",
                                  command=function() DrawPalette(is.n=TRUE))
    frame7.rb.4 <- tcltk::ttkradiobutton(frame7, variable=colorblind.type.var,
                                  value="protan", text="protan",
                                  command=function() DrawPalette(is.n=TRUE))
    frame7.rb.5 <- tcltk::ttkradiobutton(frame7, variable=colorblind.type.var,
                                  value="tritan", text="tritan",
                                  command=function() DrawPalette(is.n=TRUE))
    ## tritan support in dichromat starting from > 1.2-4
    if(compareVersion(getNamespaceVersion("dichromat"), "1.2-4") > 0) {
      tcltk::tkgrid(frame7.chk.1, frame7.chk.2, frame7.rb.3, frame7.rb.4, frame7.rb.5, "x",
           pady=c(2, 0), sticky="w")    
    } else {
      tcltk::tkgrid(frame7.chk.1, frame7.chk.2, frame7.rb.3, frame7.rb.4, "x",
           pady=c(2, 0), sticky="w")
    }
    tcltk::tkgrid.configure(frame7.chk.2, padx=c(7, 0))
    tcltk::tkgrid.configure(frame7.cvs, columnspan=5)
    tcltk::tkgrid.columnconfigure(frame7, 4, weight=1)
  } else {
    tcltk::tkgrid(frame7.chk.1, "x", pady=c(2, 0), sticky="w")
    tcltk::tkgrid.configure(frame7.cvs, columnspan=2)
    tcltk::tkgrid.columnconfigure(frame7, 1, weight=1)
  }
  tcltk::tkgrid.configure(frame7.chk.1, padx=c(10, 0))
  tcltk::tkpack(frame7, fill="x")
  
  # Initial commands
  ConvertPaletteToAttributes(pal)
  UpdateDataType()

  # Bind events
  tcltk::tclServiceMode(TRUE)

  tcltk::tkbind(tt, "<Control-o>", OpenPaletteFromFile)
  tcltk::tkbind(tt, "<Shift-Control-S>", SavePaletteToFile)
  
  tcltk::tkbind(frame1.box.2, "<<ComboboxSelected>>", UpdateDataType)
  tcltk::tkbind(frame6.box.2, "<<ComboboxSelected>>", ShowExample)

  tcltk::tkbind(frame2.cvs, "<ButtonPress>", function(x, y) SelectDefaultPalette(x, y))

  tcltk::tkbind(frame3.ent.1.3, "<KeyRelease>",
         function() EntryChange("h1", h.lim, h1.ent.var, h1.scl.var))
  tcltk::tkbind(frame3.ent.2.3, "<KeyRelease>",
         function() EntryChange("h2", h.lim, h2.ent.var, h2.scl.var))
  tcltk::tkbind(frame3.ent.3.3, "<KeyRelease>",
         function() EntryChange("c1", c.lim, c1.ent.var, c1.scl.var))
  tcltk::tkbind(frame3.ent.4.3, "<KeyRelease>",
         function() EntryChange("c2", c.lim, c2.ent.var, c2.scl.var))
  tcltk::tkbind(frame3.ent.5.3, "<KeyRelease>",
         function() EntryChange("l1", l.lim, l1.ent.var, l1.scl.var))
  tcltk::tkbind(frame3.ent.6.3, "<KeyRelease>",
         function() EntryChange("l2", l.lim, l2.ent.var, l2.scl.var))
  tcltk::tkbind(frame3.ent.7.3, "<KeyRelease>",
         function() EntryChange("p1", p.lim, p1.ent.var, p1.scl.var))
  tcltk::tkbind(frame3.ent.8.3, "<KeyRelease>",
         function() EntryChange("p2", p.lim, p2.ent.var, p2.scl.var))

  tcltk::tkbind(frame5.ent.3, "<KeyRelease>",
         function() EntryChange("n", n.lim, n.ent.var, n.scl.var))

  tcltk::tkbind(tt, "<Destroy>", function() tcltk::tclvalue(tt.done.var) <- 1)

  # GUI control
  tcltk::tkfocus(tt)
  tcltk::tkgrab(tt)

  tcltk::tkwait.variable(tt.done.var)

  tcltk::tclServiceMode(FALSE)
  tcltk::tkgrab.release(tt)
  tcltk::tkdestroy(tt)
  tcltk::tclServiceMode(TRUE)

  if (dev.example %in% dev.list())
    dev.off(which=dev.example)

  invisible(pal.rtn)
}
