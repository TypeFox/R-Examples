# A GUI for choosing point symbols.

ChoosePch <- function(pch=NA, parent=NULL) {

  ## Additional functions

  # Save pch and quit
   SavePch <- function() {
     rtn.pch <<- TxtToPch(tclvalue(pch.var))
     tclvalue(tt.done.var) <- 1
   }

  # FFrame cell based on mouse selection
  MouseSelect <- function(x, y) {
    i <- ceiling(as.numeric(y) / dy)
    j <- ceiling(as.numeric(x) / dx)
    if (i == 0)
      i <- 1
    if (j == 0)
      j <- 1
    tcl(frame1.cvs, "delete", "browse")
    idx <- (n * (i - 1)) + j
    if (idx > length(pch.show)) {
      pch <- "NA"
    } else {
      DrawPolygon(i, j, fill="", outline="#CA0020", tag="browse")
      pch <- PchToTxt(pch.show[[idx]])
    }
    tclvalue(pch.var) <- pch
  }

  # Draw polygon
  DrawPolygon <- function(i, j, fill, outline, tag) {
    x1 <- j * dx - dx - 0.5
    y1 <- i * dy - dy - 0.5
    x2 <- j * dx - 0.5
    y2 <- i * dy - 0.5
    pts <- .Tcl.args(c(x1, y1, x2, y1, x2, y2, x1, y2))
    tkcreate(frame1.cvs, "polygon", pts, fill=fill, outline=outline, tag=tag)
  }

  # Draw image
  DrawImage <- function() {
    tkcreate(frame1.cvs, "image", center, anchor="center", image=img.var)
    for (i in seq_len(nrow(x.lines)))
      tkcreate(frame1.cvs, "line", x.lines[i, ], fill="#CCCCCC", tag="grid")
    for (i in seq_len(nrow(y.lines)))
      tkcreate(frame1.cvs, "line", y.lines[i, ], fill="#CCCCCC", tag="grid")
  }

  # Pch to text string
  PchToTxt <- function(pch) {
    if (is.na(pch)) {
      txt <- "NA"
    } else if (is.numeric(pch)) {
      txt <- as.character(as.integer(pch))
    } else {
      txt <- paste0("\"", pch, "\"")
    }
    return(txt)
  }

  # Text string to pch
  TxtToPch <- function(txt) {
    txt <- as.character(txt)
    if (txt %in% c("NA", "\"\"", "")) {
      pch <- NA
    } else if (suppressWarnings(!is.na(as.integer(txt)))) {
      pch <- as.integer(txt)
      if (!pch %in% pch.show)
        pch <- NA
    } else {
      txt.1 <- substr(txt, 1, 1)
      if (txt.1 == "\"")
        pch <- substr(txt, 2, 2)
      else
        pch <- txt.1
    }
    return(pch)
  }

  # Draw pch image template (included for development purposes)
  DrawPchImageTemplate <- function() {
    m <- 15
    n <- 15
    cex <- rep(1, length(pch.show))
    cex[1:26] <- rep(1.5, 26)
    postscript(file="pch.template.eps")
    op <- par(pty="s")
    plot(c(-1, m), c(-1, n), type="n", xlab="", ylab="", xaxs="i", yaxs="i")
    grid(2 * (m + 1), 2 * (n + 1), lty=1)
    for (i in seq_along(pch.show))
      points((i - 1) %% m, (n - 1) - ((i - 1) %/% n), pch=pch.show[i],
             col="black", bg="green", cex=cex[i])
    par(op)
    dev.off()
  }

  ## Main program

  if ("package:RSurvey" %in% search())
    image.path <- file.path(system.file("images", package="RSurvey"), "pch.gif")
  else
    image.path <- file.path(getwd(), "inst", "images", "pch.gif")

  pch.show <- c(0:25, 33:126, 161:255)

  w <- 375
  h <- 375
  m <- 15
  n <- 15

  dx <- w / n
  dy <- h / m
  center <- .Tcl.args(c(w / 2, h / 2))

  x.seq <- seq(dx, w - dx, dx)
  y.seq <- seq(dy, h - dy, dy)
  x.lines <- cbind(x1=0, y1=y.seq, x2=w + 1, y2=y.seq)
  y.lines <- cbind(x1=x.seq, y1=0, x2=x.seq, y2=h + 1)

  rtn.pch <- NULL

  # Assign variables linked to Tk widgets
  pch.var <- tclVar(PchToTxt(pch))
  img.var <- tclVar()
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
  tktitle(tt) <- "Choose A Graphic Symbol"

  # Create image
  tkimage.create("photo", img.var, format="GIF", file=image.path)

  # Frame 0 contains ok and cancel buttons

  frame0 <- ttkframe(tt, relief="flat")

  frame0.lab.1 <- ttklabel(frame0, text="pch =")
  frame0.ent.2 <- ttkentry(frame0, textvariable=pch.var, width=4)
  frame0.but.3 <- ttkbutton(frame0, width=12, text="OK", command=SavePch)
  frame0.but.4 <- ttkbutton(frame0, width=12, text="Cancel",
                            command=function() {
                              pch <<- NULL
                              tclvalue(tt.done.var) <- 1
                            })

  tkgrid(frame0.lab.1, frame0.ent.2, "x", frame0.but.3, frame0.but.4,
         pady=c(0, 10))
  tkgrid.columnconfigure(frame0, 2, weight=1)

  tkgrid.configure(frame0.lab.1, sticky="w", padx=c(10, 0))
  tkgrid.configure(frame0.ent.2, padx=c(2, 0))


  tkgrid.configure(frame0.but.3, sticky="e", padx=c(0, 4))
  tkgrid.configure(frame0.but.4, sticky="w", padx=c(0, 10))

  tkpack(frame0, fill="x", side="bottom", anchor="e")

  # Canvas
  frame1 <- ttkframe(tt, relief="flat")
  frame1.cvs <- tkcanvas(frame1, relief="flat", width=w + 1, height=h + 1,
                         background="white", confine=TRUE, closeenough=0,
                         borderwidth=0, highlightthickness=0)
  tkgrid(frame1.cvs, padx=10, pady=10)
  tkpack(frame1)

  # Draw image and intial selection polyon
  DrawImage()
  if (pch %in% pch.show) {
    idx <- which(pch.show %in% pch)
    i <- ceiling(idx / n)
    j <- idx - n * (i - 1L)
    DrawPolygon(i, j, fill="", outline="#CA0020", tag="browse")
  }

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame0.ent.2, "<Return>", SavePch)

  tkbind(frame1.cvs, "<ButtonPress>", function(x, y) MouseSelect(x, y))

  # GUI control

  tkfocus(tt)
  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  invisible(rtn.pch)
}
